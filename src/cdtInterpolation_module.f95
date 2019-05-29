
  !!!!!!!!!!!!!!!!! Ordiranry kriging !!!!!!!!!!!!!!!!!

  subroutine kriging_interp(xyStn, valStn, nStn, xyGrd, nGrd, nmin, nmax, spheric,&
                            vgm, nug, sill, rg, outM)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: nStn, nGrd, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyStn(nStn, 2), valStn(nStn), xyGrd(nGrd, 2)
    integer, intent(in) :: vgm
    real(rp), intent(in) :: nug, sill, rg
    real(rp), intent(out) :: outM(nGrd, 4)
    real(rp) :: out(2), outNaN(2), Fillval, x
    integer :: i, ret

    x = -1.0_rp
    Fillval = sqrt(x)
    outNaN(1:2) = Fillval
    outM(:, 1:2) = xyGrd
    do i = 1, nGrd
      call okriging_pixel(xyGrd(i, :), xyStn, valStn, nStn, nmin, nmax, spheric,&
                          vgm, nug, sill, rg, out, ret)
      if(ret == 0) then
        outM(i, 3:4) = out
      else
        outM(i, 3:4) = outNaN
      end if
    end do
  end subroutine kriging_interp

  subroutine okriging_pixel(xyp, xym, values, n, nmin, nmax, spheric,&
                            vgm, nug, sill, rg, out, ret)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2), values(n)
    integer, intent(in) :: vgm
    real(rp), intent(in) :: nug, sill, rg
    real(rp), intent(out) :: out(2)
    integer, intent(out) :: ret
    real(rp) :: dst(n), rmax
    integer :: order(n), npts
    real(rp) :: zs0(n), xys0(n, 2)
    real(rp), dimension(:), allocatable :: zs, xdst, dG, vG, predv, Wk
    real(rp), dimension(:, :), allocatable :: xys1, dS, vS, ivS, predm

    call distance_pixel(xyp, xym, n, nmin, nmax, spheric, dst, rmax, npts, order)

    allocate(zs(npts), xdst(npts))

    zs0 = values(order)
    zs = zs0(1:npts)
    xdst = dst(1:npts)

    if(all(xdst < 1.e-10)) then
      out(1) = sum(zs)/npts
      out(2) = 0
      deallocate(zs, xdst)
      return
    end if

    if(abs(maxval(zs) - minval(zs)) < 1.e-10) then
      out(1) = zs(1)
      out(2) = 0
      deallocate(zs, xdst)
      return
    end if

    allocate(xys1(npts, 2), dS(npts, npts), dG(npts))
    allocate(predm(npts, npts), vS(npts + 1, npts + 1), ivS(npts + 1, npts + 1))
    allocate(predv(npts), vG(npts + 1), Wk(npts + 1))

    xys0 = xym(order, :)
    xys1 = xys0(1:npts, :)

    call distance_matrix(xys1, xys1, npts, npts, spheric, dS)
    call distance_vector(xyp, xys1, npts, spheric, dG)
    call predict_vgm_matrix(dS, npts, vgm, nug, sill, rg, predm)
    vS(1:npts, 1:npts) = sill - predm
    vS(:, npts + 1) = 1.0_rp
    vS(npts + 1, :) = 1.0_rp
    vS(npts + 1, npts + 1) = 0.0_rp
    call inverse_matrix(vS, npts + 1, ivS, ret)

    if(ret > 0) return

    call predict_vgm_vector(dG, npts, vgm, nug, sill, rg, predv)
    vG(1:npts) = sill - predv
    vG(npts + 1) = 1.0_rp

    Wk = matmul(ivS, vG)
    out(1) = sum(zs * Wk(1:npts))
    out(2) = sill - sum(Wk * vG)
    deallocate(zs, xdst, xys1, dS, dG, predm, vS, ivS, predv, vG, Wk)
  end subroutine okriging_pixel

  subroutine predict_vgm_vector(x, n, vgm, nug, sill, rg, pred)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    integer, intent(in) :: vgm
    real(rp), intent(in) :: x(n), nug, sill, rg
    real(rp), intent(out) :: pred(n)

    if(vgm == 1) then
      ! "Gau"
      pred = nug + (sill - nug) * (1.0_rp - exp(-3.0_rp * (x/rg)**2))
    else if(vgm == 2) then
      ! "Exp"
      pred = nug + (sill - nug) * (1.0_rp - exp(-3.0_rp * x/rg))
    else if(vgm == 3) then
      ! "Sph"
      where(x <= rg)
        pred = nug + (sill - nug) * (1.5_rp * (x/rg) - 0.5_rp * (x/rg)**3)
      elsewhere
        pred = sill
      end where
    else if(vgm == 4) then
      ! "Pen"
      where(x <= rg)
        pred = nug + (sill - nug) * (1.875_rp * (x/rg) - 1.25_rp * (x/rg)**3 + 0.375_rp * (x/rg)**5)
      elsewhere
        pred = sill
      end where
    else
      stop 1
    end if
  end subroutine predict_vgm_vector

  subroutine predict_vgm_matrix(x, n, vgm, nug, sill, rg, pred)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    integer, intent(in) :: vgm
    real(rp), intent(in) :: x(n, n), nug, sill, rg
    real(rp), intent(out) :: pred(n, n)

    if(vgm ==1) then
      pred = nug + (sill - nug) * (1.0_rp - exp(-3.0_rp * (x/rg)**2))
    else if(vgm == 2) then
      pred = nug + (sill - nug) * (1.0_rp - exp(-3.0_rp * x/rg))
    else if(vgm == 3) then
      where(x <= rg)
        pred = nug + (sill - nug) * (1.5_rp * (x/rg) - 0.5_rp * (x/rg)**3)
      elsewhere
        pred = sill
      end where
    else if(vgm == 4) then
      where(x <= rg)
        pred = nug + (sill - nug) * (1.875_rp * (x/rg) - 1.25_rp * (x/rg)**3 + 0.375_rp * (x/rg)**5)
      elsewhere
        pred = sill
      end where
    else
      stop 1
    end if
  end subroutine predict_vgm_matrix

  !!!!!!!!!!!!!!!!! Spheremap Interpolation !!!!!!!!!!!!!!!!!

  subroutine spheremap_interp(xyStn, valStn, nStn, xyGrd, nGrd, nmin, nmax,&
                              spheric, outM)
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: nStn, nGrd, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyStn(nStn, 2), valStn(nStn), xyGrd(nGrd, 2)
    real(rp), intent(out) :: outM(nGrd, 3)
    real(rp) :: out
    integer :: i

    outM(:, 1:2) = xyGrd
    do i = 1, nGrd
      call spheremap_pixel(xyGrd(i, :), xyStn, valStn, nStn, nmin, nmax, spheric, out)
      outM(i, 3) = out
    end do
  end subroutine spheremap_interp

  subroutine spheremap_pixel(xyp, xym, values, n, nmin, nmax, spheric, out)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2), values(n)
    real(rp), intent(out) :: out

    real(rp) :: xdst(n), rmax
    integer :: order(n), npts, i, k, ninf
    real(rp) :: zs0(n), xys0(n, 2), v, cosd_jk, sind_jk
    integer, allocatable :: Skinf(:), mk(:)
    real(rp), dimension(:), allocatable :: zs, dst, Sk, sinf, Wk, Tk, dlon_jk, dlat_jk
    real(rp), dimension(:), allocatable :: dz_lonk, dz_latk, dst_kl, dlon_lk, dlat_lk
    real(rp), dimension(:), allocatable :: cos_theta_kl, cosd_kl, cosd_jl, sind_jl
    real(rp), dimension(:), allocatable :: xx_kl, yy_kl, deltaZ_k
    real(rp), dimension(:, :), allocatable :: xys1
    real(rp), parameter :: infini = huge(0.0_rp)
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

    call distance_pixel(xyp, xym, n, nmin, nmax, spheric, xdst, rmax, npts, order)

    allocate(zs(npts), dst(npts))

    dst = xdst(1:npts)
    zs0 = values(order)
    zs = zs0(1:npts)

    if(all(dst < 1.e-10)) then
      out = sum(zs)/npts
      deallocate(zs, dst)
      return
    end if

    if(abs(maxval(zs) - minval(zs)) < 1.e-10) then
      out = zs(1)
      deallocate(zs, dst)
      return
    end if

    allocate(xys1(npts, 2), Sk(npts), Skinf(npts))
    allocate(Wk(npts), Tk(npts), mk(npts - 1), cos_theta_kl(npts - 1))
    allocate(xx_kl(npts - 1), yy_kl(npts - 1), dlat_jk(npts), dlon_jk(npts), deltaZ_k(npts))
    allocate(dz_lonk(npts), dz_latk(npts), dst_kl(npts - 1), dlat_lk(npts - 1), dlon_lk(npts - 1))
    allocate(cosd_kl(npts - 1), cosd_jl(npts - 1), sind_jl(npts - 1))

    xys0 = xym(order, :)
    xys1 = xys0(1:npts, :)

    Sk(1:npts) = 0.0_rp
    where(dst <= rmax/3.0_rp) Sk = 1.0_rp/dst
    where((dst > rmax/3.0_rp) .and. (dst <= rmax))
      Sk = ((dst/rmax) - 1.0_rp)**2 * (27.0_rp/(4.0_rp * rmax))
    end where

    if(any(Sk > infini)) then
      Skinf(1:npts) = 0
      where(.not. Sk > infini) Skinf = 1
      ninf = sum(Skinf)
      allocate(sinf(ninf))
      sinf = pack(Sk, .not. Sk > infini)
      where(.not. Sk > infini)
        Sk = 0.1_rp * (sinf - minval(sinf))/(maxval(sinf) - minval(sinf))
      elsewhere
        Sk = 1.0_rp
      end where
      deallocate(sinf, Skinf)
    end if

    Tk(1:npts) = 0.0_rp
    Wk(1:npts) = 0.0_rp

    do k = 1, npts
      if(k == 1) then
        mk = (/ (i, i = 2, npts) /)
      else if(k == npts) then
        mk = (/ (i, i = 1, npts - 1) /)
      else
        mk = (/ (i, i = 1, k - 1), (i, i = k + 1, npts) /)
      end if

      if (spheric == 1) then
        call cos_spheric_distance(xys1(k, :), xys1(mk, :), npts - 1, cosd_kl)
        cosd_jk = cos(dst(k))
        cosd_jl = cos(dst(mk))
        sind_jk = sin(dst(k))
        sind_jl = sin(dst(mk))
        cos_theta_kl = (cosd_kl - cosd_jk * cosd_jl)/(sind_jk * sind_jl)
      else
        xx_kl = (xys1(k, 1) - xyp(1)) * (xys1(mk, 1) - xyp(1))
        yy_kl = (xys1(k, 2) - xyp(2)) * (xys1(mk, 2) - xyp(2))
        cos_theta_kl = (xx_kl + yy_kl)/(dst(k) * dst(mk))
      end if
      Tk(k) = sum(Sk(mk) * (1 - cos_theta_kl))
      Wk(k) = (Sk(k)**2) * (1 + Tk(k) /sum(Sk(mk)))
    end do

    if (spheric == 1) then
      dlon_jk = (xyp(1) - xys1(:, 1)) * cos(xyp(2) * pi/180.0_rp)
      dlat_jk = xyp(2) - xys1(:, 2)
    else
      dlon_jk = xyp(1) - xys1(:, 1)
      dlat_jk = xyp(2) - xys1(:, 2)
    end if

    dz_lonk(1:npts) = 0.0_rp
    dz_latk(1:npts) = 0.0_rp

    do k = 1, npts
      if(k == 1) then
        mk = (/ (i, i = 2, npts) /)
      else if(k == npts) then
        mk = (/ (i, i = 1, npts - 1) /)
      else
        mk = (/ (i, i = 1, k - 1), (i, i = k + 1, npts) /)
      end if

      call distance_vector(xys1(k, :), xys1(mk, :), npts - 1, spheric, dst_kl)
      if (spheric == 1) then
        dlon_lk = (xys1(mk, 1) - xys1(k, 1)) * cos(xys1(mk, 2) * pi/180.0_rp)
        dlat_lk = xys1(mk, 2) - xys1(k, 2)
        dst_kl = 6378.388_rp * dst_kl
      else
        dlon_lk = xys1(mk, 1) - xys1(k, 1)
        dlat_lk = xys1(mk, 2) - xys1(k, 2)
      end if

      dz_lonk(k) = sum(Wk(mk) * (zs(mk) - zs(k)) * dlon_lk * (1/dst_kl**2))/sum(Wk(mk))
      dz_latk(k) = sum(Wk(mk) * (zs(mk) - zs(k)) * dlat_lk * (1/dst_kl**2))/sum(Wk(mk))
    end do

    v = 0.1 * (maxval(zs0) - minval(zs0)) / maxval(sqrt(dz_lonk**2 + dz_latk**2))
    deltaZ_k = (dz_lonk * dlon_jk + dz_latk * dlat_jk) * (v / (v + dst))

    out = sum(Wk * (zs + deltaZ_k))/sum(Wk)

    deallocate(zs, dst, Sk, Wk, Tk, dlon_jk, dlat_jk, mk)
    deallocate(dz_lonk, dz_latk, dst_kl, dlon_lk, dlat_lk)
    deallocate(cos_theta_kl, cosd_kl, cosd_jl, sind_jl)
    deallocate(xx_kl, yy_kl, deltaZ_k, xys1)
  end subroutine spheremap_pixel

  !!!!!!!!!!!!!!!!! Inverse distance weighted (IDW)  !!!!!!!!!!!!!!!!!

  subroutine idw_interp(xyStn, valStn, nStn, xyGrd, nGrd, nmin, nmax,&
                        spheric, p, outM)
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: nStn, nGrd, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyStn(nStn, 2), valStn(nStn), xyGrd(nGrd, 2), p
    real(rp), intent(out) :: outM(nGrd, 3)
    real(rp) :: out
    integer :: i

    outM(:, 1:2) = xyGrd
    do i = 1, nGrd
      call idw_pixel(xyGrd(i, :), xyStn, valStn, nStn, nmin, nmax, spheric, p, out)
      outM(i, 3) = out
    end do
  end subroutine idw_interp

  subroutine idw_pixel(xyp, xym, values, n, nmin, nmax, spheric, p, out)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2), values(n), p
    real(rp), intent(out) :: out

    real(rp) :: xdst(n), rmax
    integer :: order(n), npts, ninf
    real(rp) :: zs0(n)
    real(rp), dimension(:), allocatable :: zs, dst, Wk, wi
    integer, allocatable :: winf(:)
    real(rp), parameter :: infini = huge(0.0_rp)

    call distance_pixel(xyp, xym, n, nmin, nmax, spheric, xdst, rmax, npts, order)

    allocate(zs(npts), dst(npts))

    dst = xdst(1:npts)
    zs0 = values(order)
    zs = zs0(1:npts)

    if(all(dst < 1.e-10)) then
      out = sum(zs)/npts
      deallocate(zs, dst)
      return
    end if

    if(abs(maxval(zs) - minval(zs)) < 1.e-10) then
      out = zs(1)
      deallocate(zs, dst)
      return
    end if

    allocate(Wk(npts), winf(npts))

    Wk = 1/dst**p

    if(any(Wk > infini)) then
      winf(1:npts) = 0
      where(.not. Wk > infini) winf = 1
      ninf = sum(winf)
      allocate(wi(ninf))
      wi = pack(Wk, .not. Wk > infini)
      where(.not. Wk > infini)
        Wk = (1.0_rp - ((maxval(wi) - wi)/(maxval(wi) - minval(wi))))/2.0_rp
      elsewhere
        Wk = 1.0_rp
      end where
      deallocate(wi, winf)
    end if

    out = sum(Wk * zs) / sum(Wk)

    deallocate(zs, dst, Wk)
  end subroutine idw_pixel

  !!!!!!!!!!!!!!!!!  Modified Shepard interpolation  !!!!!!!!!!!!!!!!!

  subroutine shepard_interp(xyStn, valStn, nStn, xyGrd, nGrd, nmin, nmax,&
                            spheric, p, outM)
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: nStn, nGrd, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyStn(nStn, 2), valStn(nStn), xyGrd(nGrd, 2), p
    real(rp), intent(out) :: outM(nGrd, 3)
    real(rp) :: out
    integer :: i

    outM(:, 1:2) = xyGrd
    do i = 1, nGrd
      call shepard_pixel(xyGrd(i, :), xyStn, valStn, nStn, nmin, nmax, spheric, p, out)
      outM(i, 3) = out
    end do
  end subroutine shepard_interp

  subroutine shepard_pixel(xyp, xym, values, n, nmin, nmax, spheric, p, out)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2), values(n), p
    real(rp), intent(out) :: out

    real(rp) :: xdst(n), rmax
    integer :: order(n), npts, ninf
    real(rp) :: zs0(n)
    real(rp), dimension(:), allocatable :: zs, dst, Wk, wi
    integer, allocatable :: winf(:)
    real(rp), parameter :: infini = huge(0.0_rp)

    call distance_pixel(xyp, xym, n, nmin, nmax, spheric, xdst, rmax, npts, order)

    allocate(zs(npts), dst(npts))

    dst = xdst(1:npts)
    zs0 = values(order)
    zs = zs0(1:npts)

    if(all(dst < 1.e-10)) then
      out = sum(zs)/npts
      deallocate(zs, dst)
      return
    end if

    if(abs(maxval(zs) - minval(zs)) < 1.e-10) then
      out = zs(1)
      deallocate(zs, dst)
      return
    end if

    allocate(Wk(npts), winf(npts))

    Wk = ((rmax - dst) / (rmax * dst))**p

    if(any(Wk > infini)) then
      winf(1:npts) = 0
      where(.not. Wk > infini) winf = 1
      ninf = sum(winf)
      allocate(wi(ninf))
      wi = pack(Wk, .not. Wk > infini)
      where(.not. Wk > infini)
        Wk = (1.0_rp - ((maxval(wi) - wi)/(maxval(wi) - minval(wi))))/2.0_rp
      elsewhere
        Wk = 1.0_rp
      end where
      deallocate(wi, winf)
    end if

    out = sum(Wk * zs) / sum(Wk)

    deallocate(zs, dst, Wk)
  end subroutine shepard_pixel

  !!!!!!!!!!!!!!!!!  Cressman weighting method  !!!!!!!!!!!!!!!!!

  subroutine cressman_interp(xyStn, valStn, nStn, xyGrd, nGrd, nmin, nmax,&
                             spheric, outM)
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: nStn, nGrd, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyStn(nStn, 2), valStn(nStn), xyGrd(nGrd, 2)
    real(rp), intent(out) :: outM(nGrd, 3)
    real(rp) :: out
    integer :: i

    outM(:, 1:2) = xyGrd
    do i = 1, nGrd
      call cressman_pixel(xyGrd(i, :), xyStn, valStn, nStn, nmin, nmax, spheric, out)
      outM(i, 3) = out
    end do
  end subroutine cressman_interp

  subroutine cressman_pixel(xyp, xym, values, n, nmin, nmax, spheric, out)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2), values(n)
    real(rp), intent(out) :: out

    real(rp) :: xdst(n), rmax
    integer :: order(n), npts
    real(rp) :: zs0(n)
    real(rp), dimension(:), allocatable :: zs, dst, Wk

    call distance_pixel(xyp, xym, n, nmin, nmax, spheric, xdst, rmax, npts, order)

    allocate(zs(npts), dst(npts))

    dst = xdst(1:npts)
    zs0 = values(order)
    zs = zs0(1:npts)

    if(all(dst < 1.e-10)) then
      out = sum(zs)/npts
      deallocate(zs, dst)
      return
    end if

    if(abs(maxval(zs) - minval(zs)) < 1.e-10) then
      out = zs(1)
      deallocate(zs, dst)
      return
    end if

    allocate(Wk(npts))

    Wk = (rmax**2 - dst**2)/(rmax**2 + dst**2)
    out = sum(Wk * zs) / sum(Wk)

    deallocate(zs, dst, Wk)
  end subroutine cressman_pixel

  !!!!!!!!!!!!!!!!!  Barnes weighting method  !!!!!!!!!!!!!!!!!

  subroutine barnes_interp(xyStn, valStn, nStn, xyGrd, nGrd, nmin, nmax,&
                           spheric, p, outM)
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: nStn, nGrd, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyStn(nStn, 2), valStn(nStn), xyGrd(nGrd, 2), p
    real(rp), intent(out) :: outM(nGrd, 3)
    real(rp) :: out
    integer :: i

    outM(:, 1:2) = xyGrd
    do i = 1, nGrd
      call barnes_pixel(xyGrd(i, :), xyStn, valStn, nStn, nmin, nmax, spheric, p, out)
      outM(i, 3) = out
    end do
  end subroutine barnes_interp

  subroutine barnes_pixel(xyp, xym, values, n, nmin, nmax, spheric, p, out)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2), values(n), p
    real(rp), intent(out) :: out

    real(rp) :: xdst(n), rmax
    integer :: order(n), npts
    real(rp) :: zs0(n)
    real(rp), dimension(:), allocatable :: zs, dst, Wk

    call distance_pixel(xyp, xym, n, nmin, nmax, spheric, xdst, rmax, npts, order)

    allocate(zs(npts), dst(npts))

    dst = xdst(1:npts)
    zs0 = values(order)
    zs = zs0(1:npts)

    if(all(dst < 1.e-10)) then
      out = sum(zs)/npts
      deallocate(zs, dst)
      return
    end if

    if(abs(maxval(zs) - minval(zs)) < 1.e-10) then
      out = zs(1)
      deallocate(zs, dst)
      return
    end if

    allocate(Wk(npts))

    Wk = exp(-(dst/(p * rmax))**2)
    out = sum(Wk * zs) / sum(Wk)

    deallocate(zs, dst, Wk)
  end subroutine barnes_pixel

  !!************************************************************************!!

  !!!!!!!!!!!!!!!!! cartesian_distance !!!!!!!!!!!!!!!!!

  subroutine cartesian_distance(xyp, xym, n, dst)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    real(rp), intent(in) :: xyp(2), xym(n, 2)
    real(rp), intent(out) :: dst(n)

    dst = sqrt((xyp(1) - xym(:, 1))**2 + (xyp(2) - xym(:, 2))**2)
  end subroutine cartesian_distance

  !!!!!!!!!!!!!!!!! cos_spheric_distance !!!!!!!!!!!!!!!!!

  subroutine cos_spheric_distance(xyp, xym, n, dst)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    real(rp), intent(in) :: xyp(2), xym(n, 2)
    real(rp), intent(out) :: dst(n)
    real(rp) :: sin_yp, cos_yp, const
    real(rp) :: sin_ym(n), cos_ym(n), d_cos_xpm(n), sin_y(n), cos_y(n)
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

    const = pi/180.0_rp
    sin_yp = sin(xyp(2) * const)
    sin_ym = sin(xym(:, 2) * const)
    cos_yp = cos(xyp(2) * const)
    cos_ym = cos(xym(:, 2) * const)
    d_cos_xpm = cos((xyp(1) - xym(:, 1)) * const)
    sin_y = sin_yp * sin_ym
    cos_y = cos_yp * cos_ym
    dst = sin_y + cos_y * d_cos_xpm
  end subroutine cos_spheric_distance

  !!!!!!!!!!!!!!!!! spheric_distance !!!!!!!!!!!!!!!!!

  subroutine spheric_distance(xyp, xym, n, dst)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    real(rp), intent(in) :: xyp(2), xym(n, 2)
    real(rp), intent(out) :: dst(n)
    real(rp) :: res(n)

    call cos_spheric_distance(xyp, xym, n, res)
    dst =  acos(res)
  end subroutine spheric_distance

  !!!!!!!!!!!!!!!!! distance_vector !!!!!!!!!!!!!!!!!

  subroutine distance_vector(xyp, xym, n, spheric, dst)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2)
    real(rp), intent(out) :: dst(n)

    if (spheric == 1) then
        call spheric_distance(xyp, xym, n, dst)
    else
        call cartesian_distance(xyp, xym, n, dst)
    end if
  end subroutine distance_vector

  !!!!!!!!!!!!!!!!! distance_matrix !!!!!!!!!!!!!!!!!

  subroutine distance_matrix(xym1, xym2, n1, n2, spheric, dst)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n1, n2
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xym1(n1, 2), xym2(n2, 2)
    real(rp), intent(out) :: dst(n2, n1)
    real(rp) :: res(n2)
    integer :: i

    do i = 1, n1
        call distance_vector(xym1(i, :), xym2, n2, spheric, res)
        dst(:, i) = res
    end do
  end subroutine distance_matrix

  !!!!!!!!!!!!!!!!! distance_pixel !!!!!!!!!!!!!!!!!

  subroutine distance_pixel(xyp, xym, n, nmin, nmax, spheric, dst, rmax, npts, order)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n, nmin, nmax
    integer, intent(in) :: spheric
    real(rp), intent(in) :: xyp(2), xym(n, 2)
    real(rp), intent(out) :: dst(n), rmax
    integer, intent(out) :: order(n), npts
    real(rp) :: rsconst
    integer :: nmean, npl(n)

    call distance_vector(xyp, xym, n, spheric, dst)
    call order_vector(dst, n, order)

    nmean = int((nmin + nmax)/2)
    if(nmean > n) nmean = n
    rsconst = sum(dst(1:nmean))/nmean

    where(dst <= rsconst)
      npl = 1
    elsewhere
      npl = 0
    endwhere
    npts = sum(npl)

    if(npts < nmin) then
      rmax = dst(nmin)
    else if(npts > nmax) then
      rmax = dst(nmax)
    else
      rmax = rsconst
    end if

    where(dst <= rmax)
      npl = 1
    elsewhere
      npl = 0
    endwhere
    npts = sum(npl)
  end subroutine distance_pixel

  !!!!!!!!!!!!!!!!! order small array !!!!!!!!!!!!!!!!!

  subroutine order_vector(vec, n, order)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    real(rp), intent(inout) :: vec(n)
    integer, intent(out) :: order(n)
    integer :: i, k, m
    real(rp) :: start, pivot
    integer :: ns, np

    order = (/ (i, i = 1, n) /)

    do i = 1, n - 1
      start = vec(i)
      ns = order(i)
      m = i + 1
      do k = m, n
        if (start <= vec(k)) cycle
        pivot = start; start = vec(k); vec(k) = pivot
        np = ns; ns = order(k); order(k) = np
      end do
      vec(i) = start
      order(i) = ns
    end do
  end subroutine order_vector

  !!!!!!!!!!!!!!!!! matrix pseudoinverse !!!!!!!!!!!!!!!!!

  subroutine inverse_matrix(a, n, inv, ret)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: n
    real(rp), intent(in) :: a(n, n)
    real(rp), intent(out) :: inv(n, n)
    integer, intent(out) :: ret
    real(rp) :: u(n, n), w(n), v(n, n), diag(n, n)
    integer :: j

    u = a
    call svdcmp(u, n, n, w, v, ret)

    if(ret > 0) return

    ! diag(1:n, 1:n) = 0.0d+00
    ! do j = 1, n
    !   if(w(j) /= 0.0d+00) then
    !     diag(j, j) = 1.0d+00/w(j)
    !   end if
    ! end do

    diag(1:n, 1:n) = 0.0_rp
    do j = 1, n
      if(w(j) > 1.e-13) then
        diag(j, j) = 1.0_rp/w(j)
      end if
    end do

    inv = matmul(v, matmul(diag, transpose(u)))
  end subroutine inverse_matrix

  !!************************************************************************!!

  !!!!!!!!!!!!!!!!! singular value decomposition !!!!!!!!!!!!!!!!!

  subroutine svdcmp(a, m, n, w, v, ret)
  !* --------------------------------------------------------------------- *
  !* Author: Jean-Pierre Moreau
  !* F90 Release 2.0 By J-P Moreau, Paris
  !* www.jpmoreau.fr
  !* --------------------------------------------------------------------- *
    implicit none

    integer, parameter :: rp = selected_real_kind(15)
    integer, intent(in) :: m, n
    real(rp), intent(inout) :: a(m,n), v(n,n), w(n)
    integer, intent(out) :: ret
    integer :: nmax
    parameter (nmax = 500)  !Maximum anticipated value of n.
    integer :: i, its, j, jj, k, l, nm 
    real(rp) :: anorm, c, f, g, h, s, scale, x, y, z, rv1(nmax)   !, pythag

    g = 0.d0  !Householder reduction to bidiagonal form.
    scale = 0.d0
    anorm = 0.d0
    ret = 0

    do i=1,n
      l=i+1
      rv1(i)=scale*g
      g=0.d0
      s=0.d0
      scale=0.d0

      if(i.le.m)then
        do k=i,m
          scale=scale+abs(a(k,i))
        end do

        if(scale.ne.0.d0)then
          do k=i,m
            a(k,i)=a(k,i)/scale
            s=s+a(k,i)*a(k,i)
          end do

          f=a(i,i)
          g=-dsign(dsqrt(s),f)
          h=f*g-s
          a(i,i)=f-g

          do j=l,n
            s=0.d0
            do k=i,m
              s=s+a(k,i)*a(k,j)
            end do

            f=s/h

            do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            end do
          end do
          do k=i,m
            a(k,i)=scale*a(k,i)
          end do
        endif
      endif

      w(i)=scale *g
      g=0.d0
      s=0.d0
      scale=0.d0

      if((i.le.m).and.(i.ne.n))then
        do k=l,n
          scale=scale+abs(a(i,k))
        end do

        if(scale.ne.0.d0)then
          do k=l,n
            a(i,k)=a(i,k)/scale
            s=s+a(i,k)*a(i,k)
          end do

          f=a(i,l)
          g=-sign(sqrt(s),f)
          h=f*g-s
          a(i,l)=f-g

          do k=l,n
            rv1(k)=a(i,k)/h
          end do

          do j=l,m
            s=0.d0
            do k=l,n
              s=s+a(j,k)*a(i,k)
            end do

            do k=l,n
              a(j,k)=a(j,k)+s*rv1(k)
            end do
          end do
          do k=l,n
            a(i,k)=scale*a(i,k)
          end do
        endif
      endif
      anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
    end do !do i=1,n


    !Accumulation of right-hand transformations.
    do i=n,1,-1
      if(i.lt.n)then
        if(g.ne.0.d0)then
          !Double division to avoid possible underflow.
          do j=l,n
            v(j,i)=(a(i,j)/a(i,l))/g
          end do

          do j=l,n
            s=0.d0
            do k=l,n
              s=s+a(i,k)*v(k,j)
            end do

            do k=l,n
              v(k,j)=v(k,j)+s*v(k,i)
            end do
          end do
        endif
        do j=l,n
          v(i,j)=0.d0
          v(j,i)=0.d0
        end do
      endif
      v(i,i)=1.d0
      g=rv1(i)
      l=i
    end do


    !Accumulation of left-hand transformations.
    do i=min(m,n),1,-1
      l=i+1
      g=w(i)

      do j=l,n
        a(i,j)=0.d0
      end do

      if(g.ne.0.d0)then
        g=1.d0/g

        do j=l,n
          s=0.d0
          do k=l,m
            s=s+a(k,i)*a(k,j)
          end do

          f=(s/a(i,i))*g

          do k=i,m
            a(k,j)=a(k,j)+f*a(k,i)
          end do
        end do

        do j=i,m
          a(j,i)=a(j,i)*g
        end do
      else
        do j= i,m
          a(j,i)=0.d0
        end do
      endif

      a(i,i)=a(i,i)+1.d0
    end do

    !Diagonalization of the bidiagonal form: Loop over
    !singular values, and over allowed iterations.
    do k=n,1,-1
      do its=1,30
        !Test for splitting.
        do l=k,1,-1
          nm=l-1 !Note that rv1(1) is always zero.
          if((abs(rv1(l))+anorm).eq.anorm) goto 2
          if((abs(w(nm))+anorm).eq.anorm) goto 1
        end do

      1 c=0.d0 !Cancellation of rv1(l), if l > 1.
        s=1.d0

        do i=l,k
          f=s*rv1(i)
          rv1(i)=c*rv1(i)
          if((abs(f)+anorm).eq.anorm) goto 2
          g=w(i)
          ! h=pythag(f,g)
          call pythag(f,g, h)
          w(i)=h
          h=1.d0/h
          c= (g*h)
          s=-(f*h)
          do j=1,m
            y=a(j,nm)
            z=a(j,i)
            a(j,nm)=(y*c)+(z*s)
            a(j,i)=-(y*s)+(z*c)
          end do
        end do

      2 z=w(k)

        !Convergence.
        if(l.eq.k)then
          !Singular value is made nonnegative.
          if(z.lt.0.d0)then
            w(k)=-z
            do j=1,n
              v(j,k)=-v(j,k)
            end do
          endif
          goto 3
        endif

        !!! if(its.eq.30) pause 'no convergence in svdcmp'
        if(its.eq.30) then
          !! write (*, *) 'no convergence in svdcmp'
          !! stop 100
          ! call intpr("No convergence in svdcmp", -1, its, 0)
          ret = 1
          return
        end if

        x=w(l) !Shift from bottom 2-by-2 minor.
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
        ! g=pythag(f,1.d0)
        call pythag(f,1.d0, g)
        f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
        c=1.d0 !Next QR transformation:
        s=1.d0
        do j=l,nm
          i=j+1
          g=rv1(i)
          y=w(i)
          h=s*g
          g=c*g
          ! z=pythag(f,h)
          call pythag(f,h, z)
          rv1(j)=z
          c=f/z
          s=h/z
          f= (x*c)+(g*s)
          g=-(x*s)+(g*c)
          h=y*s
          y=y*c
          do jj=1,n
            x=v(jj,j)
            z=v(jj,i)
            v(jj,j)= (x*c)+(z*s)
            v(jj,i)=-(x*s)+(z*c)
          end do
          ! z=pythag(f,h)
          call pythag(f,h, z)
          w(j)=z !Rotation can be arbitrary if z = 0.
          if(z.ne.0.d0)then
            z=1.d0/z
            c=f*z
            s=h*z
          endif
          f= (c*g)+(s*y)
          x=-(s*g)+(c*y)
          do jj=1,m
            y=a(jj,j)
            z=a(jj,i)
            a(jj,j)= (y*c)+(z*s)
            a(jj,i)=-(y*s)+(z*c)
          end do
        end do !j=l;nm

        rv1(l)=0.d0
        rv1(k)=f
        w(k)=x
      end do !its=1,30

    3 continue

    end do !k=n,1,-1

    return

  end subroutine svdcmp

  subroutine pythag(a, b, z)
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    real(rp), intent(in) :: a, b
    real(rp), intent(out) :: z
    !Computes sqrt(a**2 + b**2) without destructive underflow or overflow.
    real(rp) :: absa, absb

    absa = abs(a)
    absb = abs(b)

    if(absa.gt.absb) then
      z = absa * sqrt(1. + (absb/absa)**2)
    else
      if(absb.eq.0.)then
        z = 0.
      else
        z = absb * sqrt(1. + (absa/absb)**2)
      endif
    endif

    return
  end subroutine pythag
