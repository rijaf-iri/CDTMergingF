#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void F77_NAME(idw_interp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(shepard_interp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(barnes_interp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spheremap_interp)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cressman_interp)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kriging_interp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(distance_vector)(void *, void *, void *, void *, void *);
extern void F77_NAME(distance_matrix)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"idw_interp",        (DL_FUNC) &F77_NAME(idw_interp),        10},
  {"shepard_interp",    (DL_FUNC) &F77_NAME(shepard_interp),    10},
  {"barnes_interp",     (DL_FUNC) &F77_NAME(barnes_interp),     10},
  {"spheremap_interp",  (DL_FUNC) &F77_NAME(spheremap_interp),   9},
  {"cressman_interp",   (DL_FUNC) &F77_NAME(cressman_interp),    9},
  {"kriging_interp",    (DL_FUNC) &F77_NAME(kriging_interp),    13},
  {"distance_vector",   (DL_FUNC) &F77_NAME(distance_vector),    5},
  {"distance_matrix",   (DL_FUNC) &F77_NAME(distance_matrix),    6},
  {NULL,  NULL,   0}
};

void R_init_CDTMergingF (DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
