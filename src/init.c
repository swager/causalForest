#include "causalTree.h"
#include "R_ext/Rdynload.h"
#include "node.h"
#include "causalTreeproto.h"

SEXP init_ctcallback(SEXP rhox, SEXP ny, SEXP nr, SEXP expr1x, SEXP expr2x);
SEXP causalTreeexp2(SEXP dtimes, SEXP seps);
SEXP pred_causalTree(SEXP dimx, SEXP nnode, SEXP nsplit, SEXP dimc,
		SEXP nnum, SEXP nodes2, SEXP vnum, SEXP split2,
		SEXP csplit2, SEXP usesur, SEXP xdata2, SEXP xmiss2);
SEXP estimate_causalTree(SEXP dimx, SEXP nnode, SEXP nsplit, SEXP dimc,
  	SEXP nnum, SEXP nodes2, SEXP vnum, SEXP split2,
		SEXP csplit2, SEXP usesur, SEXP xdata2, SEXP xmiss2);

static const R_CallMethodDef CallEntries[] = {
    {"init_ctcallback", (DL_FUNC) &init_ctcallback, 5},
    //{"causalTree", (DL_FUNC) &causalTree, 11},
    //{"causalTree", (DL_FUNC) &causalTree, 12},
    {"causalTree", (DL_FUNC) &causalTree, 16},
    {"xpred", (DL_FUNC) &xpred, 18},
    {"causalTreeexp2", (DL_FUNC) &causalTreeexp2, 2},
    {"pred_causalTree", (DL_FUNC) &pred_causalTree, 12},
    {"estimate_causalTree", (DL_FUNC) &estimate_causalTree, 12},
    {NULL, NULL, 0}
};

#include <Rversion.h>
void
R_init_causalForest(DllInfo * dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
