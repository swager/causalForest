/*
 * The main entry point for recursive partitioning routines.
 *
 * Input variables:
 *      ncat    = # categories for each var, 0 for continuous variables.
 *      method  = 1 - anova
 *                2 - exponential survival
                  3 - anova2
 *		  3 - classification
 *	          4 - user defright_wt_sum / right_wt- (right_sum - right_wt_sum) / (right_n - right_wt);ined callback
 *      opt     = vector of options.  Same order as causalTree.control, as a vector
 *                   of doubles.
 *      parms   = extra parameters for the split function, e.g. poissoninit
 *      xvals    = number of cross-validations to do
 *      xgrp     = indices for the cross-validations
 *      ymat    = vector of response variables
 *      xmat    = matrix of continuous variables
 *      ny      = number of columns of the y matrix (it is passed in as a
 *                  vector)
 *      wt       = vector of case weights
 *
 * Returned: a list with elements
 *      which = vector of final node numbers for each input obs
 *      cptable = the complexity table
 *      dsplit = for each split, numeric variables (doubles)
 *      isplit = for each split, integer variables
 *      dnode =  for each node, numeric variables
 *      inode =  for each node, integer variables
 *
 * Naming convention: ncat = pointer to an integer vector, ncat2 = the
 *   input R object (SEXP) containing that vector, ncat3 = an output S object
 *   containing that vector.
 */
#define MAINRP
#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "func_table.h"
#include "causalTreeproto.h"

SEXP
causalTree(SEXP ncat2, SEXP method2, SEXP opt2,
      SEXP parms2, SEXP minsize2, SEXP p2, SEXP xvals2, SEXP xgrp2,
      SEXP ymat2, SEXP xmat2, SEXP wt2, SEXP treatment2, SEXP ny2, SEXP cost2, SEXP xvar2, SEXP alpha2)
{

    pNode tree;          /* top node of the tree */
    char *errmsg;
    int i, j, k, n;
    int maxcat;
    double temp, temp2;
    int *savesort = NULL /* -Wall */ ;
    double *dptr;               /* temp */
    int *iptr;
    int method;
    /*
     * pointers to R objects
     */
    int *ncat, *xgrp;
    int xvals;
    //double *wt, *parms;
    double *wt;
    double *treatment;
    double *parms;
    int minsize;
    // add propensity score:
    double p;
    double alpha;
    /*
     * Return objects for R -- end in "3" to avoid overlap with internal names
     */
    SEXP which3, cptable3, dsplit3, isplit3, csplit3 = R_NilValue, /* -Wall */
	dnode3, inode3;

    /* work arrays for the return process */
    int nodecount, catcount, splitcount;
    double **ddnode, *ddsplit[3];
    int *iinode[6], *iisplit[3];
    int **ccsplit;
    double scale;
    CpTable cp;

    ncat = INTEGER(ncat2);
    xgrp = INTEGER(xgrp2);
    xvals = asInteger(xvals2);
    wt = REAL(wt2);
    treatment = REAL(treatment2);
    parms = REAL(parms2);
    //parms = asInteger(parms2);
    minsize = asInteger(minsize2);
    p = asReal(p2);
    alpha = asReal(alpha2);
    method = asInteger(method2);
    
    
    /*
     * initialize the splitting functions from the function table
     */
    if (asInteger(method2) <= NUM_METHODS) {
      i = asInteger(method2) - 1;
	    ct_init = func_table[i].init_split;
	    ct_choose = func_table[i].choose_split;
	    ct_eval = func_table[i].eval;
	    ct_error = func_table[i].error;
	    ct.num_y = asInteger(ny2);
    } else
      error(_("Invalid value for 'method'"));

    /*
     * set some other parameters
     */
    dptr = REAL(opt2);
    ct.min_node = (int) dptr[1];
    ct.min_split = (int) dptr[0];
    ct.complexity = dptr[2];
    ct.maxpri = (int) dptr[3] + 1;      /* max primary splits =
					   max competitors + 1 */
    if (ct.maxpri < 1)
      ct.maxpri = 1;
    ct.maxsur = (int) dptr[4];
    ct.usesurrogate = (int) dptr[5];
    ct.sur_agree = (int) dptr[6];
    ct.maxnode = (int) pow((double) 2.0, (double) dptr[7]) - 1;
    ct.n = nrows(xmat2);
    n = ct.n;                   
    /* I get tired of typing "ct.n" 100 times 
    * below */
    //Rprintf("ct.n = %d\n", n);
    
    ct.nvar = ncols(xmat2);
    ct.numcat = INTEGER(ncat2);
    ct.wt = wt;
    ct.treatment = treatment;
    ct.iscale = 0.0;
    ct.vcost = REAL(cost2);
    
    // pass the variance of predictors into ct:
    ct.xvar = REAL(xvar2);

    /*
     * create the "ragged array" pointers to the matrix
     *   x and missmat are in column major order
     *   y is in row major order
     */
    dptr = REAL(xmat2);
    ct.xdata = (double **) ALLOC(ct.nvar, sizeof(double *));
    for (i = 0; i < ct.nvar; i++) {
      ct.xdata[i] = dptr;
	    dptr += n;
    }
    
    ct.ydata = (double **) ALLOC(n, sizeof(double *));

    dptr = REAL(ymat2);
    temp2 = 0;
    for (i = 0; i < n; i++) {
	    ct.ydata[i] = dptr;
      for (j = 0; j <  ct.num_y; j++) {
        if (fabs(ct.ydata[i][j]) > temp2) temp2 = fabs(ct.ydata[i][j]);        
      }
      	dptr += ct.num_y;
    }
    ct.max_y = temp2;
    /*
     * allocate some scratch
     */
    ct.tempvec = (int *) ALLOC(n, sizeof(int));
    ct.xtemp = (double *) ALLOC(n, sizeof(double));
    ct.ytemp = (double **) ALLOC(n, sizeof(double *));
    ct.wtemp = (double *) ALLOC(n, sizeof(double));
    ct.trtemp = (double *) ALLOC(n, sizeof(double));

    /*
     * create a matrix of sort indices, one for each continuous variable
     *   This sort is "once and for all".
     * I don't have to sort the categoricals.
     */
    ct.sorts = (int **) ALLOC(ct.nvar, sizeof(int *));
    ct.sorts[0] = (int *) ALLOC(n * ct.nvar, sizeof(int));
    maxcat = 0;
    for (i = 0; i < ct.nvar; i++) {
      ct.sorts[i] = ct.sorts[0] + i * n;
	    for (k = 0; k < n; k++) {
        if (!R_FINITE(ct.xdata[i][k])) {
		      ct.tempvec[k] = -(k + 1);       /* this variable is missing */
		      ct.xtemp[k] = 0;        /* avoid weird numerics in S's NA */
	       } else {
		      ct.tempvec[k] = k;
		      ct.xtemp[k] = ct.xdata[i][k];
	      }
	    }
	    if (ncat[i] == 0)
	      mysort(0, n - 1, ct.xtemp, ct.tempvec);
	    else if (ncat[i] > maxcat)
	      maxcat = ncat[i];
	    for (k = 0; k < n; k++)
	      ct.sorts[i][k] = ct.tempvec[k];
    }

    /*
     * save away a copy of the ct.sorts, if needed for xval
     */
    if (xvals > 1) {
	    savesort = (int *) ALLOC(n * ct.nvar, sizeof(int));
	    memcpy(savesort, ct.sorts[0], n * ct.nvar * sizeof(int));
    }

    /*
     * And now the last of my scratch space
     */
    if (maxcat > 0) {
	    ct.csplit = (int *) ALLOC(3 * maxcat, sizeof(int));
	    ct.left = ct.csplit + maxcat;
	    ct.right = ct.left + maxcat;
      ct.lwt = (double *) ALLOC(2 * maxcat, sizeof(double));
	    ct.rwt = ct.lwt + maxcat;
      ct.ltr = (double *) ALLOC(2 * maxcat, sizeof(double));
      ct.rtr = ct.ltr + maxcat;
    } else
	    ct.csplit = (int *) ALLOC(1, sizeof(int));

    /*
     * initialize the top node of the tree
     */
    errmsg = _("unknown error");
    which3 = PROTECT(allocVector(INTSXP, n));
    ct.which = INTEGER(which3);
    temp = 0;
    temp2 = 0;
    
    for (i = 0; i < n; i++) {
	    ct.which[i] = 1;
	    temp += wt[i];
      temp2 += treatment[i];
    }

    // question: do I need to put treatment in init function？
    i = (*ct_init) (n, ct.ydata, maxcat, &errmsg, parms, &ct.num_resp, 1, wt, treatment);
    if (i > 0)
	    error(errmsg);
    
    nodesize = sizeof(Node) + (ct.num_resp - 20) * sizeof(double);
    tree = (pNode) ALLOC(1, nodesize);
    memset(tree, 0, nodesize);
    tree->num_obs = n;
    tree->sum_wt = temp;
    tree->sum_tr = temp2;

    // add ct.maxx_y inside the evaluation funciton
    // (*ct_eval) (n, ct.ydata, tree->response_est, &(tree->risk), wt)
    // here I am now >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (method == 5) // anova2:
      (*ct_eval) (n, ct.ydata, tree->response_est, &(tree->risk), wt, treatment, ct.max_y, alpha);
    else 
      (*ct_eval) (n, ct.ydata, tree->response_est, &(tree->risk), wt, treatment, ct.max_y);
    tree->complexity = tree->risk;
    ct.alpha = ct.complexity * tree->risk;
    // for debug only:
    //Rprintf("The risk is %f\n", tree->complexity);
    //Rprintf("ct.alpha = %f\n", ct.alpha);

    /*
     * Do the basic tree
     */
    //partition(1, tree, &temp, 0, n);
    partition(1, tree, &temp, 0, n, minsize, method, alpha);
    
    CpTable cptable = (CpTable) ALLOC(1, sizeof(cpTable));
    cptable->cp = tree->complexity;
    cptable->risk = tree->risk;
    cptable->nsplit = 0;
    cptable->forward = 0;
    cptable->xrisk = 0;
    cptable->xstd = 0;
    ct.num_unique_cp = 1;

    if (tree->rightson) {
      // need to look at...................................................
      make_cp_list(tree, tree->complexity, cptable);
	    make_cp_table(tree, tree->complexity, 0);
	if (xvals > 1) {
      xval(xvals, cptable, xgrp, maxcat, &errmsg, parms, minsize, p, savesort, method, alpha);
	}
    }
    /*
     * all done, create the return list for R
     * first the cp table
     */
    scale = 1 / tree->risk;
    //Rprintf("the scale for ct is %f\n", scale);
    i = 0;
    cptable3 = PROTECT(allocMatrix(REALSXP, xvals > 1 ? 5 : 3,
				   ct.num_unique_cp));
    dptr = REAL(cptable3);
    for (cp = cptable; cp; cp = cp->forward) {
	dptr[i++] = cp->cp * scale;
	dptr[i++] = cp->nsplit;
	dptr[i++] = cp->risk * scale;
 
	if (xvals > 1) {
      
      //Rprintf("cp->xrisk = %f, cp->std = %f\n", cp->xrisk, cp->xstd);
	    dptr[i++] = cp->xrisk * scale ;
	    dptr[i++] = cp->xstd * scale ;
	}
    }

    /*
     * Return the body of the tree
     *  For each component we first create a vector to hold the
     *  result, then a ragged array index into the vector.
     * The ctmatrix routine then fills everything in.
     */
    ctcountup(tree, &nodecount, &splitcount, &catcount);
    dnode3 = PROTECT(allocMatrix(REALSXP, nodecount, (3 + ct.num_resp)));
    ddnode = (double **) ALLOC(3 + ct.num_resp, sizeof(double *));
    dptr = REAL(dnode3);
    for (i = 0; i < 3 + ct.num_resp; i++) {
	ddnode[i] = dptr;
	dptr += nodecount;
    }

    dsplit3 = PROTECT(allocMatrix(REALSXP, splitcount, 3));
    dptr = REAL(dsplit3);
    for (i = 0; i < 3; i++) {
	ddsplit[i] = dptr;
	dptr += splitcount;
	for (j = 0; j < splitcount; j++)
	    ddsplit[i][j] = 0.0;
    }

    inode3 = PROTECT(allocMatrix(INTSXP, nodecount, 6));
    iptr = INTEGER(inode3);
    for (i = 0; i < 6; i++) {
	iinode[i] = iptr;
	iptr += nodecount;
    }

    isplit3 = PROTECT(allocMatrix(INTSXP, splitcount, 3));
    iptr = INTEGER(isplit3);
    for (i = 0; i < 3; i++) {
	iisplit[i] = iptr;
	iptr += splitcount;
    }

    if (catcount > 0) {
	csplit3 = PROTECT(allocMatrix(INTSXP, catcount, maxcat));
	ccsplit = (int **) ALLOC(maxcat, sizeof(int *));
	iptr = INTEGER(csplit3);
	for (i = 0; i < maxcat; i++) {
	    ccsplit[i] = iptr;
	    iptr += catcount;
	    for (j = 0; j < catcount; j++)
		ccsplit[i][j] = 0;      /* zero it out */
	}
    } else
	ccsplit = NULL;

    ctmatrix(tree, ct.numcat, ddsplit, iisplit, ccsplit, ddnode, iinode, 1);
    free_tree(tree, 0);         /* let the memory go */

    /*
     * Fix up the 'which' array
     *  Nodes are sometimes trimmed during the
     *  tree building, and 'which' is not updated in that case
     */
    for (i = 0; i < n; i++) {
	k = ct.which[i];
	do {
	    for (j = 0; j < nodecount; j++)
		if (iinode[0][j] == k) {
		    ct.which[i] = j + 1;
		    break;
		}
	    k /= 2;
	} while (j >= nodecount);
    }

    /* Create the output list */
    int nout = catcount > 0 ? 7 : 6;
    SEXP rlist = PROTECT(allocVector(VECSXP, nout));
    SEXP rname = allocVector(STRSXP, nout);
    setAttrib(rlist, R_NamesSymbol, rname);
    SET_VECTOR_ELT(rlist, 0, which3);
    SET_STRING_ELT(rname, 0, mkChar("which"));
    SET_VECTOR_ELT(rlist, 1, cptable3);
    SET_STRING_ELT(rname, 1, mkChar("cptable"));
    SET_VECTOR_ELT(rlist, 2, dsplit3);
    SET_STRING_ELT(rname, 2, mkChar("dsplit"));
    SET_VECTOR_ELT(rlist, 3, isplit3);
    SET_STRING_ELT(rname, 3, mkChar("isplit"));
    SET_VECTOR_ELT(rlist, 4, dnode3);
    SET_STRING_ELT(rname, 4, mkChar("dnode"));
    SET_VECTOR_ELT(rlist, 5, inode3);
    SET_STRING_ELT(rname, 5, mkChar("inode"));
    if (catcount > 0) {
	SET_VECTOR_ELT(rlist, 6, csplit3);
	SET_STRING_ELT(rname, 6, mkChar("csplit"));
    }

    UNPROTECT(1 + nout);
    return rlist;
}
