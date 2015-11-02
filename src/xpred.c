/*
 *  This has almost exactly the same inputs are causalTree, but returns
 *   cross-validated predictions instead of the fitted tree.
 *
 * Input variables:
 *      ncat    = # categories for each var, 0 for continuous variables.
 *      method  = 1 - anova
 *                2 - exponential survival
 *		  3 - classification
 *	          4 - user defined callback
 *      opt     = vector of options.  Same order as causalTree.control, as a vector
 *                   of doubles.
 *      parms   = min_node_size
 *      xvals    = number of cross-validations to do
 *      xgrp     = indices for the cross-validations
 *      ymat    = vector of response variables
 *      xmat    = matrix of continuous variables
 *      ny      = number of columns of the y matrix (it is passed in as a
 *                  vector)
 *      wt       = vector of case weights
 *      all     = if 1 return all the predictions, otherwise just the
 *                 first element
 *      cp      = vector of cp values to use as cut points
 *      toprisk = risk for the top node of the tree
 *      nresp   = number of response values
 */
#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "func_table.h"
#include "causalTreeproto.h"

SEXP
xpred(SEXP ncat2, SEXP method2, SEXP opt2,
      SEXP parms2, SEXP minsize2, SEXP xvals2, SEXP xgrp2,
      SEXP ymat2, SEXP xmat2, SEXP wt2, SEXP treatment2,
      SEXP ny2, SEXP cost2, SEXP all2, SEXP cp2, SEXP toprisk2, SEXP nresp2, SEXP alpha2)
{
    char *errmsg;
    int i, j, k, n;
    int last, ii;
    int maxcat, ncp;
    int xgroup;
    double temp, total_wt, old_wt;
    int *savesort;
    double *dptr;               /* temp */
    int nresp;
    double toprisk;

    pNode xtree;
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
    double *predict;
    double *cp;
    double alpha = asReal(alpha2);
    int method = asInteger(method2);
    

    /*
     *        Return objects for R
     */
    SEXP predict2;

    /*
     *  the first half of the routine is almost identical to causalTree.c
     * first get copies of some input variables
     */
    ncat = INTEGER(ncat2);
    xgrp = INTEGER(xgrp2);
    xvals = asInteger(xvals2);
    wt = REAL(wt2);
    treatment = REAL(treatment2);
    //parms = REAL(parms2);
    parms = REAL(parms2);
    minsize = asInteger(minsize2);
    ncp = LENGTH(cp2);
    cp = REAL(cp2);
    toprisk = asReal(toprisk2);

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
    ct.maxpri = (int) dptr[3] + 1;      /* max primary splits = max
					 * competitors + 1 */
    if (ct.maxpri < 1)
	    ct.maxpri = 1;
    ct.maxsur = (int) dptr[4];
    ct.usesurrogate = (int) dptr[5];
    ct.sur_agree = (int) dptr[6];
    ct.maxnode = (int) pow((double) 2.0, (double) dptr[7]) - 1;
    ct.n = nrows(xmat2);
    n = ct.n;                   
    /* I get tired of typing "ct.n" 100 times * below */
    ct.nvar = ncols(xmat2);
    ct.numcat = INTEGER(ncat2);
    ct.wt = wt;
    ct.treatment = treatment;
    ct.iscale = 0.0;
    ct.vcost = REAL(cost2);
    ct.num_resp = asInteger(nresp2);

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
    for (i = 0; i < n; i++) {
	    ct.ydata[i] = dptr;
	    dptr += ct.num_y;
    }
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
     * save away a copy of the ct.sorts
     */
    savesort = (int *) ALLOC(n * ct.nvar, sizeof(int));
    memcpy(savesort, ct.sorts[0], n * ct.nvar * sizeof(int));

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
     * Initialize the top node
     */

    ct.which = (int *) ALLOC(n, sizeof(int));
    xtree = (pNode) ALLOC(1, nodesize);
    
    (*ct_init) (n, ct.ydata, maxcat, &errmsg, parms, &ct.num_resp, 1, wt, treatment);

    /*
     * From this point on we look much more like xval.c
     */
    ct.alpha = ct.complexity * toprisk;
    for (i = 0; i < ncp; i++) {
      cp[i] *= toprisk;       /* scale to internal units */
      // rescale the cp:
      //cp[i] *= (xvals - 1) * 1.0 / xvals; 
      //Rprintf("cp[%d] = %f\n", i, cp[i]);   
    }
    //ct.alpha *= (xvals - 1) * 1.0 / xvals;
    
    //Rprintf("ct.alpha = %f\n", ct.alpha);

    /*
     *        allocate the output vector
     */
    if (asInteger(all2) == 1)
	    nresp = ct.num_resp;    /* number returned */
    else
    	nresp = 1;
    predict2 = PROTECT(allocVector(REALSXP, n * ncp * nresp));
    predict = REAL(predict2);

    /*
     * do the validations
     */
    total_wt = 0;
    for (i = 0; i < ct.n; i++)
	    total_wt += ct.wt[i];
    old_wt = total_wt;

    k = 0;                      /* -Wall */
    for (xgroup = 0; xgroup < xvals; xgroup++) {
      
	    /*
	     * restore ct.sorts, with the data for this run at the top
	     * this requires one pass per variable
	     */
	    for (j = 0; j < ct.nvar; j++) {
	      k = 0;
	      for (i = 0; i < ct.n; i++) {
		      ii = savesort[j * n + i];       /* walk through the variables
						 * in order */
		      if (ii < 0)
		        ii = -(1 + ii);     /* missings move too */
		      if (xgrp[ii] != xgroup + 1) {      
           /*
		        * this obs is left in --
		        * copy to the front half of ct.sorts
		        */
            ct.sorts[j][k] = savesort[j * n + i];
		        k++;
	        }
	      }
	    }

	    /*
	     * Fix up the y vector, and save a list of "left out" obs in
	     * the tail, unused end of ct.sorts[0][i];
	     */
	    last = k;
	    k = 0;
	    temp = 0;
	    for (i = 0; i < n; i++) {
        ct.which[i] = 1;    /* everyone starts in the top node */
	      if (xgrp[i] == xgroup + 1) {
		      ct.sorts[0][last] = i;
		      last++;
	      } else {
		      ct.ytemp[k] = ct.ydata[i];
		      ct.wtemp[k] = ct.wt[i];
          ct.trtemp[k] = ct.treatment[i];
		      temp += ct.wt[i];
		      k++;
	      }
	    }

	    /* at this point k = #obs in the prediction group */
	    /* rescale the cp */
	    for (j = 0; j < ct.num_unique_cp; j++)
	       cp[j] *= temp / old_wt;
	    ct.alpha *= temp / old_wt;
	    old_wt = temp;

	  /*
	   * partition the new tree
	   */
	    xtree->num_obs = k;
	    (*ct_init) (k, ct.ytemp, maxcat, &errmsg, parms, &ii, 2, ct.wtemp, ct.trtemp);
	    //(*ct_eval) (k, ct.ytemp, xtree->response_est, &(xtree->risk), ct.wtemp);
      if (method == 5)
        (*ct_eval) (k, ct.ytemp, xtree->response_est, &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha);
      else
        (*ct_eval) (k, ct.ytemp, xtree->response_est, &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y);
	    xtree->complexity = xtree->risk;
	    //partition(1, xtree, &temp, 0, k);
      partition(1, xtree, &temp, 0, k, minsize, method, alpha);
	    fix_cp(xtree, xtree->complexity);
     // Rprintf("%d th xtree\n", xgroup + 1);
	   // print_tree(xtree, 5);        /* debug line */
      //Rprintf("end.\n");
      
	/*
	 * run the extra data down the new tree
	 */
	  for (i = k; i < ct.n; i++) {
      j = ct.sorts[0][i];
      // need to be tested:
      //Rprintf("%d obs:\n", j + 1);
	    rundown2(xtree, j, cp, (predict + j * ncp * nresp), nresp);
	  }

	  free_tree(xtree, 0);
	  R_CheckUserInterrupt();
   }

    UNPROTECT(1);
    return predict2;
}
