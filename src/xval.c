#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

#ifndef DEBUG
# define DEBUG 0
#endif
#if DEBUG > 1
static int debug = 0;           
/* if it is odd, print out every tree */
/* if >= 2, print out every risk value we see */
#endif

	void
//xval(int n_xval, CpTable cptable_head, int *x_grp,
//		int maxcat, char **errmsg, double *parms, int *savesort)
xval(int n_xval, CpTable cptable_head, int *x_grp,
   	int maxcat, char **errmsg, double *parms, int minsize, double p, int *savesort, int method, double alpha)
{
	int i, j, k, ii, jj;
	int last;
	int xgroup;
	double *xtemp, *xpred, *xpred2;
	int *savew;
	double *cp;
	double alphasave;
	pNode xtree;
	CpTable cplist;
	double temp;
	double old_wt, total_wt;
  int neighbor; // nearest neighbor number

	alphasave = ct.alpha;
  
  // only for debugging
  
  int cv_count = 0;
  //Rprintf("n_xval = %d\n", n_xval);

	/*
	 * Allocate a set of temporary arrays
	 */
	xtemp = (double *) CALLOC(4 * ct.num_unique_cp, sizeof(double));
	xpred = xtemp + ct.num_unique_cp;
  xpred2 = xpred + ct.num_unique_cp;
	cp = xpred2 + ct.num_unique_cp;
	savew = (int *) CALLOC(ct.n, sizeof(int));
	for (i = 0; i < ct.n; i++)
		savew[i] = ct.which[i]; /* restore at the end */

	/*
	 * Make the list of CPs that I will compare against
	 */
   // test for 
	// cp[0] = 10 * cptable_head->cp;      /* close enough to infinity */
  cp[0] = 10000 * cptable_head->cp;
	for (cplist = cptable_head, i = 1; i < ct.num_unique_cp;cplist = cplist->forward, i++) {  
    //Rprintf("old cp[%d] = %f\n", i, cplist->cp);
		cp[i] = sqrt(cplist->cp * (cplist->forward)->cp);
    //Rprintf("geometric cp[%d] = %f\n", i, cp[i]);
    //rescale the cp here:
    // cp[i] = (n_xval - 1) * 1.0 / n_xval * cp[i];
    //Rprintf("scaled cp[%d] = %f\n", i, cp[i]);
	}
   //rescale alpha:
  // ct.alpha *= (n_xval - 1) * 1.0 / n_xval;
 // Rprintf("ct.alpha = %f\n", ct.alpha);
  

	/* why we need to concern about wt> */
	total_wt = 0;
	for (i = 0; i < ct.n; i++)
		total_wt += ct.wt[i];
	old_wt = total_wt;

	/*
	 * do the validations
	 */
	k = 0;                      /* -Wall */
	for (xgroup = 0; xgroup < n_xval; xgroup++) {
		/*
		 * restore ct.sorts, with the data for this run at the top
		 * this requires one pass per variable
		 */
		for (j = 0; j < ct.nvar; j++) {
			k = 0;
			for (i = 0; i < ct.n; i++) {
				ii = savesort[j * ct.n + i];
				if (ii < 0)
					ii = -(1 + ii);     /* missings move too */
				if (x_grp[ii] != xgroup + 1) { 
					// samples not belong to the test fold:
					/*
					 * this obs is left in --
					 *  copy to the front half of ct.sorts
					 */
					ct.sorts[j][k] = savesort[j * ct.n + i]; // the reason to store in savesort
					k++;
				}
			} 
		}

		/*
		 *  Fix up the y vector, and save a list of "left out" obs *   in
		 * the tail, unused end of ct.sorts[0][i];
		 */
		last = k;

		k = 0;
		temp = 0;
		for (i = 0; i < ct.n; i++) {
			ct.which[i] = 1;    /* everyone starts in group 1 */
			if (x_grp[i] == xgroup + 1) {
        //Rprintf("validation data is %d\n", i + 1);
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
  
    /* rescale the cp */
    for (j = 0; j < ct.num_unique_cp; j++) 
      cp[j] *= temp / old_wt;
    ct.alpha *= temp / old_wt;
    old_wt = temp;


		/*
		 * partition the new tree
		 */
		xtree = (pNode) CALLOC(1, nodesize);
		xtree->num_obs = k;
		(*ct_init) (k, ct.ytemp, maxcat, errmsg, parms, &temp, 2, ct.wtemp, ct.trtemp);
		//(*ct_eval) (k, ct.ytemp, xtree->response_est, &(xtree->risk), ct.wtemp);
    if (method ==5)
      (*ct_eval) (k, ct.ytemp, xtree->response_est, &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha);
    else 
      (*ct_eval) (k, ct.ytemp, xtree->response_est, &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y);
		xtree->complexity = xtree->risk;
    //Rprintf("xtree->complexity = %f\n", xtree->complexity);
		//partition(1, xtree, &temp, 0, k);

    partition(1, xtree, &temp, 0, k, minsize, method, alpha);
    
    //Rprintf("now, xtree->complexity = %f\n", xtree->complexity);
		//the complexity should be min(me, any-node-above-me). This routine fixes that.
		fix_cp(xtree, xtree->complexity);
    //Rprintf("after fixation, xtree->complexity = %f\n", xtree->complexity);
    //Rprintf("cv_count = %d\n\n", ++cv_count);
		
    /*
		 * run the extra data down the new tree
		 */
     
		for (i = k; i < ct.n; i++) {
      j = ct.sorts[0][i]; // left-out samples for testing
			//rundown(xtree, j, cp, xpred, xtemp); 
      // for testing only
      //Rprintf("validation %d ", j+1);
     // Rprintf("x1 variable %f ", ct.xdata[0][j]);
      
      if (p < 0) {
       // Rprintf("--matching: ");
        //matching method:
        neighbor = findNeighbor(j, k); 
        //Rprintf("k = %d.\n", k);
        
        //Rprintf("its neighbor %d\n", neighbor+1);
        rundown3(xtree, j, neighbor, cp, xpred, xpred2, xtemp);
        // Rprintf("the error is %f\n", xtemp[]);
        
       
      } else {
        // TOT method:
        //Rprintf("--TOT: ");
        rundown(xtree, j, cp, xpred, xtemp, p);
      }
     
#if DEBUG > 1
			if (debug > 1) {
				jj = j + 1;
				Rprintf("\nObs %d, y=%f \n", jj, ct.ydata[j][0]);
			}
#endif
			/* add it in to the risk */
			cplist = cptable_head;
			for (jj = 0; jj < ct.num_unique_cp; jj++) {
				cplist->xrisk += xtemp[jj] * ct.wt[j];
				//cplist->xrisk += xtemp[jj];
				cplist->xstd += xtemp[jj] * xtemp[jj] * ct.wt[j];
        //cplist->xstd += xtemp[jj] * xtemp[jj];
        
#if DEBUG > 1
			//	if (debug > 1)
					Rprintf("  cp=%f, pred=%f, xtemp=%f\n",
							cp[jj] , xpred[jj], xtemp[jj]);
#endif
				cplist = cplist->forward;
			}
      // debug only:
      //round ++;
		}
    //Rprintf("%d cv round!\n", round);
    //Rprintf("count = %d\n", count);
		free_tree(xtree, 1);    // Calloc-ed
		R_CheckUserInterrupt();
	}

	for (cplist = cptable_head; cplist; cplist = cplist->forward) {
		cplist->xstd = sqrt(cplist->xstd -
				cplist->xrisk * cplist->xrisk / total_wt);
    //cplist->xstd = sqrt(cplist->xstd -
  		//	cplist->xrisk * cplist->xrisk / ct.n);
	}
	ct.alpha = alphasave;
	for (i = 0; i < ct.n; i++)
		ct.which[i] = savew[i];
	Free(savew);
	Free(xtemp);
}
