/*
 * Run an observation down the tree, and return the prediction error,
 *    for several CP values at once.
 *
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
//rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp)
rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, double p)
{
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    pNode otree =  tree;
    
    // for debug only:
    int opnumber = 0;

    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
      while (cp[i] < tree->complexity) {
	    tree = branch(tree, obs);
	    if (tree == 0)
		goto oops;
	    otree = tree;
	}
	xpred[i] = tree->response_est[0];
  
  //Rprintf("cp = %f, %d's prediction value: %f, ",cp[i], obs2, xpred[i]);
  // error functions where we need to change:
 // xtemp[i] = (*ct_error) (ct.ydata[obs2], tree->response_est);
 
	xtemp[i] = (*ct_error) (ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tree->response_est, p);
  //Rprintf("error: %f\n", xtemp[i]);
    }

    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	for (; i < ct.num_unique_cp; i++)
	    xpred[i] = otree->response_est[0];
  // xtemp[i] = (*ct_error) (ct.ydata[obs2], otree->response_est);
	xtemp[i] = (*ct_error) (ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], otree->response_est, p);
	Rprintf("oops number %d.\n", opnumber++);
  return;
    }
    /*
     * I never really expect to get to this code.  It can only happen if
     *  the last cp on my list is smaller than the terminal cp of the
     *  xval tree just built.  This is impossible (I think).  But just in
     *  case I put a message here.
     */
    warning("Warning message--see rundown.c");
}
