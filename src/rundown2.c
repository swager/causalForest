/*
 * Run an observation down the tree, and return the predicted value,
 *    for several CP values at once.
 * (A subset of rundown.c, which also returns the prediction error).
 */
#include "node.h"
#include "causalTree.h"
#include "causalTreeproto.h"

void
rundown2(pNode tree, int obs, double *cp, double *xpred, int nresp)
{
    int i, j, k = 0;
    pNode otree = tree;

    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
      //Rprintf("cp[%d] = %f, old_tree->compleixity = %f\n", i, cp[i], tree->complexity);
      while (cp[i] < tree->complexity) {
        tree = branch(tree, obs);
        //Rprintf("tree->complexity = %f, tree->response = %f\n", tree->complexity,tree->response_est[0]);
	      if (tree == 0)
		      goto oops;
	      otree = tree;
      }
	    for (j = 0; j < nresp; j++)
	      xpred[k++] = tree->response_est[j];
    }

    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	for (; i < ct.num_unique_cp; i++)
	    for (j = 0; j < nresp; j++)
		xpred[k++] = otree->response_est[j];
	return;
    }
    /*
     * I never really expect to get to this code.  It can only happen if
     *  the last cp on my list is smaller than the terminal cp of the
     *  xval tree just built.  This is impossible (I think).  But just in
     *  case I put a message here.
     */
    warning("Warning message--see rundown2.c");
}
