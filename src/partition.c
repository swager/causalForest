/*
 * The main workhorse of the recursive partitioning module.  When called
 *   with a node, it partitions it and then calls itself to partition the
 *   children it has created.
 * If the node is not splittable (too few people, or complexity is too small)
 *   it simply returns.  The routine may not be able to discover that the
 *   complexity is too small until after the children have been partitioned,
 *   so it needs to check this at the end.
 * The vector who[n] indexes which observations are in this node, to speed
 *   up the routine.
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

int
//partition(int nodenum, pNode splitnode, double *sumrisk, int n1, int n2)
partition(int nodenum, pNode splitnode, double *sumrisk, int n1, int n2, int minsize, int method, double alpha)
{
    pNode me;
    double tempcp;
    int i, j, k;
    double tempcp2;
    double left_risk, right_risk;
    int left_split, right_split;
    double twt, ttr;
    int nleft, nright;
    int n;
    int min_node_size = minsize;

    me = splitnode;
    n = n2 - n1;                /* total number of observations */

    if (nodenum > 1) {
      twt = 0;
      ttr = 0;
	    k = 0;
	    for (i = n1; i < n2; i++) {
	      j = ct.sorts[0][i]; /* any variable would do, use first */
	      if (j < 0)
		      j = -(1 + j);   /* if missing, value = -(1+ true index) */
	      ct.wtemp[k] = ct.wt[j];
        ct.trtemp[k] = ct.treatment[j];
        //Rprintf("trtemp[%d] = %f\n", k, ct.trtemp[k]);
	      ct.ytemp[k] = ct.ydata[j];
	      twt += ct.wt[j];
        ttr += ct.treatment[j];
	      k++;
	    }
	    //(*ct_eval) (n, ct.ytemp, me->response_est, &(me->risk), ct.wtemp);
      if (method == 5)
        (*ct_eval) (n, ct.ytemp, me->response_est, &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha);
      else
        (*ct_eval) (n, ct.ytemp, me->response_est, &(me->risk), ct.wtemp, ct.trtemp, ct.max_y);
	    me->num_obs = n;
	    me->sum_wt = twt;
      me->sum_tr = ttr;
	    tempcp = me->risk;
	    if (tempcp > me->complexity)
	      tempcp = me->complexity;
    } else
	    tempcp = me->risk; // nodenum == 1 i.e.: root node
   // Rprintf("me->num_obs= %d\n", me->num_obs); 
    //Rprintf("temcp = %f\n", tempcp);

    /*
     * Can I quit now ?
     */
    // Rprintf("min_split = %d\n", ct.min_split);
    //Rprintf("ct.maxnode = %d, nodenum = %d\n", ct.maxnode, nodenum);
    //Rprintf("me->num_obs= %d, ct.min_split =%d\n", me->num_obs , ct.min_split);
    //Rprintf("tempcp = %f\n", tempcp);
    if (me->num_obs < ct.min_split || tempcp <= ct.alpha || nodenum > ct.maxnode) {
      me->complexity = ct.alpha;
  	  *sumrisk = me->risk;
      //Rprintf("quit now!\n");
	/*
	 * make sure the split doesn't have random pointers to somewhere
	 * i.e., don't trust that whoever allocated memory set it to zero
	 */
	    me->leftson = (pNode)  NULL;
	    me->rightson = (pNode) NULL;
	    me->primary = (pSplit) NULL;
	    me->surrogate = (pSplit) NULL;
	    return 0;
    }
    /*
     * Guess I have to do the split
     */
    //bsplit(me, n1, n2);
    bsplit(me, n1, n2, min_node_size, method, alpha);
    if (!me->primary) {
     // Rprintf("stop here!\n");
	/*
	 * This is rather rare -- but I couldn't find a split worth doing
	 */
	    me->complexity = ct.alpha;
	    me->leftson = (pNode) NULL;
	    me->rightson = (pNode) NULL;
	    me->primary = (pSplit) NULL;
	    me->surrogate = (pSplit) NULL;
	    *sumrisk = me->risk;
	    return 0;
    }
#ifdef DEBUG
    print_tree(me, 4);
#endif
    if (ct.maxsur > 0)
	surrogate(me, n1, n2);
    else
	me->surrogate = (pSplit) NULL;
    nodesplit(me, nodenum, n1, n2, &nleft, &nright);

    /*
     * split the leftson
     */
    me->leftson = (pNode) CALLOC(1, nodesize);
    (me->leftson)->complexity = tempcp - ct.alpha;
    //Rprintf("me->leftson->complexity = %f\n", (me->leftson)->complexity);
    left_split =
	//partition(2 * nodenum, me->leftson, &left_risk, n1, n1 + nleft);
    partition(2 * nodenum, me->leftson, &left_risk, n1, n1 + nleft, min_node_size, method, alpha);

    /*
     * Update my estimate of cp, and split the right son.
     */
    tempcp = (me->risk - left_risk) / (left_split + 1);
    //Rprintf("left_risk = %f, left_split = %d\n", left_risk, left_split);
    tempcp2 = (me->risk - (me->leftson)->risk);
    if (tempcp < tempcp2)
	tempcp = tempcp2;
   //Rprintf("now tempcp = %f\n", tempcp);
    if (tempcp > me->complexity)
	tempcp = me->complexity;

    me->rightson = (pNode) CALLOC(1, nodesize);
    (me->rightson)->complexity = tempcp - ct.alpha;
    //Rprintf("me->rightson->complexity = %f\n", (me->rightson)->complexity);
    //right_split = partition(1 + 2 * nodenum, me->rightson, &right_risk,
		//	    n1 + nleft, n1 + nleft + nright);
    right_split = partition(1 + 2 * nodenum, me->rightson, &right_risk,
  		    n1 + nleft, n1 + nleft + nright, min_node_size, method, alpha);


    /*
     * Now calculate my actual C.P., which depends on children nodes, and
     *  on grandchildren who do not collapse before the children.
     * The calculation is done assuming that I am the top node of the
     *  whole tree, an assumption to be fixed up later.
     */
    tempcp = (me->risk - (left_risk + right_risk)) /
	  (left_split + right_split + 1);
    //Rprintf("the final tempcp = %f\n", tempcp);

    /* Who goes first -- minimum of tempcp, leftson, and rightson */
    if ((me->rightson)->complexity > (me->leftson)->complexity) {
      if (tempcp > (me->leftson)->complexity) {
	    /* leftson collapses first */
	      left_risk = (me->leftson)->risk;
	      left_split = 0;

	      tempcp = (me->risk - (left_risk + right_risk)) /
        (left_split + right_split + 1);
	      if (tempcp > (me->rightson)->complexity) {
		/* right one goes too */
		      right_risk = (me->rightson)->risk;
		      right_split = 0;
	      }
	    }
    } else if (tempcp > (me->rightson)->complexity) {
	/* right hand child goes first */
	  right_split = 0;
	  right_risk = (me->rightson)->risk;

	  tempcp = (me->risk - (left_risk + right_risk)) /
	    (left_split + right_split + 1);
	  if (tempcp > (me->leftson)->complexity) {
	    /* left one goes too */
	    left_risk = (me->leftson)->risk;
	    left_split = 0;
      }
    }
    
    me->complexity = (me->risk - (left_risk + right_risk)) /
	(left_split + right_split + 1);
 //Rprintf("the final cp = %f\n", me->complexity);
  

    if (me->complexity <= ct.alpha) {
	/*
	 * All was in vain!  This node doesn't split after all.
	 */
	free_tree(me, 0);
	*sumrisk = me->risk;
	for (i = n1; i < n2; i++) {
	    j = ct.sorts[0][i];
	    if (j < 0)
		j = -(1 + j);
	    ct.which[j] = nodenum;      /* revert to the old nodenumber */
	}
	return 0;               /* return # of splits */
    } else {
	*sumrisk = left_risk + right_risk;
	return left_split + right_split + 1;
    }
}
