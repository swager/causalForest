
/*
 * The routine which will find the best split for a node
 *
 * Input :      node
 *              node number
 *
 * Output:      Fills in the node's
 *                      primary splits
 *                      competitor splits
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
//bsplit(pNode me, int n1, int n2)
bsplit(pNode me, int n1, int n2, int minsize, int method, double alpha)
{
	int i, j, k;
	int kk;
	int nc;
	double improve;
	double split = 0.0;
	pSplit tsplit;
	int *index;
	double *xtemp;              /* these 3 because I got tired of typeing
								 * "ct.xtemp", etc */
	double **ytemp;
	double *wtemp;
  double *trtemp;

	xtemp = ct.xtemp;
	ytemp = ct.ytemp;
	wtemp = ct.wtemp;
  trtemp = ct.trtemp;

	/*
	 * test out the variables 1 at at time
	 */
	me->primary = (pSplit) NULL;
	for (i = 0; i < ct.nvar; i++) {
		index = ct.sorts[i];
		nc = ct.numcat[i];
		/* extract x and y data */
		k = 0;
		for (j = n1; j < n2; j++) {
			kk = index[j];
      
			//if (kk >= 0 && ct.wt[kk] > 0) {
			//}
			/* x data not missing and wt > 0 */
			if(kk >= 0 && ct.wt[kk] > 0) { // here we can assgin the weight to be 0, kk >= 0 means x data not missing
				xtemp[k] = ct.xdata[i][kk];
				ytemp[k] = ct.ydata[kk];
				wtemp[k] = ct.wt[kk];
        trtemp[k] = ct.treatment[kk];
        //Rprintf("trtemp[%d] = %f\n", k, trtemp[k]);
				k++;
			}
		}
    

		if (k == 0 || (nc == 0 && xtemp[0] == xtemp[k - 1]))
			continue;           /* no place to split */

		//(*ct_choose) (k, ytemp, xtemp, nc, ct.min_node, &improve,
		//	      &split, ct.csplit, me->risk, wtemp);
    if (method == 5)
      (*ct_choose) (k, ytemp, xtemp, nc, ct.min_node, &improve,
				&split, ct.csplit, me->risk, wtemp, trtemp, minsize, alpha);
	  else 	
      (*ct_choose) (k, ytemp, xtemp, nc, ct.min_node, &improve,
				&split, ct.csplit, me->risk, wtemp, trtemp, minsize);
        
        //Rprintf("%d predictor has improve = %f\n", i, improve);


		/*
		 * Originally, this just said "if (improve > 0)", but rounding
		 * error will sometimes create a non zero that should be 0.  Yet we
		 * want to retain invariance to the scale of "improve".
		 */
		if (improve > ct.iscale)
			ct.iscale = improve;  /* largest seen so far */
		//Rprintf("improve = %f\n", improve);
		if (improve > (ct.iscale * 1e-10)) {
			improve /= ct.vcost[i];     /* scale the improvement */
			tsplit = insert_split(&(me->primary), nc, improve, ct.maxpri);
			if (tsplit) {
				tsplit->improve = improve;
				tsplit->var_num = i;
				tsplit->spoint = split;
				tsplit->count = k;
				if (nc == 0) {
					tsplit->spoint = split;
					tsplit->csplit[0] = ct.csplit[0];
				} else
					for (k = 0; k < nc; k++)
						tsplit->csplit[k] = ct.csplit[k];
			}
		}
	}
}
