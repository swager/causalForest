/*
* The four routines for anova + sample variance splitting:
*/
#include "causalTree.h"
#include "causalTreeproto.h"

static double *mean, *sums, *wtsums;
static double *wts, *trs, *trsums;
static double *wtsqrsums, *wttrsqrsums;
static int *countn;
static int *tsplit;

// upper bound for y_max, only for debugging:
// static int MAX = 61;
 // only for debugging

//static int min_node_size = 2;


int
tstatsinit(int n, double *y[], int maxcat, char **error,
	  double *parm, int *size, int who, double *wt, double *treatment)
{
  if (who == 1 && maxcat > 0) {
	graycode_init0(maxcat);
	countn = (int *) ALLOC(2 * maxcat, sizeof(int));
	tsplit = countn + maxcat;
  // mean = (double *) ALLOC(3 * maxcat, sizeof(double));
	mean = (double *) ALLOC(6 * maxcat, sizeof(double));
	wts = mean + maxcat;
  trs = wts + maxcat;
	sums = trs + maxcat;
  // change part
  wtsums = sums + maxcat;
  trsums = wtsums + maxcat;
  
    }
    *size = 1;
    return 0;
}

/*
* The anova evaluation function.  Return the mean and the ss.
*/
void
tstatsass(int n, double *y[], double *value, double *risk, double *wt, double *treatment, double max_y)
{
    int i;
    double temp = 0., temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
    double ttreat = 0.;
    double effect;
    // double ss;

    for (i = 0; i < n; i++) {
      temp1 += *y[i] * wt[i] * treatment[i];
      //temp0 += *y[i] * (1 - wt[i]);
      temp0 += *y[i] * wt[i] * (1 - treatment[i]);
      twt += wt[i];
      ttreat += wt[i] * treatment[i];
    }
    //mean = temp / twt;
    //effect = temp1 / twt - temp0 / (n - twt);
    effect = temp1 / ttreat - temp0 / (twt - ttreat);

    /* 
    ss = 0;
    for (i = 0; i < n; i++) {
      temp = *y[i] - mean;
	    ss += temp * temp * wt[i];
    } 
    */

    *value = effect;
    //*risk = 4 * n * MAX * MAX - effect * effect * n;
    //*risk = n * MAX * MAX - effect * effect * n;
    //max_y = MAX;
    *risk = 4 * n * max_y * max_y - n * effect * effect ;
}

/*
 * The anova splitting function.  Find that split point in x such that
 *  the sum of squares of y within the two groups is decreased as much
 *  as possible.  It is not necessary to actually calculate the SS, the
 *  improvement involves only means in the two groups.
 */
 
void
//anova(int n, double *y[], double *x, int nclass,
//     int edge, double *improve, double *split, int *csplit,
//      double myrisk, double *wt)
// the ct_choose function:

tstats(int n, double *y[], double *x, int nclass,
    int edge, double *improve, double *split, int *csplit,
     double myrisk, double *wt, double *treatment, int minsize)
{
    int i, j;
    double temp;
    double left_sum, right_sum;
    double left_tr_sum, right_tr_sum;
    // add squared sum:
    double left_tr_sqr_sum, right_tr_sqr_sum;
    double left_sqr_sum, right_sqr_sum;
    double tr_var, con_var;
    double left_tr_var, left_con_var, right_tr_var, right_con_var;
    double left_tr, right_tr;
    double left_wt, right_wt;
    int left_n, right_n;
    double best;
    int direction = LEFT;
    int where = 0;
    double node_effect, left_effect, right_effect;
    double left_temp, right_temp;
   // double min_node_size = parm[0];
    int min_node_size = minsize;
    
    right_wt = 0;
    right_tr = 0;
    right_sum = 0;
    right_tr_sum = 0;
    right_sqr_sum = 0;
    right_tr_sqr_sum = 0;
    right_n = n;
    for (i = 0; i < n; i++) {
      right_wt += wt[i];
      right_tr += wt[i] * treatment[i];
      right_sum += *y[i] * wt[i];
      right_tr_sum += *y[i] * wt[i] * treatment[i];
      right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
      right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
     // Rprintf("treatment[%d] = %f,", i, treatment[i] );
    }
    
    temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
    tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
    con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
              - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) / ((right_wt - right_tr) * (right_wt - right_tr));
    //node_effect = temp * temp * n
    node_effect = temp * temp * n - 2 * (tr_var + con_var) * n;
    
    //Rprintf("n = %d, node_effect = %f\n", n, node_effect);
    
    
    if (nclass == 0) {
      /* continuous predictor */
      left_wt = 0;
      left_tr = 0;
      left_n = 0;
      left_sum = 0;
      left_tr_sum = 0;
      left_sqr_sum = 0;
      left_tr_sqr_sum = 0;
    
	    best = 0;
      
	    for (i = 0; right_n > edge; i++) {
	      left_wt += wt[i];
	      right_wt -= wt[i];
        
        left_tr += wt[i] * treatment[i];
        right_tr -= wt[i] * treatment[i];
        
        left_n++;
        right_n--;
        
	      temp = *y[i] * wt[i] * treatment[i];
	      left_tr_sum += temp;
	      right_tr_sum -= temp;
        
        left_sum += *y[i] * wt[i];
        right_sum -= *y[i] * wt[i];
        
        temp = (*y[i]) * (*y[i]) * wt[i] * treatment[i];
        left_tr_sqr_sum += temp;
        right_tr_sqr_sum -= temp;
        
        temp = (*y[i]) * (*y[i]) * wt[i];
        left_sqr_sum += temp;
        right_sqr_sum -= temp;
        
        //Rprintf("left_n = %d, edge = %d\n", left_n, edge);
        //Rprintf("left_tr = %d, left_wt = %d, right_tr = %d, right_wt = %d\n", left_tr, left_wt, right_tr, right_wt);
        //Rprintf("minsize = %d\n", min_node_size);
        // I have a question about weights here in report:
	      if (x[i + 1] != x[i] && left_n >= edge &&
            left_tr >= min_node_size &&
            left_wt - left_tr >= min_node_size &&
            right_tr >= min_node_size &&
            right_wt - right_tr >= min_node_size) {
              //Rprintf("get in!\n");
              
              //taumean = left_tr / left_wt;
              //temp = (left_tr_sum - left_sum * taumean) /
              //((1 - taumean) * left_tr);
              left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
              left_tr_var = left_tr_sqr_sum / left_tr - left_tr_sum  * left_tr_sum / (left_tr * left_tr);
              left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr) 
                             - (left_sum - left_tr_sum) * (left_sum - left_tr_sum)/ ((left_wt - left_tr) * (left_wt - left_tr));        
              
              left_effect = left_temp * left_temp * left_n
              - 2 * (left_tr_var + left_con_var) * left_n;
              
              //taumean = right_tr / right_wt;
              //temp = (right_tr_sum - right_sum * taumean) /
              //((1 - taumean) * right_tr);
              right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
              right_tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
              right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
                             - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) / ((right_wt - right_tr) * (right_wt - right_tr));
              right_effect = right_temp * right_temp * right_n
              - 2 * (right_tr_var + right_con_var) * right_n;    
              
              temp = left_effect + right_effect - node_effect;
              //Rprintf("at %f,leftn: %d, lefteffect: %f, rightn: %d, righteffect: %f\n", x[i], left_n, left_effect,right_n, right_effect, node_effect);
              //Rprintf("current best is %f, and temp improv = %f.\n", best, temp);
              
		          if (temp > best) {
		              best = temp;
		              where = i;
                  //left_temp = left_tr_sum / left_wt - (left_sum - left_wt_sum) / (left_n - left_wt);
                  //right_temp = right_wt_sum / right_wt- (right_sum - right_wt_sum) / (right_n - right_wt);
                  
		              if (left_temp < right_temp)
			              direction = LEFT;
		              else
			              direction = RIGHT;
		          }             
	    }
	   }
     // debug here:
     
     //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	  //*improve = best / myrisk;
    
    *improve = best;
	  if (best > 0) {         /* found something */
//    Rprintf("best = %f, split = %f\n", best, (x[where] + x[where + 1]) / 2 );
	      csplit[0] = direction;
        *split = (x[where] + x[where + 1]) / 2; /* where to split!!!!!!!!! */ 
        }
    }
  
    /*
     * Categorical predictor
     */
    else {
      for (i = 0; i < nclass; i++) {
	      countn[i] = 0;
	      wts[i] = 0;
        trs[i] = 0;
        sums[i] = 0;
        wtsums[i] = 0;
        trsums[i] = 0;
        wtsqrsums[i] = 0;
        wttrsqrsums[i] = 0;
	    }

       /* rank the classes by their mean y value */
       /* RANK THE CLASSES BY THEI */
       // Rprintf("nclass = %d ", nclass);
	    for (i = 0; i < n; i++) {
	        j = (int) x[i] - 1;
          // Rprintf("%d cat, ", j);
	        countn[j]++;
	        wts[j] += wt[i];
          trs[j] += wt[i] * treatment[i];
	        sums[j] += *y[i];
          // adding part
          wtsums[j] += *y[i] * wt[i];
          trsums[j] += *y[i] * wt[i] * treatment[i];
          wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
          wttrsqrsums[j] += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
	    }
      
    	for (i = 0; i < nclass; i++) {
	        if (countn[i] > 0) {
            tsplit[i] = RIGHT;
		        mean[i] = sums[i] / wts[i];
            // mean[i] = sums[i] / countn[i];
            //Rprintf("countn[%d] = %d, mean[%d] = %f\n", i, countn[i], i, mean[i]);
	        } else
            tsplit[i] = 0;
	    }
	    graycode_init2(nclass, countn, mean);

	/*
	 * Now find the split that we want
	 */
   
	left_wt = 0;
  left_tr = 0;
  left_n = 0;
	left_sum = 0;
  left_tr_sum = 0;
  left_sqr_sum = 0;
  left_tr_sqr_sum = 0;
  
	best = 0;
	where = 0;
	while ((j = graycode()) < nclass) {
    //Rprintf("graycode()= %d\n", j);
	    tsplit[j] = LEFT;
	    left_n += countn[j];
	    right_n -= countn[j];
      
	    left_wt += wts[j];
	    right_wt -= wts[j];
      
      left_tr += trs[j];
      right_tr -= trs[j];
      
	    left_sum += wtsums[j];
	    right_sum -= wtsums[j];
      
      left_tr_sum += trsums[j];
      right_tr_sum -= trsums[j];
      
      left_sqr_sum += wtsqrsums[j];
      right_sqr_sum -= wtsqrsums[j];
      
      left_tr_sqr_sum += wttrsqrsums[j];
      right_tr_sqr_sum -= wttrsqrsums[j];

      
	    if (left_n >= edge && right_n >= edge &&
          left_tr >= min_node_size &&
          left_wt - left_tr >= min_node_size &&
          right_tr >= min_node_size &&
          right_wt - right_tr >= min_node_size) {
            
            left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
            left_tr_var = left_tr_sqr_sum / left_tr - left_tr_sum  * left_tr_sum / (left_tr * left_tr);
            left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr) 
                             - (left_sum - left_tr_sum) * (left_sum - left_tr_sum)/ ((left_wt - left_tr) * (left_wt - left_tr));        
            left_effect = left_temp * left_temp * left_n
              - 2 * (left_tr_var + left_con_var) * left_n;
            
             //Rprintf("left_sum = %f, left_wt_sum = %f, left_wt = %f, left_n = %d\n", left_sum, left_wt_sum, left_wt, left_n);             
             right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
             right_tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
             right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
                             - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) / ((right_wt - right_tr) * (right_wt - right_tr));
             right_effect = right_temp * right_temp * right_n
              - 2 * (right_tr_var + right_con_var) * right_n; 
          
             temp = left_effect + right_effect - node_effect;
            //Rprintf("left_n= %d, lefteffect = %f, right_n = %d, righteffect = %f\n", left_n, left_effect, right_n, right_effect);    
            
            //temp = left_sum * left_sum / left_wt +
		        //right_sum * right_sum / right_wt;
		        if (temp > best) {
		            best = temp;
                //left_temp = left_wt_sum / left_wt - (left_sum - left_wt_sum) / (left_n - left_wt);
                //right_temp = right_wt_sum / right_wt- (right_sum - right_wt_sum) / (right_n - right_wt);
		            
                if (left_temp > right_temp)
                  for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
                else
                  for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
		        }
	      }
	 }
   *improve = best;
   // Rprintf("for %f variable, improv = %f\n", x[0], *improve);

	//*improve = best / myrisk;       /* % improvement */
  }
  //Rprintf("for %f variable, improv = %f\n", x[0], *improve);
}
