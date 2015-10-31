/*
 * Run an observation down the tree, and return the prediction error,
 *    for several CP values at once.
 *
 */
#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"
static double INFTY = 9999999999; // infinity

int findNeighbor(int obs, int k) { 
  // k is the starting index of the validation set
  int i, j, temp, neighbor;
  int obs2 = (obs < 0)? -(1 + obs) : obs;
  double dist, min = INFTY;
 // this found is to test whether I can find the neighbor in the validation set!
 // int found = 0;
 
 neighbor = 0;
  
  for (i = k; i < ct.n; i++) {
    j = ct.sorts[0][i];
    temp = (j < 0)? -(1 + j) : j;
    if (ct.treatment[temp] != ct.treatment[obs2]) {
      //found = 1;
      //Rprintf("found one!\n");
      dist = measureDistance(obs2, temp);
      if (dist < min) {
        neighbor = j;
        min = dist;       
      } 
    }       
  }
  // for dubgging only:
  //if(found == 0) Rprintf("There is only one group in validation set!");
  //else Rprintf("two groups in validation! ");
  return neighbor;
}


double measureDistance(int i, int j) {
  int k;
  double distance = 0;
  for (k = 0; k < ct.nvar; k++) {
    distance += (ct.xdata[k][i] - ct.xdata[k][j]) * (ct.xdata[k][i] - ct.xdata[k][j]) / ct.xvar[k];   
  }
  return distance;
}
