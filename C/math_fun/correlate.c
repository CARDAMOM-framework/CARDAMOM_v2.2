#pragma once
#include "math.h"
#include "stdlib.h"
#include "malloc.h"
/*************************************************************************/
/*                                                                       */
/* This function computes the correlation between vectors x and y.       */
/* The algorithm is quite fast and is considered numerically stable.     */
/*                                                                       */
/*************************************************************************/
/*                                                                       */
/*            Author: Francesco Pozzi                                    */
/*            Adapted by: AA Bloom 18 Jul 2012
/*            E-Mail: francesco.pozzi@anu.edu.au                         */
/*            Date: 4 December 2008                                      */
/*                                                                       */
/*************************************************************************/

double correlate(double const *x,double const *y, int const dimx, int const dimy)
{

int i;
double sum_sq_x, sum_sq_y, sum_coproduct, mean_x, mean_y, sweep, delta_x, delta_y, pop_sd_x, pop_sd_y, cov_x_y, correlation;





if((dimx == dimy) && (dimx > 1)){

  sum_sq_x = 0;
  sum_sq_y = 0;
  sum_coproduct = 0;
  mean_x = *x;
  mean_y = *y;
  for (i = 1; i < dimx; i++){


    sweep = (double) (i + 0.0) / (i + 1.0);
    delta_x = *(x + i) - mean_x;
    delta_y = *(y + i) - mean_y;
    sum_sq_x = sum_sq_x + delta_x * delta_x * sweep;
    sum_sq_y = sum_sq_y + delta_y * delta_y * sweep;
    sum_coproduct = sum_coproduct + delta_x * delta_y * sweep;
    mean_x = mean_x + delta_x / (i + 1.0);
    mean_y = mean_y + delta_y / (i + 1.0);
  }
  pop_sd_x = sqrt(sum_sq_x / dimx);
  pop_sd_y = sqrt(sum_sq_y / dimx);
  cov_x_y = sum_coproduct / dimx;
  correlation = cov_x_y / (pop_sd_x * pop_sd_y);
}

/***************************************************************************************/
/***************************************************************************************/
return correlation;

}
