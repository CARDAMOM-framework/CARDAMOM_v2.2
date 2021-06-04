#pragma once
#include "math.h"


double covariance(double a[], double b[], int n) {
    if(n == 0){return 0.0;}
    double sum_a = 0, sum_b=0;
	int i;
    for(i = 0; i < n; i++){
       sum_a += a[i];
       sum_b += b[i];}

    double mean_a = sum_a /(double)n;
    double mean_b = sum_b /(double)n;
    double sq_diff_sum = 0;
    for(i=0;i<n;i++) {
       sq_diff_sum = sq_diff_sum + (a[i] - mean_a) * (b[i] - mean_b);
    }
    return sq_diff_sum / ((double)n-1);
}
