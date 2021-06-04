#pragma once
#include "math.h"


double std(double a[], int n) {
    if(n == 0){return 0.0;}
    double sum = 0;
	int i;
    for(i = 0; i < n; i++)
       sum += a[i];
    double mean = sum / n;
    double sq_diff_sum = 0;
    for(i=0;i<n;i++) {
       double diff = a[i] - mean;
       sq_diff_sum = sq_diff_sum + diff * diff;
    }
    double variance = sq_diff_sum / (n-1);
    return sqrt(variance);
}
