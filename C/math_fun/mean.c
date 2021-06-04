#pragma once
#include "math.h"

double mean(double a[], int n) {
    if(n == 0){return 0.0;}
    double sum = 0;
	int i;
    for(i = 0; i < n; i++)
       sum += a[i];
    return sum / n;
}
