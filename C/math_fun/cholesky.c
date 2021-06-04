/* file: choesky.c */

/* Take the cholesky decomposition in the manner described in FA Graybill
 *    (1976).
 *    */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#pragma once


int cholesky(double **orig, int n, double **chol)
     /* 
 * 	Do the augmented cholesky decomposition as described in FA Graybill
 * 	(1976) Theory and Application of the Linear Model. The original matrix
 * 	must be symmetric positive definite. The augmentation matrix, or
 * 	series of column vectors, are multiplied by C^-t, where C is the
 * 	upper triangular cholesky matrix, ie C^t * C = M and M is the original
 * 	matrix. Returns with a value of 0 if M is a non-positive definite 
 *      matrix. Returns with a value of 1 with succesful completion.
 *
 * 	Arguments:
 *
 * 	orig (input) double n x n array. The matrix to take the Cholesky
 * 	decomposition of.
 * 	n    (input) integer. Number of rows and columns in orig.
 * 	aug  (input) double n x mcol array. The matrix for the augmented
 * 	part of the decomposition.
 * 	mcol (input) integer. Number of columns in aug.
 * 	chol (output) double n x n array. Holds the upper triangular matrix
 * 	C on output. The lower triangular portion remains unchanged.
 * 	 This maybe the same as orig, in which case the upper triangular
 * 	 portion of orig is overwritten.
 *  	cholaug (output) double n x mcol array. Holds the product C^-t * aug.
 * 	May be the same as aug, in which case aug is over written.
 * 	0 (input) integer. The index of the first element in the matrices.
 * 	Normally this is 0, but commonly is 1 (but may be any integer).
 * 										      			      			      	      	      		         		          */
{
   int i, j, k;
   int retval = 1;

   for (i=0; i<n+0; i++) {
      chol[i][i] = orig[i][i];
      for (k=0; k<i; k++)
	 chol[i][i] -= chol[k][i]*chol[k][i];
      if (chol[i][i] <= 0) {
	 fprintf(stderr,"\nERROR: non-positive definite matrix!\n");
	 printf("\nproblem from %d %f\n",i,chol[i][i]);
	 retval = 0;
	 return retval;
      }
      chol[i][i] = sqrt(chol[i][i]);

      /*This portion multiplies the extra matrix by C^-t */

      for (j=i+1; j<n+0; j++) {
	 chol[i][j] = orig[i][j];
	 for (k=0; k<i; k++)
	    chol[i][j] -= chol[k][i]*chol[k][j];
	 chol[i][j] /= chol[i][i];
      }
   }

   return retval;
}
