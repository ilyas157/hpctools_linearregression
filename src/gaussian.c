#include "gaussian.h"
#include <stdlib.h>   
#include <math.h> 

/* -------------------------------------------------------------------------
 * TODO (STUDENT): gaussian_elimination_solve
 *
 * Solve the p x p system:
 *
 *   XtX * beta = Xty
 *
 * using Gaussian elimination with partial pivoting, followed by back
 * substitution:
 *
 *   1. Build an augmented p x (p+1) matrix [XtX | Xty] (work on a local
 *      copy — do not modify XtX/Xty in place, you may want to keep them
 *      for the report).
 *   2. Forward elimination: for each pivot column k = 0..p-1,
 *        a. partial pivoting: find the row r >= k with the largest
 *           absolute value in column k, and swap rows k and r if r != k
 *           (this avoids dividing by a very small/zero pivot).
 *        b. eliminate column k from all rows below k by subtracting an
 *           appropriate multiple of row k.
 *   3. Back substitution: once the augmented matrix is in upper
 *      triangular form, solve for beta[p-1], beta[p-2], ..., beta[0]
 *      from the bottom row upward.
 *
 * XtX  : p x p, row-major (read-only)
 * Xty  : p (right-hand side, read-only)
 * beta : p (output, caller-allocated)
 * ---------------------------------------------------------------------- */
void gaussian_elimination_solve(const double *XtX, const double *Xty,
                                double *beta, int p) {

   /* TODO: implement Gaussian elimination with partial pivoting +
    * back substitution here. A scratch p x (p+1) augmented matrix can be
    * allocated locally with malloc and freed before returning.
    */
   int w = p + 1;
   double *A = malloc((size_t)p * w * sizeof(double));

   // fill A = [XtX | Xty]
   for (int i = 0; i < p; i++)
   {
      for (int j = 0; j < p; j++)
      {
         A[i * (p + 1) + j] = XtX[i * p + j];
      }
      A[i * (p + 1) + p] = Xty[i];
   }

   for (int k = 0; k < p; k++)
   {
      // find the pivot
      int pivot = k;
      for (int i = k + 1; i < p; i++)
      {
         if (fabs(A[i * (p + 1) + k]) > fabs(A[pivot * (p + 1) + k]))
         {
            pivot = i;
         }
      }

      // swap with row0
      if (pivot != k)
      {
         double temp;
         for (int j = 0; j < p + 1; j++)
         {
            temp = A[k * (p + 1) + j];
            A[k * (p + 1) + j] = A[pivot * (p + 1) + j];
            A[pivot * (p + 1) + j] = temp;
         }
      }
      for (int i = k + 1; i < p; i++)
      {
         double factor = A[i * (p + 1) + k] / A[k * (p + 1) + k];
         for (int j = k; j < p + 1; j++)
         {
            A[i * (p + 1) + j] = A[i * (p + 1) + j] - factor * A[(p + 1) * k + j];
         }
      }
   }

   for (int i = p - 1; i >= 0; i--){
      double sum = A[i * (p+1) + p];
      for (int j = i + 1; j < p; j++)
      {
         sum -= A[i* (p+1) + j] * beta[j];
      }
      beta[i] = sum / A[i* (p+1) + i];
   }


   free(A);
}
