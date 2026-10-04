#include "gaussjordan.h"
#include <stdlib.h>   
#include <math.h> 

void gauss_jordan_solve(const double *XtX, const double *Xty,
                                double *beta, int p)
{
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

      for (int i = 0; i < p; i++)
      {
         if (i == k){
            continue;
         }
         double factor = A[i * (p + 1) + k] / A[k * (p + 1) + k];
         for (int j = k; j < p + 1; j++)
         {
            A[i * (p + 1) + j] = A[i * (p + 1) + j] - factor * A[(p + 1) * k + j];
         }
      }
   }

   for (int i = 0; i < p; i++) {
      beta[i] = A[i * w + p] / A[i * w + i];
   }


   free(A);
}
