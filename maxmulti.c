//Matrix Multiplication 
//COMS 480
//Michael Davis

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void matMult(double **A, double **B, double**C, int n)
{
  
  for( int r = 0; r < n; r++)
      for(int c = 0; c < n; c++)
	for(int k = 0; k < n; k++)
	  C[r][c] += A[r][k] * B[k][c];       
}

void blockMatMulti(double **A, double **B, double**C, int n, int blockSize)
{
  
  double sum;
  int br, bc, r, c, k;
  int myMin(int a, int b)
  {
    return (a > b) ? b : a;
  }
  
  for (br = 0; br < n; br += blockSize){
      for(bc = 0; bc < n; bc += blockSize){
	  for(r = 0; r < n; ++r){
	      for(c = br; c < myMin(br + blockSize, n); ++c){
		  sum = 0.0;
		  for(k = bc; k < myMin(bc + blockSize, n); ++k)
		    sum += A[r][k] * B[k][c];
		    C[r][c] += sum;
		}
	    }
	}
    }
}

int main()
{
  int  r ,c , set, i, size, block;
  double val;
  clock_t t;
 
  printf("Enter size of the matrix: ");
  scanf("%d", &size);

  printf("\nEnter size of Block: ");
  scanf("%d", &block);
   
  int nRows = size;
  int nCols = size;
 
      
  //Allocate Matrix 1
      double** mat1 = (double**)malloc(nRows * sizeof(double*));
      mat1[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat1[r] = &(mat1[0][r * nCols]);
      }

  //Allocate Matrix 2
      double** mat2 = (double**)malloc(nRows * sizeof(double*));
      mat2[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat2[r] = &(mat2[0][r * nCols]);
      }

  //Allocate Matrix 3
      double** mat3 = (double**)malloc(nRows * sizeof(double*));
      mat3[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat3[r] = &(mat3[0][r * nCols]);
      }
 
  //Allocate Matrix 4
      double** mat4 = (double**)malloc(nRows * sizeof(double*));
      mat4[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat4[r] = &(mat4[0][r * nCols]);
      }

  //Allocate Matrix 5
      double** mat5 = (double**)malloc(nRows * sizeof(double*));
      mat5[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat5[r] = &(mat5[0][r * nCols]);
      }
  
  // Fill Matrix
      for(int r = 0; r < nRows; ++r)
	for(int c = 0; c < nCols; ++c){
	  mat1[r][c] = r * nCols + c + 1;
	}

  // Printing Values of Matrix 1
  /*
  printf("\nThe data of Matrix 1:\n");
  for( int r = 0; r < nRows; ++r){
    for(int c = 0; c < nCols; ++c)
      printf("%lf\t", mat1[r][c]);
    printf("\n");
  }
  */

      for(int r = nRows -1; r >= 0; r--)
	for(int c = nCols - 1; c >= 0; c--)
	  mat2[nRows - 1 - r][nCols - 1 - c] = r * nCols + c + 1;
  
  /*
  //Priting values of Matrix 2
  printf("\nThe data of Matrix 2:\n");
  for( int r = 0; r < nRows; ++r){
    for( int c = 0; c < nCols; ++c)
      printf("%lf\t", mat2[r][c]);
    printf("\n");
  }
  */

      t = clock();
  
  // Matrix Multiplication between Mat1 and Mat2 then placed values into Mat
      for(int r = 0; r < nRows; ++r)
	for(int c = 0; c < nCols; ++c){
	  mat3[r][c] = 0.0;
	  for(int i = 0; i < nRows; ++i)
	    mat3[r][c] += mat1[r][i] * mat2[i][c];
	}
      t = clock() - t;
      double time_taken = ((double)t)/CLOCKS_PER_SEC;
      printf("\nMM size %d took %f seconds to complete\n",size, time_taken);
     
      
  
  /*
  //Priting out Matrix 3 with multiplication results 
  printf("\nResults of Multiplication...printing Matrix 3:\n");
  for( int r = 0; r < nRows; ++r){
    
    for(int c = 0; c < nCols; ++c)
      printf("%lf\t", mat3[r][c]);
  
    printf("\n");
  }
  */
      t = clock();
  
  //Transpose Matrix 2
      for (int r = 0; r < nRows; ++r)
	for(int c = r + 1; c < nCols; ++c){
	  val = mat2[r][c];
	  mat2[r][c] = mat2[c][r];
	  mat2[c][r] = val;

	}

      t = clock() - t;
      double transpose_time = ((double)t)/CLOCKS_PER_SEC;
      printf("\nTranspose of size %d took %f seconds to complete\n", size, transpose_time);
      
  /*
  //Printing results of Tranpose
  printf("\nResult of Transpose:\n");
  for(int r = 0; r < nRows; ++r){
    for(int c = 0; c < nCols; ++c)
      printf("%lf\t", mat2[r][c]);
    printf("\n");
  }
  */

      t = clock();
  
  //Transpose Multiplication
      for(int r = 0; r < nRows; ++r)
	for(int c = 0; c < nCols; ++c){
	  mat4[r][c] = 0.0;
	  for(int i = 0; i < nRows; ++i)
	    mat4[r][c] += mat2[r][c] * mat3[c][r];
	}
      t = clock() - t;
      double tranposemulti_time = ((double)t)/CLOCKS_PER_SEC;
      printf("\nTranspose Multi  of size %d took %f seconds to complete\n",size,  tranposemulti_time);
     


 /*
 //Priting out Matrix 4
  printf("\nResult of Transpose Multiplication...Printing out Matrix 4:\n");
  for(int r = 0; r < nRows; ++r){
    for(int c = 0; c < nCols; ++c)
      printf("%lf\t", mat4[r][c]);
    printf("\n");
  }
  */
 
  // Testing Block Multiplication
      t = clock();
      blockMatMulti(mat1, mat2, mat5, size, block);
      t = clock() - t;
      double block_time = ((double)t)/CLOCKS_PER_SEC;
      printf("\nBlock Multiplication of block size %d took %f seconds to complete\n", block,  block_time);
     

  /*
  printf("\n Printing out Matrix 5:\n");
  for(int r = 0; r < nRows; ++r){
    for(int c = 0; c < nCols; ++c)
      printf("%lf\t", mat5[r][c]);
    printf("\n");
  }
  */
  
      free(mat1);
      free(mat2);
      free(mat3);
      free(mat4);
      free(mat5);
    
  return 0;
 
  
  
}
