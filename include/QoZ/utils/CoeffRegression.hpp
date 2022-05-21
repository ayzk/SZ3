
 
//#include <stdio.h>
//#include <stdlib.h>

#ifndef SZ_MATRIX_OPERATION_HPP
#define SZ_MATRIX_OPERATION_HPP
namespace QoZ {

    inline void matrixTranspose(double *a, double *b, size_t row, size_t column)
    {//row and column are for the matrix a
        size_t i,j;
        for(i = 0; i<column; i++){
            for(j = 0; j<row; j++){
                b[i*row+j] = a[j*column+i];
            }
        }
    } 

    inline void matrixMul(double *A, double *B, double *C, size_t rowA, size_t columnB, size_t columnA)
    {
        for (size_t i=0;i<rowA;i++){
            for (size_t j=0; j<columnB;j++){
                C[i*columnB+j] = 0;
                for (size_t k=0;k<columnA;k++){
                    C[i*columnB+j] +=A[i*columnA+k]*B[k*columnB+j];
                }
             }
         }
    }

    inline void matrixVecMul(double *A, double *b, double *c, size_t rowA, size_t columnA)
    {
        for (size_t i=0;i<rowA;i++){
    		c[i] = 0;
    		for (size_t k=0;k<columnA;k++){
    			c[i] += A[i*columnA+k]*b[k];
    		}
         }
    }

    int Gauss(double *a, double *b, size_t size, double** result)
    {
    	*result = new double[size];
        double A[12][12];//the matrix A storing elements
        double B[12];//the matrix sotring B 

        int i, j, k; 
        for(i=0;i<size;i++)
        {
            for(j=0;j<size;j++)
            {
                A[i][j]=a[i*size+j];
            }
        }

        for(i=0; i<size;i++)
        {
            B[i]=b[i];
        } 

        //judge whether we can use Guass(i.e., whether the diagonal elements are all 0.
        for(i=0;i<size;i++)
        {
            if(A[i][i]==0)
            {
                //printf("This matrix cannot be solved by Gauss\n");
                return 1; //error
            }
        } 

        float C[size];//Store coefficients

        //Elimitation method
        for(k=0; k<size-1;k++)
        {
            for(i=k+1;i<size;i++)
            {
                C[i]=A[i][k]/A[k][k];
            }

            for(i=k+1;i<size;i++)
            {
                for(j=0;j<size;j++)
                {
                    A[i][j]=A[i][j]-C[i]*A[k][j];   
                }
                B[i]=B[i]-C[i]*B[k];
            }
        } 

        //Use X to store the results
        float X[size];
        X[size-1]=B[size-1]/A[size-1][size-1];

        //Get the values of unknowns
        for(i=size-2;i>=0;i--)
        {
            double Sum=0;
            for(j=i+1;j<size;j++)
            {
                Sum+=A[i][j]*X[j];
            }
            X[i]=(B[i]-Sum)/A[i][i];
        }

    	for(i=0;i<size;i++)
    		(*result)[i] = X[i];

    	return 0;
    } 
     /*
    void printMatrix(double* matrix, int m, int n)
    {
    	for(int i = 0;i<m;i++)
    	{
    		for(int j = 0;j<n;j++)
    		{
    			printf("%f\t", matrix[i*n+j]);
    		}
    		printf("\n");
    	}
    }

    void printVector(double* vector, int n)
    {
    	for(int i = 0;i<n;i++)
    		printf("%f\t", vector[i]);
    	printf("\n");		
    }
    */
    double* Regression(double * A,size_t numPoints,size_t numFeatures,double * b)
    {
        
        double* AT = new double[numPoints*numFeatures]; //transpose
        
        matrixTranspose(A, AT, numPoints, numFeatures);
        
        double* c = new double [numFeatures];
        matrixVecMul(AT, b, c, numFeatures, numPoints);
        
        double* ATA = new double [numFeatures*numFeatures];
        matrixMul(AT,A, ATA, numFeatures, numFeatures, numPoints);
        /*
        printMatrix(ATA, interpSize, interpSize);
        printf("-----------------------------------\n");
        printVector(c, interpSize);
        */
        double* result = NULL;
        Gauss(ATA, c, numFeatures, &result);
        
        //delete A(extractedMatrix);
        delete []AT;
        delete []ATA;  
        //free(b);
        delete []c;
        
        return result;
    }
}

#endif