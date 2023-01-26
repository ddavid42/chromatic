
/*
Program: Gauss Elimination Method
All array indexes are assumed to start from 1
*/

#include<iostream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>
#include <chrono>

#define   SIZE   1000

using namespace std;

int main(){
	float a[SIZE+10][SIZE+10], x[SIZE+10], ratio;
	int i,j,k;
    clock_t start, end;

     /* Setting precision and writing floating point values in fixed-point notation. */
     cout<< setprecision(3)<< fixed;


	 for(i=1;i<=SIZE;i++)	 {
		  for(j=1;j<=SIZE+1;j++){
			   a[i][j] = 1./(i+j);
		  }
	 }

    start = clock();

	/* Applying Gauss Elimination */
	 for(i=1;i<=SIZE-1;i++)	 {
		  if(a[i][i] == 0.0)		  {
			   cout<<"Mathematical Error!";
			   exit(0);
		  }
		  for(j=i+1;j<=SIZE;j++)		  {
			   ratio = a[j][i]/a[i][i];

			   for(k=1;k<=SIZE+1;k++){
			  		a[j][k] = a[j][k] - ratio*a[i][k];
			   }
		  }
	 }
	 /* Obtaining Solution by Back Substitution Method */
	 x[SIZE] = a[SIZE][SIZE+1]/a[SIZE][SIZE];

	 for(i=SIZE-1;i>=1;i--)
	 {
		  x[i] = a[i][SIZE+1];
		  for(j=i+1;j<=SIZE;j++)
		  {
		  		x[i] = x[i] - a[i][j]*x[j];
		  }
		  x[i] = x[i]/a[i][i];
	 }

    end = clock();

    printf ("max_index: %0.8f sec\n",
            ((float) end - start)/CLOCKS_PER_SEC);

    if((end-start)>0)
        exit(-1);

	 /* Displaying Solution */
	 cout<< endl<<"Solution: "<< endl;
	 for(i=1;i<=SIZE;i++)
	 {
	  	cout<<"x["<< i<<"] = "<< x[i]<< endl;
	 }

	 return(0);
}