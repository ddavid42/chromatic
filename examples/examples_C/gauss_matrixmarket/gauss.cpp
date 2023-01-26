#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>
#include "origins.hpp"
#include "mmio.h"

std::atomic<uint64_t> BaseOriginVector::atomic_id_counts = {0};

static const int n = 4;
void
print_np_mat_ChNbr(OriginDouble a[n][n+1]) {
  printf("INDEX\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n+1; ++j)
      printf("%5ld ", a[i][j].contributions[0].symbol_id);
    printf("\n");
  }
  printf("VALUE\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n+1; ++j)
      printf("%e\t", a[i][j].value);
    printf("\n");
  }
}


void read_matrix(char *input, std::vector< std::vector<OriginDouble> >& a){
  FILE *in;
  MM_typecode matcode;
  int ret_code;
  int M, N;   
  int i, j, nz;
  double v;

  if ((in = fopen(input, "r")) == NULL) 
      exit(1);


  if (mm_read_banner(in, &matcode) != 0)
  {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
  }

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
          mm_is_sparse(matcode) )
  {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(2);
  }
  /* find out size of sparse matrix .... */
  if ((ret_code = mm_read_mtx_crd_size(in, &M, &N, &nz)) !=0)
      exit(3);

  if(M != N){
    printf("Sorry, this application support only square matrices and your is (%d x %d) \n",M,N);
    exit(4);  
  }

  /* reserve memory for matrices a => (M)x(M+1)*/
  a.resize(M);
  for(int i = 0 ; i < M ; ++i)
    a[i].resize(M+1);


  for (int k=0; k<nz; k++){
    // Reading augmented matrix coefficients
    fscanf(in, "%d %d %lg\n", &i, &j, &v);
    /* adjust from 1-based to 0-based */
    i-=1; j-=1;
    // printf("%d %d %d %e \n",k,i,j,v);
    a[i][j]=OriginDouble(v, true);
  }

  if (in !=stdin) fclose(in);

}

int main(int argc, char** argv) {
  std::ofstream out;
  int N;

  // Read MatrixMarket file passed as an input argument 
  if (argc < 3){
		fprintf(stderr, "Usage: %s [matrix-market-input-filename] [chromatic-output-data]\n", argv[0]);
		exit(1);
	} 
  
  std::vector<std::vector<OriginDouble>> a;

  read_matrix(argv[1], a);
  N = a.size();
  OriginDouble x[N];
  memset(x, 0, sizeof(x));

  
  // Set an arbitrary vector results
  for (int i = 0; i < N; ++i)
    a[i][N] = OriginDouble(i, true);

    
  // Applying Gauss Elimination
  for (int i = 0; i < N; ++i) {
    printf("Gauss, line :%d\n",i);
    if (a[i][i] == 0.0){
      printf("Gaussian elemination not possible a[i][i] == 0.0 for i= %d\n",i);
      exit(1);
    }
    for (int j = i+1; j < N; ++j) {
      OriginDouble ratio = a[j][i]/a[i][i];
      for (int k = 0; k < N+1; ++k)
        a[j][k] = a[j][k] - ratio * a[i][k];
    }
  }

  // Back Substitution
  x[N-1] = a[N-1][N]/a[N-1][N-1];

  for (int i = N-2; i >= 0; --i) {
    x[i] = a[i][n];
    for (int j = i+1; j < N; ++j)
       x[i] = x[i] - a[i][j]*x[j];
    x[i] = x[i]/a[i][i];
  }

  // Displaying solution
  printf("\nRequired solution is: \n");
  for (int i = 0; i < N; ++i) {
    printf("X%d value: %e\n", i, x[i].value);
    printf("Idx : [-1, %e]\n", x[i].contribution_without_origin);
    for (int origin_index = 0; origin_index < x[i].contributions_size; ++origin_index)
      printf("[%ld, %e]\n", x[i].contributions[origin_index].symbol_id, x[i].contributions[origin_index].coefficient);
    printf("\n");
  }

  // Plot the Norm1 of each input values on the resulting system (ie. Detailed Condition Number)
  double heatmap1[N][N+1]{};
  for (int i = 0; i < N; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        float v = x[i].contributions[origin].coefficient;
        heatmap1[(k-1)/(N+1)][(k-1)%(N+1)] += v;
     }
  }

  if (argc > 1) {
     std::ofstream out(argv[2]);
     if ( out.is_open() )
        std::cout<<"File " << argv[2] <<" opened"<<std::endl;
     for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            out << heatmap1[i][j] << ", ";
        out << heatmap1[i][n] << '\n';
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[2] <<" closed"<<std::endl;
        out.close();
     }
  }
}

