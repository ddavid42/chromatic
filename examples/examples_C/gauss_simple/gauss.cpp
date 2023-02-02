#include <cstdio>
#include <cstring>
#include <fstream>
#include "origins.hpp"

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

int main(int argc, char** argv) {
  double b[n][n+1] = {{ 1/1.,1/2.,1/3.,1/4., 10. },
                      { 1/2.,1/3.,1/4.,1/5., 20. },
                      { 1/3.,1/4.,1/5.,1/6., 30. },
                      { 1/4.*(1+0.001),1/5.,1/6.,1/7., 40. }};

  OriginDouble a[n][n+1];
  memset(a, 0, sizeof(a));
  OriginDouble x[n];
  memset(x, 0, sizeof(x));

  // Reading augmented matrix coefficients
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n+1; ++j)
       a[i][j] = OriginDouble(b[i][j], true);

  printf("Initial System\n");
  print_np_mat_ChNbr(a);
    
  // Applying Gauss Elimination
  for (int i = 0; i < n; ++i) {
    if (a[i][i] == 0.0)
      exit(1);
    for (int j = i+1; j < n; ++j) {
      OriginDouble ratio = a[j][i]/a[i][i];
      for (int k = 0; k < n+1; ++k)
        a[j][k] = a[j][k] - ratio * a[i][k];
    }
  }

  // Back Substitution
  x[n-1] = a[n-1][n]/a[n-1][n-1];

  for (int i = n-2; i >= 0; --i) {
    x[i] = a[i][n];
    for (int j = i+1; j < n; ++j)
       x[i] = x[i] - a[i][j]*x[j];
    x[i] = x[i]/a[i][i];
  }

  // Displaying solution
  printf("\nRequired solution is: \n");
  for (int i = 0; i < n; ++i) {
    printf("X%d value: %e\n", i, x[i].value);
    printf("Idx : [-1, %e]\n", x[i].contribution_without_origin);
    for (int origin_index = 0; origin_index < x[i].contributions_size; ++origin_index)
      printf("[%ld, %e]\n", x[i].contributions[origin_index].symbol_id, x[i].contributions[origin_index].coefficient);
    printf("\n");
  }

  // Plot the Norm1 of each input values on the resulting system (ie. Detailed Condition Number)
  double heatmap1[n][n+1]{};
  for (int i = 0; i < n; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        float v = x[i].contributions[origin].coefficient;
        heatmap1[(k-1)/(n+1)][(k-1)%(n+1)] += v;
     }
  }

  if (argc > 1) {
     std::ofstream out(argv[1]);
     if ( out.is_open() )
        std::cout<<"File " << argv[1] <<" opened"<<std::endl;
     for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            out << heatmap1[i][j] << ", ";
        out << heatmap1[i][n] << '\n';
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[1] <<" closed"<<std::endl;
        out.close();
     }
  }
}

