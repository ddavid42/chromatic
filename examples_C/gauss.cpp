#include <cstdio>
#include <cstring>
#include "origins.hpp"

std::atomic<uint64_t> BaseOriginVector::atomic_id_counts = 0;

static const int n = 4;
void
print_np_mat_ChNbr(OriginDouble a[n][n+1]) {
  printf("INDEX\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n+1; ++j)
      printf("%5ld ", a[i][j].origins[0].symbol_id);
    printf("\n");
  }
  printf("VALUE\n");
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n+1; ++j)
      printf("%e\t", a[i][j].value);
    printf("\n");
  }
}

int main() {
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

  printf("Initial System");
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
    printf("Idx\n");
    for (int origin_index = 0; origin_index < x[i].origins_size; ++origin_index)
      printf("[%ld, %e]\n", x[i].origins[origin_index].symbol_id, x[i].origins[origin_index].coefficient*pow(2.0, 53));
    printf("\n");
  }
}
