#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <cassert>
#include "mmio.h"

bool isLessIndex(const std::pair<int, double>& first, const std::pair<int, double>& second)
   { return first.first < second.first; }

class SparseMatrix {
  public:
   int lines = 0, columns = 0;
  private:
   std::vector< std::vector<std::pair<int, double>> > content;

  public:
   SparseMatrix() = default;
   SparseMatrix(const SparseMatrix&) = default;
   SparseMatrix(SparseMatrix&&) = default;
   SparseMatrix& operator=(const SparseMatrix&) = default;
   SparseMatrix& operator=(SparseMatrix&&) = default;

   void setSize(int line, int column) { content.resize(line); lines = line; columns = column; }
   int getLines() const { return lines; }
   int getColumns() const { return columns; }
   void set(int i, int j, double value)
      {  auto& line = content[i];
         if (line.empty())
            line.push_back(std::make_pair(j, value));
         else {
            auto columnFound = std::lower_bound(line.begin(), line.end(), std::make_pair(j, double(0.0)), &isLessIndex);
            if (columnFound == line.end() || columnFound->first > j)
               line.insert(columnFound, std::make_pair(j, value));
            else
               columnFound->second = value;
         }
      }
   class ColumnIterator {
     private:
      std::vector<std::pair<int, double>>::const_iterator iter, iterEnd;
      int beforeIter;
      int columns;

     public:
      ColumnIterator(const std::vector<std::pair<int, double>>& ref, int acolumns, bool isFirstValue=false)
         :  iter(ref.begin()), iterEnd(ref.end()),
            beforeIter((ref.empty() || isFirstValue) ? 0 : ref[0].first), columns(acolumns) {}
      ColumnIterator(const std::vector<std::pair<int, double>>& ref, int j, int acolumns)
         :  iter(ref.begin()), iterEnd(ref.end()),
            beforeIter(ref.empty() ? 0 : ref[0].first), columns(acolumns)
         { iter = std::lower_bound(iter, iterEnd, std::make_pair(j, double(0.0)), &isLessIndex);
           if (iter == iterEnd)
              beforeIter = columns-j;
           else
              beforeIter = iter->first-j;
         }
      ColumnIterator(const ColumnIterator&) = default;
      ColumnIterator& operator=(const ColumnIterator&) = default;

      double getValue() const
         {  if (beforeIter) return 0.0;
            return iter->second;
         }
      bool hasValue() const { return iter != iterEnd && !beforeIter; }
      const double& getSValue() const { assert(!beforeIter); return iter->second; }
      bool isValid() const { return iter != iterEnd || beforeIter; }
      ColumnIterator& operator++()
         {  assert(isValid());
            if (beforeIter) --beforeIter;
            else {
               auto oldIter = iter;
               ++iter;
               if (iter == iterEnd)
                  beforeIter = columns - oldIter->first - 1;
               else
                  beforeIter = iter->first - oldIter->first - 1;
            }
            return *this;
         }
      int getIndex() const { assert(isValid()); return ((iter != iterEnd) ? iter->first : columns) - beforeIter; }
      ColumnIterator& setToNextValue()
         {  assert(isValid());
            if (beforeIter) beforeIter = 0;
            else ++iter;
            return *this;
         }
   };
   class RefColumnIterator {
     private:
      std::vector<std::pair<int, double>>& ref;
      std::vector<std::pair<int, double>>::iterator iter, iterEnd;
      int beforeIter;
      int columns;

     public:
      RefColumnIterator(std::vector<std::pair<int, double>>& aref, int acolumns)
         :  ref(aref), iter(ref.begin()), iterEnd(ref.end()) {}
      RefColumnIterator(std::vector<std::pair<int, double>>& aref, int j, int acolumns)
         :  ref(aref), iter(aref.begin()), iterEnd(aref.end()),
            beforeIter(aref.empty() ? 0 : aref[0].first), columns(acolumns)
         { iter = std::lower_bound(iter, iterEnd, std::make_pair(j, double(0.0)), &isLessIndex);
           if (iter == iterEnd)
              beforeIter = columns-j;
           else
              beforeIter = iter->first-j;
         }
      RefColumnIterator(const RefColumnIterator&) = default;
      RefColumnIterator& operator=(const RefColumnIterator&) = default;

      double& getSValue()
         {  int index = (iter == iterEnd) ? columns : iter->first;
            if (beforeIter) {
               index -= beforeIter;
               iter = ref.insert(iter, std::make_pair(index, double(0.0)));
               iterEnd = ref.end();
               beforeIter = 0;
            }
            return iter->second;
         }
      bool hasValue() const { return iter != iterEnd && !beforeIter; }
      bool isValid() const { return iter != iterEnd || beforeIter; }
      RefColumnIterator& operator++()
         {  assert(isValid());
            if (beforeIter) --beforeIter;
            else {
               auto oldIter = iter;
               ++iter;
               if (iter == iterEnd)
                  beforeIter = columns - oldIter->first - 1;
               else
                  beforeIter = iter->first - oldIter->first - 1;
            }
            return *this;
         }
      RefColumnIterator& setToNextIndex(int index)
         {  assert(isValid());
            int currentIndex = ((iter != iterEnd) ? iter->first : columns) - beforeIter;
            int shift = index - currentIndex;
            assert(shift > 0);
            while (beforeIter < shift && iter != iterEnd) {
               shift -= (beforeIter+1);
               int oldIndex = iter->first;
               ++iter;
               beforeIter = ((iter != iterEnd) ? iter->first : columns) - oldIndex - 1;
            }
            if (beforeIter >= shift)
               beforeIter -= shift;
            else
               beforeIter = 0;
            return *this;
         }

   };
   class ReverseColumnIterator {
     private:
      std::vector<std::pair<int, double>>::const_reverse_iterator iter, iterEnd;
      int afterIter;
      int columns;

     public:
      ReverseColumnIterator(const std::vector<std::pair<int, double>>& ref, int acolumns)
         :  iter(ref.rbegin()), iterEnd(ref.rend()),
            afterIter(ref.empty() ? 0 : (acolumns-ref.back().first-1)), columns(acolumns) {}
      ReverseColumnIterator(const ReverseColumnIterator&) = default;
      ReverseColumnIterator& operator=(const ReverseColumnIterator&) = default;

      double getValue() const
         {  if (afterIter) return 0.0;
            return iter->second;
         }
      bool hasValue() const { return iter != iterEnd && !afterIter; }
      const double& getSValue() const { return iter->second; }
      bool isValid() const { return iter != iterEnd || afterIter; }
      ReverseColumnIterator& operator++()
         {  assert(isValid());
            if (afterIter) --afterIter;
            else {
               auto oldIter = iter;
               ++iter;
               if (iter == iterEnd)
                  afterIter = oldIter->first - 1;
               else
                  afterIter = oldIter->first - iter->first - 1;
            }
            return *this;
         }
      int getIndex() const { assert(isValid()); return ((iter != iterEnd) ? iter->first : 0) + afterIter; }
      ReverseColumnIterator& setToNextValue()
         {  assert(isValid());
            if (afterIter) afterIter = 0;
            else ++iter;
            return *this;
         }
   };

   ColumnIterator getColumnsStartingAt(int i, int j) const
     {  auto& line = content[i];
        if (!line.empty() && line.front().first == j)
          return ColumnIterator(line, columns, true);
        return ColumnIterator(line, j, columns);
     }
   RefColumnIterator getRefColumnsStartingAt(int i, int j)
     {  auto& line = content[i];
        return RefColumnIterator(line, j, columns);
     }
   ReverseColumnIterator getReverseColumns(int i) const
     {  auto& line = content[i];
        return ReverseColumnIterator(line, columns);
     }
};

void read_matrix(char *input, SparseMatrix& a){
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

  a.setSize(M, N+1);
  /* reserve memory for matrices a => (M)x(M+1)*/
//a.resize(M);
//for(int i = 0 ; i < M ; ++i)
//  a[i].resize(M+1);

  for (int k=0; k<nz; k++){
    // Reading augmented matrix coefficients
    fscanf(in, "%d %d %lg\n", &i, &j, &v);
    /* adjust from 1-based to 0-based */
    i-=1; j-=1;
    // printf("%d %d %d %e \n",k,i,j,v);
    a.set(i, j, v);
    a.set(j, i, v);
  }

  if (in !=stdin) fclose(in);

}

int main(int argc, char** argv) {
  std::ofstream out;
  int N;
  clock_t start, end;

  // Read MatrixMarket file passed as an input argument 
  if (argc < 2){
		fprintf(stderr, "Usage: %s [matrix-market-input-filename]\n", argv[0]);
		exit(1);
	} 
  
  SparseMatrix a;

  read_matrix(argv[1], a);
  N = a.getLines();
  double x[N];
  memset(x, 0, sizeof(x));

  
  // Set an arbitrary vector results
  for (int i = 0; i < N; ++i)
    a.set(i, N, double(i));

  start = clock(); /* Lancement de la mesure */

  // Applying Gauss Elimination
  for (int i = 0; i < N; ++i) {
    printf("Gauss, line :%d\n",i);
    SparseMatrix::ColumnIterator firstColumnIterator = a.getColumnsStartingAt(i, i);
    if (!firstColumnIterator.hasValue()) {
      printf("Gaussian elemination not possible a[i][i] == 0.0 for i= %d\n",i);
      exit(1);
    }
    for (int j = i+1; j < N; ++j) {
      SparseMatrix::RefColumnIterator secondColumnIterator = a.getRefColumnsStartingAt(j, i);
      if (!secondColumnIterator.hasValue())
        continue;
      double ratio = secondColumnIterator.getSValue() / firstColumnIterator.getValue(); // a[j][i]/a[i][i];
      if (ratio == 0.0)
         continue;
      auto columnIterator = firstColumnIterator;
      while (columnIterator.setToNextValue().isValid()) {
        secondColumnIterator.setToNextIndex(columnIterator.getIndex());
        secondColumnIterator.getSValue() -= ratio * columnIterator.getValue();
      }
    }
  }

  // Back Substitution
  SparseMatrix::ReverseColumnIterator reverseColumnIterator = a.getReverseColumns(N-1);
  double res = reverseColumnIterator.getValue();
  ++reverseColumnIterator;
  
  x[N-1] = res/reverseColumnIterator.getValue();

  for (int i = N-2; i >= 0; --i) {
    SparseMatrix::ReverseColumnIterator reverseColumnIterator = a.getReverseColumns(i);
    x[i] = reverseColumnIterator.getValue();
    int index;
    while (reverseColumnIterator.setToNextValue().isValid()
         && (index = reverseColumnIterator.getIndex()) > i)
//  for (int j = N-1; j > i; --j)
//     ++reverseColumnIterator;
//     if (reverseColumnIterator.hasValue())
       x[i] -= x[reverseColumnIterator.getIndex()]*reverseColumnIterator.getValue();
//  }
//  ++reverseColumnIterator;
    x[i] /= (index == i) ? reverseColumnIterator.getValue() : 0.0;
  }

  end = clock();  /* Arret de la mesure */
  printf("\nTIME: %lf",((double)end - start) / CLOCKS_PER_SEC);
 
  // Displaying solution
  printf("\nRequired solution is: \n");
  for (int i = 0; i < N; ++i) {
    printf("X%d value: %e\n", i, x[i]);
  }
}

