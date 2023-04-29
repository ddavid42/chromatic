#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <algorithm>
#include <cassert>
#include "mmio.h"
#include "workplan.h"

#ifdef _ORIGINS
void
print_np_mat_ChNbr(double** a/*[n][n+1]*/, int n) {
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
#endif

bool
setToSymbolId(uint64_t symbol, WorkPlan::Cursor& cursor) {
   if (!cursor.empty()) {
      QuadTreeElement* current = cursor.back();
      while (current->hasContent()) {
         if (!cursor.setToNext()) return false;
         current = cursor.back();
      }
      if (current->sharedReference().contributions[0].symbol_id <= symbol) {
         while (current->sharedReference().contributions[0].symbol_id < symbol) {
            if (!cursor.setToNext()) return false;
            current = cursor.back();
            while (current->hasContent()) {
               if (!cursor.setToNext()) return false;
               current = cursor.back();
            }
         }
         if (current->sharedReference().contributions[0].symbol_id != symbol) return false;
         return true;
      }
   }
   if (!cursor.setToFirst()) return false;
   QuadTreeElement* current = cursor.back();
   while (current->hasContent()) {
      if (!cursor.setToNext()) return false;
      current = cursor.back();
   }
   while (current->sharedReference().contributions[0].symbol_id < symbol) {
      if (!cursor.setToNext()) return false;
      current = cursor.back();
      while (current->hasContent()) {
         if (!cursor.setToNext()) return false;
         current = cursor.back();
      }
   }
   if (current->sharedReference().contributions[0].symbol_id != symbol) return false;
   return true;
}

typedef std::vector<std::pair<int, int>> IdDictionary;

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
   double& set(int i, int j, double value)
      {  auto& line = content[i];
         double* result;
         if (line.empty()) {
            line.push_back(std::make_pair(j, value));
            result = &line.back().second;
         }
         else {
            auto columnFound = std::lower_bound(line.begin(), line.end(), std::make_pair(j, double(0.0)), &isLessIndex);
            if (columnFound == line.end() || columnFound->first > j) {
               columnFound = line.insert(columnFound, std::make_pair(j, value));
            else
               columnFound->second = value;
            result = &columnFound->second;
         }
         return result;
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

void read_matrix(char *input, SparseMatrix& a, WorkPlan* plan, IdDictionary* dictionary){
  FILE *in;
  MM_typecode matcode;
  int ret_code;
  int M, N;   
  int i, j, nz;

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

if (plan) {

#if defined(_ORIGINS) && defined(_MATRIX)
  WorkPlan::Cursor cursor(*plan);
  if (cursor.setToFirst()) {
    do {
      if (!cursor.elementAt()->hasContent() && !cursor.elementAt()->hasContribution()) {
        cursor.elementAt()->sharedReference() = double(1.0, true);
        cursor.elementAt()->sharedReference().contributions[0].is_atomic_var
          = cursor.elementAt()->getWidth() == 1 && cursor.elementAt()->getLength() == 1;
      }
    } while (cursor.setToNext());
  }
//double ref_val[4][4];
//for (int i = 0; i < 4; ++i)
//  for (int j = 0; j < 4; ++j)
//    ref_val[i][j] = double(1.0, true);
#endif

  for (int k=0; k<nz; k++){
#if defined(_ORIGINS)
    old_double v;
#else
    double v;
#endif
    // Reading augmented matrix coefficients
    fscanf(in, "%d %d %lg\n", &i, &j, &v);
    /* adjust from 1-based to 0-based */
    i-=1; j-=1;
    
    // printf("%d %d %d %e \n",k,i,j,v);
#if defined(_ORIGINS) && defined(_MATRIX)
    if (!cursor.setToPosition(i, j))
      throw 0;
    if (cursor.elementAt()->hasContent())
      throw 0;
    double* val;
    if (cursor.elementAt()->hasContribution())
      val = &a.set(i, j, v);
    else {
      val = &a.set(i, j, cursor.elementAt()->sharedReference();
      val->value = v;
      val->contributions[0].coefficient = v;
      val->contributions[0].over_coefficient = std::fabs(v);
    }
#else
    val = &a.set(i, j, v);
#endif
    a.set(j, i, *val);
  }

}
else {
  for (int k=0; k<nz; k++){
#if defined(_ORIGINS)
    old_double v;
#else
    double v;
#endif
    // Reading augmented matrix coefficients
    fscanf(in, "%d %d %lg\n", &i, &j, &v);
    /* adjust from 1-based to 0-based */
    i-=1; j-=1;
    
    // printf("%d %d %d %e \n",k,i,j,v);
#if defined(_ORIGINS) && defined(_MATRIX)
    double* val = &a.set(i, j, double(v, true));
    val->contributions[0].is_atomic_var = true;
    dictionary->push_back(std::make_pair(i, j));
#else
    double* val = &a.set(i, j, v);
#endif
    a.set(j, i, *val);
  }
}

  if (in !=stdin) fclose(in);

}

int main(int argc, char** argv) {
  std::ofstream out;
  int N;
  clock_t start, end;

  // Read MatrixMarket file passed as an input argument 
  if (argc < 3){
		fprintf(stderr, "Usage: %s init_map.txt out_map.txt matrix-market-input-filename [chromatic-output-data] [chromatic-output-over-data] [chromatic-erroutput-data]\n", argv[0]);
		fprintf(stderr, "       %s -x matrix-market-input-filename [chromatic-output-data] [chromatic-output-over-data] [chromatic-erroutput-data]\n", argv[0]);
		exit(1);
	} 
  
  SparseMatrix a;
  WorkPlan theplan;
  WorkPlan* plan = &theplan;
  if (strcmp(argv[1], "-x") == 0) {
    plan = nullptr;
    --argv;
    ++argc;
  }
  else
  { std::ifstream in(argv[1]);
    plan->read(in);
  }

  IdDictionary dictionary;
  read_matrix(argv[3], a, plan, &dictionary);
  N = a.getLines();
  double x[N];
  memset(x, 0, sizeof(x));

  
  // Set an arbitrary vector results
  for (int seg = 0; seg <= 10; ++seg) {
#if defined(_ORIGINS) && defined(_INPUT_VECTOR)
    double ref_val(1.0, true, false);
#endif
    for (int j = 0; j < N/10; ++j) {
      int i = (N/10)*seg + j;
      if (i >= N)
        break;
#if defined(_ORIGINS) && defined(_INPUT_VECTOR)
      a.set(i, N, ref_val)->value = i;
#else
      a.set(i, N, i);
#endif
    }
  }
  start = clock(); /* Lancement de la mesure */

  // Applying Gauss Elimination
  for (int i = 0; i < N; ++i) {
    printf("Gauss, line :%d\n",i);
    fflush(stdout);
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
       x[i] -= x[reverseColumnIterator.getIndex()]*reverseColumnIterator.getValue();
    x[i] /= (index == i) ? reverseColumnIterator.getValue() : 0.0;
  }

  end = clock();  /* Arret de la mesure */
#ifdef _ORIGINS
  printf("\nTIME: %lf",((old_double)end - start) / CLOCKS_PER_SEC);
#else
  printf("\nTIME: %lf",((double)end - start) / CLOCKS_PER_SEC);
#endif
  // Displaying solution
  printf("\nRequired solution for %d is: \n", N);
  old_double total_over_contribution = 0.0;
  old_double total_over_approximation= 0.0;
  old_double total_value = 0.0;
  for (int i = 0; i < N; ++i) {
#ifdef _ORIGINS
    old_double over_contribution = 0.0;
    printf("X%d value: %e\n", i, x[i].value);
    printf("Idx : [-2, %e]\n", x[i].contribution_with_unknown_origin);
    printf("      [-1, %e, %e]\n", x[i].contribution_without_origin, x[i].over_contribution_without_origin);
    for (int origin_index = 0; origin_index < x[i].contributions_size; ++origin_index) {
      printf("[%ld, %e, %e]\n", x[i].contributions[origin_index].symbol_id, x[i].contributions[origin_index].coefficient, x[i].contributions[origin_index].over_coefficient);
      over_contribution += x[i].contributions[origin_index].over_coefficient
          - std::fabs(x[i].contributions[origin_index].coefficient);
    }
    over_contribution += x[i].contribution_with_unknown_origin;
    printf("over contribution: %e, %e\n", over_contribution, over_contribution/std::fabs(x[i].value)*100);
    total_over_approximation += over_contribution;
    total_over_contribution += over_contribution
       + x[i].over_contribution_without_origin - std::fabs(x[i].contribution_without_origin);
    total_value += std::fabs(x[i].value);
#else
    printf("X%d value: %e\n", i, x[i]);
#endif
    printf("\n");
    fflush(stdout);
  }

#ifdef _ORIGINS_ERROR
  old_double error = 0.0;
  for (int i = 0; i < N; ++i) {
    std::cout << "X" << i << " value: absolute = ";
    x[i].double_st::print(std::cout);
    error += std::fabs(x[i].accumulated_error);
    std::cout << '\n';
    std::cout << "            relative = ";
    x[i].printRelative(std::cout);
    std::cout << std::endl;
  }
  std::cout << "The total error on the vector is " << error << std::endl;
  std::cout << "The approximation on contribution of values is " << total_over_approximation << ", relatively on values = " << (total_over_approximation/total_value)*100 << '%' << std::endl;
  std::cout << "The imprecision on contribution of values is " << total_over_contribution << ", relatively on values = "  << (total_over_contribution/total_value)*100 << '%' << std::endl;
#endif

#ifdef _ORIGINS
#ifdef _MATRIX

// #define M1 4
// #define M2 4
struct MapResult {
   old_double contribution_value = 0.0;
   old_double contribution_over_value = 0.0;
#ifdef _ORIGINS_ERROR
   old_double contribution_error = 0.0;
   old_double contribution_over_error = 0.0;
#endif
};
std::map<uint64_t, MapResult> heatmap;
if (plan) {
  for (int i = 0; i < N; ++i) {
     WorkPlan::Cursor cursor(*plan);
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        if (!setToSymbolId(k, cursor))
            return 3;
        cursor.elementAt()->getSContribution() += x[i].contributions[origin].coefficient;
        cursor.elementAt()->getSOverContribution() += x[i].contributions[origin].over_coefficient;
#ifdef _ORIGINS_ERROR
        cursor.elementAt()->getSErrorContribution() += std::fabs(x[i].contributions[origin].coefficient*x[i].accumulated_error);
        cursor.elementAt()->getSErrorOverContribution() += x[i].contributions[origin].over_coefficient*std::fabs(x[i].accumulated_error);
#endif
     }
  }
  {  WorkPlan::Cursor cursor(*plan);
     if (cursor.setToFirst()) {
        QuadTreeElement* current;
        do {
           current = cursor.back();
           while (current->hasContent()) {
              if (!cursor.setToNext()) return 1;
              current = cursor.back();
           }
           if (!current->hasContribution())
              current->setContribution();
        } while (cursor.setToNext());
     }
  }
  {  std::ofstream out(argv[2]);
     plan->write(out);
  }
}
else {
  for (int i = 0; i < N; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        // printf("%lu, %d, %d, %d\n", k, i, (int) ((k-1)/M2), (int)((k-1)%M2));
        // fflush(stdout);
        heatmap[k-1].contribution_value += x[i].contributions[origin].coefficient;
        heatmap[k-1].contribution_over_value += x[i].contributions[origin].over_coefficient;
        heatmap[k-1].contribution_error += x[i].contributions[origin].coefficient*std::fabs(x[i].accumulated_error);
        heatmap[k-1].contribution_over_error += x[i].contributions[origin].over_coefficient*std::fabs(x[i].accumulated_error);
     }
  }
}
  if (argc > 4) {
     std::ofstream out(argv[4]);
     if ( out.is_open() )
        std::cout<<"File " << argv[4] <<" opened"<<std::endl;
     if (plan) {
     out << plan->getWidth() << ' ' << plan->getLength() << ' ' << plan->getNonZeroCount() << '\n';
     WorkPlan::Cursor cursor(*plan);
     if (cursor.setToFirst()) {
       QuadTreeElement* current;
       do {
         current = cursor.back();
         while (current->hasContent()) {
            if (!cursor.setToNext()) return 1;
            current = cursor.back();
         }
         if (current->getContribution() != 0.0 || current->getOverContribution() != 0.0)
            out << current->getX() << ' ' << current->getY() << ' '
                << current->getWidth() << ' ' << current->getLength() << ' '
                << current->getContribution() << '\n';
       } while (cursor.setToNext());
     }
     }
     else {
     out << a.size() << ' ' << a[0].size() << ' ' << heatmap.size() << '\n';
     for (const auto& element : heatmap) {
         if (element.second.contribution_value != 0.0 || element.second.contribution_over_value != 0.0)
            out << dictionary[element.first].first << ' ' << dictionary[element.first].second
                << " 1 1 " << element.second.contribution_value << '\n';
     }
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[4] <<" closed"<<std::endl;
        out.close();
     }
  }

  if (argc > 5) {
     std::ofstream out(argv[5]);
     if ( out.is_open() )
        std::cout<<"File " << argv[5] <<" opened"<<std::endl;
     if (plan) {
     out << plan->getWidth() << ' ' << plan->getLength() << ' ' << plan->getCount() << '\n';
     WorkPlan::Cursor cursor(*plan);
     if (cursor.setToFirst()) {
       QuadTreeElement* current;
       do {
         current = cursor.back();
         while (current->hasContent()) {
            if (!cursor.setToNext()) return 1;
            current = cursor.back();
         }
         if (current->getContribution() != 0.0 || current->getOverContribution() != 0.0)
            out << current->getX() << ' ' << current->getY() << ' '
                << current->getWidth() << ' ' << current->getLength() << ' '
                << current->getOverContribution() << '\n';
       } while (cursor.setToNext());
     }
     }
     else {
     out << a.size() << ' ' << a[0].size() << ' ' << heatmap.size() << '\n';
     for (const auto& element : heatmap) {
         if (element.second.contribution_value != 0.0 || element.second.contribution_over_value != 0.0)
            out << dictionary[element.first].first << ' ' << dictionary[element.first].second
                << " 1 1 " << element.second.contribution_over_value << '\n';
     }
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[5] <<" closed"<<std::endl;
        out.close();
     }
  }

#ifdef _ORIGINS_ERROR
  if (argc > 6) {
     std::ofstream out(argv[6]);
     if ( out.is_open() )
        std::cout<<"File " << argv[6] <<" opened"<<std::endl;
     if (plan) {
     out << plan->getWidth() << ' ' << plan->getLength() << ' ' << plan->getCount() << '\n';
     WorkPlan::Cursor cursor(*plan);
     if (cursor.setToFirst()) {
       QuadTreeElement* current;
       do {
         current = cursor.back();
         while (current->hasContent()) {
            if (!cursor.setToNext()) return 1;
            current = cursor.back();
         }
         if (current->getContribution() != 0.0 || current->getOverContribution() != 0.0)
            out << current->getX() << ' ' << current->getY() << ' '
                << current->getWidth() << ' ' << current->getLength() << ' '
                << current->getErrorContribution() << '\n';
       } while (cursor.setToNext());
     }
     }
     else {
     out << a.size() << ' ' << a[0].size() << ' ' << heatmap.size() << '\n';
     for (const auto& element : heatmap) {
         if (element.second.contribution_error != 0.0)
            out << dictionary[element.first].first << ' ' << dictionary[element.first].second
                << " 1 1 " << element.second.contribution_error << '\n';
     }
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[6] <<" closed"<<std::endl;
        out.close();
     }
  }
#endif
     
#else
#define M1 N
#define M2 (N+1)

  // Plot the Norm1 of each input values on the resulting system (ie. Detailed Condition Number)
  old_double* heatmap1 = new old_double[M1*M2]{};
  for (int i = 0; i < N; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        // printf("%lu, %d, %d, %d\n", k, i, (int) ((k-1)/M2), (int)((k-1)%M2));
        // fflush(stdout);
        old_float v = x[i].contributions[origin].coefficient;
        heatmap1[(k-1)/M2*M2 + (k-1)%M2] += v;
     }
  }

  old_double* heatmap2 = new old_double[M1*M2]{};
  for (int i = 0; i < N; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        // printf("%lu, %d, %d, %d\n", k, i, (int) ((k-1)/M2), (int)((k-1)%M2));
        // fflush(stdout);
        old_float v = x[i].contributions[origin].over_coefficient;
        heatmap2[(k-1)/M2*M2 + (k-1)%M2] += v;
     }
  }

#ifdef _ORIGINS_ERROR
  old_double* heatmap3 = new old_double[M1*M2]{};
  for (int i = 0; i < N; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        // printf("%lu, %d, %d, %d\n", k, i, (int) ((k-1)/M2), (int)((k-1)%M2));
        // fflush(stdout);
        old_float v = std::fabs(x[i].contributions[origin].coefficient*x[i].error);
        heatmap3[(k-1)/M2*M2 + (k-1)%M2] += v;
     }
  }
#endif

#ifdef _ORIGINS_DOMAIN
  old_double* heatmap4 = new old_double[M1*M2]{};
  for (int i = 0; i < N; ++i) {
     for (int origin=0; origin < x[i].contributions_size; ++origin) {
        uint64_t k = x[i].contributions[origin].symbol_id;
        // printf("%lu, %d, %d, %d\n", k, i, (int) ((k-1)/M2), (int)((k-1)%M2));
        // fflush(stdout);
        old_float v = x[i].contributions[origin].coefficient
            * (x[i].max_domain - x[i].min_domain)/x[i].value;
        heatmap4[(k-1)/M2*M2 + (k-1)%M2] += v;
     }
  }
#endif

  if (argc > 4) {
     std::ofstream out(argv[4]);
     if ( out.is_open() )
        std::cout<<"File " << argv[4] <<" opened"<<std::endl;
     for (int i = 0; i < M1; ++i) {
        for (int j = 0; j < M2-1; ++j)
            out << heatmap1[i*M2+j] << ", ";
        out << heatmap1[i*M2+M2-1] << '\n';
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[4] <<" closed"<<std::endl;
        out.close();
     }
  }
  delete [] heatmap1;

  if (argc > 5) {
     std::ofstream out(argv[5]);
     if ( out.is_open() )
        std::cout<<"File " << argv[5] <<" opened"<<std::endl;
     for (int i = 0; i < M1; ++i) {
        for (int j = 0; j < M2-1; ++j)
            out << heatmap2[i*M2+j] << ", ";
        out << heatmap2[i*M2+M2-1] << '\n';
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[5] <<" closed"<<std::endl;
        out.close();
     }
  }
  delete [] heatmap2;

#ifdef _ORIGINS_ERROR
  if (argc > 6) {
     std::ofstream out(argv[6]);
     if ( out.is_open() )
        std::cout<<"File " << argv[6] <<" opened"<<std::endl;
     for (int i = 0; i < M1; ++i) {
        for (int j = 0; j < M2-1; ++j)
            out << heatmap3[i*M2+j] << ", ";
        out << heatmap3[i*M2+M2-1] << '\n';
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[6] <<" closed"<<std::endl;
        out.close();
     }
  }
  delete [] heatmap3;
#endif
     
#ifdef _ORIGINS_DOMAIN
  if (argc > 7) {
     std::ofstream out(argv[7]);
     if ( out.is_open() )
        std::cout<<"File " << argv[7] <<" opened"<<std::endl;
     for (int i = 0; i < M1; ++i) {
        for (int j = 0; j < M2-1; ++j)
            out << heatmap4[i*M2+j] << ", ";
        out << heatmap4[i*M2+M2-1] << '\n';
     }
     if ( out.is_open() ){
        std::cout<<"File " << argv[7] <<" closed"<<std::endl;
        out.close();
     }
  }
  delete [] heatmap4;
#endif
#endif
#endif
}

