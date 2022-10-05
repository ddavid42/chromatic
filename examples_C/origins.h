#include <atomic>
#include <algorithm>
#include <memory>
#include <cassert>
#include <cmath>
#include <limits>
// #include <iosfwd>
#include <iostream>

#pragma once

// sorted array or hash table ?

struct BaseOriginVector {
   static std::atomic<uint64_t> atomic_id_counts;
};

struct Origin {
   uint64_t symbol_id = 0; // 0 if a mixture of symbols
   float coefficient = 0;  // absolute contribution = relative contribution*value

   Origin& operator=(const Origin&) = default;

   bool operator<(const Origin& source) const
      {  return (symbol_id < source.symbol_id); }
   bool operator<=(const Origin& source) const
      {  return (symbol_id <= source.symbol_id); }
   bool operator==(const Origin& source) const
      {  return symbol_id == source.symbol_id; }
   bool operator!=(const Origin& source) const
      {  return symbol_id != source.symbol_id; }
   bool operator>(const Origin& source) const
      {  return symbol_id > source.symbol_id; }
   bool operator>=(const Origin& source) const
      {  return (symbol_id >= source.symbol_id); }
   int compare(const Origin& source) const 
      {  return (symbol_id < source.symbol_id) ? -1 : ((symbol_id > source.symbol_id) ? 1 : 0); }
   void print(std::ostream& out, bool& isFirst) const
      {  if (coefficient != 0) {
            if (!isFirst)
               out << " + ";
            else
               isFirst = false;
            out << coefficient << "*id_" << symbol_id;
         }
      }
};

template <typename Type, int N>
struct OriginVector : public BaseOriginVector {
   Type value = 0.0;
   std::array<Origin, N> origins = std::array<Origin, N>{};
   float coeff_without_origin = 0.0;
   // Origin origins[N];
   int origins_size = 0;

   static const int MaxSize = N;

   OriginVector() = default;
   OriginVector(const Type& avalue, bool doesFollow=false);
   OriginVector(const OriginVector<Type, N>& source) = default;
   template <typename OtherType>
   OriginVector(const OriginVector<OtherType, N>& source);

   void print(std::ostream& out) const
      {  bool isFirst = false;
         out << value << ", ";
         std::for_each(origins.begin(), origins.begin() + origins_size,
               [&out, &isFirst](const Origin& origin) { origin.print(out, isFirst); });
         if (coeff_without_origin) {
            if (!isFirst)
               out << " + ";
            out << coeff_without_origin << "*id_unknown";
         }
      }
   OriginVector& operator=(const OriginVector<Type, N>& source) = default;
   template <typename OtherType>
   void assign(const OriginVector<OtherType, N>& source);

   template <int P>
   static void insertSorted(std::array<Origin, P>& destinations, int& destination_size,
         const Origin& origin)
      {  auto insertionIter = std::lower_bound(destinations.begin(),
               destinations.begin() + destination_size, origin);
         if (insertionIter < destinations.begin() + destination_size)
            std::move(insertionIter, destinations.begin()+destination_size, insertionIter+1);
         destination_size++;
         *insertionIter = origin;
      }
   template <int P>
   static void copy(std::array<Origin, N>& destinations, int& destination_size,
         const std::array<Origin, P>& origins, int origins_size,
         float& coeffWithoutOrigin);
   Type add(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         bool is_plus, float& newCoeffWithoutOrigin) const;
   void addAssign(const OriginVector<Type, N>& source, bool is_plus);

   void multAssign(const Type& source);
   void divAssign(const Type& source);
   Type mult(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin) const;
   void multAssign(const OriginVector<Type, N>& source);
   Type div(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin) const;
   void divAssign(const OriginVector<Type, N>& source);
   void fabsAssign();
   void inverseAssign();
   void logAssign();
   void expAssign();
   void sqrtAssign();
   void sinAssign();
   void cosAssign();
   void tanAssign();
   void asinAssign();
   void acosAssign();
   void atanAssign();
   Type atan2(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin) const;
   void atan2Assign(const OriginVector<Type, N>& source);
   Type pow(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin) const;
   void powAssign(const OriginVector<Type, N>& source);
   void oppositeAssign() { value = -value; }
   void maxAssign(const OriginVector<Type, N>& source)
      {  if (value < source.value)
            operator=(source);
      }
   void minAssign(const OriginVector<Type, N>& source)
      {  if (value > source.value)
            operator=(source);
      }
   void absAssign()
      {  if (value < 0)
            oppositeAssign();
      }
   void floorAssign()
      {  std::for_each(origins.begin(), origins.begin() + origins_size,
            [](Origin& origin) { origin = Origin {}; });
         coeff_without_origin = 0.0;
         origins_size = 0;
         value = std::floor(value);
      }
   void fmodAssign(const OriginVector<Type, N>& source)
      {  value = std::fmod(value, source.value); }
   void modfAssign(OriginVector<Type, N>& source)
      {  value = std::modf(value, &source.value);
         source.coeff_without_origin = 0.0;
         source.origins_size = 0;
      }
};

template <typename Type, int N>
inline
OriginVector<Type, N>::OriginVector(const Type& avalue, bool doesFollow)
   :  value(avalue)
{  if (doesFollow) {
      origins[0] = Origin { ++atomic_id_counts, (float) 1.0 };
      origins_size = 1;
   }
   else
      coeff_without_origin = 1;
}

template <typename Type, int N>
template <typename OtherType>
inline
OriginVector<Type, N>::OriginVector(const OriginVector<OtherType, N>& source)
   :  origins(source.origins), coeff_without_origin(source.coeff_without_origin),
      origins_size(source.origins_size), value(source.value) {}

template <typename Type, int N>
template <typename OtherType>
inline void
OriginVector<Type, N>::assign(const OriginVector<OtherType, N>& source) {
   origins = source.origins;
   coeff_without_origin = source.coeff_without_origin;
   origins_size = source.origins_size;
   value = source.value;
}

template <typename Type, int N>
template <int P>
void
OriginVector<Type, N>::copy(std::array<Origin, N>& destinations, int& destinations_size,
         const std::array<Origin, P>& origins, int origins_size, float& coeffWithoutOrigin) {
   float amplitude = 0;
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [&amplitude](const Origin& origin)
         {  amplitude += std::fabs(origin.coefficient); });
   float error_limit = std::fabs(amplitude * 1e-6);
   std::array<const Origin*, P-N> smallestCoefficients;
   int smallestCoefficientsIndex = 0;
   int smallestCoefficientsSize = origins_size - N;
   auto compareCoefficient = [](const Origin* fst, const Origin* snd)
      {  return std::fabs(fst->coefficient) < std::fabs(snd->coefficient); };
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [error_limit, &smallestCoefficients, &smallestCoefficientsIndex,
            &smallestCoefficientsSize, compareCoefficient,
            &destinations, &destinations_size]
      (const Origin& origin)
         {  if (std::fabs(origin.coefficient) < error_limit) {
               if (smallestCoefficientsIndex > 0) {
                  insertSorted<N>(destinations, destinations_size,
                        *smallestCoefficients[smallestCoefficientsIndex-1]);
                  --smallestCoefficientsIndex;
               }
               --smallestCoefficientsSize;
               return;
            }
            if (smallestCoefficientsSize > 0) {
               if (smallestCoefficientsIndex < smallestCoefficientsSize) {
                  smallestCoefficients[smallestCoefficientsIndex++] = &origin;
                  if (smallestCoefficientsIndex == smallestCoefficientsSize)
                     std::sort(smallestCoefficients.begin(), smallestCoefficients.begin()+smallestCoefficientsSize,
                        compareCoefficient);
               }
               else {
                  auto iter = std::lower_bound(smallestCoefficients.begin(),
                        smallestCoefficients.begin()+smallestCoefficientsSize, &origin,
                        compareCoefficient);
                  if (iter < smallestCoefficients.begin()+smallestCoefficientsSize) {
                     insertSorted<N>(destinations, destinations_size,
                           *smallestCoefficients[smallestCoefficientsSize-1]);
                     std::move(iter, smallestCoefficients.begin()+smallestCoefficientsSize-1,
                           iter+1);
                     *iter = &origin;
                  }
                  else
                     destinations[destinations_size++] = origin;
               }
            }
            else
               destinations[destinations_size++] = origin;
         });

   // we could introduce a definition for tne new origin to be stored
   //   in a global linked list of definitions
   // operation_id is then false
   auto compareOrigin = [](const Origin* fst, const Origin* snd)
      {  return *fst < *snd; };
   if (smallestCoefficientsSize > 1)
      std::sort(smallestCoefficients.begin(), smallestCoefficients.begin()+smallestCoefficientsSize,
            compareOrigin);
   Origin newOrigin;
   const Origin* previousOrigin = nullptr;
   for (int smallIndex = 0; smallIndex < smallestCoefficientsSize; ++smallIndex) {
      if (!previousOrigin) {
         previousOrigin = smallestCoefficients[smallIndex];
         newOrigin = *smallestCoefficients[smallIndex];
         newOrigin.symbol_id = 0;
      }
      else
         newOrigin.coefficient += smallestCoefficients[smallIndex]->coefficient;
   }
   if (previousOrigin) {
      auto insertionIter = std::lower_bound(destinations.begin(),
            destinations.begin() + destinations_size, newOrigin);
      if (insertionIter != destinations.end() && !(newOrigin < *insertionIter))
         insertionIter->coefficient += newOrigin.coefficient;
      else
         coeffWithoutOrigin += std::fabs(newOrigin.coefficient);
   }
}

template <typename Type, int N>
Type
OriginVector<Type, N>::add(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      bool is_plus, float& newCoeffWithoutOrigin) const {
   Type newValue = is_plus ? (value + source.value) : (value - source.value);
   Type totalValue = std::fabs(value) + std::fabs(source.value);
   newCoeffWithoutOrigin = coeff_without_origin*std::fabs(value)/totalValue
      + source.coeff_without_origin*std::fabs(source.value)/totalValue;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient *= std::fabs(value)/totalValue;
         new_origins[new_origins_size].coefficient += source.origins[sourceIndex].coefficient*std::fabs(source.value)/totalValue;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient *= std::fabs(value)/totalValue;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.origins[sourceIndex];
         new_origins[new_origins_size].coefficient *= std::fabs(source.value)/totalValue;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size] = origins[thisIndex];
      new_origins[new_origins_size].coefficient *= std::fabs(value)/totalValue;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size] = source.origins[sourceIndex];
      new_origins[new_origins_size].coefficient *= std::fabs(source.value)/totalValue;
      ++new_origins_size;
      ++sourceIndex;
   }
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::addAssign(const OriginVector<Type, N>& source, bool is_plus) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = add(source, new_origins, new_origins_size, is_plus, coeff_without_origin);
   origins_size = 0;
   OriginVector<Type, N>::copy<2*N>(origins, origins_size, new_origins, new_origins_size,
         coeff_without_origin);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::multAssign(const Type& source) {
   value *= source;
   coeff_without_origin = 0.5*coeff_without_origin + 0.5;
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [source](Origin& origin) { origin.coefficient *= 0.5; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::divAssign(const Type& source) {
   value /= source;
   coeff_without_origin = 0.5*coeff_without_origin + 0.5;
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [source](Origin& origin) { origin.coefficient *= 0.5; });
}

template <typename Type, int N>
Type
OriginVector<Type, N>::mult(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin) const {
   Type newValue = value * source.value;
   newCoeffWithoutOrigin = 0.5*coeff_without_origin + 0.5*source.coeff_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*origins[thisIndex].coefficient + 0.5*source.origins[sourceIndex].coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.origins[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size] = origins[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size] = source.origins[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++sourceIndex;
   }
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::multAssign(const OriginVector<Type, N>& source) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = mult(source, new_origins, new_origins_size, coeff_without_origin);
   origins_size = 0;
   OriginVector<Type, N>::copy<2*N>(origins, origins_size, new_origins, new_origins_size,
         coeff_without_origin);
}

template <typename Type, int N>
Type
OriginVector<Type, N>::div(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin) const {
   Type newValue = value / source.value;
   newCoeffWithoutOrigin = 0.5*coeff_without_origin + 0.5*source.coeff_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*origins[thisIndex].coefficient + 0.5*source.origins[sourceIndex].coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.origins[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size] = origins[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size] = source.origins[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++sourceIndex;
   }
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::divAssign(const OriginVector<Type, N>& source) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = div(source, new_origins, new_origins_size, coeff_without_origin);
   origins_size = 0;
   copy<2*N>(origins, origins_size, new_origins, new_origins_size, coeff_without_origin);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::fabsAssign() {
   value = ::fabs(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::inverseAssign() {
   value = 1.0 / value;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::logAssign() {
   value = std::log(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::expAssign() {
   value = std::exp(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::sqrtAssign() {
   value = std::sqrt(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::sinAssign() {
   value = std::sin(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::cosAssign() {
   value = std::cos(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::tanAssign() {
   value = std::tan(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::asinAssign() {
   value = std::asin(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::acosAssign() {
   value = std::acos(value);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::atanAssign() {
   value = std::atan(value);
}

template <typename Type, int N>
Type
OriginVector<Type, N>::atan2(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin) const {
   Type newValue = std::atan2(value, source.value);
   newCoeffWithoutOrigin = 0.5*coeff_without_origin + 0.5*source.coeff_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*origins[thisIndex].coefficient + 0.5*source.origins[sourceIndex].coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.origins[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size] = origins[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size] = source.origins[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++sourceIndex;
   }
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::atan2Assign(const OriginVector<Type, N>& source) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = atan2(source, new_origins, new_origins_size, coeff_without_origin);
   origins_size = 0;
   copy<2*N>(origins, origins_size, new_origins, new_origins_size, coeff_without_origin);
}

template <typename Type, int N>
Type
OriginVector<Type, N>::pow(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin) const {
   Type newValue = std::pow(value, source.value);
   newCoeffWithoutOrigin = 0.5*coeff_without_origin + 0.5*source.coeff_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*origins[thisIndex].coefficient + 0.5*source.origins[sourceIndex].coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = origins[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.origins[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size] = origins[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size] = source.origins[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      ++new_origins_size;
      ++sourceIndex;
   }
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::powAssign(const OriginVector<Type, N>& source) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = pow(source, new_origins, new_origins_size, coeff_without_origin);
   origins_size = 0;
   copy<2*N>(origins, origins_size, new_origins, new_origins_size, coeff_without_origin);
}

struct double_st;
struct float_st : public OriginVector<float, 20> {
  private:
   typedef OriginVector<float, 20> inherited;

  public:
   float_st() = default;
   float_st(const float& avalue, bool)
      :  inherited(avalue, true) {}
   float_st(const float_st& source) = default;
   float_st(const double_st& source);
};

struct double_st : public OriginVector<double, 20> {
  private:
   typedef OriginVector<double, 20> inherited;

  public:
   double_st() = default;
   double_st(const double& avalue, bool)
      :  inherited(avalue, true) {}
   double_st(const double_st& source) = default;
   double_st(const float_st& source);
};

inline
float_st::float_st(const double_st& source) : inherited(source) {}

inline
double_st::double_st(const float_st& source) : inherited(source) {}

