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

   Origin& operator+=(const Origin& source)
      {  if (symbol_id && symbol_id == source.symbol_id)
            coefficient += source.coefficient;
         else {
            symbol_id = 0;
            coefficient = std::fabs(coefficient) + std::fabs(source.coefficient);
         }
         return *this;
      }
   Origin operator+(const Origin& source) const
      {  if (symbol_id && symbol_id == source.symbol_id)
            return Origin { symbol_id, coefficient + source.coefficient };
         return Origin { 0, std::fabs(coefficient) + std::fabs(source.coefficient) };
      }
   Origin& operator-=(const Origin& source)
      {  if (symbol_id && symbol_id == source.symbol_id)
            coefficient -= source.coefficient;
         else {
            symbol_id = 0;
            coefficient = std::fabs(coefficient) + std::fabs(source.coefficient);
         }
         return *this;
      }
   Origin operator-(const Origin& source) const
      {  if (symbol_id && symbol_id == source.symbol_id)
            return Origin { symbol_id, coefficient - source.coefficient };
         return Origin { 0, std::fabs(coefficient) + std::fabs(source.coefficient) };
      }
   template <typename Type>
   Origin& operator*=(const Type& factor)
      {  coefficient *= factor;
         return *this;
      }
   template <typename Type>
   Origin& operator/=(const Type& factor)
      {  coefficient /= factor;
         return *this;
      }
   template <typename Type>
   Origin operator*(const Type& factor) const
      {  return Origin { symbol_id, (float) (coefficient * factor) }; }
   template <typename Type>
   Origin operator/(const Type& factor) const
      {  return Origin { symbol_id, (float) (coefficient / factor) }; }
   Origin& addAssign(const Origin& source, bool is_plus=true)
      {  assert(symbol_id == source.symbol_id);
         if (!symbol_id)
            coefficient = std::fabs(coefficient) + std::fabs(source.coefficient);
         else if (is_plus)
            coefficient += source.coefficient;
         else
            coefficient -= source.coefficient;
         return *this;
      }
   Origin add(const Origin& source, bool is_plus=true) const
      {  assert(symbol_id == source.symbol_id);
         if (!symbol_id)
            return Origin { 0, std::fabs(coefficient) + std::fabs(source.coefficient) };
         if (is_plus)
            return Origin { symbol_id, coefficient + source.coefficient };
         else
            return Origin { symbol_id, coefficient - source.coefficient };
      }
   Origin& oppositeAssign() { if (symbol_id) coefficient = -coefficient; return *this; }
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
   OriginVector(const Type& avalue, bool);
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
   void oppositeAssign()
      {  std::for_each(origins.begin(), origins.begin() + origins_size,
            [](Origin& origin)
               {  origin.coefficient = -origin.coefficient; });
         coeff_without_origin = -coeff_without_origin;
         value = -value;
      }
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
OriginVector<Type, N>::OriginVector(const Type& avalue, bool)
   :  value(avalue) {
   int exp;
   std::frexp(std::fabs(avalue * 256.0), &exp);
   Type intValue = avalue / std::ldexp((Type) 1.0, exp-9);
   int ceilIntValue = intValue;
   if ((Type) ceilIntValue != intValue) {
      Type ulp = std::ldexp(std::numeric_limits<Type>::epsilon(), exp - 10);
      origins[0] = Origin { ++atomic_id_counts, (float) (ulp/4.0) /* average error */ };
      origins_size = 1;
   }
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
         newOrigin += *smallestCoefficients[smallIndex];
   }
   if (previousOrigin) {
      auto insertionIter = std::lower_bound(destinations.begin(),
            destinations.begin() + destinations_size, newOrigin);
      if (insertionIter != destinations.end() && !(newOrigin < *insertionIter))
         *insertionIter += newOrigin;
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
   newCoeffWithoutOrigin = coeff_without_origin + source.coeff_without_origin; // fabs is not useful here since both should be >= 0
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size++] = origins[thisIndex].add(
               source.origins[sourceIndex], is_plus);
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size++] = origins[thisIndex];
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size++] = source.origins[sourceIndex];
         if (!is_plus)
            new_origins[new_origins_size-1].oppositeAssign();
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size++] = origins[thisIndex];
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size++] = source.origins[sourceIndex];
      if (!is_plus)
         new_origins[new_origins_size-1].oppositeAssign();
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
   coeff_without_origin *= source;
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [source](Origin& origin) { origin *= source; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::divAssign(const Type& source) {
   value /= source;
   coeff_without_origin /= source;
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [source](Origin& origin) { origin /= source; });
}

template <typename Type, int N>
Type
OriginVector<Type, N>::mult(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin) const {
   Type newValue = value * source.value;
   newCoeffWithoutOrigin = std::fabs(coeff_without_origin * source.value)
      + std::fabs(source.coeff_without_origin * value);
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size++] = (origins[thisIndex] * source.value)
            .add(source.origins[sourceIndex] * value);
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size++] = origins[thisIndex] * source.value;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size++] = source.origins[sourceIndex] * value;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size++] = origins[thisIndex] * source.value;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size++] = source.origins[sourceIndex] * value;
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
   float sourceSquare = source.value * source.value;
   newCoeffWithoutOrigin = (std::fabs(coeff_without_origin * source.value)
      + std::fabs(source.coeff_without_origin * value)) / sourceSquare;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size++] = (origins[thisIndex] * source.value)
            .add(source.origins[sourceIndex] * value, false /* isPlus */);
         new_origins[new_origins_size-1] /= sourceSquare;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size++] = origins[thisIndex] * source.value;
         new_origins[new_origins_size-1] /= sourceSquare;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size++] = source.origins[sourceIndex] * value;
         new_origins[new_origins_size-1].oppositeAssign();
         new_origins[new_origins_size-1] /= sourceSquare;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size++] = origins[thisIndex] * source.value;
      new_origins[new_origins_size-1] /= sourceSquare;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size++] = source.origins[sourceIndex] * value;
      new_origins[new_origins_size-1].oppositeAssign();
      new_origins[new_origins_size-1] /= sourceSquare;
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
   bool changeSign = value < 0;
   value = ::fabs(value);
   coeff_without_origin = ::fabs(coeff_without_origin);
   if (changeSign)
      std::for_each(origins.begin(), origins.begin() + origins_size,
         [](Origin& origin) { origin.oppositeAssign(); });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::inverseAssign() {
   float valueSquare = value*value;
   value = 1.0 / value;
   coeff_without_origin = -coeff_without_origin / valueSquare;
   std::for_each(origins.begin(), origins.begin() + origins_size,
      [valueSquare](Origin& origin) { origin /= valueSquare; origin.oppositeAssign(); });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::logAssign() {
   Type oldValue = value;
   value = std::log(value);
   coeff_without_origin /= oldValue;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [oldValue](Origin& origin) { origin /= oldValue; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::expAssign() {
   value = std::exp(value);
   coeff_without_origin /= value;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [this](Origin& origin) { origin *= value; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::sqrtAssign() {
   value = std::sqrt(value);
   Type factor = 1.0/(2.0*value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::sinAssign() {
   Type factor = std::cos(value);
   value = std::sin(value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::cosAssign() {
   Type factor = -std::sin(value);
   value = std::cos(value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::tanAssign() {
   Type factor = 1.0/std::cos(value);
   value = std::tan(value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::asinAssign() {
   Type factor = 1.0/(std::sqrt(1.0 - value*value));
   value = std::asin(value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::acosAssign() {
   Type factor = -1.0/(std::sqrt(1.0 - value*value));
   value = std::acos(value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::atanAssign() {
   Type factor = 1.0/(1.0 + value*value);
   value = std::atan(value);
   coeff_without_origin *= factor;

   std::for_each(origins.begin(), origins.begin() + origins_size,
      [factor](Origin& origin) { origin *= factor; });
}

template <typename Type, int N>
Type
OriginVector<Type, N>::atan2(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin) const {
   Type newValue = std::atan2(value, source.value);
   float thisFactor = -source.value/(value*value + source.value*source.value);
   float sourceFactor = value/(value*value + source.value*source.value);
   newCoeffWithoutOrigin = (std::fabs(coeff_without_origin * thisFactor)
      + std::fabs(source.coeff_without_origin * sourceFactor));
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size++] = (origins[thisIndex] * thisFactor)
            .add(source.origins[sourceIndex] * sourceFactor, false /* isPlus */);
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size++] = origins[thisIndex] * thisFactor;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size++] = source.origins[sourceIndex] * sourceFactor;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size++] = origins[thisIndex] * thisFactor;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size++] = source.origins[sourceIndex] * sourceFactor;
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
   float thisFactor = source.value*std::pow(value, source.value-1.0);
   float sourceFactor = std::log(value)*newValue;
   newCoeffWithoutOrigin = (std::fabs(coeff_without_origin * thisFactor)
      + std::fabs(source.coeff_without_origin * sourceFactor));
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < origins_size && sourceIndex < source.origins_size) {
      auto compare = origins[thisIndex].compare(source.origins[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size++] = (origins[thisIndex] * thisFactor)
            .add(source.origins[sourceIndex] * sourceFactor, false /* isPlus */);
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size++] = origins[thisIndex] * thisFactor;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size++] = source.origins[sourceIndex] * sourceFactor;
         ++sourceIndex;
      }
   }
   while (thisIndex < origins_size) {
      new_origins[new_origins_size++] = origins[thisIndex] * thisFactor;
      ++thisIndex;
   }
   while (sourceIndex < source.origins_size) {
      new_origins[new_origins_size++] = source.origins[sourceIndex] * sourceFactor;
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
struct float_st : public OriginVector<float, 12> {
  private:
   typedef OriginVector<float, 12> inherited;

  public:
   float_st() = default;
   float_st(const float& avalue, bool)
      :  inherited(avalue, true) {}
   float_st(const float_st& source) = default;
   float_st(const double_st& source);
};

struct double_st : public OriginVector<double, 12> {
  private:
   typedef OriginVector<double, 12> inherited;

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

