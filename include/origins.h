#include <atomic>
#include <algorithm>
#include <memory>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctype.h>
#include "eft.h"


#pragma once

// sorted array or hash table ?

struct BaseOriginVector {
   static std::atomic<uint64_t> atomic_id_counts;
};

struct Origin {
   uint64_t symbol_id = 0; // 0 if a mixture of symbols
   float coefficient = 0;  // absolute contribution = relative contribution*value
   float over_coefficient = 0;  // absolute contribution = relative contribution*value
   bool is_atomic_var = false; // true iff a mixture of variables

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
            out << '[' << coefficient << ", " << over_coefficient << ']' << "*id_" << symbol_id;
         }
      }
   template <typename Type>
   void printRelative(std::ostream& out, bool& isFirst, const Type& value) const
      {  if (coefficient != 0) {
            if (!isFirst)
               out << " + ";
            else
               isFirst = false;
            out << '[' << coefficient/value << ", " << over_coefficient/std::fabs(value) << ']' << "*id_" << symbol_id;
         }
      }
};

template <typename Type, int N>
struct OriginVector : public BaseOriginVector {
   Type value = 0.0;
   std::array<Origin, N> contributions = std::array<Origin, N>{};
   float contribution_without_origin = 0.0;
   float over_contribution_without_origin = 0.0;
   int contributions_size = 0;
#ifdef _ORIGINS_DOMAIN
   std::array<Origin, N> domain = std::array<Origin, N>{};
   int domain_size = 0;
   float domain_additional_interval = 0.0;
   Type domain_center = 0.0;
   float min_domain = 0.0, max_domain = 0.0;
#endif
#ifdef _ORIGINS_ERROR
   double accumulated_error = 0.0;
#endif

   static const int MaxSize = N;

   OriginVector() = default;
   OriginVector(const Type& avalue, bool doesFollow=false, bool isAtomic=true);
   OriginVector(const OriginVector<Type, N>& source) = default;
   template <typename OtherType>
   OriginVector(const OriginVector<OtherType, N>& source);

   void clear_origin() { contribution_without_origin = 0.0; over_contribution_without_origin = 0.0; contributions_size = 0; }
   void print(std::ostream& out) const
      {  bool isFirst = false;
         out << value << ", ";
         std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [&out, &isFirst](const Origin& origin) { origin.print(out, isFirst); });
         if (contribution_without_origin || over_contribution_without_origin) {
            if (!isFirst)
               out << " + ";
            out << '[' << contribution_without_origin << ", " << over_contribution_without_origin << ']' << "*id_unknown";
         }
#ifdef _ORIGINS_DOMAIN
         out << ", [" << value + min_domain << ", " << value + max_domain << ']';
#endif
#ifdef _ORIGINS_ERROR
         out << ", error=" << accumulated_error;
#endif
      }
   void printRelative(std::ostream& out) const
      {  bool isFirst = false;
         out << value << ", ";
         std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [&out, &isFirst, this](const Origin& origin) { origin.printRelative(out, isFirst, value); });
         if (contribution_without_origin || over_contribution_without_origin) {
            if (!isFirst)
               out << " + ";
            out << '[' << contribution_without_origin/value << ", " << over_contribution_without_origin/std::fabs(value) << ']' << "*id_unknown";
         }
#ifdef _ORIGINS_DOMAIN
         out << ", [" << (value + min_domain)/value << ", " << (value + max_domain)/value << ']';
#endif
#ifdef _ORIGINS_ERROR
         out << ", relative error=" << accumulated_error/value;
#endif
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
         const std::array<Origin, P>& contributions, int contributions_size,
         float& coeffWithoutOrigin, float* overCoeffWithoutOrigin=nullptr);
   void addContribution(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         bool is_plus, float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const;
#ifdef _ORIGINS_DOMAIN
   void addDomain(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         bool is_plus, float& newCoeffWithoutOrigin, Type& newDomainCenter,
         float& newMinDomain, float& newMaxDomain) const;
#endif
   void addAssign(const OriginVector<Type, N>& source, bool is_plus);

   void multAssign(const Type& source);
   void divAssign(const Type& source);
   void multContribution(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const;
#ifdef _ORIGINS_DOMAIN
   void multDomain(const Origin* sourceBegin, int sizeSource, Type sourceValue,
         float sourceAdditionalInterval, Type sourceCenter, float sourceMin, float sourceMax,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, Type& newDomainCenter, float& newMinDomain,
         float& newMaxDomain) const;
   void multDomain(const std::array<Origin, (2*N)>& source, int sizeSource,
         Type sourceValue, float sourceAdditionalInterval, Type sourceCenter,
         float sourceMin, float sourceMax,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, Type& newDomainCenter, float& newMinDomain,
         float& newMaxDomain) const
      {  multDomain(&*source.begin(), sizeSource,
            sourceValue, sourceAdditionalInterval, sourceCenter, sourceMin, sourceMax,
            new_origins, new_origins_size, newCoeffWithoutOrigin, newDomainCenter,
            newMinDomain, newMaxDomain);
      }
   void multDomain(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, Type& newDomainCenter,
         float& newMinDomain, float& newMaxDomain) const
      {  multDomain(&*source.domain.begin(), source.domain_size,
            source.value, source.domain_additional_interval, source.domain_center,
            source.min_domain, source.max_domain,
            new_origins, new_origins_size, newCoeffWithoutOrigin, newDomainCenter,
            newMinDomain, newMaxDomain);
      }
#endif
   void multAssign(const OriginVector<Type, N>& source);
   void divContribution(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const;
   void divDomain(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, Type& newDomainCenter, float& newMinDomain,
         float& newMaxDomain) const;
   void divAssign(const OriginVector<Type, N>& source);
   void fabsAssign() { absAssign(); }
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
         float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const;
   void atan2Assign(const OriginVector<Type, N>& source);
   Type pow(const OriginVector<Type, N>& source,
         std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
         float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const;
   void powAssign(const OriginVector<Type, N>& source);
   void oppositeAssign()
      {  value = -value;
         contribution_without_origin = -contribution_without_origin;
         std::for_each(contributions.begin(), contributions.begin() + contributions_size,
            [](Origin& origin) { origin.coefficient = -origin.coefficient; });
#ifdef _ORIGINS_DOMAIN
         std::for_each(domain.begin(), domain.begin() + domain_size,
            [](Origin& origin)
               {  origin.coefficient = -origin.coefficient; });
         domain_center = -domain_center;
         min_domain = -min_domain;
         max_domain = -max_domain;
         std::swap(min_domain, max_domain);
#endif
#ifdef _ORIGINS_ERROR
         accumulated_error = -accumulated_error;
#endif
      }
   void inverseAssign()
      {  contribution_without_origin = contribution_without_origin/(value*value);
         over_contribution_without_origin = over_contribution_without_origin/(value*value);
         std::for_each(contributions.begin(), contributions.begin() + contributions_size,
            [this](Origin& origin) { origin.coefficient /= (value*value); origin.over_coefficient /= (value*value); });
#ifdef _ORIGINS_DOMAIN
         float sum = domain_additional_interval;
         std::for_each(domain.begin(), domain.begin() + domain_size,
            [&sum](const Origin& origin)
               {  sum += std::fabs(origin.coefficient); });
         if (min_domain > 0 || max_domain < 0) {
            float newMinDomain = 1.0/max_domain;
            max_domain = 1.0/min_domain;
            min_domain = newMinDomain;
         }
         else {
            min_domain = -std::numeric_limits<Type>::infinity();
            max_domain = std::numeric_limits<Type>::infinity();
            domain_additional_interval = std::numeric_limits<Type>::infinity();
            value = 1.0/value;
            return;
         }
         Type affineCenter = value + domain_center;
         if (sum >= std::fabs(affineCenter)) {
            domain_additional_interval = std::numeric_limits<Type>::infinity();
            domain_size = 0;
            value = 1.0/value;
            return;
         }
         domain_center = 1.0/affineCenter - 1.0/value;
         std::for_each(domain.begin(), domain.begin() + domain_size,
            [affineCenter](Origin& origin)
               {  origin.coefficient /= -affineCenter*affineCenter; });
         float newAdditionalInterval = sum*sum/(std::fabs(affineCenter) - sum);
         domain_additional_interval = domain_additional_interval/(affineCenter*affineCenter)
            + newAdditionalInterval;
#endif
#ifdef _ORIGINS_ERROR
         Type result = 1.0/value;
         double new_error = -EFT::RemainderDiv((Type) 1.0, value, result);
         accumulated_error = (new_error - result*accumulated_error)/(value + accumulated_error);
         value = result;
#else
         value = 1.0/value;
#endif
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
      {  
#ifdef _ORIGINS_DOMAIN
         if (max_domain < 0) {
            value = -value;
            std::for_each(domain.begin(), domain.begin() + domain_size,
               [](Origin& origin)
                  {  origin.coefficient = -origin.coefficient; });
            domain_center = -domain_center;
            min_domain = -min_domain;
            max_domain = -max_domain;
            std::swap(min_domain, max_domain);
         }
         else if (min_domain < 0) {
            if (value < 0) {
               std::for_each(domain.begin(), domain.begin() + domain_size,
                  [](Origin& origin)
                     {  origin.coefficient = -origin.coefficient; });
               domain_center = -domain_center;
            }
            if (-min_domain > max_domain)
               max_domain = min_domain;
            min_domain = 0.0;
         }
#endif
         if (value < 0) {
            contribution_without_origin = -contribution_without_origin;
            std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [](Origin& origin) { origin.coefficient = -origin.coefficient; });
            value = -value;
#ifdef _ORIGINS_ERROR
            accumulated_error = -accumulated_error;
#endif
         }
      }
   void floorAssign()
      {  
         Type result = std::floor(value);
         if (value != 0) {
            contribution_without_origin *= (result / value);
            over_contribution_without_origin *= std::fabs(result / value);
            std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [this, result](Origin& origin) { origin.coefficient *= (result / value); origin.over_coefficient *= std::fabs(result / value); });
         }
         else
            clear_origin();
#ifdef _ORIGINS_DOMAIN
         min_domain = std::floor(value + min_domain) - value;
         max_domain = std::floor(value + max_domain) - value;
#endif
#ifdef _ORIGINS_ERROR
         accumulated_error = std::floor(value - accumulated_error) - result;
#endif
         value = result;
      }
   void fmodAssign(const OriginVector<Type, N>& source)
      {  
#ifdef _ORIGINS_DOMAIN
         if (value >= 0) {
            if (max_domain + value > source.value + source.max_domain) {
               max_domain = source.max_domain + source.value - value;
               min_domain = -value;
            }
            else if (min_domain < -value)
               min_domain = -value;
         }
#endif
         Type result = std::fmod(value, source.value);
         if (value != 0) {
            contribution_without_origin *= (result / value);
            over_contribution_without_origin *= std::fabs(result / value);
            std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [this, result](Origin& origin) { origin.coefficient *= (result / value); origin.over_coefficient *= std::fabs(result / value); });
         }
         else
            clear_origin();
#ifdef _ORIGINS_ERROR
         accumulated_error = std::fmod(value - accumulated_error, source.value - source.accumulated_error)
               - result;
#endif
         value = result;
      }
   void modfAssign(OriginVector<Type, N>& source)
      {  
#ifdef _ORIGINS_DOMAIN
         source.min_domain = std::floor(value + min_domain) - value;
         source.max_domain = std::floor(value + max_domain) - value;
         source.domain_center = domain_center - 0.5;
         source.domain_additional_interval = domain_additional_interval + 0.5;
         if (value >= 0) {
            if (max_domain + value > 1.0)
               max_domain = 1.0 - value;
            if (min_domain + value < 0.0)
               min_domain = -value;
         }
         else {
            if (max_domain + value > 0.0)
               max_domain = -value;
            if (min_domain + value < -1.0)
               min_domain = -1.0 - value;
         }
#endif
         Type result, sourceResult = 0;
         result = std::modf(value, &sourceResult);
#ifdef _ORIGINS_ERROR
         double sourceWithError = 0;
         accumulated_error = std::modf(value - accumulated_error, &sourceWithError) - result;
         source.accumulated_error = sourceResult - sourceWithError;
#endif
         if (value != 0) {
            contribution_without_origin *= (result / value);
            over_contribution_without_origin *= std::fabs(result / value);
            std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [this, result](Origin& origin) { origin.coefficient *= (result / value); origin.over_coefficient *= std::fabs(result / value); });
         }
         else
            clear_origin();
         value = result;
         source.contributions = contributions;
         source.contribution_without_origin = contribution_without_origin;
         source.over_contribution_without_origin = over_contribution_without_origin;
         source.contributions_size = contributions_size;
      }
   void restrictContributions()
      {  
#ifdef _ORIGINS_DOMAIN
         float amplitude = contribution_without_origin;
         std::for_each(contributions.begin(), contributions.begin() + contributions_size,
            [&amplitude](const Origin& origin)
               {  amplitude += std::fabs(origin.coefficient); });
         if (std::max(std::fabs(min_domain), std::fabs(max_domain)) < amplitude*std::fabs(value)) {
            float newAmplitude = std::max(std::fabs(min_domain), std::fabs(max_domain)) /std::fabs(value);
            float factor = newAmplitude / amplitude;
            std::for_each(contributions.begin(), contributions.begin() + contributions_size,
               [factor](Origin& origin)
                  {  origin.coefficient *= factor; });
            contribution_without_origin *= factor;
         }
#endif
      }
};

template <typename Type, int N>
inline
OriginVector<Type, N>::OriginVector(const Type& avalue, bool doesFollow, bool isAtomic)
   :  value(avalue)
{  if (value == 0.0)
      return;
   if (doesFollow) {
      contributions[0] = Origin { ++atomic_id_counts, (float) value, (float) std::fabs(value), isAtomic };
      contributions_size = 1;
#ifdef _ORIGINS_DOMAIN
      domain[0] = Origin { atomic_id_counts, (float) 0.0001 };
      domain_size = 1;
      max_domain = 1.0001*value;
      min_domain = 0.9999*value;
      if (value < 0)
         std::swap(min_domain, max_domain);
#endif
   }
   else {
      contribution_without_origin = value;
      over_contribution_without_origin = std::fabs(value);
   }
}

template <typename Type, int N>
template <typename OtherType>
inline
OriginVector<Type, N>::OriginVector(const OriginVector<OtherType, N>& source)
   :  contributions(source.contributions), contribution_without_origin(source.contribution_without_origin),
      over_contribution_without_origin(source.over_contribution_without_origin),
      contributions_size(source.contributions_size), value(source.value)
#ifdef _ORIGINS_DOMAIN
      , domain(source.domain), domain_size(source.domain_size),
      domain_additional_interval(source.domain_additional_interval), domain_center(source.domain_center),
      min_domain(source.min_domain), max_domain(source.max_domain)
#endif
   {
#ifdef _ORIGINS_ERROR
      accumulated_error = value - source.value + source.accumulated_error;
#endif
   }

template <typename Type, int N>
template <typename OtherType>
inline void
OriginVector<Type, N>::assign(const OriginVector<OtherType, N>& source) {
   contributions = source.contributions;
   contribution_without_origin = source.contribution_without_origin;
   over_contribution_without_origin = source.over_contribution_without_origin;
   contributions_size = source.contributions_size;
#ifdef _ORIGINS_DOMAIN
   domain = source.domain;
   domain_size = source.domain_size;
   domain_additional_interval = source.domain_additional_interval;
   domain_center = source.domain_center;
   min_domain = source.min_domain;
   max_domain = source.max_domain;
#endif
   value = source.value;
#ifdef _ORIGINS_ERROR
   accumulated_error = value - source.value + source.accumulated_error;
#endif
}

template <typename Type, int N>
template <int P>
void
OriginVector<Type, N>::copy(std::array<Origin, N>& destinations, int& destinations_size,
         const std::array<Origin, P>& contributions, int contributions_size, float& coeffWithoutOrigin,
         float* overCoeffWithoutOrigin) {
   float amplitude = 0;
   std::for_each(contributions.begin(), contributions.begin() + contributions_size,
      [&amplitude](const Origin& origin)
         {  amplitude += std::fabs(origin.coefficient); });
   float error_limit = std::fabs(amplitude * 1e-6);
   std::array<const Origin*, P-N> smallestCoefficients;
   int smallestCoefficientsIndex = 0;
   int smallestCoefficientsSize = contributions_size - N;
   auto compareCoefficient = [](const Origin* fst, const Origin* snd)
      {  return std::fabs(fst->coefficient) < std::fabs(snd->coefficient); };
   std::for_each(contributions.begin(), contributions.begin() + contributions_size,
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
      else {
         coeffWithoutOrigin += newOrigin.coefficient;
         if (overCoeffWithoutOrigin)
            *overCoeffWithoutOrigin += std::fabs(newOrigin.coefficient);
      }
   }
}

template <typename Type, int N>
void
OriginVector<Type, N>::addContribution(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      bool is_plus, float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const {
// Type totalValue = std::fabs(is_plus ? (value + source.value) : (value - source.value));
   newCoeffWithoutOrigin = is_plus ? (contribution_without_origin + source.contribution_without_origin)
         :  (contribution_without_origin - source.contribution_without_origin);
   newOverCoeffWithoutOrigin = over_contribution_without_origin + source.over_contribution_without_origin;
// if (totalValue != 0)
//    newCoeffWithoutOrigin = contribution_without_origin*std::fabs(value)/totalValue
//       + source.contribution_without_origin*std::fabs(source.value)/totalValue;
// else
//    newCoeffWithoutOrigin = 1.0;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < contributions_size && sourceIndex < source.contributions_size) {
      auto compare = contributions[thisIndex].compare(source.contributions[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         if (is_plus)
            new_origins[new_origins_size].coefficient += source.contributions[sourceIndex].coefficient;
         else
            new_origins[new_origins_size].coefficient -= source.contributions[sourceIndex].coefficient;
         new_origins[new_origins_size].over_coefficient += source.contributions[sourceIndex].over_coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.contributions[sourceIndex];
         if (!is_plus)
            new_origins[new_origins_size].coefficient = -new_origins[new_origins_size].coefficient;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < contributions_size) {
      new_origins[new_origins_size] = contributions[thisIndex];
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.contributions_size) {
      new_origins[new_origins_size] = source.contributions[sourceIndex];
      if (!is_plus)
         new_origins[new_origins_size].coefficient = -new_origins[new_origins_size].coefficient;
      ++new_origins_size;
      ++sourceIndex;
   }
}

#ifdef _ORIGINS_DOMAIN
template <typename Type, int N>
void
OriginVector<Type, N>::addDomain(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      bool is_plus, float& newCoeffWithoutOrigin, Type& newDomainCenter,
      float& newMinDomain, float& newMaxDomain) const {
   newCoeffWithoutOrigin = domain_additional_interval + source.domain_additional_interval;
   newDomainCenter = domain_center + source.domain_center;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < domain_size && sourceIndex < source.domain_size) {
      auto compare = domain[thisIndex].compare(source.domain[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = domain[thisIndex];
         if (is_plus)
            new_origins[new_origins_size].coefficient += source.domain[sourceIndex].coefficient;
         else
            new_origins[new_origins_size].coefficient -= source.domain[sourceIndex].coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = domain[thisIndex];
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.domain[sourceIndex];
         if (!is_plus)
            new_origins[new_origins_size].coefficient = -new_origins[new_origins_size].coefficient;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < domain_size) {
      new_origins[new_origins_size] = domain[thisIndex];
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.domain_size) {
      new_origins[new_origins_size] = source.domain[sourceIndex];
      if (!is_plus)
         new_origins[new_origins_size].coefficient = -new_origins[new_origins_size].coefficient;
      ++new_origins_size;
      ++sourceIndex;
   }
   newMinDomain = min_domain;
   newMaxDomain = max_domain;
   if (is_plus) {
      newMinDomain += source.min_domain;
      newMaxDomain += source.max_domain;
   }
   else {
      newMinDomain -= source.min_domain;
      newMaxDomain -= source.max_domain;
   }

   // apply domain reduction
   float sum = newCoeffWithoutOrigin;
   std::for_each(new_origins.begin(), new_origins.begin() + new_origins_size,
      [&sum](const Origin& origin)
         {  sum += std::fabs(origin.coefficient); });
   if (newMinDomain < sum + value + domain_center)
      newMinDomain = sum + value + domain_center;
   if (newMaxDomain > sum + value + domain_center)
      newMaxDomain = sum + value + domain_center;
}
#endif

template <typename Type, int N>
inline void
OriginVector<Type, N>::addAssign(const OriginVector<Type, N>& source, bool is_plus) {
   Type newValue = is_plus ? (value + source.value) : (value - source.value);
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   {  addContribution(source, new_origins, new_origins_size, is_plus, contribution_without_origin, over_contribution_without_origin);
      contributions_size = 0;
      OriginVector<Type, N>::copy<2*N>(contributions, contributions_size, new_origins, new_origins_size,
            contribution_without_origin, &over_contribution_without_origin);
   }
#ifdef _ORIGINS_DOMAIN
   new_origins_size = 0;
   {  addDomain(source, new_origins, new_origins_size, is_plus, domain_additional_interval,
         domain_center, min_domain, max_domain);
      domain_size = 0;
      OriginVector<Type, N>::copy<2*N>(domain, domain_size, new_origins, new_origins_size,
            domain_additional_interval);
   }
#endif
#ifdef _ORIGINS_ERROR
   double newError = -EFT::TwoSum(value, source.value, newValue);
   accumulated_error += source.accumulated_error;
   accumulated_error += newError;
#endif
   value = newValue;
   if (value == 0.0)
      clear_origin();
   restrictContributions();
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::multAssign(const Type& source) {
   contribution_without_origin = source*(0.5*contribution_without_origin + 0.5);
   over_contribution_without_origin = std::fabs(source)*(0.5*over_contribution_without_origin + 0.5);
   std::for_each(contributions.begin(), contributions.begin() + contributions_size,
      [source](Origin& origin) { origin.coefficient *= 0.5*source; origin.over_coefficient *= 0.5*std::fabs(source); });
#ifdef _ORIGINS_DOMAIN
   domain_additional_interval *= std::fabs(source);
   domain_center *= source;
   std::for_each(domain.begin(), domain.begin() + domain_size,
      [source](Origin& origin) { origin.coefficient *= source; });
   min_domain *= source;
   max_domain *= source;
   if (source < 0)
      std::swap(min_domain, max_domain);
#endif
#ifdef _ORIGINS_ERROR
   Type result = value * source;
   double newError = -EFT::FastTwoProd(value, source, result);
   accumulated_error *= source;
   accumulated_error += newError;
   value = result;
#else
   value *= source;
#endif
   if (value == 0.0)
      clear_origin();
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::divAssign(const Type& source) {
   contribution_without_origin = (0.5*contribution_without_origin + 0.5)/source;
   over_contribution_without_origin = (0.5*contribution_without_origin + 0.5)/std::fabs(source);
   std::for_each(contributions.begin(), contributions.begin() + contributions_size,
      [source](Origin& origin) { origin.coefficient *= 0.5/source; origin.over_coefficient *= 0.5/std::fabs(source); });
#ifdef _ORIGINS_DOMAIN
   domain_additional_interval /= std::fabs(source);
   domain_center /= source;
   std::for_each(domain.begin(), domain.begin() + domain_size,
      [source](Origin& origin) { origin.coefficient /= source; });
   min_domain /= source;
   max_domain /= source;
   if (source < 0)
      std::swap(min_domain, max_domain);
#endif
#ifdef _ORIGINS_ERROR
   Type result = value / source;
   double newError = -EFT::RemainderDiv(value, source, result);
   accumulated_error = (newError + accumulated_error) / source;
   value = result;
#else
   value /= source;
#endif
   if (value == 0.0)
      clear_origin();
}

template <typename Type, int N>
void
OriginVector<Type, N>::multContribution(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const {
   newCoeffWithoutOrigin = source.value*0.5*contribution_without_origin
         + value*0.5*source.contribution_without_origin;
   newOverCoeffWithoutOrigin = std::fabs(source.value)*0.5*over_contribution_without_origin
         + std::fabs(value)*0.5*source.over_contribution_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < contributions_size && sourceIndex < source.contributions_size) {
      auto compare = contributions[thisIndex].compare(source.contributions[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient
            = source.value*0.5*contributions[thisIndex].coefficient
               + value*0.5*source.contributions[sourceIndex].coefficient;
         new_origins[new_origins_size].over_coefficient
            = std::fabs(source.value)*0.5*contributions[thisIndex].over_coefficient
               + std::fabs(value)*0.5*source.contributions[sourceIndex].over_coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5*source.value;
         new_origins[new_origins_size].over_coefficient *= 0.5*std::fabs(source.value);
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.contributions[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5*value;
         new_origins[new_origins_size].over_coefficient *= 0.5*std::fabs(value);
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < contributions_size) {
      new_origins[new_origins_size] = contributions[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5*source.value;
      new_origins[new_origins_size].over_coefficient *= 0.5*std::fabs(source.value);
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.contributions_size) {
      new_origins[new_origins_size] = source.contributions[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5*value;
      new_origins[new_origins_size].over_coefficient *= 0.5*std::fabs(value);
      ++new_origins_size;
      ++sourceIndex;
   }
}

#ifdef _ORIGINS_DOMAIN
template <typename Type, int N>
void
OriginVector<Type, N>::multDomain(const Origin* sourceBegin, int sizeSource,
      Type sourceValue, float sourceAdditionalInterval, Type sourceCenter,
      float sourceMin, float sourceMax,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin, Type& newDomainCenter, float& newMinDomain,
      float& newMaxDomain) const {
   newCoeffWithoutOrigin = std::fabs(domain_additional_interval*(sourceValue+sourceCenter))
      + std::fabs(sourceAdditionalInterval*(value+domain_center));
   newDomainCenter = domain_center*(sourceValue + sourceCenter) + value*sourceCenter;

   float newMin = std::min(std::min(min_domain*sourceMin, min_domain*sourceMax),
         std::min(max_domain*sourceMax, max_domain*sourceMin));
   newMaxDomain = std::max(std::max(min_domain*sourceMin, min_domain*sourceMax),
         std::max(max_domain*sourceMax, max_domain*sourceMin));
   newMinDomain = newMin;
   // multplication of noise symbols
   for (int thisIndex = 0; thisIndex < domain_size; ++thisIndex) {
      for (int sourceIndex = 0; sourceIndex < sizeSource; ++sourceIndex) {
         auto compare = domain[thisIndex].compare(*(sourceBegin + sourceIndex));
         if (compare < 0)
            break;
         if (compare > 0) {
            int thisFoundIndex = -1;
            for (int thisSymmetricIndex = 0; thisSymmetricIndex < thisIndex; ++thisSymmetricIndex) {
               auto otherCompare = domain[thisSymmetricIndex].compare(*(sourceBegin+sourceIndex));
               if (otherCompare >= 0) {
                  if (otherCompare == 0)
                     thisFoundIndex = thisSymmetricIndex;
                  break;
               }
            }
            int sourceFoundIndex = -1;
            if (thisFoundIndex >= 0) {
               for (int sourceSymmetricIndex = sourceIndex+1; sourceSymmetricIndex < sizeSource; ++sourceSymmetricIndex) {
                  auto otherCompare = domain[thisIndex].compare(*(sourceBegin + sourceSymmetricIndex));
                  if (otherCompare <= 0) {
                     if (otherCompare == 0)
                        sourceFoundIndex = sourceSymmetricIndex;
                     break;
                  }
               }
            }
            if (thisFoundIndex >= 0 && sourceFoundIndex >= 0)
               newCoeffWithoutOrigin += domain[thisIndex].coefficient*(sourceBegin + sourceIndex)->coefficient
                  + domain[thisFoundIndex].coefficient*(sourceBegin +sourceFoundIndex)->coefficient;
            else
               newCoeffWithoutOrigin += domain[thisIndex].coefficient*(sourceBegin + sourceIndex)->coefficient;
         }
         else { // compare == 0
            float newShift = 0.5*domain[thisIndex].coefficient*(sourceBegin + sourceIndex)->coefficient;
            newDomainCenter += newShift;
            newCoeffWithoutOrigin += std::fabs(newShift);
         }
      }
   }

   // multplication by the constant
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < domain_size && sourceIndex < sizeSource) {
      auto compare = domain[thisIndex].compare(*(sourceBegin + sourceIndex));
      if (compare == 0) {
         new_origins[new_origins_size] = domain[thisIndex];
         new_origins[new_origins_size].coefficient *= (sourceValue + sourceCenter);
         new_origins[new_origins_size].coefficient += (sourceBegin + sourceIndex)->coefficient*(value + domain_center);
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = domain[thisIndex];
         new_origins[new_origins_size].coefficient *= (sourceValue + sourceCenter);
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = *(sourceBegin + sourceIndex);
         new_origins[new_origins_size].coefficient *= (value + domain_center);
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < domain_size) {
      new_origins[new_origins_size] = domain[thisIndex];
      new_origins[new_origins_size].coefficient *= (sourceValue + sourceCenter);
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < sizeSource) {
      new_origins[new_origins_size] = *(sourceBegin + sourceIndex);
      new_origins[new_origins_size].coefficient *= (value + domain_center);
      ++new_origins_size;
      ++sourceIndex;
   }

   // apply domain reduction
   float sum = newCoeffWithoutOrigin;
   std::for_each(new_origins.begin(), new_origins.begin() + new_origins_size,
      [&sum](const Origin& origin)
         {  sum += std::fabs(origin.coefficient); });
   if (newMinDomain < -sum + value + domain_center)
      newMinDomain = -sum + value + domain_center;
   if (newMaxDomain > sum + value + domain_center)
      newMaxDomain = sum + value + domain_center;
}
#endif

template <typename Type, int N>
inline void
OriginVector<Type, N>::multAssign(const OriginVector<Type, N>& source) {
   Type newValue = value*source.value;
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   {  multContribution(source, new_origins, new_origins_size, contribution_without_origin, over_contribution_without_origin);
      contributions_size = 0;
      OriginVector<Type, N>::copy<2*N>(contributions, contributions_size, new_origins, new_origins_size,
            contribution_without_origin, &over_contribution_without_origin);
   }
#ifdef _ORIGINS_DOMAIN
   new_origins_size = 0;
   {  multDomain(source, new_origins, new_origins_size, domain_additional_interval, domain_center,
         min_domain, max_domain);
      domain_size = 0;
      OriginVector<Type, N>::copy<2*N>(domain, domain_size, new_origins, new_origins_size,
            domain_additional_interval);
   }
#endif
#ifdef _ORIGINS_ERROR
   double newError = -EFT::FastTwoProd(value, source.value, newValue);
   accumulated_error *= (source.value + source.accumulated_error);
   accumulated_error += newError;
   accumulated_error += value*source.accumulated_error;
#endif
   value = newValue;
   if (value == 0.0)
      clear_origin();
   restrictContributions();
}

template <typename Type, int N>
void
OriginVector<Type, N>::divContribution(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const {
   newCoeffWithoutOrigin = 0.5*contribution_without_origin/source.value
         + value*0.5*source.contribution_without_origin/(source.value*source.value);
   newOverCoeffWithoutOrigin = 0.5*over_contribution_without_origin/std::fabs(source.value)
         + std::fabs(value)*0.5*source.over_contribution_without_origin/(source.value*source.value);
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < contributions_size && sourceIndex < source.contributions_size) {
      auto compare = contributions[thisIndex].compare(source.contributions[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*contributions[thisIndex].coefficient/source.value
                  + value*0.5*source.contributions[sourceIndex].coefficient/(source.value*source.value);
         new_origins[new_origins_size].over_coefficient
            = 0.5*contributions[thisIndex].over_coefficient/std::fabs(source.value)
                  + std::fabs(value)*0.5*source.contributions[sourceIndex].over_coefficient/(source.value*source.value);
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5/source.value;
         new_origins[new_origins_size].over_coefficient *= 0.5/std::fabs(source.value);
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.contributions[sourceIndex];
         new_origins[new_origins_size].coefficient *= value*0.5/(source.value*source.value);
         new_origins[new_origins_size].over_coefficient *= std::fabs(value)*0.5/(source.value*source.value);
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < contributions_size) {
      new_origins[new_origins_size] = contributions[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5/source.value;
      new_origins[new_origins_size].over_coefficient *= 0.5/std::fabs(source.value);
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.contributions_size) {
      new_origins[new_origins_size] = source.contributions[sourceIndex];
      new_origins[new_origins_size].coefficient *= value*0.5/(source.value*source.value);
      new_origins[new_origins_size].over_coefficient *= std::fabs(value)*0.5/(source.value*source.value);
      ++new_origins_size;
      ++sourceIndex;
   }
}

#ifdef _ORIGINS_DOMAIN
template <typename Type, int N>
void
OriginVector<Type, N>::divDomain(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin, Type& newDomainCenter,
      float& newMinDomain, float& newMaxDomain) const {
   Type factor = 0;
   float min_derivative = 0, max_derivative = 0;
   typedef std::pair<float, std::pair<int, int>> GradientElement;
   std::array<GradientElement, N> gradients;
   int size_gradient = 0;
   int thisIndex = 0, sourceIndex = 0;
   // compute first min_derivative, max_derivative
   while (sourceIndex < source.domain_size) {
      int compare = -1;
      for (; compare < 0 && thisIndex < domain_size; ++thisIndex)
         compare = domain[thisIndex].compare(source.domain[sourceIndex]);
      if (compare < 0)
         break;
      if (source.domain[sourceIndex].coefficient != 0) {
         if (compare > 0 || domain[thisIndex].coefficient == 0) {
            auto diff = std::fabs(source.domain[sourceIndex].coefficient);
            min_derivative -= diff;
            max_derivative += diff;
         }
         else if (domain[thisIndex].coefficient > 0) {
            min_derivative -= source.domain[sourceIndex].coefficient;
            max_derivative -= source.domain[sourceIndex].coefficient;
         }
         else { // domain[sourceIndex].coefficient < 0
            min_derivative += source.domain[sourceIndex].coefficient;
            max_derivative += source.domain[sourceIndex].coefficient;
         }
         gradients[size_gradient++] =
            {  domain[thisIndex].coefficient / source.domain[sourceIndex].coefficient,
               { sourceIndex, thisIndex } };
      }
      ++sourceIndex;
      if (compare == 0)
         ++thisIndex;
   }
   std::sort(gradients.begin(), gradients.begin() + size_gradient);
   int gradient_index = 0;
   // find gradient as better coefficient that reduces the amplitude of *this - gradient*source
   float gradient = 0.0;
   if (min_derivative > 0) {
      auto iter = std::upper_bound(gradients.begin(), gradients.begin() + size_gradient,
            GradientElement( 0.0, { 0, 0 } ));
      int gradient_index = iter - gradients.begin();
      while (gradient_index >= 0 && gradients[gradient_index].first == 0.0)
         --gradient_index;
      do {
         if (gradient_index < 0)
            break;
         if (domain[gradients[gradient_index].second.second].coefficient > 0)
            min_derivative -= 2*source.domain[gradients[gradient_index].second.first].coefficient;
         else
            min_derivative += 2*source.domain[gradients[gradient_index].second.first].coefficient;
         while (gradient_index >= 0 && gradients[gradient_index].first == gradient) {
            --gradient_index;
            if (domain[gradients[gradient_index].second.second].coefficient > 0)
               min_derivative -= 2*source.domain[gradients[gradient_index].second.first].coefficient;
            else
               min_derivative += 2*source.domain[gradients[gradient_index].second.first].coefficient;
         }
         if (min_derivative > 0 && gradient_index >= 0)
            gradient = gradients[gradient_index].first;
      } while (min_derivative > 0);
   }
   else if (max_derivative < 0) {
      auto iter = std::lower_bound(gradients.begin(), gradients.begin() + size_gradient,
            GradientElement(0.0, { 0, 0 }));
      int gradient_index = iter - gradients.begin();
      while (gradient_index < size_gradient && gradients[gradient_index].first == 0.0)
         --gradient_index;
      do {
         if (gradient_index == size_gradient)
            break;
         if (domain[gradients[gradient_index].second.second].coefficient > 0)
            max_derivative -= 2*source.domain[gradients[gradient_index].second.first].coefficient;
         else
            max_derivative += 2*source.domain[gradients[gradient_index].second.first].coefficient;
         while (gradient_index < size_gradient && gradients[gradient_index].first == gradient) {
            ++gradient_index;
            if (domain[gradients[gradient_index].second.second].coefficient > 0)
               max_derivative -= 2*source.domain[gradients[gradient_index].second.first].coefficient;
            else
               max_derivative += 2*source.domain[gradients[gradient_index].second.first].coefficient;
         }
         if (max_derivative > 0 && gradient_index < size_gradient)
            gradient = gradients[gradient_index].first;
      } while (max_derivative < 0);
   }

   OriginVector<Type, N> sourceInverse(source);
   sourceInverse.inverseAssign();
   // std::array<Origin, (2*N)> inverseSource;
   // int inverseSourceSize = 0;
   // float inverseCoeffWithoutOrigin = 0, inverseCenter = 0;
   // source.inverseDomain(inverseSource, inverseSourceSize, inverseCoeffWithoutOrigin, inverseCenter);

   OriginVector<Type, N> sourceCopy(source);
   sourceCopy.multAssign(gradient);
   std::array<Origin, (2*N)> minimizeThis;
   int minimizeThisSize = 0;
   float minimizeCoeffWithoutOrigin = 0;
   Type minimizeCenter = 0;
   float minimizeMinDomain = 0, minimizeMaxDomain = 0;
   addDomain(sourceCopy, minimizeThis, minimizeThisSize, false /* is_plus */,
         minimizeCoeffWithoutOrigin, minimizeCenter, minimizeMinDomain, minimizeMaxDomain);

   sourceInverse.multDomain(&*minimizeThis.begin(), minimizeThisSize, value - gradient*source.value,
         minimizeCoeffWithoutOrigin, minimizeCenter, minimizeMinDomain, minimizeMaxDomain,
         new_origins, new_origins_size, newCoeffWithoutOrigin, newDomainCenter, newMinDomain,
         newMaxDomain);
   newDomainCenter -= gradient;

   // apply domain reduction
   float sum = newCoeffWithoutOrigin;
   std::for_each(new_origins.begin(), new_origins.begin() + new_origins_size,
      [&sum](const Origin& origin)
         {  sum += std::fabs(origin.coefficient); });
   if (newMinDomain < sum + value + domain_center)
      newMinDomain = sum + value + domain_center;
   if (newMaxDomain > sum + value + domain_center)
      newMaxDomain = sum + value + domain_center;
}
#endif

template <typename Type, int N>
inline void
OriginVector<Type, N>::divAssign(const OriginVector<Type, N>& source) {
   Type newValue = value/source.value;
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   {  divContribution(source, new_origins, new_origins_size, contribution_without_origin, over_contribution_without_origin);
      contributions_size = 0;
      OriginVector<Type, N>::copy<2*N>(contributions, contributions_size, new_origins, new_origins_size,
            contribution_without_origin, &over_contribution_without_origin);
   }
#ifdef _ORIGINS_DOMAIN
   new_origins_size = 0;
   {  divDomain(source, new_origins, new_origins_size, domain_additional_interval, domain_center,
         min_domain, max_domain);
      domain_size = 0;
      OriginVector<Type, N>::copy<2*N>(domain, domain_size, new_origins, new_origins_size,
            domain_additional_interval);
   }
#endif
#ifdef _ORIGINS_ERROR
   double newError = -EFT::RemainderDiv(value, source.value, newValue);
   accumulated_error = ((newError + accumulated_error) - newValue*source.accumulated_error)
         / (source.value + source.accumulated_error);
#endif
   value = newValue;
   if (value == 0.0)
      clear_origin();
   restrictContributions();
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::logAssign() {
   value = std::log(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::expAssign() {
   value = std::exp(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::sqrtAssign() {
   value = std::sqrt(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::sinAssign() {
   value = std::sin(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::cosAssign() {
   value = std::cos(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::tanAssign() {
   value = std::tan(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::asinAssign() {
   value = std::asin(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::acosAssign() {
   value = std::acos(value);
   assert(false);
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::atanAssign() {
   value = std::atan(value);
   assert(false);
}

template <typename Type, int N>
Type
OriginVector<Type, N>::atan2(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const {
   Type newValue = std::atan2(value, source.value);
   // TODO = use divAssign to correctly implement contributions
   newCoeffWithoutOrigin = 0.5*contribution_without_origin + 0.5*source.contribution_without_origin;
   newOverCoeffWithoutOrigin = 0.5*over_contribution_without_origin + 0.5*source.over_contribution_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < contributions_size && sourceIndex < source.contributions_size) {
      auto compare = contributions[thisIndex].compare(source.contributions[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*contributions[thisIndex].coefficient + 0.5*source.contributions[sourceIndex].coefficient;
         new_origins[new_origins_size].over_coefficient
            = 0.5*contributions[thisIndex].over_coefficient + 0.5*source.contributions[sourceIndex].over_coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         new_origins[new_origins_size].over_coefficient *= 0.5;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.contributions[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         new_origins[new_origins_size].over_coefficient *= 0.5;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < contributions_size) {
      new_origins[new_origins_size] = contributions[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      new_origins[new_origins_size].over_coefficient *= 0.5;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.contributions_size) {
      new_origins[new_origins_size] = source.contributions[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      new_origins[new_origins_size].over_coefficient *= 0.5;
      ++new_origins_size;
      ++sourceIndex;
   }
   assert(false);
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::atan2Assign(const OriginVector<Type, N>& source) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = atan2(source, new_origins, new_origins_size, contribution_without_origin, over_contribution_without_origin);
   contributions_size = 0;
   copy<2*N>(contributions, contributions_size, new_origins, new_origins_size, contribution_without_origin, &over_contribution_without_origin);
   assert(false);
}

template <typename Type, int N>
Type
OriginVector<Type, N>::pow(const OriginVector<Type, N>& source,
      std::array<Origin, (2*N)>& new_origins, int& new_origins_size,
      float& newCoeffWithoutOrigin, float& newOverCoeffWithoutOrigin) const {
   Type newValue = std::pow(value, source.value);
   newCoeffWithoutOrigin = 0.5*contribution_without_origin + 0.5*source.contribution_without_origin;
   newOverCoeffWithoutOrigin = 0.5*over_contribution_without_origin + 0.5*source.over_contribution_without_origin;
   int thisIndex = 0, sourceIndex = 0;
   while (thisIndex < contributions_size && sourceIndex < source.contributions_size) {
      auto compare = contributions[thisIndex].compare(source.contributions[sourceIndex]);
      if (compare == 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient
            = 0.5*contributions[thisIndex].coefficient + 0.5*source.contributions[sourceIndex].coefficient;
         new_origins[new_origins_size].over_coefficient
            = 0.5*contributions[thisIndex].over_coefficient + 0.5*source.contributions[sourceIndex].over_coefficient;
         ++new_origins_size;
         ++thisIndex;
         ++sourceIndex;
      }
      else if (compare < 0) {
         new_origins[new_origins_size] = contributions[thisIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         new_origins[new_origins_size].over_coefficient *= 0.5;
         ++new_origins_size;
         ++thisIndex;
      }
      else {
         new_origins[new_origins_size] = source.contributions[sourceIndex];
         new_origins[new_origins_size].coefficient *= 0.5;
         new_origins[new_origins_size].over_coefficient *= 0.5;
         ++new_origins_size;
         ++sourceIndex;
      }
   }
   while (thisIndex < contributions_size) {
      new_origins[new_origins_size] = contributions[thisIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      new_origins[new_origins_size].over_coefficient *= 0.5;
      ++new_origins_size;
      ++thisIndex;
   }
   while (sourceIndex < source.contributions_size) {
      new_origins[new_origins_size] = source.contributions[sourceIndex];
      new_origins[new_origins_size].coefficient *= 0.5;
      new_origins[new_origins_size].over_coefficient *= 0.5;
      ++new_origins_size;
      ++sourceIndex;
   }
   assert(false);
   return newValue;
}

template <typename Type, int N>
inline void
OriginVector<Type, N>::powAssign(const OriginVector<Type, N>& source) {
   std::array<Origin, (2*N)> new_origins;
   int new_origins_size = 0;
   value = pow(source, new_origins, new_origins_size, contribution_without_origin, over_contribution_without_origin);
   contributions_size = 0;
   copy<2*N>(contributions, contributions_size, new_origins, new_origins_size, contribution_without_origin, &over_contribution_without_origin);
   assert(false);
}

struct double_st;
struct float_st : public OriginVector<float, 16> {
  private:
   typedef OriginVector<float, 16> inherited;

  public:
   float_st() = default;
   float_st(const float& avalue, bool doesTrack=false)
      :  inherited(avalue, doesTrack) {}
   float_st(const float_st& source) = default;
   float_st(const double_st& source);
};

struct double_st : public OriginVector<double, 16> {
  private:
   typedef OriginVector<double, 16> inherited;

  public:
   double_st() = default;
   double_st(const double& avalue, bool doesTrack=false)
      :  inherited(avalue, doesTrack) {}
   double_st(const double_st& source) = default;
   double_st(const float_st& source);
};

inline
float_st::float_st(const double_st& source) : inherited(source) {}

inline
double_st::double_st(const float_st& source) : inherited(source) {}

