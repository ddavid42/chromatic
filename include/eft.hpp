#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include "eft.h"

template <typename Type>
struct EnhancedDouble {
   Type value = 0.0;
   double error = 0.0;

   double accumulated_error() const { return error; }

   EnhancedDouble() = default;
   EnhancedDouble(const Type& avalue);
   EnhancedDouble(const EnhancedDouble<Type>& source) = default;
   template <typename OtherType>
   EnhancedDouble(const EnhancedDouble<OtherType>& source);

   void print(std::ostream& out) const
      {  out << value << ", error=" << error; }
   void printRelative(std::ostream& out) const
      {  out << value << ", relative error=" << error/value; }
   EnhancedDouble& operator=(const EnhancedDouble<Type>& source) = default;
   template <typename OtherType>
   void assign(const EnhancedDouble<OtherType>& source);

   void addAssign(const EnhancedDouble<Type>& source, bool is_plus);
   void multAssign(const Type& source);
   void divAssign(const Type& source);
   void multAssign(const EnhancedDouble<Type>& source);
   void divAssign(const EnhancedDouble<Type>& source);
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
   void atan2Assign(const EnhancedDouble<Type>& source);
   void powAssign(const EnhancedDouble<Type>& source);
   void oppositeAssign()
      {  value = -value;
         error = -error;
      }
   void inverseAssign()
      {  Type result = 1.0/value;
         double new_error = -EFT::RemainderDiv((Type) 1.0, value, result);
         error = (new_error - result*error)/(value + error);
         value = result;
      }
   void maxAssign(const EnhancedDouble<Type>& source)
      {  if (value < source.value)
            operator=(source);
      }
   void minAssign(const EnhancedDouble<Type>& source)
      {  if (value > source.value)
            operator=(source);
      }
   void absAssign()
      {  if (value < 0) {
            value = -value;
            error = -error;
         }
      }
   void floorAssign()
      {  Type result = std::floor(value);
         error = std::floor(value - error) - result;
         value = result;
      }
   void fmodAssign(const EnhancedDouble<Type>& source)
      {  
         Type result = std::fmod(value, source.value);
         error = std::fmod(value - error, source.value - source.error)
               - result;
         value = result;
      }
   void modfAssign(EnhancedDouble<Type>& source)
      {  
         Type result, sourceResult = 0;
         result = std::modf(value, &sourceResult);
         double sourceWithError = 0;
         error = std::modf(value - error, &sourceWithError) - result;
         source.value = sourceResult;
         source.error = sourceResult - sourceWithError;
         value = result;
      }
};

template <typename Type>
inline
EnhancedDouble<Type>::EnhancedDouble(const Type& avalue) :  value(avalue) {}

template <typename Type>
template <typename OtherType>
inline
EnhancedDouble<Type>::EnhancedDouble(const EnhancedDouble<OtherType>& source)
   :  value(source.value)
   {
      error = value - source.value + source.error;
   }

template <typename Type>
template <typename OtherType>
inline void
EnhancedDouble<Type>::assign(const EnhancedDouble<OtherType>& source) {
   value = source.value;
   error = value - source.value + source.error;
}

template <typename Type>
inline void
EnhancedDouble<Type>::addAssign(const EnhancedDouble<Type>& source, bool is_plus) {
   Type newValue = is_plus ? (value + source.value) : (value - source.value);
   double newError = -EFT::TwoSum(value, is_plus ? source.value : -source.value, newValue);
   error += source.error;
   error += newError;
   value = newValue;
}

template <typename Type>
inline void
EnhancedDouble<Type>::multAssign(const Type& source) {
   Type result = value * source;
   double newError = -EFT::FastTwoProd(value, source, result);
   error *= source;
   error += newError;
   value = result;
}

template <typename Type>
inline void
EnhancedDouble<Type>::divAssign(const Type& source) {
   Type result = value / source;
   double newError = -EFT::RemainderDiv(value, source, result);
   error = (newError + error) / source;
   value = result;
}

template <typename Type>
inline void
EnhancedDouble<Type>::multAssign(const EnhancedDouble<Type>& source) {
   Type newValue = value*source.value;
   double newError = -EFT::FastTwoProd(value, source.value, newValue);
   error *= (source.value + source.error);
   error += newError;
   error += value*source.error;
   value = newValue;
}

template <typename Type>
inline void
EnhancedDouble<Type>::divAssign(const EnhancedDouble<Type>& source) {
   Type newValue = value/source.value;
   double newError = -EFT::RemainderDiv(value, source.value, newValue);
   error = ((newError + error) - newValue*source.error)
         / (source.value + source.error);
   value = newValue;
}

template <typename Type>
inline void
EnhancedDouble<Type>::logAssign() {
   value = std::log(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::expAssign() {
   value = std::exp(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::sqrtAssign() {
   value = std::sqrt(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::sinAssign() {
   value = std::sin(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::cosAssign() {
   value = std::cos(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::tanAssign() {
   value = std::tan(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::asinAssign() {
   value = std::asin(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::acosAssign() {
   value = std::acos(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::atanAssign() {
   value = std::atan(value);
   assert(false);
}

template <typename Type>
inline void
EnhancedDouble<Type>::atan2Assign(const EnhancedDouble<Type>& source) {
   Type newValue = std::atan2(value, source.value);
   value = newValue;
   assert(false);
}


template <typename Type>
inline void
EnhancedDouble<Type>::powAssign(const EnhancedDouble<Type>& source) {
   Type newValue = std::pow(value, source.value);
   value = newValue;
   assert(false);
}

struct double_st;
struct float_st : public EnhancedDouble<float> {
  private:
   typedef EnhancedDouble<float> inherited;

  public:
   float_st() = default;
   float_st(const float& avalue) : inherited(avalue) {}
   float_st(const float_st& source) = default;
   float_st(const double_st& source);
};

struct double_st : public EnhancedDouble<double> {
  private:
   typedef EnhancedDouble<double> inherited;

  public:
   double_st() = default;
   double_st(const double& avalue) : inherited(avalue) {}
   double_st(const double_st& source) = default;
   double_st(const float_st& source);
};

inline
float_st::float_st(const double_st& source) : inherited(source) {}

inline
double_st::double_st(const float_st& source) : inherited(source) {}


extern "C" {
   void eft_instrument_init(void* adiagnosis_file);
   void eft_instrument_end(void* adiagnosis_file);
   bool eft_has_custom_print();
}

// inline bool eft_has_custom_print() { return true; }

class DoubleSt;
class FloatSt : public float_st {
  private:
   typedef float_st inherited;

  public:
   FloatSt() {}
   FloatSt(float val)
      : inherited(val) {}
   FloatSt(double val)
      : inherited(val) {}
   FloatSt(const char* val)
      : inherited(::atof(val)) {}
   FloatSt(int val)
      : inherited(val) {}
   FloatSt(long long int val)
      : inherited(val) {}
   FloatSt(const FloatSt&) = default;
   FloatSt(FloatSt&& source) = default;
   FloatSt(const DoubleSt& source);
   FloatSt(DoubleSt&& source);

   FloatSt& operator=(float val) { return (FloatSt&) inherited::operator=(FloatSt(val)); }
   FloatSt& operator=(double val) { return (FloatSt&) inherited::operator=(FloatSt(val)); }
   FloatSt& operator=(int val) { return (FloatSt&) inherited::operator=(FloatSt(val)); }
   FloatSt& operator=(long long int val) { return (FloatSt&) inherited::operator=(FloatSt(val)); }
   FloatSt& operator=(const FloatSt&) = default;
   FloatSt& operator=(FloatSt&& source) = default;
   FloatSt& operator=(const DoubleSt&);
   FloatSt& operator=(DoubleSt&& source);

   void print(std::ostream& adiagnosis_file, const char* varname) const
      {  adiagnosis_file << varname << ": ";
         inherited::print(adiagnosis_file);
         adiagnosis_file << '\n';
      }
   void retrieveRelativeError(double& error) const;
   void retrieveBoundsAndAbsoluteError(float& min, float& max, double& error) const;

   operator long long int() const { return inherited::value; }
   operator int() const { return inherited::value; }

   friend FloatSt operator-(FloatSt&& source);
   FloatSt& operator+=(const FloatSt& source);
   FloatSt& operator+=(FloatSt&& source);
   FloatSt& operator+=(DoubleSt&& source);
   friend FloatSt operator+(const FloatSt& first, const FloatSt& second);
   friend FloatSt operator+(const FloatSt& first, FloatSt&& second);
   friend FloatSt operator+(FloatSt&& first, FloatSt&& second);
   friend DoubleSt operator+(const FloatSt& first, DoubleSt&& second);
   friend DoubleSt operator+(FloatSt&& first, DoubleSt&& second);
   FloatSt& operator-=(const FloatSt& source);
   FloatSt& operator-=(FloatSt&& source);
   FloatSt& operator-=(DoubleSt&& source);
   friend FloatSt operator-(const FloatSt& first, const FloatSt& second);
   friend FloatSt operator-(const FloatSt& first, FloatSt&& second);
   friend FloatSt operator-(FloatSt&& first, FloatSt&& second);
   friend DoubleSt operator-(const FloatSt& first, DoubleSt&& second);
   friend DoubleSt operator-(FloatSt&& first, DoubleSt&& second);
   FloatSt& operator*=(const FloatSt& source);
   FloatSt& operator*=(FloatSt&& source);
   FloatSt& operator*=(DoubleSt&& source);
   friend FloatSt operator*(const FloatSt& first, const FloatSt& second);
   friend FloatSt operator*(const FloatSt& first, FloatSt&& second);
   friend FloatSt operator*(FloatSt&& first, FloatSt&& second);
   friend DoubleSt operator*(const FloatSt& first, DoubleSt&& second);
   friend DoubleSt operator*(FloatSt&& first, DoubleSt&& second);
   FloatSt& operator/=(const FloatSt& source);
   FloatSt& operator/=(FloatSt&& source);
   FloatSt& operator/=(DoubleSt&& source);
   friend FloatSt operator/(const FloatSt& first, const FloatSt& second);
   friend FloatSt operator/(const FloatSt& first, FloatSt&& second);
   friend FloatSt operator/(FloatSt&& first, FloatSt&& second);
   friend DoubleSt operator/(const FloatSt& first, DoubleSt&& second);
   friend DoubleSt operator/(FloatSt&& first, DoubleSt&& second);

   friend bool operator<(const FloatSt& first, const FloatSt& second)
      {  return first.value < second.value; }
   friend bool operator<=(const FloatSt& first, const FloatSt& second)
      {  return first.value <= second.value; }
   friend bool operator==(const FloatSt& first, const FloatSt& second)
      {  return first.value == second.value; }
   friend bool operator!=(const FloatSt& first, const FloatSt& second)
      {  return first.value != second.value; }
   friend bool operator>(const FloatSt& first, const FloatSt& second)
      {  return first.value > second.value; }
   friend bool operator>=(const FloatSt& first, const FloatSt& second)
      {  return first.value >= second.value; }

   friend bool operator<(const FloatSt& first, float second)
      {  return first.value < second; }
   friend bool operator<=(const FloatSt& first, float second)
      {  return first.value <= second; }
   friend bool operator==(const FloatSt& first, float second)
      {  return first.value == second; }
   friend bool operator!=(const FloatSt& first, float second)
      {  return first.value != second; }
   friend bool operator>(const FloatSt& first, float second)
      {  return first.value > second; }
   friend bool operator>=(const FloatSt& first, float second)
      {  return first.value >= second; }
   friend bool operator<(const FloatSt& first, const DoubleSt& second);
   friend bool operator<=(const FloatSt& first, const DoubleSt& second);
   friend bool operator==(const FloatSt& first, const DoubleSt& second);
   friend bool operator!=(const FloatSt& first, const DoubleSt& second);
   friend bool operator>(const FloatSt& first, const DoubleSt& second);
   friend bool operator>=(const FloatSt& first, const DoubleSt& second);
};

inline FloatSt log(const FloatSt& source)
{  FloatSt result(source);
   result.logAssign();
   return result;
}

inline FloatSt fabs(const FloatSt& source)
{  FloatSt result(source);
   result.fabsAssign();
   return result;
}

inline FloatSt exp(const FloatSt& source)
{  FloatSt result(source);
   result.expAssign();
   return result;
}

inline FloatSt sqrt(const FloatSt& source)
{  FloatSt result(source);
   result.sqrtAssign();
   return result;
}

inline FloatSt sin(const FloatSt& source)
{  FloatSt result(source);
   result.sinAssign();
   return result;
}

inline FloatSt cos(const FloatSt& source)
{  FloatSt result(source);
   result.cosAssign();
   return result;
}

inline FloatSt tan(const FloatSt& source)
{  FloatSt result(source);
   result.tanAssign();
   return result;
}

inline FloatSt asin(const FloatSt& source)
{  FloatSt result(source);
   result.asinAssign();
   return result;
}

inline FloatSt acos(const FloatSt& source)
{  FloatSt result(source);
   result.acosAssign();
   return result;
}

inline FloatSt atan(const FloatSt& source)
{  FloatSt result(source);
   result.atanAssign();
   return result;
}

inline FloatSt atan2(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.atan2Assign(second);
   return result;
}

inline FloatSt pow(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.powAssign(second);
   return result;
}

inline FloatSt max(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.maxAssign(second);
   return result;
}

inline FloatSt min(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.minAssign(second);
   return result;
}

inline FloatSt abs(const FloatSt& source, unsigned instructionId, unsigned typeId)
{  FloatSt result(source);
   result.absAssign();
   return result;
}

inline FloatSt floor(const FloatSt& source)
{  FloatSt result(source);
   result.floorAssign();
   return result;
}

inline FloatSt fmod(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.fmodAssign(second);
   return result;
}

inline FloatSt modf(const FloatSt& first, FloatSt* second)
{  FloatSt result(first);
   result.modfAssign(*second);
   return result;
}

inline bool finite(const FloatSt& source)
{  return finite(source.value); }

inline bool isfinite(const FloatSt& source)
{  return std::isfinite(source.value); }

inline bool isnan(const FloatSt& source)
{  return std::isnan(source.value); }

inline bool isinf(const FloatSt& source)
{  return std::isinf(source.value); }

inline FloatSt operator-(FloatSt&& source)
{  FloatSt result(source);
   result.inherited::oppositeAssign();
   return result;
}

inline FloatSt& FloatSt::operator+=(const FloatSt& source)
{  inherited::addAssign(source, true /* is_plus */);
   return *this;
}

inline FloatSt& FloatSt::operator+=(FloatSt&& source)
{  inherited::addAssign(source, true /* is_plus */);
   return *this;
}

inline FloatSt operator+(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline FloatSt operator+(const FloatSt& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline FloatSt operator+(FloatSt&& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline FloatSt& FloatSt::operator-=(FloatSt&& source)
{  inherited::addAssign(source, false /* is_plus */);
   return *this;
}

inline FloatSt operator-(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline FloatSt operator-(const FloatSt& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline FloatSt operator-(FloatSt&& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline FloatSt& FloatSt::operator*=(const FloatSt& source)
{  inherited::multAssign(source);
   return *this;
}

inline FloatSt& FloatSt::operator*=(FloatSt&& source)
{  inherited::multAssign(source);
   return *this;
}

inline FloatSt operator*(const FloatSt& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::multAssign(second);
   return result;
}

inline FloatSt operator*(FloatSt&& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::multAssign(second);
   return result;
}

inline FloatSt& FloatSt::operator/=(const FloatSt& source)
{  inherited::divAssign(source);
   return *this;
}

inline FloatSt& FloatSt::operator/=(FloatSt&& source)
{  inherited::divAssign(source);
   return *this;
}

inline FloatSt operator/(const FloatSt& first, const FloatSt& second)
{  FloatSt result(first);
   result.inherited::divAssign(second);
   return result;
}

inline FloatSt operator/(const FloatSt& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::divAssign(second);
   return result;
}

inline FloatSt operator/(FloatSt&& first, FloatSt&& second)
{  FloatSt result(first);
   result.inherited::divAssign(second);
   return result;
}

class DoubleSt : public double_st {
  private:
   typedef double_st inherited;

  public:
   DoubleSt() {}
   DoubleSt(double val) : inherited(val) {}
   DoubleSt(const char* val) : inherited(::atof(val)) {}
   DoubleSt(int val) : inherited(val) {}
   DoubleSt(long long int val) : inherited(val) {}
   DoubleSt(const DoubleSt&) = default;
   DoubleSt(DoubleSt&& source) = default;
   DoubleSt(const FloatSt& source)
      {  inherited::assign(source); }
   DoubleSt(FloatSt&& source)
      {  inherited::assign(std::move(source)); }
   DoubleSt(double min, double max);
   DoubleSt(double min, double max, double err, unsigned instructionId, unsigned typeId);

   DoubleSt& operator=(double val) { return (DoubleSt&) inherited::operator=(DoubleSt(val)); }
   DoubleSt& operator=(int val) { return (DoubleSt&) inherited::operator=(DoubleSt(val)); }
   DoubleSt& operator=(long long int val) { return (DoubleSt&) inherited::operator=(DoubleSt(val)); }
   DoubleSt& operator=(const DoubleSt&) = default;
   DoubleSt& operator=(DoubleSt&& source) = default;
   DoubleSt& operator=(FloatSt&& source);

   void print(std::ostream& adiagnosis_file, const char* varname) const
      {  adiagnosis_file << varname << ": ";
         inherited::print(adiagnosis_file);
         adiagnosis_file << '\n';
      }
   void retrieveRelativeError(double& error) const;
   void retrieveBoundsAndAbsoluteError(double& min, double& max, double& error) const;

   operator long long int() const { return inherited::value; }
   operator int() const { return inherited::value; }

   friend DoubleSt operator-(DoubleSt&& source);
   DoubleSt& operator+=(const DoubleSt& source);
   DoubleSt& operator+=(DoubleSt&& source);
   DoubleSt& operator+=(FloatSt&& source);
   friend DoubleSt operator+(const DoubleSt& first, const DoubleSt& second);
   friend DoubleSt operator+(const DoubleSt& first, DoubleSt&& second);
   friend DoubleSt operator+(DoubleSt&& first, DoubleSt&& second);
   friend DoubleSt operator+(const DoubleSt& first, FloatSt&& second);
   friend DoubleSt operator+(DoubleSt&& first, FloatSt&& second);
   DoubleSt& operator-=(const DoubleSt& source);
   DoubleSt& operator-=(DoubleSt&& source);
   DoubleSt& operator-=(FloatSt&& source);
   friend DoubleSt operator-(const DoubleSt& first, const DoubleSt& second);
   friend DoubleSt operator-(const DoubleSt& first, DoubleSt&& second);
   friend DoubleSt operator-(DoubleSt&& first, DoubleSt&& second);
   friend DoubleSt operator-(const DoubleSt& first, FloatSt&& second);
   friend DoubleSt operator-(DoubleSt&& first, FloatSt&& second);
   DoubleSt& operator*=(const DoubleSt& source);
   DoubleSt& operator*=(DoubleSt&& source);
   DoubleSt& operator*=(FloatSt&& source);
   friend DoubleSt operator*(const DoubleSt& first, const DoubleSt& second);
   friend DoubleSt operator*(const DoubleSt& first, DoubleSt&& second);
   friend DoubleSt operator*(DoubleSt&& first, DoubleSt&& second);
   friend DoubleSt operator*(const DoubleSt& first, FloatSt&& second);
   friend DoubleSt operator*(DoubleSt&& first, FloatSt&& second);
   DoubleSt& operator/=(const DoubleSt& source);
   DoubleSt& operator/=(DoubleSt&& source);
   DoubleSt& operator/=(FloatSt&& source);
   friend DoubleSt operator/(const DoubleSt& first, const DoubleSt& second);
   friend DoubleSt operator/(const DoubleSt& first, DoubleSt&& second);
   friend DoubleSt operator/(DoubleSt&& first, DoubleSt&& second);
   friend DoubleSt operator/(const DoubleSt& first, FloatSt&& second);
   friend DoubleSt operator/(DoubleSt&& first, FloatSt&& second);

   friend bool operator<(const DoubleSt& first, double second)
      {  return first.value < second; }
   friend bool operator<=(const DoubleSt& first, double second)
      {  return first.value <= second; }
   friend bool operator==(const DoubleSt& first, double second)
      {  return first.value == second; }
   friend bool operator!=(const DoubleSt& first, double second)
      {  return first.value != second; }
   friend bool operator>(const DoubleSt& first, double second)
      {  return first.value > second; }
   friend bool operator>=(const DoubleSt& first, double second)
      {  return first.value >= second; }

   friend bool operator<(const DoubleSt& first, const DoubleSt& second)
      {  return first.value < second.value; }
   friend bool operator<=(const DoubleSt& first, const DoubleSt& second)
      {  return first.value <= second.value; }
   friend bool operator==(const DoubleSt& first, const DoubleSt& second)
      {  return first.value == second.value; }
   friend bool operator!=(const DoubleSt& first, const DoubleSt& second)
      {  return first.value != second.value; }
   friend bool operator>(const DoubleSt& first, const DoubleSt& second)
      {  return first.value > second.value; }
   friend bool operator>=(const DoubleSt& first, const DoubleSt& second)
      {  return first.value >= second.value; }
   friend bool operator<(const DoubleSt& first, const FloatSt& second)
      {  return first.value < second.value; }
   friend bool operator<=(const DoubleSt& first, const FloatSt& second)
      {  return first.value <= second.value; }
   friend bool operator==(const DoubleSt& first, const FloatSt& second)
      {  return first.value == second.value; }
   friend bool operator!=(const DoubleSt& first, const FloatSt& second)
      {  return first.value != second.value; }
   friend bool operator>(const DoubleSt& first, const FloatSt& second)
      {  return first.value > second.value; }
   friend bool operator>=(const DoubleSt& first, const FloatSt& second)
      {  return first.value >= second.value; }
};

inline bool operator<(const FloatSt& first, const DoubleSt& second)
   {  return first.value < second.value; }
inline bool operator<=(const FloatSt& first, const DoubleSt& second)
   {  return first.value <= second.value; }
inline bool operator==(const FloatSt& first, const DoubleSt& second)
   {  return first.value == second.value; }
inline bool operator!=(const FloatSt& first, const DoubleSt& second)
   {  return first.value != second.value; }
inline bool operator>(const FloatSt& first, const DoubleSt& second)
   {  return first.value > second.value; }
inline bool operator>=(const FloatSt& first, const DoubleSt& second)
   {  return first.value >= second.value; }

inline
FloatSt::FloatSt(const DoubleSt& source)
   {  inherited::assign(source); }

inline
FloatSt::FloatSt(DoubleSt&& source)
   {  inherited::assign(std::move(source)); }

inline DoubleSt fabs(const DoubleSt& source)
{  DoubleSt result(source);
   result.fabsAssign();
   return result;
}

inline DoubleSt log(const DoubleSt& source)
{  DoubleSt result(source);
   result.logAssign();
   return result;
}

inline DoubleSt exp(const DoubleSt& source)
{  DoubleSt result(source);
   result.expAssign();
   return result;
}

inline DoubleSt sqrt(const DoubleSt& source)
{  DoubleSt result(source);
   result.sqrtAssign();
   return result;
}

inline DoubleSt sin(const DoubleSt& source)
{  DoubleSt result(source);
   result.sinAssign();
   return result;
}

inline DoubleSt cos(const DoubleSt& source)
{  DoubleSt result(source);
   result.cosAssign();
   return result;
}

inline DoubleSt tan(const DoubleSt& source)
{  DoubleSt result(source);
   result.tanAssign();
   return result;
}

inline DoubleSt asin(const DoubleSt& source)
{  DoubleSt result(source);
   result.asinAssign();
   return result;
}

inline DoubleSt acos(const DoubleSt& source)
{  DoubleSt result(source);
   result.acosAssign();
   return result;
}

inline DoubleSt atan(const DoubleSt& source)
{  DoubleSt result(source);
   result.atanAssign();
   return result;
}

inline DoubleSt atan2(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.atan2Assign(second);
   return result;
}

inline DoubleSt pow(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.powAssign(second);
   return result;
}

inline DoubleSt max(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.maxAssign(second);
   return result;
}

inline DoubleSt min(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.minAssign(second);
   return result;
}

inline DoubleSt abs(const DoubleSt& source)
{  DoubleSt result(source);
   result.absAssign();
   return result;
}

inline DoubleSt floor(const DoubleSt& source)
{  DoubleSt result(source);
   result.floorAssign();
   return result;
}

inline DoubleSt fmod(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.fmodAssign(second);
   return result;
}

inline DoubleSt modf(const DoubleSt& first, DoubleSt* second)
{  DoubleSt result(first);
   result.modfAssign(*second);
   return result;
}

inline bool finite(const DoubleSt& source)
{  return finite(source.value); }

inline bool isfinite(const DoubleSt& source)
{  return std::isfinite(source.value); }

inline bool isnan(const DoubleSt& source)
{  return std::isnan(source.value); }

inline bool isinf(const DoubleSt& source)
{  return std::isinf(source.value); }

inline FloatSt&
FloatSt::operator=(DoubleSt&& source)
{  inherited::assign(source);
   return *this;
}

inline DoubleSt&
DoubleSt::operator=(FloatSt&& source)
{  inherited::assign(source);
   return *this;
}

inline DoubleSt operator-(DoubleSt&& source)
{  DoubleSt result(source);
   result.inherited::oppositeAssign();
   return result;
}

inline DoubleSt& DoubleSt::operator+=(const DoubleSt& source)
{  inherited::addAssign(source, true /* is_plus */);
   return *this;
}

inline DoubleSt& DoubleSt::operator+=(DoubleSt&& source)
{  inherited::addAssign(source, true /* is_plus */);
   return *this;
}

inline FloatSt& FloatSt::operator+=(DoubleSt&& source)
{  FloatSt sourceConvert(source);
   inherited::addAssign(sourceConvert, true /* is_plus */);
   return *this;
}

inline DoubleSt& DoubleSt::operator+=(FloatSt&& source)
{  DoubleSt sourceConvert(source);
   inherited::addAssign(sourceConvert, true /* is_plus */);
   return *this;
}

inline DoubleSt operator+(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline DoubleSt operator+(const DoubleSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline DoubleSt operator+(DoubleSt&& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline DoubleSt operator+(const DoubleSt& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::addAssign(secondConvert, true /* is_plus */);
   return result;
}

inline DoubleSt operator+(DoubleSt&& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::addAssign(secondConvert, true /* is_plus */);
   return result;
}

inline DoubleSt operator+(const FloatSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.addAssign(second, true /* is_plus */);
   return result;
}

inline DoubleSt operator+(FloatSt&& first, DoubleSt&& second)
{  DoubleSt result(std::move(first));
   result.addAssign(second, true /* is_plus */);
   return result;
}

inline DoubleSt& DoubleSt::operator-=(const DoubleSt& source)
{  inherited::addAssign(source, false /* is_plus */);
   return *this;
}

inline DoubleSt& DoubleSt::operator-=(DoubleSt&& source)
{  inherited::addAssign(source, false /* is_plus */);
   return *this;
}

inline FloatSt& FloatSt::operator-=(DoubleSt&& source)
{  FloatSt sourceConvert(source);
   inherited::addAssign(sourceConvert, false /* is_plus */);
   return *this;
}

inline DoubleSt& DoubleSt::operator-=(FloatSt&& source)
{  DoubleSt sourceConvert(source);
   inherited::addAssign(sourceConvert, false /* is_plus */);
   return *this;
}

inline DoubleSt operator-(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline DoubleSt operator-(const DoubleSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline DoubleSt operator-(DoubleSt&& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline DoubleSt operator-(const DoubleSt& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::addAssign(secondConvert, false /* is_plus */);
   return result;
}

inline DoubleSt operator-(DoubleSt&& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::addAssign(secondConvert, false /* is_plus */);
   return result;
}

inline DoubleSt operator-(const FloatSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.addAssign(second, false /* is_plus */);
   return result;
}

inline DoubleSt operator-(FloatSt&& first, DoubleSt&& second)
{  DoubleSt result(std::move(first));
   result.addAssign(second, false /* is_plus */);
   return result;
}

inline DoubleSt& DoubleSt::operator*=(const DoubleSt& source)
{  inherited::multAssign(source);
   return *this;
}

inline DoubleSt& DoubleSt::operator*=(DoubleSt&& source)
{  inherited::multAssign(source);
   return *this;
}

inline FloatSt& FloatSt::operator*=(DoubleSt&& source)
{  FloatSt sourceConvert(source);
   inherited::multAssign(sourceConvert);
   return *this;
}

inline DoubleSt& DoubleSt::operator*=(FloatSt&& source)
{  DoubleSt sourceConvert(source);
   inherited::multAssign(sourceConvert);
   return *this;
}

inline DoubleSt operator*(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.inherited::multAssign(second);
   return result;
}

inline DoubleSt operator*(const DoubleSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::multAssign(second);
   return result;
}

inline DoubleSt operator*(DoubleSt&& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::multAssign(second);
   return result;
}

inline DoubleSt operator*(const DoubleSt& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::multAssign(secondConvert);
   return result;
}

inline DoubleSt operator*(DoubleSt&& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::multAssign(secondConvert);
   return result;
}

inline DoubleSt operator*(const FloatSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.multAssign(second);
   return result;
}

inline DoubleSt operator*(FloatSt&& first, DoubleSt&& second)
{  DoubleSt result(std::move(first));
   result.multAssign(second);
   return result;
}

inline DoubleSt& DoubleSt::operator/=(const DoubleSt& source)
{  inherited::divAssign(source);
   return *this;
}

inline DoubleSt& DoubleSt::operator/=(DoubleSt&& source)
{  inherited::divAssign(source);
   return *this;
}

inline FloatSt& FloatSt::operator/=(DoubleSt&& source)
{  FloatSt sourceConvert(source);
   inherited::divAssign(sourceConvert);
   return *this;
}

inline DoubleSt& DoubleSt::operator/=(FloatSt&& source)
{  DoubleSt sourceConvert(source);
   inherited::divAssign(sourceConvert);
   return *this;
}

inline DoubleSt operator/(const DoubleSt& first, const DoubleSt& second)
{  DoubleSt result(first);
   result.inherited::divAssign(second);
   return result;
}

inline DoubleSt operator/(const DoubleSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::divAssign(second);
   return result;
}

inline DoubleSt operator/(DoubleSt&& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.inherited::divAssign(second);
   return result;
}

inline DoubleSt operator/(const DoubleSt& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::divAssign(secondConvert);
   return result;
}

inline DoubleSt operator/(DoubleSt&& first, FloatSt&& second)
{  DoubleSt result(first);
   DoubleSt secondConvert(second);
   result.inherited::divAssign(secondConvert);
   return result;
}

inline DoubleSt operator/(const FloatSt& first, DoubleSt&& second)
{  DoubleSt result(first);
   result.divAssign(second);
   return result;
}

inline DoubleSt operator/(FloatSt&& first, DoubleSt&& second)
{  DoubleSt result(std::move(first));
   result.divAssign(second);
   return result;
}

typedef float old_float;
typedef double old_double;

#define float FloatSt
#define double DoubleSt

