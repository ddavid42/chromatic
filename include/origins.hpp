#pragma once

#include "origins.h"

extern "C" {
   void origins_instrument_init(void* adiagnosis_file);
   void origins_instrument_end(void* adiagnosis_file);
   bool origins_has_custom_print();
}

// inline bool origins_has_custom_print() { return true; }

class OriginDouble;
class OriginFloat : public float_st {
  private:
   typedef float_st inherited;

  public:
   OriginFloat() {}
   OriginFloat(float val, bool doesTrack=false)
      : inherited(val, doesTrack) {}
   OriginFloat(double val, bool doesTrack=false)
      : inherited(val, doesTrack) {}
   OriginFloat(const char* val, bool doesTrack=false)
      : inherited(::atof(val), doesTrack) {}
   OriginFloat(int val, bool doesTrack=false)
      : inherited(val, doesTrack) {}
   OriginFloat(long long int val, bool doesTrack=false)
      : inherited(val, doesTrack) {}
   OriginFloat(const OriginFloat&) = default;
   OriginFloat(OriginFloat&& source) = default;
   OriginFloat(const OriginDouble& source);
   OriginFloat(OriginDouble&& source);
   OriginFloat(float min, float max) : inherited((min+max)/2.0, true) {}
   OriginFloat(float min, float max, double err)
      : inherited((min+max)/2.0, true) {}

   OriginFloat& operator=(float val) { return (OriginFloat&) inherited::operator=(OriginFloat(val)); }
   OriginFloat& operator=(double val) { return (OriginFloat&) inherited::operator=(OriginFloat(val)); }
   OriginFloat& operator=(int val) { return (OriginFloat&) inherited::operator=(OriginFloat(val)); }
   OriginFloat& operator=(long long int val) { return (OriginFloat&) inherited::operator=(OriginFloat(val)); }
   OriginFloat& operator=(const OriginFloat&) = default;
   OriginFloat& operator=(OriginFloat&& source) = default;
   OriginFloat& operator=(const OriginDouble&);
   OriginFloat& operator=(OriginDouble&& source);

   void print(std::ostream& adiagnosis_file, const char* varname) const
      {  adiagnosis_file << varname << ": ";
         inherited::print(adiagnosis_file);
         adiagnosis_file << '\n';
      }
   void retrieveRelativeError(double& error) const;
   void retrieveBoundsAndAbsoluteError(float& min, float& max, double& error) const;

   operator long long int() const { return inherited::value; }
   operator int() const { return inherited::value; }

   friend OriginFloat operator-(OriginFloat&& source);
   OriginFloat& operator+=(OriginFloat&& source);
   OriginFloat& operator+=(OriginDouble&& source);
   friend OriginFloat operator+(const OriginFloat& first, OriginFloat&& second);
   friend OriginFloat operator+(OriginFloat&& first, OriginFloat&& second);
   friend OriginDouble operator+(const OriginFloat& first, OriginDouble&& second);
   friend OriginDouble operator+(OriginFloat&& first, OriginDouble&& second);
   OriginFloat& operator-=(OriginFloat&& source);
   OriginFloat& operator-=(OriginDouble&& source);
   friend OriginFloat operator-(const OriginFloat& first, OriginFloat&& second);
   friend OriginFloat operator-(OriginFloat&& first, OriginFloat&& second);
   friend OriginDouble operator-(const OriginFloat& first, OriginDouble&& second);
   friend OriginDouble operator-(OriginFloat&& first, OriginDouble&& second);
   OriginFloat& operator*=(OriginFloat&& source);
   OriginFloat& operator*=(OriginDouble&& source);
   friend OriginFloat operator*(const OriginFloat& first, OriginFloat&& second);
   friend OriginFloat operator*(OriginFloat&& first, OriginFloat&& second);
   friend OriginDouble operator*(const OriginFloat& first, OriginDouble&& second);
   friend OriginDouble operator*(OriginFloat&& first, OriginDouble&& second);
   OriginFloat& operator/=(OriginFloat&& source);
   OriginFloat& operator/=(OriginDouble&& source);
   friend OriginFloat operator/(const OriginFloat& first, OriginFloat&& second);
   friend OriginFloat operator/(OriginFloat&& first, OriginFloat&& second);
   friend OriginDouble operator/(const OriginFloat& first, OriginDouble&& second);
   friend OriginDouble operator/(OriginFloat&& first, OriginDouble&& second);

   friend bool operator<(const OriginFloat& first, const OriginFloat& second)
      {  return first.value < second.value; }
   friend bool operator<=(const OriginFloat& first, const OriginFloat& second)
      {  return first.value <= second.value; }
   friend bool operator==(const OriginFloat& first, const OriginFloat& second)
      {  return first.value == second.value; }
   friend bool operator!=(const OriginFloat& first, const OriginFloat& second)
      {  return first.value != second.value; }
   friend bool operator>(const OriginFloat& first, const OriginFloat& second)
      {  return first.value > second.value; }
   friend bool operator>=(const OriginFloat& first, const OriginFloat& second)
      {  return first.value >= second.value; }

   friend bool operator<(const OriginFloat& first, float second)
      {  return first.value < second; }
   friend bool operator<=(const OriginFloat& first, float second)
      {  return first.value <= second; }
   friend bool operator==(const OriginFloat& first, float second)
      {  return first.value == second; }
   friend bool operator!=(const OriginFloat& first, float second)
      {  return first.value != second; }
   friend bool operator>(const OriginFloat& first, float second)
      {  return first.value > second; }
   friend bool operator>=(const OriginFloat& first, float second)
      {  return first.value >= second; }
   friend bool operator<(const OriginFloat& first, const OriginDouble& second);
   friend bool operator<=(const OriginFloat& first, const OriginDouble& second);
   friend bool operator==(const OriginFloat& first, const OriginDouble& second);
   friend bool operator!=(const OriginFloat& first, const OriginDouble& second);
   friend bool operator>(const OriginFloat& first, const OriginDouble& second);
   friend bool operator>=(const OriginFloat& first, const OriginDouble& second);
};

inline OriginFloat log(const OriginFloat& source)
{  OriginFloat result(source);
   result.logAssign();
   return result;
}

inline OriginFloat fabs(const OriginFloat& source)
{  OriginFloat result(source);
   result.fabsAssign();
   return result;
}

inline OriginFloat exp(const OriginFloat& source)
{  OriginFloat result(source);
   result.expAssign();
   return result;
}

inline OriginFloat sqrt(const OriginFloat& source)
{  OriginFloat result(source);
   result.sqrtAssign();
   return result;
}

inline OriginFloat sin(const OriginFloat& source)
{  OriginFloat result(source);
   result.sinAssign();
   return result;
}

inline OriginFloat cos(const OriginFloat& source)
{  OriginFloat result(source);
   result.cosAssign();
   return result;
}

inline OriginFloat tan(const OriginFloat& source)
{  OriginFloat result(source);
   result.tanAssign();
   return result;
}

inline OriginFloat asin(const OriginFloat& source)
{  OriginFloat result(source);
   result.asinAssign();
   return result;
}

inline OriginFloat acos(const OriginFloat& source)
{  OriginFloat result(source);
   result.acosAssign();
   return result;
}

inline OriginFloat atan(const OriginFloat& source)
{  OriginFloat result(source);
   result.atanAssign();
   return result;
}

inline OriginFloat atan2(const OriginFloat& first, const OriginFloat& second)
{  OriginFloat result(first);
   result.atan2Assign(second);
   return result;
}

inline OriginFloat pow(const OriginFloat& first, const OriginFloat& second)
{  OriginFloat result(first);
   result.powAssign(second);
   return result;
}

inline OriginFloat max(const OriginFloat& first, const OriginFloat& second)
{  OriginFloat result(first);
   result.maxAssign(second);
   return result;
}

inline OriginFloat min(const OriginFloat& first, const OriginFloat& second)
{  OriginFloat result(first);
   result.minAssign(second);
   return result;
}

inline OriginFloat abs(const OriginFloat& source, unsigned instructionId, unsigned typeId)
{  OriginFloat result(source);
   result.absAssign();
   return result;
}

inline OriginFloat floor(const OriginFloat& source)
{  OriginFloat result(source);
   result.floorAssign();
   return result;
}

inline OriginFloat fmod(const OriginFloat& first, const OriginFloat& second)
{  OriginFloat result(first);
   result.fmodAssign(second);
   return result;
}

inline OriginFloat modf(const OriginFloat& first, OriginFloat* second)
{  OriginFloat result(first);
   result.modfAssign(*second);
   return result;
}

inline bool finite(const OriginFloat& source)
{  return finite(source.value); }

inline bool isfinite(const OriginFloat& source)
{  return std::isfinite(source.value); }

inline bool isnan(const OriginFloat& source)
{  return std::isnan(source.value); }

inline bool isinf(const OriginFloat& source)
{  return std::isinf(source.value); }

inline OriginFloat operator-(OriginFloat&& source)
{  OriginFloat result(source);
   result.inherited::oppositeAssign();
   return result;
}

inline OriginFloat& OriginFloat::operator+=(OriginFloat&& source)
{  inherited::addAssign(source, true /* is_plus */);
   return *this;
}

inline OriginFloat operator+(const OriginFloat& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline OriginFloat operator+(OriginFloat&& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline OriginFloat& OriginFloat::operator-=(OriginFloat&& source)
{  inherited::addAssign(source, false /* is_plus */);
   return *this;
}

inline OriginFloat operator-(const OriginFloat& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline OriginFloat operator-(OriginFloat&& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline OriginFloat& OriginFloat::operator*=(OriginFloat&& source)
{  inherited::multAssign(source);
   return *this;
}

inline OriginFloat operator*(const OriginFloat& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::multAssign(second);
   return result;
}

inline OriginFloat operator*(OriginFloat&& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::multAssign(second);
   return result;
}

inline OriginFloat& OriginFloat::operator/=(OriginFloat&& source)
{  inherited::divAssign(source);
   return *this;
}

inline OriginFloat operator/(const OriginFloat& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::divAssign(second);
   return result;
}

inline OriginFloat operator/(OriginFloat&& first, OriginFloat&& second)
{  OriginFloat result(first);
   result.inherited::divAssign(second);
   return result;
}

class OriginDouble : public double_st {
  private:
   typedef double_st inherited;

  public:
   OriginDouble() {}
   OriginDouble(double val, bool doesTrack=false) : inherited(val, doesTrack) {}
   OriginDouble(const char* val, bool doesTrack=false) : inherited(::atof(val), doesTrack) {}
   OriginDouble(int val, bool doesTrack=false) : inherited(val, doesTrack) {}
   OriginDouble(long long int val, bool doesTrack=false) : inherited(val, doesTrack) {}
   OriginDouble(const OriginDouble&) = default;
   OriginDouble(OriginDouble&& source) = default;
   OriginDouble(const OriginFloat& source)
      {  inherited::assign(source); }
   OriginDouble(OriginFloat&& source)
      {  inherited::assign(std::move(source)); }
   OriginDouble(double min, double max);
   OriginDouble(double min, double max, double err, unsigned instructionId, unsigned typeId);

   OriginDouble& operator=(double val) { return (OriginDouble&) inherited::operator=(OriginDouble(val)); }
   OriginDouble& operator=(int val) { return (OriginDouble&) inherited::operator=(OriginDouble(val)); }
   OriginDouble& operator=(long long int val) { return (OriginDouble&) inherited::operator=(OriginDouble(val)); }
   OriginDouble& operator=(const OriginDouble&) = default;
   OriginDouble& operator=(OriginDouble&& source) = default;
   OriginDouble& operator=(OriginFloat&& source);

   void print(std::ostream& adiagnosis_file, const char* varname) const
      {  adiagnosis_file << varname << ": ";
         inherited::print(adiagnosis_file);
         adiagnosis_file << '\n';
      }
   void retrieveRelativeError(double& error) const;
   void retrieveBoundsAndAbsoluteError(double& min, double& max, double& error) const;

   operator long long int() const { return inherited::value; }
   operator int() const { return inherited::value; }

   friend OriginDouble operator-(OriginDouble&& source);
   OriginDouble& operator+=(OriginDouble&& source);
   OriginDouble& operator+=(OriginFloat&& source);
   friend OriginDouble operator+(const OriginDouble& first, OriginDouble&& second);
   friend OriginDouble operator+(OriginDouble&& first, OriginDouble&& second);
   friend OriginDouble operator+(const OriginDouble& first, OriginFloat&& second);
   friend OriginDouble operator+(OriginDouble&& first, OriginFloat&& second);
   OriginDouble& operator-=(OriginDouble&& source);
   OriginDouble& operator-=(OriginFloat&& source);
   friend OriginDouble operator-(const OriginDouble& first, OriginDouble&& second);
   friend OriginDouble operator-(OriginDouble&& first, OriginDouble&& second);
   friend OriginDouble operator-(const OriginDouble& first, OriginFloat&& second);
   friend OriginDouble operator-(OriginDouble&& first, OriginFloat&& second);
   OriginDouble& operator*=(OriginDouble&& source);
   OriginDouble& operator*=(OriginFloat&& source);
   friend OriginDouble operator*(const OriginDouble& first, OriginDouble&& second);
   friend OriginDouble operator*(OriginDouble&& first, OriginDouble&& second);
   friend OriginDouble operator*(const OriginDouble& first, OriginFloat&& second);
   friend OriginDouble operator*(OriginDouble&& first, OriginFloat&& second);
   OriginDouble& operator/=(OriginDouble&& source);
   OriginDouble& operator/=(OriginFloat&& source);
   friend OriginDouble operator/(const OriginDouble& first, OriginDouble&& second);
   friend OriginDouble operator/(OriginDouble&& first, OriginDouble&& second);
   friend OriginDouble operator/(const OriginDouble& first, OriginFloat&& second);
   friend OriginDouble operator/(OriginDouble&& first, OriginFloat&& second);

   friend bool operator<(const OriginDouble& first, double second)
      {  return first.value < second; }
   friend bool operator<=(const OriginDouble& first, double second)
      {  return first.value <= second; }
   friend bool operator==(const OriginDouble& first, double second)
      {  return first.value == second; }
   friend bool operator!=(const OriginDouble& first, double second)
      {  return first.value != second; }
   friend bool operator>(const OriginDouble& first, double second)
      {  return first.value > second; }
   friend bool operator>=(const OriginDouble& first, double second)
      {  return first.value >= second; }

   friend bool operator<(const OriginDouble& first, const OriginDouble& second)
      {  return first.value < second.value; }
   friend bool operator<=(const OriginDouble& first, const OriginDouble& second)
      {  return first.value <= second.value; }
   friend bool operator==(const OriginDouble& first, const OriginDouble& second)
      {  return first.value == second.value; }
   friend bool operator!=(const OriginDouble& first, const OriginDouble& second)
      {  return first.value != second.value; }
   friend bool operator>(const OriginDouble& first, const OriginDouble& second)
      {  return first.value > second.value; }
   friend bool operator>=(const OriginDouble& first, const OriginDouble& second)
      {  return first.value >= second.value; }
   friend bool operator<(const OriginDouble& first, const OriginFloat& second)
      {  return first.value < second.value; }
   friend bool operator<=(const OriginDouble& first, const OriginFloat& second)
      {  return first.value <= second.value; }
   friend bool operator==(const OriginDouble& first, const OriginFloat& second)
      {  return first.value == second.value; }
   friend bool operator!=(const OriginDouble& first, const OriginFloat& second)
      {  return first.value != second.value; }
   friend bool operator>(const OriginDouble& first, const OriginFloat& second)
      {  return first.value > second.value; }
   friend bool operator>=(const OriginDouble& first, const OriginFloat& second)
      {  return first.value >= second.value; }
};

inline bool operator<(const OriginFloat& first, const OriginDouble& second)
   {  return first.value < second.value; }
inline bool operator<=(const OriginFloat& first, const OriginDouble& second)
   {  return first.value <= second.value; }
inline bool operator==(const OriginFloat& first, const OriginDouble& second)
   {  return first.value == second.value; }
inline bool operator!=(const OriginFloat& first, const OriginDouble& second)
   {  return first.value != second.value; }
inline bool operator>(const OriginFloat& first, const OriginDouble& second)
   {  return first.value > second.value; }
inline bool operator>=(const OriginFloat& first, const OriginDouble& second)
   {  return first.value >= second.value; }

inline
OriginFloat::OriginFloat(const OriginDouble& source)
   {  inherited::assign(source); }

inline
OriginFloat::OriginFloat(OriginDouble&& source)
   {  inherited::assign(std::move(source)); }

inline OriginDouble fabs(const OriginDouble& source)
{  OriginDouble result(source);
   result.fabsAssign();
   return result;
}

inline OriginDouble log(const OriginDouble& source)
{  OriginDouble result(source);
   result.logAssign();
   return result;
}

inline OriginDouble exp(const OriginDouble& source)
{  OriginDouble result(source);
   result.expAssign();
   return result;
}

inline OriginDouble sqrt(const OriginDouble& source)
{  OriginDouble result(source);
   result.sqrtAssign();
   return result;
}

inline OriginDouble sin(const OriginDouble& source)
{  OriginDouble result(source);
   result.sinAssign();
   return result;
}

inline OriginDouble cos(const OriginDouble& source)
{  OriginDouble result(source);
   result.cosAssign();
   return result;
}

inline OriginDouble tan(const OriginDouble& source)
{  OriginDouble result(source);
   result.tanAssign();
   return result;
}

inline OriginDouble asin(const OriginDouble& source)
{  OriginDouble result(source);
   result.asinAssign();
   return result;
}

inline OriginDouble acos(const OriginDouble& source)
{  OriginDouble result(source);
   result.acosAssign();
   return result;
}

inline OriginDouble atan(const OriginDouble& source)
{  OriginDouble result(source);
   result.atanAssign();
   return result;
}

inline OriginDouble atan2(const OriginDouble& first, const OriginDouble& second)
{  OriginDouble result(first);
   result.atan2Assign(second);
   return result;
}

inline OriginDouble pow(const OriginDouble& first, const OriginDouble& second)
{  OriginDouble result(first);
   result.powAssign(second);
   return result;
}

inline OriginDouble max(const OriginDouble& first, const OriginDouble& second)
{  OriginDouble result(first);
   result.maxAssign(second);
   return result;
}

inline OriginDouble min(const OriginDouble& first, const OriginDouble& second)
{  OriginDouble result(first);
   result.minAssign(second);
   return result;
}

inline OriginDouble abs(const OriginDouble& source)
{  OriginDouble result(source);
   result.absAssign();
   return result;
}

inline OriginDouble floor(const OriginDouble& source)
{  OriginDouble result(source);
   result.floorAssign();
   return result;
}

inline OriginDouble fmod(const OriginDouble& first, const OriginDouble& second)
{  OriginDouble result(first);
   result.fmodAssign(second);
   return result;
}

inline OriginDouble modf(const OriginDouble& first, OriginDouble* second)
{  OriginDouble result(first);
   result.modfAssign(*second);
   return result;
}

inline bool finite(const OriginDouble& source)
{  return finite(source.value); }

inline bool isfinite(const OriginDouble& source)
{  return std::isfinite(source.value); }

inline bool isnan(const OriginDouble& source)
{  return std::isnan(source.value); }

inline bool isinf(const OriginDouble& source)
{  return std::isinf(source.value); }

inline OriginFloat&
OriginFloat::operator=(OriginDouble&& source)
{  inherited::assign(source);
   return *this;
}

inline OriginDouble&
OriginDouble::operator=(OriginFloat&& source)
{  inherited::assign(source);
   return *this;
}

inline OriginDouble operator-(OriginDouble&& source)
{  OriginDouble result(source);
   result.inherited::oppositeAssign();
   return result;
}

inline OriginDouble& OriginDouble::operator+=(OriginDouble&& source)
{  inherited::addAssign(source, true /* is_plus */);
   return *this;
}

inline OriginFloat& OriginFloat::operator+=(OriginDouble&& source)
{  OriginFloat sourceConvert(source);
   inherited::addAssign(sourceConvert, true /* is_plus */);
   return *this;
}

inline OriginDouble& OriginDouble::operator+=(OriginFloat&& source)
{  OriginDouble sourceConvert(source);
   inherited::addAssign(sourceConvert, true /* is_plus */);
   return *this;
}

inline OriginDouble operator+(const OriginDouble& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline OriginDouble operator+(OriginDouble&& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::addAssign(second, true /* is_plus */);
   return result;
}

inline OriginDouble operator+(const OriginDouble& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::addAssign(secondConvert, true /* is_plus */);
   return result;
}

inline OriginDouble operator+(OriginDouble&& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::addAssign(secondConvert, true /* is_plus */);
   return result;
}

inline OriginDouble operator+(const OriginFloat& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.addAssign(second, true /* is_plus */);
   return result;
}

inline OriginDouble operator+(OriginFloat&& first, OriginDouble&& second)
{  OriginDouble result(std::move(first));
   result.addAssign(second, true /* is_plus */);
   return result;
}

inline OriginDouble& OriginDouble::operator-=(OriginDouble&& source)
{  inherited::addAssign(source, false /* is_plus */);
   return *this;
}

inline OriginFloat& OriginFloat::operator-=(OriginDouble&& source)
{  OriginFloat sourceConvert(source);
   inherited::addAssign(sourceConvert, false /* is_plus */);
   return *this;
}

inline OriginDouble& OriginDouble::operator-=(OriginFloat&& source)
{  OriginDouble sourceConvert(source);
   inherited::addAssign(sourceConvert, false /* is_plus */);
   return *this;
}

inline OriginDouble operator-(const OriginDouble& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline OriginDouble operator-(OriginDouble&& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::addAssign(second, false /* is_plus */);
   return result;
}

inline OriginDouble operator-(const OriginDouble& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::addAssign(secondConvert, false /* is_plus */);
   return result;
}

inline OriginDouble operator-(OriginDouble&& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::addAssign(secondConvert, false /* is_plus */);
   return result;
}

inline OriginDouble operator-(const OriginFloat& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.addAssign(second, false /* is_plus */);
   return result;
}

inline OriginDouble operator-(OriginFloat&& first, OriginDouble&& second)
{  OriginDouble result(std::move(first));
   result.addAssign(second, false /* is_plus */);
   return result;
}

inline OriginDouble& OriginDouble::operator*=(OriginDouble&& source)
{  inherited::multAssign(source);
   return *this;
}

inline OriginFloat& OriginFloat::operator*=(OriginDouble&& source)
{  OriginFloat sourceConvert(source);
   inherited::multAssign(sourceConvert);
   return *this;
}

inline OriginDouble& OriginDouble::operator*=(OriginFloat&& source)
{  OriginDouble sourceConvert(source);
   inherited::multAssign(sourceConvert);
   return *this;
}

inline OriginDouble operator*(const OriginDouble& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::multAssign(second);
   return result;
}

inline OriginDouble operator*(OriginDouble&& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::multAssign(second);
   return result;
}

inline OriginDouble operator*(const OriginDouble& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::multAssign(secondConvert);
   return result;
}

inline OriginDouble operator*(OriginDouble&& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::multAssign(secondConvert);
   return result;
}

inline OriginDouble operator*(const OriginFloat& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.multAssign(second);
   return result;
}

inline OriginDouble operator*(OriginFloat&& first, OriginDouble&& second)
{  OriginDouble result(std::move(first));
   result.multAssign(second);
   return result;
}

inline OriginDouble& OriginDouble::operator/=(OriginDouble&& source)
{  inherited::divAssign(source);
   return *this;
}

inline OriginFloat& OriginFloat::operator/=(OriginDouble&& source)
{  OriginFloat sourceConvert(source);
   inherited::divAssign(sourceConvert);
   return *this;
}

inline OriginDouble& OriginDouble::operator/=(OriginFloat&& source)
{  OriginDouble sourceConvert(source);
   inherited::divAssign(sourceConvert);
   return *this;
}

inline OriginDouble operator/(const OriginDouble& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::divAssign(second);
   return result;
}

inline OriginDouble operator/(OriginDouble&& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.inherited::divAssign(second);
   return result;
}

inline OriginDouble operator/(const OriginDouble& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::divAssign(secondConvert);
   return result;
}

inline OriginDouble operator/(OriginDouble&& first, OriginFloat&& second)
{  OriginDouble result(first);
   OriginDouble secondConvert(second);
   result.inherited::divAssign(secondConvert);
   return result;
}

inline OriginDouble operator/(const OriginFloat& first, OriginDouble&& second)
{  OriginDouble result(first);
   result.divAssign(second);
   return result;
}

inline OriginDouble operator/(OriginFloat&& first, OriginDouble&& second)
{  OriginDouble result(std::move(first));
   result.divAssign(second);
   return result;
}

#ifdef _ORIGINS
typedef float old_float;
typedef double old_double;

#define float OriginFloat
#define double OriginDouble

#ifdef _ORIGINS_MAIN
std::atomic<uint64_t> BaseOriginVector::atomic_id_counts = {0};
#endif
#endif
