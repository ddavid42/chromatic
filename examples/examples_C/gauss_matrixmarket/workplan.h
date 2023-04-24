#pragma once
#include <iostream>
#include <memory>
#include <vector>

#ifndef double
typedef double old_double;
#endif

static const int SUBDIVISION = 1;
static const int EXP_SUBDIVISION = 1 << (SUBDIVISION+1);
static const int SQUARE_EXP_SUBDIVISION = EXP_SUBDIVISION*EXP_SUBDIVISION;

class QuadTreeElement;
struct QuadTreeSubElements {
   std::unique_ptr<std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>> content;

   QuadTreeSubElements() = default;
   QuadTreeSubElements(const QuadTreeSubElements& source);
   QuadTreeSubElements(QuadTreeSubElements&& source) = default;
   QuadTreeSubElements& operator=(const QuadTreeSubElements& source);
   QuadTreeSubElements& operator=(QuadTreeSubElements&& source) = default;
};

class QuadTreeElement {
  private:
   int uPosX=0, uPosY=0;
   int uWidth=0, uLength=0;
   QuadTreeSubElements apaqteSubTrees;
   bool fHasContributions = false;
   old_double dContributions = 0.0;
   old_double dOverContributions = 0.0;
#ifdef _ORIGINS_ERROR
   old_double dErrorContributions = 0.0;
   old_double dErrorOverContributions = 0.0;
#endif
   double dSharedReference = 0.0;

  public:
   QuadTreeElement() = default;
   QuadTreeElement(const QuadTreeElement& source) = default;
   QuadTreeElement(QuadTreeElement&& source) = default;
   QuadTreeElement& operator=(const QuadTreeElement& source) = default;
   QuadTreeElement& operator=(QuadTreeElement&& source) = default;

   QuadTreeSubElements& subTrees() { return apaqteSubTrees; }
   bool hasContent() const { return apaqteSubTrees.content.get(); }
   bool hasContribution() const { return fHasContributions; }
   const old_double& getContribution() const { return dContributions; }
   old_double& getSContribution() { return dContributions; }
   const old_double& getOverContribution() const { return dOverContributions; }
   old_double& getSOverContribution() { return dOverContributions; }
#ifdef _ORIGINS_ERROR
   const old_double& getErrorContribution() const { return dErrorContributions; }
   old_double& getSErrorContribution() { return dErrorContributions; }
   const old_double& getErrorOverContribution() const { return dErrorOverContributions; }
   old_double& getSErrorOverContribution() { return dErrorOverContributions; }
#endif
   double& sharedReference() { return dSharedReference; }
   std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>& content() { return *apaqteSubTrees.content; }
   QuadTreeElement& setPosition(int posX, int posY) { uPosX = posX; uPosY = posY; return *this; }
   QuadTreeElement& setWidth(int width, int length) { uWidth = width; uLength = length; return *this; }
   QuadTreeElement& setContribution(bool hasContribution, old_double contribution, old_double over_contribution) { fHasContributions = hasContribution; dContributions = contribution; dOverContributions = over_contribution; return *this; }
   QuadTreeElement& setErrorContribution(old_double error)
      {  dErrorContributions = error*dContributions; dErrorOverContributions = error*dOverContributions; return *this; }
   QuadTreeElement& setErrorContribution(old_double error_contribution, old_double error_over_contribution)
      {  dErrorContributions = error_contribution; dErrorOverContributions = error_over_contribution; return *this; }
   QuadTreeElement& setContribution(old_double contribution, old_double over_contribution)
      {  fHasContributions = true; dContributions = contribution; dOverContributions = over_contribution; return *this; }
   QuadTreeElement& setContribution()
      {  fHasContributions = true; return *this; }

   bool isValid() const { return uPosX >= 0 && uPosY >= 0; }
   int getX() const { return uPosX; }
   int getY() const { return uPosY; }
   int getWidth() const { return uWidth; }
   int getLength() const { return uLength; }
};

inline
QuadTreeSubElements::QuadTreeSubElements(const QuadTreeSubElements& source)
   {  if (source.content)
         content.reset(new std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>(*source.content));
   }

inline QuadTreeSubElements&
QuadTreeSubElements::operator=(const QuadTreeSubElements& source)
   {  if (source.content)
         content.reset(new std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>(*source.content));
      else
         content.reset();
      return *this;
   }

class WorkPlan {
  public:
   class Cursor;

  private:
   std::unique_ptr<std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>> apaqteRoot;

   bool setToNextBrother(std::vector<QuadTreeElement*>& stack) const;
   friend class Cursor;

  protected:
   void init(std::unique_ptr<std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>>& root,
         int shiftX, int shiftY, int lines, int columns)
      {  int width = (lines+EXP_SUBDIVISION-1) / EXP_SUBDIVISION;
         int length = (columns+EXP_SUBDIVISION-1) / EXP_SUBDIVISION;
         // if (width < 1) width = 1;
         // if (length < 1) length = 1;
         root.reset(new std::array<QuadTreeElement, SQUARE_EXP_SUBDIVISION>());
         for (int x = 0; x < EXP_SUBDIVISION; ++x) {
            for (int y = 0; y < EXP_SUBDIVISION; ++y) {
               if (x*width < lines && y*length < columns)
                  (*root)[x*EXP_SUBDIVISION + y].setPosition(shiftX + x*width, shiftY + y*length).
                     setWidth((x+1)*width <= lines ? width : (lines - x*width),
                           (y+1)*length <= columns ? length : (columns - y*length));
               else
                  (*root)[x*EXP_SUBDIVISION + y].setPosition(-1, -1);
            }
         }
      }

  public:
   WorkPlan() = default;
   WorkPlan(WorkPlan&& source) = default;
   WorkPlan& operator=(WorkPlan&& source) = default;

   class Cursor : public std::vector<QuadTreeElement*> {
     private:
      WorkPlan* wpReference;

     public:
      Cursor(WorkPlan& source) : wpReference(&source) {}
      Cursor(const WorkPlan& source) : wpReference(&const_cast<WorkPlan&>(source)) {}
      Cursor(const Cursor& source) = default;
      Cursor(Cursor&& source) = default;
      Cursor& operator=(const Cursor& source) = default;
      Cursor& operator=(Cursor&& source) = default;

      bool isValid() const { return !empty(); }
      QuadTreeElement* elementAt() const { return empty() ? nullptr : back(); }
      bool setToFirst()
         {  clear();
            bool result = wpReference->apaqteRoot.get();
            if (result) push_back(&(*wpReference->apaqteRoot)[0]);
            return result;
         }
      bool setToNext()
         {  if (!empty()) {
               if (back()->hasContent()) {
                  push_back(&back()->content()[0]);
                  return true;
               }
            }
            bool result = wpReference->setToNextBrother(*this);
            if (!result)
               return false;
            return !empty();
         }
      bool setToPosition(int x, int y)
         {  if (!empty()) {
               QuadTreeElement* current = back();
               while (current->hasContent()) {
                  if (!setToNext()) return false;
                  current = back();
               }
               if (current->getX() <= x && current->getY() <= y) {
                  if (x < current->getX() + current->getWidth() &&
                        y < current->getY() + current->getLength())
                     return true;
                  while (setToNext()) {
                     current = back();
                     while (current->hasContent()) {
                        if (!setToNext()) return false;
                        current = back();
                     }
                     if (current->getX() <= x && x < current->getX() + current->getWidth() &&
                           current->getY() <= y && y < current->getY() + current->getLength())
                        return true;
                  }
                  return false;
               }
            }
            if (setToFirst()) {
               do {
                  QuadTreeElement* current = back();
                  while (current->hasContent()) {
                     if (!setToNext()) return false;
                     current = back();
                  }
                  if (current->getX() <= x && x < current->getX() + current->getWidth() &&
                        current->getY() <= y && y < current->getY() + current->getLength())
                     return true;
               } while (setToNext());
            }
            return false;
         }
   };

   int getWidth() const { return apaqteRoot->back().getX() + apaqteRoot->back().getWidth(); }
   int getLength() const { return apaqteRoot->back().getY() + apaqteRoot->back().getLength(); }
   int getCount() const
      {  Cursor cursor(*this);
         int result = 0;
         if (cursor.setToFirst()) {
            do {
               QuadTreeElement* current = cursor.back();
               while (current->hasContent()) {
                  if (!cursor.setToNext()) return false;
                  current = cursor.back();
               }
               ++result;
            } while (cursor.setToNext());
         }
         return result;
      }
   int getNonZeroCount() const
      {  Cursor cursor(*this);
         int result = 0;
         if (cursor.setToFirst()) {
            do {
               QuadTreeElement* current = cursor.back();
               while (current->hasContent()) {
                  if (!cursor.setToNext()) return false;
                  current = cursor.back();
               }
               if (current->getContribution() != 0.0 || current->getOverContribution() != 0.0)
                  ++result;
            } while (cursor.setToNext());
         }
         return result;
      }
   
   bool isValid() const { return apaqteRoot.get(); }
   bool read(std::istream& in);
   void write(std::ostream& out) const;
};

