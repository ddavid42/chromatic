#include "workplan.h"

bool
WorkPlan::setToNextBrother(std::vector<QuadTreeElement*>& stack) const {
   QuadTreeElement* son;
   while (!stack.empty()) {
      son = stack.back();
      stack.pop_back();
      auto iter = apaqteRoot->begin(), iterEnd = apaqteRoot->end();
      if (!stack.empty())
         iter = stack.back()->content().begin(), iterEnd = stack.back()->content().end();
      bool hasFound = false;
      for (; !hasFound && iter != iterEnd; ++iter) {
         if (&*iter == son)
            hasFound = true;
      }
      if (!hasFound)
         return false;
      if (iter != iterEnd) {
         stack.push_back(&(*iter));
         return true;
      }
   };
   return true;
}

bool
WorkPlan::read(std::istream& in) {
   int lines, columns;
   in >> lines;
   in >> columns;
   init(apaqteRoot, 0, 0, lines, columns);
   int ch;
   while ((ch = in.get()) != EOF) {
      if (!std::isspace(ch)) {
         if (ch != '1')
            return false;
         break;
      }
   }
   std::vector<QuadTreeElement*> stack;
   stack.push_back(&((*apaqteRoot)[0]));

   while (!stack.empty() && (ch = in.get()) != EOF) {
      if (std::isspace(ch))
         continue;
      if (ch == '1') {
         init(stack.back()->subTrees().content, stack.back()->getX(), stack.back()->getY(),
               stack.back()->getWidth(), stack.back()->getLength());
         stack.push_back(&stack.back()->content()[0]);
      }
      else if (ch == '0') {
         old_double contribution, over_contribution;
         in >> contribution;
         in >> over_contribution;
         stack.back()->setContribution(contribution, over_contribution);
         if (!setToNextBrother(stack))
            return false;
         if (stack.empty())
            return true;
      }
      else if (ch == 'x') {
         if (!setToNextBrother(stack))
            return false;
         if (stack.empty())
            return true;
      }
   }
   return false;
}

void
WorkPlan::write(std::ostream& out) const {
   if (!apaqteRoot)
      return;
   int lines = 0, columns = 0;
   const QuadTreeElement& last = (*apaqteRoot).back();
   lines = last.getX() + last.getWidth();
   columns = last.getY() + last.getLength();
   out << lines << ' ' << columns << '\n';

   out << '1';
   std::vector<QuadTreeElement*> stack;
   stack.push_back(&((*apaqteRoot)[0]));

   while (!stack.empty()) {
      if (stack.back()->hasContent()) {
         out << " 1";
         stack.push_back(&stack.back()->content()[0]);
      }
      else {
         if (stack.back()->hasContribution())
            out << "0 " << stack.back()->getContribution() << ' ' << stack.back()->getOverContribution() << ' ';
         else
            out << "x";
         int len = stack.size();
         if (!setToNextBrother(stack))
            return;
         int newlen = stack.size();
         if (newlen < len) {
            out << '\n';
            for (;newlen > 1; --newlen)
               out << "  ";
            out << ' ';
         }
      };
   }
}
