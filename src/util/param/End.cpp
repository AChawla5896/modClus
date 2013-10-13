#ifndef UTIL_END_CPP
#define UTIL_END_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "End.h"
#include <util/global.h>

namespace Util
{

   /* 
   * Constructor.
   */
   End::End()
    : label_("}")
   {
   }

   /* 
   * Destructor.
   */
   End::~End()
   {}

   /* 
   * Read and check end bracket.
   */
   void End::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         // std::string expected;
         // expected = "}";
         // in >> Label(expected.c_str());
         in >> label_;
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
   }

   /* 
   * End::writeParam()
   */
   void End::writeParam(std::ostream &out)
   {
      out << indent() << "}" << std::endl;
   }

   /* 
   * Empty implementation of virtual resetParam() method.
   */
   void End::resetParam()
   {}

} 
#endif
