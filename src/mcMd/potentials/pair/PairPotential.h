#ifndef MCMD_PAIR_POTENTIAL_H
#define MCMD_PAIR_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <string>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   using namespace Util;

   /**
   * Interface for a Pair Potential.
   *
   * \ingroup McMd_Pair_Module
   */
   class PairPotential
   {

   public:

      /**
      * Destructor (does nothing)
      */
      virtual ~PairPotential()
      {}

      /// \name Pair Interaction Interface
      //@{ 

      /**
      * Return pair energy for a single pair.
      */
      virtual 
      double energy(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return force / separation for a single pair.
      */
      virtual 
      double forceOverR(double rsq, int iAtomType, int jAtomType) const = 0;

      /**
      * Return maximum cutoff distance.
      */
      virtual double maxPairCutoff() const = 0;

      /**
      * Return name of pair interaction class (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const = 0;

      //@}
      /// \name Global Energy, Force, and Stress Evaluators
      //@{ 

      /**
      * Calculate the total nonBonded pair energy for the associated System.
      */
      virtual double energy() = 0;

      /**
      * Compute total nonbonded pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const = 0;

      /**
      * Compute x, y, z nonbonded pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const = 0;

      /**
      * Compute stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const = 0;

      //@}
   };

} 
#endif
