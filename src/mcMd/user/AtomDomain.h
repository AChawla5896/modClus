#ifndef MCMD_ATOM_DOMAIN_H
#define MCMD_ATOM_DOMAIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd
{

   class Atom;

   /**
   * Molecule in a cluster.
   *
   * A cluster is defined by a linked list of ClusterLink objects.
   */
   struct AtomDomain
   {

   public:

      /**
      * Set to default state.
      */
      void clear();

      /**
      * Set pointer to associated molecule.
      */
      void setAtom(Atom& atom);

      void setNeighborCount(double totalNeighbors, double selectNeighbors); 

      /**
      * Get associated molecule by reference.
      */
      Atom& atom() const;

      /**
      * Get pointer to next link in linked list.
      */
      int totalNeighbors() const;

      /** 
      * Get pointer to next link in linked list.
      */
      int selectNeighbors() const;

      /**
      * Get cluster identifier.
      */
      double domainPurity() const;

   private:

      /// Pointer to the associated molecule.
      Atom* atomPtr_;

      /// Pointer to the next link.
      double totalNeighbors_;

      /// Integer id of the associated cluster.
      double selectNeighbors_;

      double domainPurity_;

   };

   // Set ClusterLink to default null state.
   inline 
   void AtomDomain::clear()
   {
      atomPtr_ = 0;  
      totalNeighbors_ = -1.0;
      selectNeighbors_ = -1.0;
      domainPurity_ = -1.0;
   }

   // Set pointer to associated Molecule.
   inline
   void AtomDomain::setAtom(Atom& atom)
   {  atomPtr_ = &atom;  }

   inline
   Atom& AtomDomain::atom() const
   {  return *atomPtr_; }

   inline
   void AtomDomain::setNeighborCount(double totalNeighbors, double selectNeighbors)
   {  
      totalNeighbors_ = totalNeighbors;
      selectNeighbors_ = selectNeighbors;
      if (totalNeighbors_ == 0 && selectNeighbors_ == 0){
         //domainPurity_ = 0;
      }
      else {
         // domainPurity_ = ((double) selectNeighbors_ / (double) totalNeighbors_) * 100.0;
         domainPurity_ = ( selectNeighbors_ /  totalNeighbors_);
      }   
   }

   // Get the id of the associated cluster.
   inline
   int AtomDomain::totalNeighbors() const
   {  return totalNeighbors_; }

   inline
   int AtomDomain::selectNeighbors() const
   {  return selectNeighbors_; }

   inline
   double AtomDomain::domainPurity() const
   {  return domainPurity_; }

}
#endif
