#ifndef MCMD_STRUCTURE_FACTOR_P_GRID_H
#define MCMD_STRUCTURE_FACTOR_P_GRID_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactorP.h"
#include <simp/crystal/LatticeSystem.h>

namespace McMd
{

   using namespace Util;

   /**
   * StructureFactorPGrid evaluates structure factors in Fourier space.
   *
   * This class evaluates the structures factors for all wavevectors
   * for all wavevectors within a grid, within a region in which all
   * of the h, k, l integer wavevector components (i.e., Miller indices) 
   * has an absolute magnitude less than or equal to a parameter hMax.
   * 
   * The class also requires a user to specify a variable lattice of
   * enum type LatticeSystem whose value (e.g., Cubic, Orthorhombic, 
   * Triclinic) determines the symmetry group that is used to group
   * group equivalent wavevectors into "stars". Use Triclinic to 
   * impose no symmetry. 
   *
   * Example: Here is parameter file for a system with nAtomType = 2, 
   * with two monomer types 0 and 1, for calculating the total, 00,
   * and 01 correlation functions on a grid with Miller indices up
   * and including 5:
   *
   * \code
   * StructureFactorPGrid{
   *    interval                     1000
   *    outputFileName     StructureFactorPGrid
   *    nAtomTypeIdPair                  3
   *    atomTypeIdPairs          -1     -1
   *                              0      0
   *                              0      1
   *    hMax                             5
        lattice                      Cubic
   * }
   * \endcode
   * At the end of a simulation, all of the structure factors are
   * output in a file with a suffix *.dat. Each line in this file
   * contains the 3 Miller indices of a wavevector, the absolute
   * magnitude of the wavevector, and a list of nAtomTypeIdPair 
   * structure factor values for the wavevector, one for each 
   * atomTypeId pair.
   * 
   * \sa \ref mcMd_analyzer_StructureFactorPGrid_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class StructureFactorPGrid : public StructureFactorP
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      StructureFactorPGrid(System &system);

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - int               nAtomTypeIdPair number of Atom typeId pairs
      *   - DArray<Pair<int>> atomTypeIdPairs integer typeId pairs
      *   - int               hMax            maximum Miller index
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Set up before a simulation.
      */
      virtual void setup();
 
      /**
      * Output structure factors, averaged over stars.
      */
      virtual void output();
 
   private:

      /// Array of ids for first wavevector in each star.
      DArray<int>  starIds_;

      /// Array of star sizes.
      DArray<int>  starSizes_;

      /// Maximum Miller index of wavevectors in grid.
      int   hMax_;

      /// Number of stars of symmetry related wavevectors.
      int   nStar_;

      /// Lattice system used to create stars.
      LatticeSystem   lattice_;

      /// Has readParam been called?
      bool    isInitialized_;

      void makeStars();

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void StructureFactorPGrid::serialize(Archive& ar, const unsigned int version)
   {
      StructureFactorP::serialize(ar, version);
      ar & hMax_;
      ar & lattice_;
      //ar & nStar_;
      //ar & starIds_;
      //ar & starSizes_;
   }

}
#endif
