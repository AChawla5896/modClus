#ifndef MCMD_ENVIRONMENT_H
#define MCMD_ENVIRONMENT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>           // base class templ
#include <mcMd/simulation/System.h>                  // class templ param
#include <mcMd/user/AtomDomain.h> // member
#include <mcMd/neighbor/CellList.h>              // member
#include <util/containers/DArray.h>              // member template
#include <util/accumulators/Distribution.h>       // member

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Identify micelle clusters in polymeric systems.
   */
   class Environment : public SystemAnalyzer<System>
   {

   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      Environment(System &system);

      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval         : sampling interval
      *   - string outputFileName   : base name for output file(s)
      *   - int    speciesId        : integer id for Species of interest
      *   - int    atomTypeId       : integer id for core atom type
      *   - double cutoff           : distance cutoff
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Clear accumulator.
      */
      virtual void setup();

      /**
      * Identify clusters in configuration.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);

      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive.
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Clusters species type
      int  speciesId_;

      /// Clusters species type
      int  atomTypeId_;

      /// Histogram distribution of domain purity
      Distribution hist_;

      /// Total number of atoms of the selected type on the selected species

      int countType_;

      /// Distance cutoff
      double cutoff_;

      DArray<AtomDomain> atomEnv_;
 
      /// Has readParam been called?
      bool  isInitialized_;

      CellList cellList_;

   };

   /**
   * Serialize to/from an archive.
   */
   template <class Archive>
   void Environment::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & hist_;
   }

}
#endif
