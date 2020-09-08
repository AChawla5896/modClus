#ifndef HOT_BEAD_CLUSTERING_H
#define HOT_BEAD_CLUSTERING_H

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
#include <mcMd/user/hotClusterIdentifier.h> // member
#include <util/containers/DArray.h>              // member template
#include <util/accumulators/Distribution.h>       // member
#include <util/accumulators/IntDistribution.h>       // member

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Identify micelle clusters in polymeric systems.
   */
   class hotBeadClustering : public SystemAnalyzer<System>
   {

   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      hotBeadClustering(System &system);

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
      *   - double purity cutoff    : cutoff for purity
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

      hotClusterIdentifier identifier_;

      /// Output file stream
      std::ofstream outputFile_;

      /// Clusters species type
      int  speciesId_;

      /// Clusters species type
      int  atomTypeId_;

      /// Histogram distribution of domain purity
      Distribution hist_;

      IntDistribution  clusHist_;

      /// Total number of atoms of the selected type on the selected species

      int countType_;

      /// Histogram minimum value
      int  clusHistMin_;

      /// Histogram maximum value.
      int  clusHistMax_;

      /// Distance cutoff
      double cutoff_;

      /// Purity cutoff
      double cutoffPurity_;

      DArray<AtomDomain> atomEnv_;

      // This is a boolean array which will tell if the bead is hot (1) or
      // not (0)
      DArray<bool > isHot_;
 
      /// Has readParam been called?
      bool  isInitialized_;

      /// Number of configurations dumped thus far (first dump is zero).
      long  nSample_;

      CellList cellList_;

   };

   /**
   * Serialize to/from an archive.
   */
   template <class Archive>
   void hotBeadClustering::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & hist_;
   }

}
#endif
