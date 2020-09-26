#ifndef HOT_BEAD_CLUSTERING_CPP
#define HOT_BEAD_CLUSTERING_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "hotBeadClustering.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>


#include <simp/boundary/Boundary.h>
#include <simp/species/Species.h>

#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <util/space/Tensor.h>
#include <util/containers/DArray.h>
#include <sstream>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

    /// Constructor.
   hotBeadClustering::hotBeadClustering(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      hist_(),
      clusHist_(),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      cutoffPurity_(),
      clusHistMin_(),
      clusHistMax_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("hotBeadClustering"); }

   /// Read parameters from file, and allocate arrays.
   void hotBeadClustering::readParameters(std::istream& in) 
   {  
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      read<int>(in, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      if (atomTypeId_ >= system().simulation().nAtomType()) {
         UTIL_THROW("nTypeId >= nAtomType");
      }

      read<double>(in, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      read<double>(in, "cutoff_purity", cutoffPurity_);
      if (cutoffPurity_ < 0) {
         UTIL_THROW("Negative cutoff");
      }  

      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);

      read<int>(in,"histMin", clusHistMin_);
      read<int>(in,"histMax", clusHistMax_);
      clusHist_.setParam(clusHistMin_, clusHistMax_);
      clusHist_.clear();

      // Generate cell list as well

      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);

      std::cout <<atomCapacity;
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void hotBeadClustering::loadParameters(Serializable::IArchive& ar)
   {
      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      loadParameter<int>(ar, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }

      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      loadParameter<double>(ar, "purity cutoff", cutoffPurity_);
      if (cutoffPurity_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      ar >> hist_;

      loadParameter<int>(ar, "histMin", clusHistMin_);
      loadParameter<int>(ar, "histMax", clusHistMax_);
      ar >> clusHist_;

      ar >> nSample_;

      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);

      // Generate cell list as well
    
      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);
  

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void hotBeadClustering::save(Serializable::OArchive& ar)
   {
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & cutoffPurity_;
      ar & hist_;
      ar & clusHistMin_;
      ar & clusHistMax_;
      ar & clusHist_;
      ar & nSample_;
   }

   /*
   * Clear accumulators.
   */
   void hotBeadClustering::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");

      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      system().begin(speciesId_, molIter);
      int n = 0;


      for ( ; molIter.notEnd(); ++molIter) {
      std::cout <<"Yahin hu main";
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            if (atomIter->typeId() == atomTypeId_) {
               n++;
               std::cout <<n<<"\n";
            }
         }
          break;
      }

      countType_ = n * system().nMolecule(speciesId_);
      atomEnv_.allocate(countType_);
      std::cout <<"countType_";

      isHot_.allocate(system().nAtom());
      for (int iAtom = 0; iAtom < system().nAtom(); iAtom++) {
         isHot_ [iAtom] = 0;
      }
      // min value = 0, max value = 1 and number of bins = 100
      // Can think of reading number of bins from parameter file 
      // as well
      hist_.setParam(-0.02, 1.02, 104); 
      hist_.clear();  

      clusHist_.clear();
      nSample_ = 0;   
   }

   /* 
   * Sample data by calling ClusterIdentifier::identifyClusters.
   */
   void hotBeadClustering::sample(long iStep) 
   { 
      if (isAtInterval(iStep)) {

         // Initialize all data structures:
         // Setup a grid of empty cells
         cellList_.setup(system().boundary(), (cutoff_ + 0.4));

         // Set variables to initial state
         for (int i = 0; i < countType_; ++i) {
            atomEnv_[i].clear();
         }

         CellList::NeighborArray neighborArray;
         Boundary& boundary = system().boundary();
         System::MoleculeIterator molIter;
         Molecule::AtomIterator atomIter;
         Atom* otherAtomPtr;
         double cutoffSq = cutoff_*cutoff_;
         double rsq;   
         int thisMolId, otherMolId;
         int neighborCount = 0;
         int selectNeighborCount = 0; 

         // Build the cellList, associate Atom with AtomDomain.
         // Iterate over molecules of species speciesId_
     
         for (int iSpecies = 0; iSpecies < system().simulation().nSpecies(); ++iSpecies) {
            system().begin(iSpecies, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  // Build the cellList     
                  system().boundary().shift(atomIter->position());
                  cellList_.addAtom(*atomIter);

               }
            }
         }   
         // Check if this is how the cell List is being generated for the MD simulation,
         // because the method you use should be the same, only the cutoff should be 
         // different
            

         int iAtom = 0;
         
         system().begin(speciesId_, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
 
            thisMolId = molIter->id();
 
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {

               if (atomIter->typeId() == atomTypeId_) {

                  // Associate this Atom with AtomDomain
                  atomEnv_[iAtom].setAtom(*atomIter);

                  cellList_.getNeighbors(atomIter->position(), neighborArray);

                     for (int i = 0; i < neighborArray.size(); i++) {
                        otherAtomPtr = neighborArray[i];
                        otherMolId = otherAtomPtr->molecule().id();
                        if (otherMolId != thisMolId) {
                           rsq = boundary.distanceSq(atomIter->position(),
                                                     otherAtomPtr->position());
                           if (rsq < cutoffSq) {
                              neighborCount++;  
                              if (otherAtomPtr->typeId() == atomTypeId_){
                                 selectNeighborCount++;
                              }
                 
                           }
                        }
                     }

                  atomEnv_[iAtom].setNeighborCount(neighborCount, selectNeighborCount);

     


                  iAtom++;

                  // Re-initializing neighborCount and selectNeighborCount
                  neighborCount = 0;
                  selectNeighborCount = 0;   

               }
            }   
         }
    
         for (iAtom = 0; iAtom < countType_; iAtom++) {
            if (atomEnv_[iAtom].domainPurity() > cutoffPurity_) {            
               isHot_[(atomEnv_[iAtom].atom()).id()] = 1;
            }
         }

         identifier_.identifyClusters(isHot_);


         for (int i = 0; i < identifier_.nCluster(); i++) {
             clusHist_.sample(identifier_.cluster(i).size());
         }
         ++nSample_;
         fileMaster().openOutputFile(outputFileName(".clusters"+toString(iStep)),outputFile_);
         //Writes all of the clusters and their component molecules
         Cluster thisCluster;
         ClusterLink* thisClusterStart;
         ClusterLink* next;
         Molecule thisMolecule;
         //Loop over each cluster
         for (int i = 0; i < identifier_.nCluster(); i++) {
             thisCluster = identifier_.cluster(i);
             thisClusterStart = thisCluster.head();
             outputFile_ << i<<"  ("<<identifier_.cluster(i).size()<<")  "<< "       " ;
             //List out every molecule in that cluster
             while (thisClusterStart) {
                next = thisClusterStart->next();
                thisMolecule = thisClusterStart->molecule();
                outputFile_ << thisMolecule.id() << "  ";
                thisClusterStart = next;
             }   
             outputFile_ << "\n";
         }   
         outputFile_.close();


         fileMaster().openOutputFile(outputFileName(".COMs"+toString(iStep)),outputFile_);
         //comArray;
         Vector clusterCOM;
         Vector r0; 
         Vector dr; 
         Tensor moment;
         Tensor rgDyad;
         DArray<Vector> allCOMs;
         DArray<Tensor> allMoments;
         allCOMs.allocate(identifier_.nCluster());
         allMoments.allocate(identifier_.nCluster());
         for (int i = 0; i < identifier_.nCluster(); i++) {
             thisCluster = identifier_.cluster(i);
             outputFile_ << i<<"  ("<<identifier_.cluster(i).size()<<")  "<< "       " ;
             //For that cluster, calculate the center of mass
             clusterCOM = thisCluster.clusterCOM(atomTypeId_, system().boundary());
             outputFile_ << clusterCOM;
             outputFile_ << "\n";
             allCOMs[i] = clusterCOM;
             //Calculate Rg
             moment = thisCluster.momentTensor(atomTypeId_, system().boundary());
             allMoments[i] = moment;
         }
         outputFile_.close();
         fileMaster().openOutputFile(outputFileName(".momentTensors"+toString(iStep)),outputFile_);
         for (int i = 0; i < identifier_.nCluster(); i++) {
             outputFile_ << i<<"  ("<<identifier_.cluster(i).size()<<")  "<< "       "<< allMoments[i] << "\n";
         }
         outputFile_.close();



        // fileMaster().openOutputFile(outputFileName(".env"+toString(iStep)),outputFile_);
         //Writes all of the clusters and their component molecules
        // for (iAtom = 0; iAtom < countType_; iAtom++) {
        //     outputFile_ << atomEnv_[iAtom].atom().id() << "	" << atomEnv_[iAtom].domainPurity()<< "    "  << int( (atomEnv_[iAtom].domainPurity() + 0.02)/0.01 );
        //     outputFile_ << "\n";
        //     hist_.sample(atomEnv_[iAtom].domainPurity());
        // }
         //outputFile_.close();

      for (int iAtom = 0; iAtom < system().nAtom(); iAtom++) {
         isHot_ [iAtom] = 0;
      } 

      }
      // Reinitialize isHot to 0
   }

   /*
   * Output results to file after simulation is completed.
   */
   void hotBeadClustering::output() 
   {
      // Write parameter file
       fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
       writeParam(outputFile_);
       outputFile_.close();

      // Write histogram output
       fileMaster().openOutputFile(outputFileName(".domain"), outputFile_);
       hist_.output(outputFile_);
      // double min = hist_.min();
      // double binWidth = hist_.binWidth();
      // outputFile_ << Dbl(min);
      // int nBin = hist_.nBin();
      // for (int i = 0; i < nBin; ++i) {
      //    outputFile_ <<  Dbl((i * binWidth) + min) << "  " 
      //                <<  Dbl(hist_.data()[i]/nSample_) << "\n";
      // }
      fileMaster().openOutputFile(outputFileName(".hist"), outputFile_);
      int min = clusHist_.min();
      int nBin = clusHist_.nBin();
      for (int i = 0; i < nBin; ++i) {
         outputFile_ << Int(i + min) << "  "
                     <<  Dbl(double(clusHist_.data()[i])/double(nSample_)) << "\n";
      }
      outputFile_.close();
   }

}
#endif 
