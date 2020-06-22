#ifndef MCMD_ENVIRONMENT_CPP
#define MCMD_ENVIRONMENT_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Environment.h"
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
   Environment::Environment(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("Environment"); }

   /// Read parameters from file, and allocate arrays.
   void Environment::readParameters(std::istream& in) 
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

      // Generate cell list as well

      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);

      std::cout <<atomCapacity;
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void Environment::loadParameters(Serializable::IArchive& ar)
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

      ar >> nSample_;

      // Generate cell list as well
    
      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);
  

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void Environment::save(Serializable::OArchive& ar)
   {
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & nSample_;
   }

   /*
   * Clear accumulators.
   */
   void Environment::setup() 
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
      
      nSample_ = 0;
   }

   /* 
   * Sample data by calling ClusterIdentifier::identifyClusters.
   */
   void Environment::sample(long iStep) 
   { 
      if (isAtInterval(iStep)) {

         // Initialize all data structures:
         // Setup a grid of empty cells
         cellList_.setup(system().boundary(), cutoff_);

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



         fileMaster().openOutputFile(outputFileName(".env"+toString(iStep)),outputFile_);
         //Writes all of the clusters and their component molecules
         for (iAtom = 0; iAtom < countType_; iAtom++) {
             outputFile_ << atomEnv_[iAtom].atom().id() << "	" << atomEnv_[iAtom].domainPurity();
             outputFile_ << "\n";
         }
         outputFile_.close();

      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void Environment::output() 
   {
      // Write parameter file
      // fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      // writeParam(outputFile_);
      // outputFile_.close();

      // Write histogram output
      // fileMaster().openOutputFile(outputFileName(".hist"), outputFile_);
      // hist_.output(outputFile_);
      // int min = hist_.min();
      // int nBin = hist_.nBin();
      // for (int i = 0; i < nBin; ++i) {
      //    outputFile_ << Int(i + min) << "  " 
      //                <<  Dbl(double(hist_.data()[i])/double(nSample_)) << "\n";
      // }
      // outputFile_.close();
   }

}
#endif 
