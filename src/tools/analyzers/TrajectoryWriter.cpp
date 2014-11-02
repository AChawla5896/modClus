/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"
#include <tools/processor/Processor.h>
#include <tools/storage/AtomStorage.h>
#include <tools/storage/GroupStorage.h>
#include <util/mpi/MpiLoader.h>

#include <fstream>
#include <ios>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   TrajectoryWriter::TrajectoryWriter(Processor& processor, 
                                      bool isBinary)
    : Analyzer(processor),
      nSample_(0),
      isInitialized_(false),
      isBinary_(isBinary),
      boundaryPtr_(0),
      atomStoragePtr_(0)
      #ifdef INTER_BOND
      , bondStoragePtr_(0)
      #endif
      #ifdef INTER_ANGLE
      , angleStoragePtr_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralStoragePtr_(0)
      #endif
   {
      setClassName("TrajectoryWriter");
      boundaryPtr_ = &processor.boundary();

      atomStoragePtr_ = &processor.atoms();
      #ifdef INTER_BOND
      bondStoragePtr_ = &processor.bonds();
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_ = &processor.angles();
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_ = &processor.dihedrals();
      #endif
   }

   /*
   * Read interval and outputFileName.
   */
   void TrajectoryWriter::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   #if 0
   /*
   * Load internal state from an archive.
   */
   void TrajectoryWriter::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to output archive.
   */
   void TrajectoryWriter::save(Serializable::OArchive& ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }
   #endif

   /*
   * Setup - open the trajectory file.
   */
   void TrajectoryWriter::setup()
   {  
      FileMaster& fileMaster = processor().fileMaster();
      if (isIoProcessor()) {
         if (isBinary()) {
            fileMaster.openOutputFile(outputFileName(), outputFile_, 
                                      std::ios::out | std::ios::binary);
         } else {
            fileMaster.openOutputFile(outputFileName(), outputFile_);
         }
      }
      writeHeader(outputFile_);
   }

   /*
   * Write a frame to file.
   */
   void TrajectoryWriter::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         writeFrame(outputFile_, iStep);
         ++nSample_;
      }
   }

   /*
   * Clear sample counter and close file.
   */
   void TrajectoryWriter::clear()
   {
      nSample_ = 0;
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
   }

   /*
   * Clear sample counter and close output file.
   */
   void TrajectoryWriter::output()
   {  clear(); }

}
