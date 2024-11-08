#ifndef PATH_MANAGER__H
#define PATH_MANAGER__H

// Includes are at the bottom of this header to improve ease-of-access.

//___________________________________________________________________________________________________

// Relevant Values
#define EXPT_DATE "063123"                  // Which experimental directory holds the data
#define CALIB_RUN 81                        // What run contains data for blank DADL (stretching and E calib)
#define FILE_NAME_FORMAT "R%s%03d.root"     // What naming convention ReducedFiles follow 
// #define MASK_RUN

// External directories
#define ROOT_FILE_DIR "/data/sjygroup/sjy25/han61940/TargetTesting/"
#define REDUCED_FILE_DIR ROOT_FILE_DIR "reduced/"
#define CALIB_FILE_DIR ROOT_FILE_DIR "calib/"
#define OUTPUT_FILE_DIR ROOT_FILE_DIR "output/"

//___________________________________________________________________________________________________


// Common C++ includes. Uncommon to do this way, but I don't feel like writing this thousands of
// times in all scripts.
#include <iostream>
#include <memory>

// Common ROOT includes. See above.
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"

// Other includes; It is uncommon to place header files in this way, but it makes life easier if
// one value needs to change for all executables (e.g., the event type)
#include <BrAppOptionManager.h>
#include <CycSrim.h>
#include "T063123Event.h"

#endif