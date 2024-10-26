# Package for Analysis of Target Thicknesses and Yields by Method of Energy Loss and Transfer (PATTYMELT)
### v. 0.4.4; Author: T. Hankins, Date Last Modified: 241025

## Quick Access

- [Purpose, Limitations, and Requirements](#purpose-limitations-and-requirements)
- [Package Contents](#package-contents)
- [Setup](#setup)
- [Setup (Alternate)](#setup-alternate)
- [Usage](#usage)
    - [`perform_histogram_measurement` File Output](#perform_thickness_measurement-file-output)
- [FAQ](#faq)
- [Changelog](#changelog)
- [To Do](#to-do)





## Purpose, Limitations, and Requirements

- This package was written to streamline working with the TAMU Cyclotron Institute Target Testing Station, which uses a position-sensitive dual-axis duo-lateral (DADL) silicon detector.
- To use this package, data in the form of **Reduced Events** is needed. This can be obtained from raw data using the tools associated with the `063123/` directory that can be checked out from `CVS`. As the production of Reduced Events is beyond the scope of this package, the details of this procedure are not offered here.
- If the experiment directory is changed from `063123/`, read below to learn how to update the package.
- The choice of alpha source is important - using <sup>229</sup>Th is recommended (especially for targets that are damaged). The package currently supports using <sup>228,229</sup>Th.
- The alpha source energy loss method itself cannot be used for targets that have significant thickness variations, such as salt-like targets (e.g., Gd<sub>2</sub>O<sub>3</sub>). At the moment, the package also cannot handle targets with backings.





## Package Contents

- The most up-to-date information about the package file structure can be obtained by running `tree -I 'build'` in the root directory of the Target Testing package. In short, this package is organized similar to many other CMake-based builds:
    - `bin/` - containes the binaries (i.e., executables) that are compiled when the package is `make install`ed. Executables for analysis will be used inside of this directory.
    - `build/` - contains the CMake build files for the package. Depending on how this package is acquired, this directory may or may not exist.
    - `src/` - contains the source files from which the package is built. Source codes are modified in this directory.
- In addition to these directories, there are several other relevant files, all of which can be found in the `src/` directory.
    - `PathManager.h` - holds the relevant paths for accessing files and directories. Setting these here will set them everywhere.
    - `AnalysisTools.h` - holds analysis functions that are used in various scripts, as well as the definition for the radioactive source energies.
    - `hdbscan.h` - has an implementation of Hierarchical Density Based Spatial Clustering with Applications for Noise (HDBSCAN) (original version written by Bryan Harvey) to handle automatic clustering of source data for `front_vs_front` and `back_vs_back`.





## Setup

- This section is dedicated to the initialization of the package, including pre-builds and file management systems.
- If the directory being used for analysis is not `063123/`, go to [Setup (Alternate)](#setup-alternate).

### Initializing the Package and Directories

0. The library for the target testing acquisition needs to exist for the user. At the time of writing, the directory being used is `063123/`, and each unique account using this module needs to have the corresponding library. This is created when the `063123/` directory is checked out from `CVS` and `make install`ed. To check if you've got the library already:
    ```
    ls ~/lib/lib063123.so
    ```
    If you get a "No such file or directory" result, address this problem before continuing. The fastest way to address this:
    ```
    cd ~/
    cvs checkout 063123
    cd 063123
    ./autogen
    make
    make install
    ```
    Wait for each operation to complete before moving on in case there are errors, compile or otherwise.
1. Return to the package directory and begin by defining your desired root file directory (`ROOT_FILE_DIR`) in the package's `src/PathManager.h` (e.g., `/data/sjygroup/sjy25/han61940/TargetTesting/`). Note that this location needs to exist. Absolute paths are recommended.
2. If the `bin/` directory does not already exist, return to the root of the package, then:
    ```
    mkdir bin
    ```
    Similarly, if the `build/` directory does not already exist, return to the root of the package, then:
    ``` 
    mkdir build
    cd build
    cmake ..
    ```
    This creates a `build/` directory alongside the `bin/` and `src/` directories then initializes CMake. If the `build/` directory does exist, delete it and recreate it; the presence of pre-existing build files *will* cause compilation errors. It is not necessary to delete and recreate the `bin/` directory if it does not already exist.
    - Depending on the computer, the CMake configuration and build file generation can take as long as 90 seconds, so please be patient. Generally, it takes somewhere between 15-30 seconds.
3. If the build files generate succesfully, run `make` from inside the `build/` directory. If several cores are available, run the variant `make -jN`, where *N* is the number of cores to allocate to the process (e.g., 4). If this is successful, run `make -jN install`.
4. This package comes with a built-in executable which automatically creates relevant directories for reading and writing files. It is recommended, though not required, that this approach be followed, otherwise additional steps will have to be taken to modify the source codes. To continue, move to the `bin/` directory and run `initialize_directories`. From inside `build/`:
    ```
    cd ../bin
    ./initialize_directories
    ```
    The program will create the directories in the `ROOT_FILE_DIR` that was defined earlier and print a message to the terminal when complete. Check that the directories were created.

### Managing Data and Defining Additional Directories

- The current naming convention for Reduced Files is "R + experiment date + file number". However, the naming "variable" `FILE_NAME_FORMAT` is preprocessor defined in `PathManager.h`. The format of this can be changed, if desired. At the time of writing, it is currently "R%s%03d.root", which takes a string (%s) representing the experiment date (from `PathManager.h`) and a signed integer (%d) representing the run number (offered as an executable runtime option).

1. If following the file management practice established above, copy the Reduced Files to be analyzed into the `reduced/` folder under the `ROOT_FILE_DIR`. This location should already exist following `initialize_directories`. Otherwise, modify the definition of `REDUCED_FILE_DIR` in `PathManager.h`, then `make install`. The other directories can also be changed, if desired.
2. If the calibration run number is already known (alpha source on a DADL with no mask or target), change the run number in `PathManager.h`, then `make install`.
3. Proceed to [Usage](#usage).





## Setup (Alternate)

- This section is dedicated to the initialization of the package, including pre-builds and file management systems specifically when using experimental directories other than `063123/`.
- For the sake of this description, let's assume the name of the experimental directory is `TargTest/`. This means that the following (or something similar) should apply for the directory in comparison to `063123/`:
    - `063123/` -> `TargTest/` (directory)
    - `lib063123.so` -> `libTargTest.so` (shared library)
    - `T063123Event.h` -> `TTargTestEvent.h` (header)
    - `T063123Event` -> `TTargTestEvent` (object)
    It is assumed for the following description that this applies. Adjust accordingly if there is a difference.

0. The library for the target testing acquisition needs to exist for the user. If the directory `TargTest/` is being used, each unique account using this module needs to have the corresponding library. This is created when the `TargTest/` directory is checked out from `CVS` and `make install`ed. To check if you've got the library already:
    ```
    ls ~/lib/libTargTest.so
    ```
    If you get a "No such file or directory" result, address this problem before continuing.
1. Writing the module with dependencies on `063123/` was inevitable (generally, the main problem is a `T063123Event` object (or pointer) in each code). To change this dependence, navigate to `src/` and run the following:
    ```
    sed -i 's/063123/TargTest/g' *
    ```
    The stream editor (`sed`) command can edit files in place without the need for backups (`-i`). It works by substituting (`s`) strings matching the first argument (`063123`) with those of the second argument (`TargTest`) (separated by `/`) globally (`g`) within a file. This action is then applied to all files present in the directory (`*`). Verify that the action worked successfully by then running
    ```
    grep -r 063123
    ```
    If nothing is found, the process was successful.
2. Return to step 1 of normal [Setup](#initializing-the-package-and-directories).





## Usage

- This package was designed to proceed stepwise from Reduced Events with no calibrations to a reliable target thickness measurement. This section outlines the procedure and methods for the package.
- The DADL detector requires a series of calibrations to use correctly. The executables can be ran in the order defined below to perform the calibration without knowing the finer details of the process, but interested users are directed to [NIM A 1050, 168130 (2023)](https://doi.org/10.1016/j.nima.2023.168130) to learn more.
- All executables are located in the `bin/` directory after successfully `make install`ing the package.
- In general, the executables come with several options that can be specified:
    - `r` - the number of the run to be analyzed. By default, it is the `CALIB_RUN` number provided in `PathManager.h`. If the file doesn't exist, the program will crash. Assuming `CALIB_RUN` is set properly, this option doesn't really need to be used until `perform_thickness_measurement`.
    - `e` - the number of events to be analyzed from the file. The default value is `std::numeric_limits<int>::max()` (2<sup>31</sup> - 1), and the number that is analyzed is the lesser between the max number of events in the event tree and the provided argument.
    - `d` - if there is an option to do so, display graphs/data. The default is `true`.
    - `c` - exclusive to `perform_thickness_measurement`; determines coarseness of the position map (in mm).
    - `m` - exclusive to `perform_thickness_measurement`; determines CycSrim material to use for the thickness calculation.
- The analysis order is as follows:
    - `front_vs_back` - determines relationship between front and back energy sums. The slope and intercept of the relationship are written to a calib file.
    - `front_vs_front` - performs gain matching between the front signals. The slope and intercept of the relationship are written to a calib file.
        - This executable uses the clustering algorithm HDBSCAN to cluster points and extract calibration parameters. To achieve a reasonable fit, this requires a considerable number of events (~ 100k within the MST, this is printed to the terminal during execution. A reasonable number of total events is about 250k). As a result, the algorithm may take a few minutes to run. The implementation is at worst *O(N<sup>2</sup>)*, which can become tedious with a significant number of entries.
        - Properties of the HDBSCAN clusterer can be edited in the source code around line 100 if the fit is to be adjusted. See `hdbscan.h` for discussion of properties.
    - `back_vs_back` - performs gain matching between the back signals. The slope and intercept of the relationship are written to a calib file.
        - See description for `front_vs_front`.
    - `check_gain_match` - checks gain matching parameters. Not required, but recommended.
    - `stretching_parameters` - determines the stretching parameters to scale the raw position data to physical limits. The parameters are written to a calib file.
    - `check_stretching_parameters` - checks stretching parameters. Not required, but recommended.
    - `energy_calibration` - determines the energy conversion between channels and MeV. The slope and intercept of the conversion are written to a calib file.
        - This executable *may* require user interaction; see [README_interactables.md](./README_interactables.md) in the root for more information.
    - `perform_thickness_measurement` - performs the actual calculation of a target's thickness. Energy spectra, hit maps, and thickness maps are written to an output file in the specified `OUTPUT_FILE_DIR` in `PathManager.h`.
        - This executable requires user interaction; see [README_interactables.md](./README_interactables.md) in the root for more information.

### `perform_thickness_measurement` File Output

- The coarseness of position maps in the file output depends on the coarseness option `c` originally specified when running `perform_thickness_measurement`.
- `h_ThicknessMap` - position map of thicknesses calculated using the consolidated energy spectra.
- `Diagnostics` Directory
	- `h_LowerError` and `h_UpperError` - complements to `h_ThicknessMap` that hold the corresponding bin error. `
	- `h_HitMap_Uncorrected` - position map of incident source particles without having applied stretching parameters to the DADL. The limits in X and Y should be roughly +/- 0.3 and +/- 0.5, respectively, excluding spurious entries spanning beyond the physical dimensions of the face.
	- `h_HitMap` - position map of incident source particles with all DADL calibrations and event requirements (i.e., front versus back) applied. The number of entries in any given bin is the number of entries in the source spectra that is fit to calculate the corresponding thickness. The bins between `h_HitMap` and `h_ThicknessMap` are 1:1 and can be directly compared.
	- `h_Energy` - total energy spectrum. Included for easy access to a calibrated energy spectrum when analyzing a calibration run with the package, but can also be used to show peak shifting due to energy loss.
	- `h_ThicknessMap_Filtered_N` - to recover thicknesses where fits failed, a recovery subroutine was designed to iteratively re-fit bad bins. A series of these "filtered" histograms exist where each subsequent histogram has more re-fit bins than the last. The production of the first filtered histogram is determined by a threshold set in the code.
- `Source Spectra` Directory
	- `h_.2f_.2f` - list of source spectra used to calculate thicknesses. The format is h_X_Y, where "X" or "Y" is the local position bin center in the corresponding dimension *measured in millimeters*. As an example, a bin in `h_ThicknessMap` with coordinates X = [0.2, 0.4] and Y = [-0.6, -0.4] would correlate with the energy spectrum `h_3.00_-5.00`, as the bin centers measured in millimeters are 3.00 and -5.00, respectively.





## FAQ

- "When I run `energy_calibration` or `perform_thickness_measurement`, not all of the peaks in my spectrum are identified. Why is that?"
    - To identify peaks in the spectra, a ROOT CERN class called `TSpectrum` is used. One of the flexibilities of this class is the ability to set the threshold at which the algorithm considers "peaks" peaks. This is a relative value from 0 -> 1, in which 0 is baseline and 1 is the maximum bin in the histogram. The user has to manually choose the threshold for identification, and I generally set the threshold somewhere between 0.1-0.4 for 228Th and 229Th spectra. This, of course, varies on a case-by-case basis and becomes significant when dealing with low statistics. One of the peaks may be missing because this threshold is too high. To fix this, find within the script the threshold declaration, adjust it, and rebuild the project as necessary, but keep in mind that lowering the threshold can make other peaks become identified by the algorithm under certain circumstances.
- "Nothing displays when running `perform_thickness_measurement`. Did something go wrong?
    - The success (or lack thereof) of running this script is entirely determined by terminal output. This executable was originally designed with a display feature, but an oversight in `TSpectrum` means that `TCanvas`es are generated on each spectrum search. This significantly slows down the procedure, even with manually taking control over the auto-generated canvases. To mitigate this, the display feature was omitted (but the source code is still present in the file). Instead, view the generated ROOT file when the executable is complete to determine the outcome.
- "My output when running `energy_calibration`/`perform_thickness_measurement` is all kinds of wrong. What's going on?
    - See corresponding section in [README_interactables.md](./README_interactables.md) in the root.





## Changelog

- 0.4.4 (241025) : more minor reorganization for preliminary upload to GitHub; generalized parameter file read in into set of methods in `AnalysisTools.h` which are now called in scripts rather than copy-pasted everywhere.
- 0.4.3 (241017) : corrected energy loss error in `perform_thickness_measurement`; added `CycSrimHandler.h` to be able to pass SRIM material as argument in executable so that it doesn't have to be recompiled each time the material is changed. Added a rudimentary generalization to `perform_thickness_measurement` so that sources that don't have radon leakage (i.e., <sup>229</sup>Th) can be analyzed with the same script without iterative correction passes. Cleaned up files in `src/` in preparation for upload to GitHub.
- 0.4.2 (241003) : refactored `perform_thickness_measurement` output to include upper and lower error TH2F objects rather than a TGraphAsyErrors. Reorganized outputs for ease of access, organizing the spectra into the final thickness map and groups of the diagnostic and source energy spectra plots.
- 0.4.1 (240912) : adjusted `perform_thickness_measurement` to account for position bin edge cases when calculating thicknesses. Due to radon emission in <sup>228</sup>Th, energy spectra beyond the active area of the target are still populated. An arbitrary threshold was added to bypass evaluation of bins where the number of entries is less than the threshold percent of the number in the most populated bin. Adjusted `perform_thickness_measurement` to recalculate thicknesses in bad bins with two peaks rather than one, then to recursively perform the recalculation based on available references around the bin being recalculated until all to be recalculated are.
- 0.4.0 (240820) : created helper file `hdbscan.h` based on Bryan Harvey's version of HDBSCAN to handle clustering for `front_vs_front` and `back_vs_back` and edited the corresponding files to implement this. Further edited these files to perform principal component analysis by means of TPrincipal for extraction of calibration parameters rather than using linear fits.
- 0.3.1 (240808) : added new relevant value to `PathManager.h` which controls the naming convention for ReducedFile read-in. Updated paths in all executables; changing the file name format in `PathManager.h` should now permit reading files with different conventions. Also corrected bug in `perform_thickness_measurement` that was causing truncation of the measured peak centroids from the source fits. Added additional code in `perform_thickness_measurement` to address misfits in the thickness map.
- 0.3.0 (240720) : created `back_vs_back` and updated dependent executables. This is currently an untested executable (as the test DADL data doesn't require back gain correction), although it's effectively a carbon copy of `front_vs_front`. Should work with minimal changes needed, if any at all. Updated `README`s based on this addition. Additionally, fixed a position calculation error and accounted for another numerical rounding error bug in `perform_thickness_measurement`.
- 0.2.1 (240718) : `std::fmod()` operation exchanged for analog `std::remainder()` due to floating point evaluation inaccuracy in `perform_thickness_measurement`.
- 0.2.0 (240718) : successful preliminary test measuring thickness. Updated `README`; created [README_interactables.md](./bin/README/README_interactables.md) in `bin/README`.
- 0.1.0 (240716) : package created. Added `generate_raw_summary_file`, `initialize_directories`, `front_vs_back`, `front_vs_front`, `check_gain_match`, `stretching_parameters`, `check_stretching_parameters`, `energy_calibration`, and `perform_thickness_measurement`, as well as `AnalysisTools.h` and `PathManager.h`.





## To Do
- Once the system works, remove unnecessary includes/dependencies in the primary CMakeLists.txt
- Test the signal orientation and math for the main executable with a mask run.

- Option to include a distance measurement between the source and DADL, then correct calibrations and measurements based on 1/cos effect (reaches >=5% when the source is placed < 3.3 cm (~1.3 in) from the target/DADL).
- Option for backing to be included?
