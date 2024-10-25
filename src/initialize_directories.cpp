// This executable creates the subdirectories that will hold the various files needed for analyzing
// target thickness data.

// C++
#include <iostream>
#include <string>
#include <cstring> // strerror_s
#include <cerrno> // errno

// Other
#include "PathManager.h"

// The cluster is not on the C++17 or newer standard, so we don't have access to the <filesystem>
// library. This *should* be a workaround that considers all possible operating systems.
#if defined(_WIN32)
    #include <direct.h>
    #define mkdir _mkdir
#elif defined(__unix__) || defined(__APPLE__)
    #include <sys/stat.h>
    #include <sys/types.h>
#endif

bool create_directory(const std::string& path) {
    #if defined(_WIN32)
        int result = _mkdir(path.c_str());
    #elif defined(__unix__) || defined(__APPLE__)
        mode_t mode = 0755; // permissions: rwxr-xr-x
        int result = mkdir(path.c_str(), mode);
    #endif

    if (result == 0) {
        return true; // Directory created successfully
    } else {
        std::cerr << "Error creating directory: " << strerror(errno) << std::endl;
        return false;
    }
}

//____________________________________________________________________________________________________
int main(int argc, char **argv) {

    // Using the ROOT_FILE_DIR provided in PathManager.h, create the directories that will hold the
    // various files.

    create_directory(static_cast<std::string>(ROOT_FILE_DIR) + "reduced/");
    create_directory(static_cast<std::string>(ROOT_FILE_DIR) + "trees/");
    create_directory(static_cast<std::string>(ROOT_FILE_DIR) + "output/");
    create_directory(static_cast<std::string>(ROOT_FILE_DIR) + "calib/");

    std::cout << "initializing_directories : finished" << std::endl;

    return 0;
}