# relative path to the header files
set(HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}) 

# Configure Path Manager
configure_file(PathManager.h PathManager.h)
# include("${CMAKE_CURRENT_SOURCE_DIR}/PathManager.h")
# add_include("PathManager.h")

# Generate the root dictionary
root_generate_dictionary("G__${LIB_NAME}" ${LIB_HEADERS} LINKDEF LinkDef.h)

# create the shared library
add_library(${LIB_NAME} SHARED ${LIB_SOURCES} "G__${LIB_NAME}.cxx")
target_include_directories(${LIB_NAME} PRIVATE ${CYCAPP_INCLUDE_DIRECTORY} ${ROOT_INCLUDE_DIRS} ${HEADER_DIR})
target_link_libraries(${LIB_NAME} PRIVATE ROOT::Physics ROOT::RIO ROOT::GX11 ROOT::RGL)
target_compile_options(${LIB_NAME} PRIVATE -Wall)


# install path for the shared library
install(TARGETS ${LIB_NAME} DESTINATION $ENV{HOME}/lib)

# install path for the generated library
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}_rdict.pcm DESTINATION $ENV{HOME}/lib)
# only use the following if you want to always have this library loaded
# install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}.rootmap DESTINATION $ENV{HOME}/lib) 

# Add executables

# initialize_directories
add_executable(initialize_directories initialize_directories.cpp)
target_include_directories(initialize_directories PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(initialize_directories PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS initialize_directories DESTINATION ${INSTALL_EXEC_DIRECTORY})

# generate_raw_summary_file
add_executable(generate_raw_summary_file generate_raw_summary_file.cpp)
target_include_directories(generate_raw_summary_file PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(generate_raw_summary_file PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS generate_raw_summary_file DESTINATION ${INSTALL_EXEC_DIRECTORY})

# front_vs_back
add_executable(front_vs_back front_vs_back.cpp)
target_include_directories(front_vs_back PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(front_vs_back PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS front_vs_back DESTINATION ${INSTALL_EXEC_DIRECTORY})

# front_vs_front
add_executable(front_vs_front front_vs_front.cpp)
target_include_directories(front_vs_front PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(front_vs_front PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS front_vs_front DESTINATION ${INSTALL_EXEC_DIRECTORY})

# back_vs_back
add_executable(back_vs_back back_vs_back.cpp)
target_include_directories(back_vs_back PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(back_vs_back PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS back_vs_back DESTINATION ${INSTALL_EXEC_DIRECTORY})

# check_gain_match
add_executable(check_gain_match check_gain_match.cpp)
target_include_directories(check_gain_match PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(check_gain_match PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS check_gain_match DESTINATION ${INSTALL_EXEC_DIRECTORY})

# stretching_parameters
add_executable(stretching_parameters stretching_parameters.cpp)
target_include_directories(stretching_parameters PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(stretching_parameters PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS stretching_parameters DESTINATION ${INSTALL_EXEC_DIRECTORY})

# check_stretching_parameters
add_executable(check_stretching_parameters check_stretching_parameters.cpp)
target_include_directories(check_stretching_parameters PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(check_stretching_parameters PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT})
install(TARGETS check_stretching_parameters DESTINATION ${INSTALL_EXEC_DIRECTORY})

# energy_calibration
add_executable(energy_calibration energy_calibration.cpp)
target_include_directories(energy_calibration PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(energy_calibration PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT} Spectrum)
install(TARGETS energy_calibration DESTINATION ${INSTALL_EXEC_DIRECTORY})

# perform_thickness_measurement
add_executable(perform_thickness_measurement perform_thickness_measurement.cpp)
target_include_directories(perform_thickness_measurement PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(perform_thickness_measurement PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT} Spectrum)
install(TARGETS perform_thickness_measurement DESTINATION ${INSTALL_EXEC_DIRECTORY})

# calculate_position_resolution
add_executable(calculate_position_resolution calculate_position_resolution.cpp)
target_include_directories(calculate_position_resolution PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
target_link_libraries(calculate_position_resolution PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT} Spectrum)
install(TARGETS calculate_position_resolution DESTINATION ${INSTALL_EXEC_DIRECTORY})

# test_executable
# add_executable(test_executable test_executable.cpp)
# target_include_directories(test_executable PRIVATE ${ROOT_INCLUDE_DIRS} ${CYCAPP_INCLUDE_DIRECTORY} ${INCLUDE_EXPT})
# target_link_libraries(test_executable PRIVATE ROOT::Core ROOT::Gui ROOT::Physics ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::RIO ${CYCAPP_LIBS} ${LIB_EXPT} Spectrum)
# install(TARGETS test_executable DESTINATION ${INSTALL_EXEC_DIRECTORY})