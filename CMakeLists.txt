#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(RPS2023)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(RPS2023  RPS2023.cc ${sources} ${headers})
target_link_libraries(RPS2023 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build neutronBall. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(RPS2023_SCRIPTS
  vis.mac
  run.mac
  recoilEmanationGPS.mac
  recoilEmanationParticleGun.mac
  uniformEmanation.mac
  plotPNames.c
  plotHistXY.c
  findRecoilRange.mac
  oneRun.mac

  plotEmanationVsH2OContent.c
  plotEmanationVsH2OContentMultipleGraphs.c
  plotLocationwiseRadonCount.c
  plotRadonCountAtTwoLocationsVsH2OContent.c 
  plotRadonCountInSurroundingGrainsVsH2OContent.c
  )

foreach(_script ${RPS2023_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS RPS2023 DESTINATION bin)
