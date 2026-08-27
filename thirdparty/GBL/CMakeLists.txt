cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(
    STATUS
    "Building General Broken Lines (GBL) from Mille as part of the Acts project"
)

# Declare where to get GBL from.
FetchContent_Declare(GBL SYSTEM ${ACTS_GBL_SOURCE} SOURCE_SUBDIR cpp)

# Now set up its build.
FetchContent_MakeAvailable(GBL)

# GBL declares an export set for its library target but never installs it,
# and the config file it ships is variable-based and defines no imported
# target. Either way, anything that exports a target linking `GBL` (such as
# ActsPluginGBL) would end up with a dangling reference. Install the export
# set alongside the Acts package config instead and let ActsConfig pull it in
# directly, mirroring what the EDM4hep plugin does for ActsPodioEdmTargets.
install(
    EXPORT GBL-targets
    FILE ActsGBLTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Acts
)
