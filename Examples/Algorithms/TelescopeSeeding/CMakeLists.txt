add_library(
  ActsExamplesTelescopeSeeding SHARED
  src/TelescopeSeedingAlgorithm.cpp)
target_include_directories(
  ActsExamplesTelescopeSeeding
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTelescopeSeeding
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField
    ActsExamplesIoCsv ActsExamplesIoRoot
    Boost::program_options Eigen3::Eigen)

install(
  TARGETS ActsExamplesTelescopeSeeding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
