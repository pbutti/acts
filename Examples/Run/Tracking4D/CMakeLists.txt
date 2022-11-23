set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesGenerators
  ActsExamplesGeneratorsPythia8
  ActsExamplesGeneratorsPythia8Options
  ActsExamplesTrackJets
  ActsExamplesIoRoot
  Boost::program_options)


add_executable(
  ActsExampleTrackJets
  TrackJetsExample.cpp)
target_link_libraries(
  ActsExampleTrackJets
  PRIVATE ${_common_libraries} )


install(
  TARGETS
  ActsExampleTrackJets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
