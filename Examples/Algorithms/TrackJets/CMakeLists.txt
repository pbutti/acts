add_library(
  ActsExamplesTrackJets SHARED
  src/TrackJetsAlgorithm.cpp)
target_include_directories(
  ActsExamplesTrackJets
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PUBLIC ${PROJECT_SOURCE_DIR}/Examples/Framework/include
  PUBLIC ${PROJECT_SOURCE_DIR}/Examples/EventData/include
  PUBLIC ${FASTJET_DIR}/build/include)
target_link_libraries(
  ActsExamplesTrackJets
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ${FASTJET_DIR}/build/lib/libfastjet.so)

install(
  TARGETS ActsExamplesTrackJets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
