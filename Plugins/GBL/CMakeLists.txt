acts_add_library(
    PluginGBL
    src/ActsToGbl.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginGBL
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginGBL PUBLIC ActsAlignment GBL)

acts_compile_headers(PluginGBL GLOB include/**/*.hpp)
