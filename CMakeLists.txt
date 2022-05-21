cmake_minimum_required(VERSION 3.13)
project(QoZ VERSION 0.1.0.0)
include(GNUInstallDirs)
include(CTest)

option(BUILD_SHARED_LIBS "build shared libraries by default" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
  set(CMAKE_BUILD_TYPE ${default_build_type} CACHE
          STRING "Choose the type to build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
          "MinSizeRel" "RelWithDebInfo")
endif()

find_package(PkgConfig)
find_package(OpenMP)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/QoZ/version.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/QoZ/version.hpp)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
target_compile_features(${PROJECT_NAME}
    INTERFACE cxx_std_17
  )

target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)

option(QoZ_USE_BUNDLED_ZSTD "prefer the bundled version of Zstd" OFF)
option(QoZ_DEBUG_TIMINGS "print debug timing information" ON)

if(QoZ_DEBUG_TIMINGS)
  target_compile_definitions(${PROJECT_NAME} INTERFACE QoZ_DEBUG_TIMINGS=1)
else()
  target_compile_definitions(${PROJECT_NAME} INTERFACE QoZ_DEBUG_TIMINGS=0)
endif()

pkg_search_module(ZSTD IMPORTED_TARGET libzstd)
if(ZSTD_FOUND AND NOT QoZ_USE_BUNDLED_ZSTD)
  message("Zstd found.")
  target_link_libraries(${PROJECT_NAME} INTERFACE PkgConfig::ZSTD)
else()
  message("Zstd not found.")
  add_subdirectory(zstd)
  target_link_libraries(${PROJECT_NAME} INTERFACE zstd)
endif()

find_package(GSL)
if (GSL_FOUND)
  target_compile_definitions(${PROJECT_NAME} INTERFACE -DENABLE_GSL="TRUE")
endif()

add_subdirectory(test)

install(TARGETS ${PROJECT_NAME} 
  EXPORT QoZTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT QoZTargets NAMESPACE QoZ:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QoZ)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/QoZConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/QoZConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QoZ
  )
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/QoZConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/QoZConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/QoZConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QoZ
  )
