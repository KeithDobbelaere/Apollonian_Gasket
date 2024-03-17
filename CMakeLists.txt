cmake_minimum_required(VERSION 3.0.0)
# Set for static linking
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

project(ApollonianProj VERSION 0.1.0 LANGUAGES C CXX)

# Set the output directory for the build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/)

# Add the include directories for the project
include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(Apollonian main.cpp)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find and add the SFML static library
add_definitions(-DSFML_STATIC)
set(SFML_STATIC_LIBRARIES TRUE)
find_package(SFML 3.0 COMPONENTS Graphics REQUIRED)

target_link_libraries(Apollonian SFML::Graphics)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
