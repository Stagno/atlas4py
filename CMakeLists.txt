cmake_minimum_required(VERSION 3.14.5)
project(atlas4py LANGUAGES CXX)

include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
FetchContent_Declare(
    pybind11_fetch
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.4
)
FetchContent_GetProperties(pybind11_fetch)
if(NOT pybind11_fetch_POPULATED)
    message(STATUS "Fetch pybind11 for python-binding")
    FetchContent_Populate(pybind11_fetch)
    add_subdirectory(${pybind11_fetch_SOURCE_DIR} ${pybind11_fetch_BINARY_DIR})
endif()

find_package(atlas REQUIRED)

pybind11_add_module(_atlas4py MODULE _atlas4py.cpp)
target_link_libraries(_atlas4py PUBLIC atlas eckit)
target_compile_features(_atlas4py PUBLIC cxx_std_17)

pybind11_add_module(computation MODULE computation.cpp)
target_link_libraries(computation PUBLIC atlas eckit)
target_include_directories(computation PUBLIC /home/lukas/documents/work/dawn/dawn/prototype/)
target_include_directories(computation PUBLIC /home/lukas/documents/work/dawn/dawn/src/)
target_compile_features(computation PUBLIC cxx_std_17)
