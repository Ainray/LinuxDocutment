path::
    CMAKE_SOURCE_DIR 
    CMAKE_BINARY_DIR
    PROJECT_SOURCE
    EXECUTABLE_OUTPUT_PATH
    CMAKE_RUNTIME_OUTPUT_PATH RUNTIME_OUTPUT_PATH
    CMAKE_RUNTIME_OUTPUT_DIRECTORY RUNTIME_OUTPUT_DIRECTORY
    CMAKE_LIBRARY_OUTPUT_DIRECTORY LIBRARY_OUTPUT_DIRECTORY
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY ARCHIVE_OUTPUT_DIRECTORY

build pipelines::
    CMake generates build pipelines, maybe .sln (VS), .xcodeproj (Xcode)
or a Unix-style Makefile and so on.

source::
binary::
    Cmake needs to know the "source" and "binary" folders to generate a
build pipeline.
    "Source" contains CMakeLists.txt, "Binary" will contain the build
pipeline. A common pracitce is to create a subdirectory "build" beneath 
CMakeLists.txt. Different "binary" folders can be created for different
build system or configuration options.
    ${CMAKE_SOURCE_DIRECTORY}
    ${CMAKE_BINARY_DIRECTORY}
    ${PROJECT_SOURCE_DIR} top full path of current project, i.e. refer to
the folder of the CMakeLists.txt containing the most recent project()
command (if exists).

cache::
    "cache" is a single text file in the "binary" named CMakeCache.txt
where "cache variables" are stored. Project-defined cache variables of which
can be checked by,

            cmake -L -N
    
    ".gitignore" file usually contains to the rule "*build*".

CMakeLists.txt::
configure::
    CMakeLists.txt is responsible for defining targets(e.g. executable,
library).

generator::
    To get lists of available generators, by 

            cmake --help

    To specify a desired generator using -G option,

            cmake -G "Visual Studio 15 2017"
            cmake -DCMAKE_C_COMPILER=bcc64x.exe -DCMAKE_CXX_COMPILER=bcc64.exe -G Ninja <cmakelists.txt-path>

options::
option::
            -DCMAKE_BUILD_TYPE=Debug/MinSizeRel/RelWithDebInfo/Release

cmake_minimum_required::
            cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project::
        no need c++ compiler, 
            project(example_Prj C)

add_executable::
        
false::
        Undefined variables, defined variables are empty or contain
    0,N,NO,OFF,FALSE, NOTFOUND, <var>-NOTFOUND

set::

            set(varname value ... CACHE type doctring [FORCE])

        + initialize an variable
            set(VAR)

        + type: FILEPATH, PATH, STRING, BOOL, INTERNAL
        + docstring: 
        + FORCE(optional), set cached variable

    for example,
        if(NOT CONFIGURED_ONCE) # as built-in variable, which is cached
                                # even before the first time our project is configured.
            set(CMAKE_CXX_FLAGS "{warnings}" 
                CACHE STRING "Flags used by the compiler during all types." FORCE)
        endif()

        set(CONFIGURED_ONCE TRUE CACHE INTERNAL "A flag showing that CMake has configured at least once")

include_directories::

CMAKE_CURRENT_SOURCE_DIR::
        Full path to the source directory the CMake is currently processing.

            include(CMAKE_CURRENT_SOURCE_DIR)

add_subdirectory::
        Directory must has a CMakeLists.txt

add_library::
        add_library(target [STATIC |SHARED |MODULE] sources...)

    + MODULE are plug-ins that arenot linked against but can be loaded
dynamically at runtime.
    + add_library(blas ${SRC})

variable::
    Two types, normal or cache, normal is first looked by CMake when both 
normal and cache names have the same names.

list::
    set(SRC)
    list(APPEND SRC ${SRC1} ${SRC2})
    list(REMOVE_DUPLICATES SRC)

macro::
    macro(lapack_install_library lib)
      install(TARGETS ${lib}
        EXPORT ${LAPACK_INSTALL_EXPORT_NAME}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
    endmacro()
sources::
    + file(GLOB *) and aux_source_directory
CMake do not know when a source is added or removed.
    + add_library(... OBJECT ...)
Not recursive.

target::
        There are three kinds of target files that may be built:
    + archive
    + library
    + runtime
        . Executables are always treated as runtime targets.
        . Static libraries are always treated as archive targets
        . Module libraries are always treated as library targets
        . Shared libraries(for non-dll platforms e.g. Linux) are treated
    as library targets.
        . Shared library (for dll platforms, e.g. Windows, Cygin) the 
    DLL part of a shared library is treated as a runtime target and 
    the corresponding import library is treated as an archive targets. 
