cmake_minimum_required (VERSION 3.8) 
cmake_policy(SET CMP0069 NEW)
project(usher)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()
set(usher_VERSION "?.?.?")
find_library(ISAL_LIB isal)
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if ("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(usher_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
message(STATUS "${usher_VERSION}")

configure_file(src/version.hpp.in version.hpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -std=c++17 -Wall -Wno-unused-function  -Wno-deprecated-declarations") 
#-DDEBUG_PARSIMONY_SCORE_CHANGE_CORRECT  -DCHECK_PAR_MAIN  -DEASY_DEBUG -DSTOP_ON_ERROR    -fsanitize=address
set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -O0 -ggdb3 -fno-eliminate-unused-debug-symbols")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

if(DEFINED Protobuf_PATH) 
    find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
    find_package(Protobuf REQUIRED) 
endif()

# Print version messages
if(Protobuf_FOUND)
    message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${Protobuf_INCLUDE_DIRS})

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

#find_package(OpenMP REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
file(GLOB WBE_SRCS "src/WBE/*.cpp" "src/WBE/*.hpp")

set_source_files_properties(src/usher_mapper.cpp PROPERTIES COMPILE_FLAGS -O3)
find_package(ZLIB)
if(DEFINED Protobuf_PATH)
    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher_common.cpp
        src/usher.cpp
        )
    if(USHER_SERVER)
        add_executable(usher_server
            src/mutation_annotated_tree.cpp
            src/usher_mapper.cpp
            src/usher_common.cpp
            src/usher_server.cpp
            )
    endif()

    protobuf_generate(
        LANGUAGE cpp
        TARGET ripples 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET usher 
        PROTOS parsimony.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET wbe
        PROTOS parsimony.proto)
    
    protobuf_generate(
        LANGUAGE cpp
        TARGET wbe
        PROTOS sam.proto)

    if(USHER_SERVER)
        protobuf_generate(
            LANGUAGE cpp
            TARGET usher_server
            PROTOS parsimony.proto)
    endif()

else()
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS
        parsimony.proto)

    protobuf_generate_cpp(
        SAM_SRCS SAM_HDRS
        sam.proto)

    protobuf_generate_cpp(
        DETAILED_MUTATIONS_PROTO_SRCS
        DETAILED_MUTATIONS_PROTO_HDRS
        mutation_detailed.proto)

    protobuf_generate_cpp(
        TAXO_SRCS TAXO_HDRS
        taxodium.proto)
    
    add_executable(usher
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        src/usher_common.cpp
        src/usher.cpp
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )

    add_executable(wbe
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
	    ${WBE_SRCS}
        ${PROTO_SRCS} ${SAM_SRCS}
        ${PROTO_HDRS} ${SAM_HDRS}
        )
    
    add_executable(closest_peak_clustering
        src/WBE/Freyja/closest_peak_clustering.cpp
        )

    if(USHER_SERVER)
        add_executable(usher_server
            src/mutation_annotated_tree.cpp
            src/usher_mapper.cpp
            src/usher_common.cpp
            src/usher_server.cpp
            ${PROTO_SRCS}
            ${PROTO_HDRS}
            )
    endif()
endif()

if(SAVE_PROFILE)
    if(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DSAVE_PROFILE=1 -DDEBUG=1)
    else(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DSAVE_PROFILE=1)
    endif(DEBUG)
else(SAVE_PROFILE)
    if(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DDEBUG=1)
    else(DEBUG)
        TARGET_COMPILE_OPTIONS(usher PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    endif(DEBUG)
endif(SAVE_PROFILE)

TARGET_LINK_LIBRARIES(usher PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
target_include_directories(usher PUBLIC "${PROJECT_BINARY_DIR}")

TARGET_LINK_LIBRARIES(wbe PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(closest_peak_clustering PRIVATE stdc++ ${TBB_IMPORTED_TARGETS}) # OpenMP::OpenMP_CXX)

if(USHER_SERVER)
    TARGET_COMPILE_OPTIONS(usher_server PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
    TARGET_LINK_LIBRARIES(usher_server PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
endif()