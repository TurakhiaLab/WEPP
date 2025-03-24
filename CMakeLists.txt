cmake_minimum_required (VERSION 3.8) 
cmake_policy(SET CMP0069 NEW)
project(wbe)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

configure_file(src/version.hpp.in version.hpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -Wno-unused-function  -Wno-deprecated-declarations") 
set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -O0 -ggdb3 -fno-eliminate-unused-debug-symbols")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

# otherwise it sometimes goes to system protobuf instead
set(Protobuf_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")

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

find_package(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

#find_package(OpenMP REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
file(GLOB WBE_SRCS "src/WBE/*.cpp" "src/WBE/*.hpp")

if(DEFINED Protobuf_PATH)
    protobuf_generate(
        LANGUAGE cpp
        TARGET wbe
        PROTOS sam.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET wbe
        PROTOS parsimony.proto)

    add_executable(wbe
	    ${WBE_SRCS} ${PANMAT_SRCS}
    )
    
    add_executable(closest_peak_clustering 
        src/closest_peak_clustering.cpp
    )
else()
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS
        parsimony.proto)

    protobuf_generate_cpp(
        SAM_SRCS SAM_HDRS
        sam.proto)

    add_executable(wbe
	    ${WBE_SRCS}
        src/mutation_annotated_tree.cpp
        src/usher_mapper.cpp
        ${PROTO_SRCS} ${SAM_SRCS}
        ${PROTO_HDRS} ${SAM_HDRS}
        )

    add_executable(closest_peak_clustering 
        src/closest_peak_clustering.cpp
    )
endif()

target_link_libraries(wbe PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES})
target_link_libraries(closest_peak_clustering PRIVATE stdc++ ${TBB_IMPORTED_TARGETS})
