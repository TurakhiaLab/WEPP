cmake_minimum_required (VERSION 3.8) 

project(wepp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -O0 -ggdb3 -fno-eliminate-unused-debug-symbols")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# otherwise it sometimes goes to system protobuf instead
set(Protobuf_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")

if(DEFINED Protobuf_PATH) 
    find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
    find_package(Protobuf REQUIRED) 
endif()

include_directories(${Protobuf_INCLUDE_DIRS})

message("USING ${Protobuf_VERSION}")

if(DEFINED CapnProto_PATH) 
	find_package(CapnProto REQUIRED HINTS ${CapnProto_PATH}) 
else()
	find_package(CapnProto REQUIRED) 
endif()

include_directories(${CAPNP_INCLUDE_DIRS})
add_definitions(${CAPNP_DEFINITIONS})

file(GLOB WEPP_SRCS "src/*.cpp" "src/*.hpp")
    
# Include JSONCPP
include(${CMAKE_TOOLCHAIN_FILE})
find_package(jsoncpp CONFIG REQUIRED)

file(GLOB PANMAT_SRCS "src/panman/panman.cpp" "src/panman/panmanUtils.hpp" "src/panman/panman.hpp")

if(DEFINED CapnProto_PATH)
	capnp_generate(
		LANGUAGE cpp
		TARGET wepp  
		PROTOS panman.capnp
    )
else()
	capnp_generate_cpp(
		CAPNP_SRCS CAPNP_HDRS 
		panman.capnp
    )
endif()

if(DEFINED Protobuf_PATH)
    protobuf_generate(
        LANGUAGE cpp
        TARGET wepp
        PROTOS sam.proto
    )
    
    add_executable(wepp
	    ${WEPP_SRCS} ${PANMAT_SRCS}
    )
else()
    protobuf_generate_cpp(
        SAM_SRCS SAM_HDRS 
        sam.proto
    )
   
    add_executable(wepp
	    ${WEPP_SRCS} ${PANMAT_SRCS}
        ${SAM_SRCS} ${PROTO_SRCS} ${CAPNP_SRCS}
		${SAM_HDRS} ${PROTO_HDRS} ${CAPNP_HDRS} 
    )

    add_executable(closest_peak_clustering
        src/panman/closest_peak_clustering.cpp
    )
endif()

TARGET_LINK_LIBRARIES(wepp PRIVATE stdc++ JsonCpp::JsonCpp ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ${CAPNP_LIBRARIES}) 
TARGET_LINK_LIBRARIES(closest_peak_clustering PRIVATE stdc++ ${TBB_IMPORTED_TARGETS})
