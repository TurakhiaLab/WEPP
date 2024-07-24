cmake_minimum_required (VERSION 3.8) 

project(wbe)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -O0 -ggdb3 -fno-eliminate-unused-debug-symbols")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

#find_package(OpenMP REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

if(DEFINED Protobuf_PATH) 
    find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
    find_package(Protobuf REQUIRED) 
endif()

include_directories(${Protobuf_INCLUDE_DIRS})
file(GLOB WBE_SRCS "src/WBE/*.cpp" "src/WBE/*.hpp")
    
# Include JSONCPP
include(${CMAKE_TOOLCHAIN_FILE})
find_package(jsoncpp CONFIG REQUIRED)

file(GLOB PANMAT_SRCS "panman/panmanUtils.cpp" "panman/panman.cpp" "panman/panmanUtils.hpp" "panman/panman.hpp")

if(DEFINED Protobuf_PATH)
	protobuf_generate(
		LANGUAGE cpp
		TARGET panmanUtils 
		PROTOS protobuf/panman.proto)

    protobuf_generate(
        LANGUAGE cpp
        TARGET wbe
        PROTOS protobuf/sam.proto)
    
    add_executable(wbe
	    ${WBE_SRCS}
        ${PROTO_SRCS} 
        ${PROTO_HDRS}
        )
else()
	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		panman.proto)
	
    protobuf_generate_cpp(
        PROTO_SRCS PROTO_HDRS 
        sam.proto)
    
    add_executable(wbe
	    ${WBE_SRCS}
        ${PROTO_SRCS}
        ${PROTO_HDRS}
        )
endif()

TARGET_LINK_LIBRARIES(wbe PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) # OpenMP::OpenMP_CXX)