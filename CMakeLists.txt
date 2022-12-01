cmake_minimum_required(VERSION 2.8.9)
project (ewald)

include_directories(include/common)

# FLAGS FOR ALL CASES 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_STANDARD 17)

## CHECK FOR DEBUG MODE
if (DEBUG)
  message ("DEBUG MODE COMPILATION")
  set(CMAKE_BUILD_TYPE Debug)
else ()
  message ("RELEASE MODE COMPILATION")
  set(CMAKE_BUILD_TYPE Release)
endif (DEBUG)

## CHECK FOR VERBOSE EXECUTION MODE
if (VERBOSE_MODE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE_MODE=${VERBOSE_MODE}")
endif (VERBOSE_MODE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -std=c++17")

if (CMAKE_BUILD_TYPE MATCHES Debug)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    if (APPLE)
      message("Adding sanitizer for Clang in Apple")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    elseif(UNIX)
      message("Adding sanitizer for Clang in Linux")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    endif()
  endif()
endif()

add_subdirectory (src)
add_subdirectory (tests) 

