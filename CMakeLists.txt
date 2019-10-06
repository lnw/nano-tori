# Copyright (c) 2019, Lukas Wirz
# All rights reserved.

# This file is part of 'nano-tori' which is released under the BSD-2-clause license.
# See file LICENSE in this project.

cmake_minimum_required(VERSION 3.2.0)
project(nano-tori-gen)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(libs)

include_directories(AFTER SYSTEM .)

find_package(GSL)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  list(APPEND libs ${GSL_LIBRARIES})
  add_definitions(-DHAS_GSL)
endif(GSL_FOUND)

set(PROGRAMS1
  torus-simple
  torus-full
  torus-grid
  torus-extra
)

set(PROGRAMS2
  botch_optimise-torus-full
)

add_library(geom geometry2.cc geometry3.cc)

foreach(program ${PROGRAMS1})
  set(program_source ${program}.cc)
  add_executable(${program} ${program_source})
  target_link_libraries(${program} geom)
endforeach(program ${PROGRAMS})

foreach(program ${PROGRAMS2})
  set(program_source ${program}.cc)
  add_executable(${program} ${program_source})
  target_link_libraries(${program} geom ${libs})
endforeach(program ${PROGRAMS})

