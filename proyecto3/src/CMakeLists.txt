# Copyright (C) 2018
# Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
#
# This file is part of the numerical analysis lecture CE3102 at TEC
#
# Author: Pablo Alvarado
# Date  : 28.12.2017

find_package (Boost COMPONENTS system filesystem program_options REQUIRED)
include_directories (${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})

file(GLOB SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
file(GLOB HEADERS RELATIVE ${CMAKE_SOURCE_DIR}/include *.?pp *.h)

list(REMOVE_ITEM SRCS "main.cpp")

include(CheckIncludeFiles)

add_library(anpi STATIC ${SRCS} ${HEADERS})
add_executable(proyecto3 main.cpp)
target_link_libraries(proyecto3 anpi ${Boost_INCLUDE_DIRS} python2.7 boost_program_options)

