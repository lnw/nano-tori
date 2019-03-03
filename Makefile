# Copyright (c) 2019, Lukas Wirz
# All rights reserved.

# This file is part of 'nano-tori' which is released under the BSD-2-clause license.
# See file LICENSE in this project.


#CXX=clang++
CXX=g++

GSL_INCLUDES_L=-lgsl -lgslcblas
FLAGS=-O2 -std=c++14 -Wunused -Wshadow -g # -Wall

HEADERS=geometry2.hh geometry3.hh auxiliary.hh rational.hh
OBJECTS=geometry2.o geometry3.o
OBJECTS_P=$(patsubst %.o, build/%.o, $(OBJECTS))

build/%.o: %.cc $(HEADERS) Makefile
	$(CXX) $(FLAGS) -c $< -o $@

all: botch_optimise-torus-full torus-extra torus-full torus-grid torus-simple

torus-simple: Makefile torus-simple.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) torus-simple.cc $(OBJECTS_P) -o $@

torus-full: Makefile torus-full.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) torus-full.cc $(OBJECTS_P) -o $@

torus-grid: Makefile torus-grid.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) torus-grid.cc $(OBJECTS_P) -o $@

torus-extra: Makefile torus-extra.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) torus-extra.cc $(OBJECTS_P) -o $@

botch_optimise-torus-full: Makefile botch_optimise-torus-full.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) botch_optimise-torus-full.cc $(OBJECTS_P) -o $@ $(GSL_INCLUDES_L)


test: Makefile test.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test.cc $(OBJECTS_P) $(GSL_INCLUDES_L) -o $@
	./test

clean:
	rm -f output/*.xyz output/*.coord test debug final
	rm -f torus*.coord torus*.xyz

distclean: clean
	rm -f torus-simple torus-full torus-extra torus-grid
	rm -f botch_optimise-torus-full
	rm -f ccb.*

