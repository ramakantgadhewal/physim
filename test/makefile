CXXFLAGS = -g -Wall --pedantic -std=c++17
DIR = -I ../  

compile: compile.cpp
	g++ $< -o compile $(CXXFLAGS) $(DIR)

execute: compile.cpp	
	g++ $(CXXFLAGS) $< -o compile
	./compile
	rm compile
