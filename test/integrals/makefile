CXXFLAGS = -g -Wall --pedantic -std=c++2a
DIR = -I ../  -I ../../

compile: main.cpp
	g++ main.cpp -o main $(CXXFLAGS) $(DIR)

execute: compile
	./main
	rm main

test: test.cpp
	g++ test.cpp -o test $(CXXFLAGS) $(DIR)
	./test 
	rm test
