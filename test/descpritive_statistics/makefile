CXXFLAGS = -g -Wall --pedantic -std=c++2a
DIR = -I ../  -I ../../

compile: main.cpp
	g++ main.cpp -o main $(CXXFLAGS) $(DIR)

clear: 	
	rm main
	
execute: main
	./main
	rm main