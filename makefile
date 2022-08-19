CXXFLAGS = -g -Wall --pedantic -std=c++17

compile: main.cpp
	g++ $(CXXFLAGS) $< -o main
	rm main

execute: main.cpp	
	g++ $(CXXFLAGS) $< -o main
	./main
	rm main
