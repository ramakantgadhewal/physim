CXXFLAGS = -g -Wall --pedantic -std=c++17

compile: compile.cpp
	g++ $(CXXFLAGS) $< -o compile
	rm compile

execute: main.cpp	
	g++ $(CXXFLAGS) $< -o main
	./main
	rm main
