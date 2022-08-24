CXXFLAGS = -g -Wall --pedantic -std=c++17
DIR = -I ../  -I ../../

harmonic: harmonic.cpp
	g++ harmonic.cpp -o harmonic $(CXXFLAGS) $(DIR)

clear: 	
	rm harmonic
	
execute_harmonic: harmonic
	./harmonic
