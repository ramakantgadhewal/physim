CXXFLAGS = -g -Wall --pedantic -std=c++2a
DIR = -I ../  -I ../../

harmonic: harmonic.cpp
	g++ harmonic.cpp -o harmonic $(CXXFLAGS) $(DIR)

forced: forced.cpp
	g++ forced.cpp -o forced $(CXXFLAGS) $(DIR)

damped: damped.cpp
	g++ damped.cpp -o damped $(CXXFLAGS) $(DIR)

forced_damped: forced_damped.cpp
	g++ forced_damped.cpp -o forced_damped $(CXXFLAGS) $(DIR)


execute_harmonic: harmonic
	./harmonic

execute_forced: forced
	./forced

execute_damped: damped
	./damped

execute_forced_damped: forced_damped
	./forced_damped


clear: 	
	rm harmonic
	rm forced
	rm damped
	rm forced_damped
	
	