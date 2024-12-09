CXX = g++
CXXFLAGS = -O3 -fopenmp -I../ACMHystBrown/src -I../external/FittingAlgorithms/src

# Target to build the main executable directly
fitUncoatedParticles: fitUncoatedParticles.cpp
	$(CXX) $(CXXFLAGS) fitUncoatedParticles.cpp -o fitUncoatedParticles

# Clean target to remove the generated executable
clean:
	rm -f fitUncoatedParticles
