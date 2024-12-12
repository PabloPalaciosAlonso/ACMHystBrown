# Compiler and flags
CXX = g++
CXXFLAGS = -O3 -fopenmp -I../src -I../external/FittingAlgorithms/src
LDFLAGS = -lgsl -lgslcblas -lm

# Sources and executables
SOURCES = fitUncoatedParticles.cpp fitDistributionOfUncoatedParticles.cpp
EXECUTABLES = fitUncoatedParticles fitDistributionOfUncoatedParticles

# Default target
all: $(EXECUTABLES)

# Pattern rule for building executables
%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Clean target to remove the generated executables
clean:
	@rm -f $(EXECUTABLES)
