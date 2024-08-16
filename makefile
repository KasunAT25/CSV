# Makefile

# Define targets
TARGETS = sali_csv lipp_csv alex_csv data_prep

# Define source files for each target
SOURCES_X = sali_csv.cpp
SOURCES_Y = lipp_csv.cpp
SOURCES_Z = alex_csv.cpp
SOURCES_A = data_prep.cpp

# Define compiler and flags for sali
CXX_X = g++
CXXFLAGS_X = -fopenmp -std=c++17 -march=native -mpopcnt

#include the tbb location here
# ============================
LDFLAGS_X = -L/opt/intel/oneapi/tbb/2021.12/lib -ltbb -Wl,-rpath,/opt/intel/oneapi/tbb/2021.12/lib -Wl,--enable-new-dtags
CPPFLAGS_X = -I/opt/intel/oneapi/tbb/2021.12/include

# Define compiler and flags for lipp and alex
CXXFLAGS_Y = -std=c++17 -march=native -mpopcnt

# # Define compiler and flags for alex
# CXXFLAGS_Z = -std=c++17 -march=native -mpopcnt

CXXFLAGS_A = -std=c++17

# Build all
all: $(TARGETS)

sali_csv: $(SOURCES_X)
	$(CXX_X) $(CPPFLAGS_X) $(CXXFLAGS_X) $(SOURCES_X) -o $@ $(LDFLAGS_X)

lipp_csv: $(SOURCES_Y)
	$(CXX_X) $(CXXFLAGS_Y) $(SOURCES_Y) -o $@

alex_csv: $(SOURCES_Z)
	$(CXX_X) $(CXXFLAGS_Y) $(SOURCES_Z) -o $@

# Build dataprep
data_prep: $(SOURCES_A)
	$(CXX_X) $(CXXFLAGS_A) $(SOURCES_A) -o $@


# Phony targets
.PHONY: all clean

# Clean target
clean:
	rm -f $(TARGETS)


