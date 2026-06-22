# Build plotMET inside the mucoll singularity environment (see ../README.md).
# Needs ROOT (root-config on PATH) and LCIO. The spack env sets $LCIO; if it's
# not set, try:  export LCIO=$(dirname $(dirname $(which anajob)))

CXX      := g++
CXXFLAGS := -O2 -std=c++17 $(shell root-config --cflags) -I$(LCIO)/include
LDFLAGS  := $(shell root-config --libs) -L$(LCIO)/lib -L$(LCIO)/lib64 -llcio

plotMET: plotMET.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f plotMET

.PHONY: clean
