# Build DiGi's C++ tools into bin/
#   make         KMC k-mer set sampler (no Gurobi needed)
#   make iqp     IQP solver; uses the Gurobi install in GUROBI_HOME if set,
#                otherwise downloads Gurobi into extra/ (Linux x86-64)
#   make clean   remove bin/

CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17

GUROBI_VERSION ?= 13.0.0
GUROBI_URL = https://packages.gurobi.com/$(basename $(GUROBI_VERSION))/gurobi$(GUROBI_VERSION)_linux64.tar.gz
GUROBI_HOME ?= $(CURDIR)/extra/gurobi$(subst .,,$(GUROBI_VERSION))/linux64

all: bin/kmer_set_sampling

iqp: bin/dip_iqp

bin/kmer_set_sampling: src/kmer_sampling/kmer_set_sampling.cpp | bin
	$(CXX) $(CXXFLAGS) $< -o $@

bin/dip_iqp: src/iqp/main.cpp src/graph/graph_file.hpp | bin
	@[ -d "$(GUROBI_HOME)" ] || [ "$(origin GUROBI_HOME)" != file ] \
	  || { mkdir -p extra && curl -fL $(GUROBI_URL) | tar -xz -C extra; }
	grb_lib=$$(basename "$$(ls "$(GUROBI_HOME)"/lib/libgurobi[0-9]*.so | grep -v _light | head -1)" .so); \
	test -n "$$grb_lib" && \
	$(CXX) $(CXXFLAGS) $< -o $@ -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib \
	  -Wl,-rpath,$(GUROBI_HOME)/lib -lgurobi_c++ -l$${grb_lib#lib} -lstdc++fs

bin:
	mkdir -p bin

clean:
	rm -rf bin

.PHONY: all iqp clean
