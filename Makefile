# Build DiGi's C++ tools into bin/
#   make         k-mer set sampler, IQP solver and KMC. Needs Linux x86-64 and
#                GCC 8 or newer. Uses the Gurobi install in GUROBI_HOME if set,
#                otherwise downloads Gurobi into extra/. Downloads KMC into bin/.
#   make clean   remove bin/

CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17

GUROBI_VERSION ?= 13.0.0
GUROBI_URL = https://packages.gurobi.com/$(basename $(GUROBI_VERSION))/gurobi$(GUROBI_VERSION)_linux64.tar.gz
GUROBI_HOME ?= $(CURDIR)/extra/gurobi$(subst .,,$(GUROBI_VERSION))/linux64

KMC_VERSION ?= 3.2.4
KMC_URL = https://github.com/refresh-bio/KMC/releases/download/v$(KMC_VERSION)/KMC$(KMC_VERSION).linux.x64.tar.gz

all: bin/kmer_set_sampling bin/dip_iqp bin/kmc

bin/kmc: | bin
	curl -fL $(KMC_URL) | tar -xz -C bin --strip-components=1 bin/kmc bin/kmc_tools

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

.PHONY: all clean
