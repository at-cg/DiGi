# Install the tools scripts/liftover.sh needs into construct_graph/bin/, without conda,
# and build junction_filter, gfa_const and node_kmer_sampler (src/) there with CXX (GCC 8 or newer)
#   make graph         from the DiGi root (or plain `make` inside construct_graph/)
#   make clean         remove construct_graph/bin and construct_graph/extra
#
# Downloads static Linux x86-64 binaries: the UCSC kent tools nf-LO chains with,
# minimap2 (with k8 and paftools.js) and LAST's maf-convert. Installs CrossMap with
# pip into a venv in construct_graph/extra/crossmap. Installs a Temurin JDK into construct_graph/extra/java
# and the standalone Nextflow launcher into construct_graph/bin, then fetches the pinned nf-LO
# release into NXF_HOME=construct_graph/extra/nextflow. Defaults are the tested pair
# (Nextflow 25.10.4 on Java 17); override with JDK_VERSION / NEXTFLOW_VERSION.
#
# Needs wget, tar, gzip, bzip2, network access, and Python 3.8 or newer
# for CrossMap. The default is `python3`; point at another one with
#   make graph PYTHON=/path/to/python3
# If wget cannot verify a server's certificate (old system CA bundle), pass a newer one:
#   make graph CA_BUNDLE=/path/to/cacert.pem

CXX              ?= g++
CXXFLAGS         ?= -O3 -std=c++17
PYTHON           ?= python3
CROSSMAP_VERSION ?= 0.7.4
MINIMAP2_VERSION ?= 2.28
LAST_VERSION     ?= 1654
NF_LO_VERSION    ?= 1.8.8
KENT_VERSION     ?= v479
JDK_VERSION      ?= 17.0.21.1+2
NEXTFLOW_VERSION ?= 25.10.4
CA_BUNDLE        ?=

WGET = wget -nv $(if $(CA_BUNDLE),--ca-certificate=$(CA_BUNDLE))

# Kent tools build: the current linux.x86_64 build needs glibc 2.34; v479 runs on
# older systems (e.g. glibc 2.28 on RHEL 8)
UCSC_URL        = https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.$(KENT_VERSION)
MINIMAP2_DIR    = minimap2-$(MINIMAP2_VERSION)_x64-linux
MINIMAP2_URL    = https://github.com/lh3/minimap2/releases/download/v$(MINIMAP2_VERSION)/$(MINIMAP2_DIR).tar.bz2
MAF_CONVERT_URL = https://gitlab.com/mcfrith/last/-/raw/$(LAST_VERSION)/bin/maf-convert

# Eclipse Temurin JDK release, e.g. 17.0.21.1+2 -> temurin17-binaries, tag
# jdk-17.0.21.1%2B2, file OpenJDK17U-jdk_x64_linux_hotspot_17.0.21.1_2.tar.gz
JDK_MAJOR = $(firstword $(subst ., ,$(JDK_VERSION)))
JDK_URL   = https://github.com/adoptium/temurin$(JDK_MAJOR)-binaries/releases/download/jdk-$(subst +,%2B,$(JDK_VERSION))/OpenJDK$(JDK_MAJOR)U-jdk_x64_linux_hotspot_$(subst +,_,$(JDK_VERSION)).tar.gz
NEXTFLOW_URL = https://github.com/nextflow-io/nextflow/releases/download/v$(NEXTFLOW_VERSION)/nextflow-$(NEXTFLOW_VERSION)-dist

# Run Nextflow with the local JDK and keep its files inside construct_graph/
NXF_ENV = JAVA_HOME=$(CURDIR)/extra/java NXF_HOME=$(CURDIR)/extra/nextflow

# Kent tools on nf-LO's minimap2 path (liftover.sh runs it with --no_maf)
KENT_TOOLS = faSize faToTwoBit twoBitInfo faSplit liftUp axtChain chainAntiRepeat \
             chainMergeSort chainSplit chainPreNet chainNet netSyntenic netChainSubset \
             chainStitchId

all: $(addprefix bin/,$(KENT_TOOLS)) bin/minimap2 bin/maf-convert bin/CrossMap nf-lo \
     bin/junction_filter bin/gfa_const bin/node_kmer_sampler

# The junction filter (step 4 of liftover.sh), then graph.gfa and the node k-mers
# of a graph folder (the last step)
bin/junction_filter bin/gfa_const bin/node_kmer_sampler: bin/%: src/%.cpp src/graph_file.hpp | bin
	$(CXX) $(CXXFLAGS) $< -o $@

$(addprefix bin/,$(KENT_TOOLS)): | bin
	$(WGET) -O $@.part $(UCSC_URL)/$(@F) && chmod +x $@.part && mv $@.part $@

bin/minimap2: | bin
	$(WGET) -O - $(MINIMAP2_URL) | tar -xj -C bin --strip-components=1 \
	  $(MINIMAP2_DIR)/minimap2 $(MINIMAP2_DIR)/k8 $(MINIMAP2_DIR)/paftools.js

bin/maf-convert: | bin
	$(WGET) -O $@.part $(MAF_CONVERT_URL) && chmod +x $@.part && mv $@.part $@

# bin/python and bin/python3 point at the same interpreter, for maf-convert and
# nf-LO's inline scripts (standard library only)
bin/CrossMap: | bin
	@$(PYTHON) -c 'import sys; sys.exit(sys.version_info < (3, 8))' 2>/dev/null || { \
	  echo "CrossMap needs Python 3.8 or newer; '$(PYTHON)' is $$($(PYTHON) --version 2>&1)." >&2; \
	  echo "Rerun with: make graph PYTHON=/path/to/python3.8+" >&2; exit 1; }
	$(PYTHON) -m venv extra/crossmap
	extra/crossmap/bin/pip install -q --upgrade pip
	extra/crossmap/bin/pip install -q --only-binary=:all: CrossMap==$(CROSSMAP_VERSION)
	ln -sfn ../extra/crossmap/bin/python3 bin/python3
	ln -sfn python3 bin/python
	ln -sfn ../extra/crossmap/bin/CrossMap $@

extra/java/bin/java:
	rm -rf extra/java extra/java.part && mkdir -p extra/java.part
	$(WGET) -O - $(JDK_URL) | tar -xz -C extra/java.part --strip-components=1
	mv extra/java.part extra/java

# The -dist launcher bundles Nextflow's runtime, so nothing is fetched on first run
bin/nextflow: | bin
	$(WGET) -O $@.part $(NEXTFLOW_URL) && chmod +x $@.part && mv $@.part $@

nf-lo: bin/nextflow extra/java/bin/java
	$(NXF_ENV) bin/nextflow pull evotools/nf-LO -r $(NF_LO_VERSION)

bin:
	mkdir -p bin

clean:
	rm -rf bin extra

.PHONY: all nf-lo clean
