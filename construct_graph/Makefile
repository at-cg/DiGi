# Installs the tools that scripts/liftover.sh needs, without conda,
# into construct_graph/bin/ and construct_graph/extra/.
#
#   make graph    install (from the DiGi root; or `make` inside construct_graph/)
#   make clean    remove bin/ and extra/
#
# Requires Linux x86-64, GCC 8 or newer, Python 3.8 or newer, wget, tar, gzip,
# bzip2 and network access. Options:
#   make graph PYTHON=/path/to/python3
#   make graph CA_BUNDLE=/path/to/cacert.pem    if a certificate cannot be verified

# =========================================================
# Settings (each can be overridden on the command line)
# =========================================================

CXX       ?= g++
CXXFLAGS  ?= -O3 -std=c++17
PYTHON    ?= python3
CA_BUNDLE ?=

# Tool versions. Nextflow 26.10.7 on Java 17 is the tested pair.
# Kent tools: the current build needs glibc 2.34; v479 also runs on older systems.
CROSSMAP_VERSION ?= 0.7.4
MINIMAP2_VERSION ?= 2.28
LAST_VERSION     ?= 1654
NF_LO_VERSION    ?= 1.8.8
KENT_VERSION     ?= v479
JDK_VERSION      ?= 17.0.21.1+2
NEXTFLOW_VERSION ?= 26.10.7

# UCSC Kent tools used by nf-LO (liftover.sh runs it with --no_maf)
KENT_TOOLS = faSize faToTwoBit twoBitInfo faSplit liftUp axtChain chainAntiRepeat \
             chainMergeSort chainSplit chainPreNet chainNet netSyntenic netChainSubset \
             chainStitchId

# ======================
# Download resources
# ======================

WGET = wget -nv $(if $(CA_BUNDLE),--ca-certificate=$(CA_BUNDLE))

UCSC_URL        = https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.$(KENT_VERSION)
MINIMAP2_DIR    = minimap2-$(MINIMAP2_VERSION)_x64-linux
MINIMAP2_URL    = https://github.com/lh3/minimap2/releases/download/v$(MINIMAP2_VERSION)/$(MINIMAP2_DIR).tar.bz2
MAF_CONVERT_URL = https://gitlab.com/mcfrith/last/-/raw/$(LAST_VERSION)/bin/maf-convert
NEXTFLOW_URL    = https://github.com/nextflow-io/nextflow/releases/download/v$(NEXTFLOW_VERSION)/nextflow-$(NEXTFLOW_VERSION)-dist

# Java (Eclipse Temurin). The address is built from JDK_VERSION: 17.0.21.1+2 gives
# temurin17-binaries/.../jdk-17.0.21.1%2B2/OpenJDK17U-jdk_x64_linux_hotspot_17.0.21.1_2.tar.gz
JDK_MAJOR = $(firstword $(subst ., ,$(JDK_VERSION)))
JDK_URL   = https://github.com/adoptium/temurin$(JDK_MAJOR)-binaries/releases/download/jdk-$(subst +,%2B,$(JDK_VERSION))/OpenJDK$(JDK_MAJOR)U-jdk_x64_linux_hotspot_$(subst +,_,$(JDK_VERSION)).tar.gz

# ==================================
# Everything that `make` installs
# ==================================

all: $(addprefix bin/,$(KENT_TOOLS)) bin/minimap2 bin/maf-convert bin/CrossMap nf-lo \
     bin/junction_filter bin/gfa_const bin/node_kmer_sampler

# =============================================================================
# DiGi's C++ tools, compiled from src/: the junction filter, the GFA writer
# and the node k-mer sampler
# =============================================================================

bin/junction_filter bin/gfa_const bin/node_kmer_sampler: bin/%: src/%.cpp src/graph_file.hpp | bin
	$(CXX) $(CXXFLAGS) $< -o $@

# =============================================================================
# Downloaded programs. Single files and Java are saved under a temporary name
# first, so that a failed download is not mistaken for an installed tool.
# =============================================================================

# UCSC Kent tools
$(addprefix bin/,$(KENT_TOOLS)): | bin
	$(WGET) -O $@.part $(UCSC_URL)/$(@F) && chmod +x $@.part && mv $@.part $@

# minimap2, with k8 and paftools.js
bin/minimap2: | bin
	$(WGET) -O - $(MINIMAP2_URL) | tar -xj -C bin --strip-components=1 \
	  $(MINIMAP2_DIR)/minimap2 $(MINIMAP2_DIR)/k8 $(MINIMAP2_DIR)/paftools.js

# maf-convert, from LAST
bin/maf-convert: | bin
	$(WGET) -O $@.part $(MAF_CONVERT_URL) && chmod +x $@.part && mv $@.part $@

# =============================================================================
# CrossMap, installed with pip into a Python environment in extra/crossmap.
# bin/python and bin/python3 point to it, for maf-convert and nf-LO's scripts.
# =============================================================================

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

# =============================================================================
# Java, Nextflow and the nf-LO pipeline
# =============================================================================

# Java, into extra/java
extra/java/bin/java:
	rm -rf extra/java extra/java.part && mkdir -p extra/java.part
	$(WGET) -O - $(JDK_URL) | tar -xz -C extra/java.part --strip-components=1
	mv extra/java.part extra/java

# Nextflow. The "-dist" file includes its runtime, so nothing is fetched on first run
bin/nextflow: | bin
	$(WGET) -O $@.part $(NEXTFLOW_URL) && chmod +x $@.part && mv $@.part $@

# nf-LO, downloaded by Nextflow into extra/nextflow using the Java above
NXF_ENV = JAVA_HOME=$(CURDIR)/extra/java NXF_HOME=$(CURDIR)/extra/nextflow

nf-lo: bin/nextflow extra/java/bin/java
	$(NXF_ENV) bin/nextflow pull evotools/nf-LO -r $(NF_LO_VERSION)

# =============================================================================
# Folders and cleanup
# =============================================================================

bin:
	mkdir -p bin

clean:
	rm -rf bin extra

.PHONY: all nf-lo clean
