# ---------------------------------------------------------------------
# File: src/mdPp/patterns.mk
#
# This makefile fragment contains the pattern rules used to compile
# all C++ sources files in the src/mdPp directory tree. The src/mdPp 
# directory contains all the source code for the MdPp C++ namespace. 
#
# This file should be included by all makefile files in the mdPp/ 
# directory. This file should be included in other makefiles after 
# inclusion of the main configuration file, src/config.mk.
#-----------------------------------------------------------------------
# Makefile fragment includes

# Build configuration files
# Define *_DEFS macro definitions, paths to libraries and executables
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mdPp/config.mk

# Source file lists
# Include recipes to build library targets
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/mdPp/sources.mk

#-----------------------------------------------------------------------
# Makefile variable definitions

# Lists of all required C preprocessor macro definitions 
UTIL_ADEF=$(UTIL_DEFS) 
SIMP_ADEF=$(UTIL_DEFS) $(SIMP_DEFS)
MDPP_ADEF=$(UTIL_DEFS) $(SIMP_DEFS) $(MDPP_DEFS)

# Lists of dependencies on configuration files 
# These are added to lists of dependencies generated by $(MAKEDEP)
UTIL_CFGS= -A$(BLD_DIR)/config.mk
UTIL_CFGS+= -A$(BLD_DIR)/util/config.mk
SIMP_CFGS= $(UTIL_CFGS)
SIMP_CFGS+= -A$(BLD_DIR)/simp/config.mk
MDPP_CFGS= $(SIMP_CFGS)
MDPP_CFGS+= -A$(BLD_DIR)/mdPp/config.mk

# All libraries needed by files in src/mdPp
LIBS=$(mdPp_LIB) $(simp_LIB) $(util_LIB)

#-----------------------------------------------------------------------
# Pattern rules

# Rule to compile all *.cpp class source files in src/mdPp
$(BLD_DIR)/mdPp/%.o: $(SRC_DIR)/mdPp/%.cpp
	$(CXX) $(INCLUDES) $(MDPP_ADEF) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(MDPP_ADEF) $(CXXFLAGS) $(MDPP_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Rule to compile all *.cpp class source files in src/simp
$(BLD_DIR)/simp/%.o: $(SRC_DIR)/simp/%.cpp
	$(CXX) $(INCLUDES) $(SIMP_ADEF) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(SIMP_ADEF) $(CXXFLAGS) $(SIMP_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Rule to compile all *.cpp class source files in src/util
$(BLD_DIR)/util/%.o: $(SRC_DIR)/util/%.cpp
	$(CXX) $(INCLUDES) $(UTIL_ADEF) $(CXXFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(UTIL_ADEF) $(CXXFLAGS) $(UTIL_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Rule to compile all *.cc test programs in src/mdPp
$(BLD_DIR)/mdPp/tests/%.o: $(SRC_DIR)/mdPp/tests/%.cc
	$(CXX) $(INCLUDES) $(MDPP_ADEF) $(TESTFLAGS) -c -o $@ $<
ifdef MAKEDEP
	$(MAKEDEP) $(INCLUDES) $(MDPP_ADEF) $(TESTFLAGS) $(MDPP_CFGS) -S$(SRC_DIR) -B$(BLD_DIR) $<
endif

# Rule to link all *.cc test program executables in src/mdPp
$(BLD_DIR)/mdPp/tests/%: $(BLD_DIR)/mdPp/tests/%.o $(LIBS)
	$(CXX) -o $@ $< $(LIBS) $(LDFLAGS)
