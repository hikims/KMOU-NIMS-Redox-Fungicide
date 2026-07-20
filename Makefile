
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3
INCLUDES ?= -Isrc

# TBB_PREFIX is optional.  When set, its include/ and lib/ are added to the
# search paths, e.g.
#     make all TBB_PREFIX=/opt/homebrew/opt/tbb
TBB_PREFIX ?=
ifneq ($(strip $(TBB_PREFIX)),)
INCLUDES += -I$(TBB_PREFIX)/include
LDFLAGS  += -L$(TBB_PREFIX)/lib
endif
TBB_LIBS ?= -ltbb
OUTPUT_ROOT ?= ./results/
MACHINE_LABEL ?= local
COMMON_SRC = src/functions.cpp src/initial.cpp
BIN_DIR = bin

TARGETS = \
	main_fig2_table3_4_dose_time_maps \
	main_fig3_fig4_verification \
	main_fig5_spatial_dynamics \
	main_fig6_front_speed \
	main_supp_figS1_threshold \
	main_supp_figS2_S3_sensitivity_curves \
	main_table_S1_S2_convergence

.PHONY: all clean check-src
all: check-src $(TARGETS)

check-src:
	@test -f src/functions.cpp || (echo "Missing src/functions.cpp. Copy the shared solver source files into src/." && exit 1)
	@test -f src/initial.cpp || (echo "Missing src/initial.cpp. Copy the shared solver source files into src/." && exit 1)
	@test -f src/functions.hpp || (echo "Missing src/functions.hpp. Copy the shared solver header files into src/." && exit 1)
	@test -f src/initial.hpp || (echo "Missing src/initial.hpp. Copy the shared solver header files into src/." && exit 1)
	@test -f src/parameters.h || (echo "Missing src/parameters.h. Copy parameters.h into src/." && exit 1)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

%: drivers/%.cpp $(COMMON_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DOUTPUT_ROOT=\"$(OUTPUT_ROOT)\" -DMACHINE_LABEL=\"$(MACHINE_LABEL)\" $< $(COMMON_SRC) $(LDFLAGS) $(TBB_LIBS) -o $(BIN_DIR)/$@

clean:
	rm -rf $(BIN_DIR)
