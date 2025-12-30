# Compiler
FC := gfortran

# Build directory for .o and .mod
MODDIR := build

# Executable directory
EXEDIR := executables

# Fortran flags:
# -J writes module files to MODDIR
# -I tells compiler where to search for modules
FFLAGS := -J$(MODDIR) -I$(MODDIR)
LIBS = -llapacke -llapack -lblas

# Machine code HSX --> Human readable HSX
PROGRAM_hsx := make_hsx_readable
SRC_hsx := src/hsx_m.f90 src/hsx2hsx.f90
OBJS_hsx := $(SRC_hsx:src/%.f90=$(MODDIR)/%.o)

# Human readable HSX --> tight-binding parameters
PROGRAM_hsx2tb := hsx2tb
SRC_hsx2tb := src/hsx2tb.f90
OBJS_hsx2tb := $(SRC_hsx2tb:src/%.f90=$(MODDIR)/%.o)

# Default target
all: $(PROGRAM_hsx) $(PROGRAM_hsx2tb) $(PROGRAM_tb2bands)

# Ensure executable directory exists
$(EXEDIR):
	mkdir -p $@

$(PROGRAM_hsx): $(OBJS_hsx) | $(EXEDIR)
	$(FC) $(OBJS_hsx) -o $(EXEDIR)/$@

$(PROGRAM_hsx2tb): $(OBJS_hsx2tb) | $(EXEDIR)
	$(FC) $(OBJS_hsx2tb) $(LIBS) -o $(EXEDIR)/$@

# Ensure build directory exists
$(MODDIR):
	mkdir -p $@

# Compile rule: src/.../*.f90 --> build/*.o
$(MODDIR)/%.o: src/%.f90 | $(MODDIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(MODDIR) $(EXEDIR)

.PHONY: all clean
