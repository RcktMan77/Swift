# Compiler
FC = gfortran


# Compiler flags
FFLAGS = -O2 \
		 -fopenmp \
		 -g \
		 -fbacktrace \
		 -fcheck=all \
		 -ffpe-trap=invalid,zero,overflow,underflow


# Linker flags
LDFLAGS = -fopenmp


# TecIO library paths
TECIO_DIR=/home/daviszs/deps/tecio/2024R1
TECIO_INC=$(TECIO_DIR)/include
TECIO_LIB=$(TECIO_DIR)/lib
TECIO_LIBS=-L$(TECIO_LIB) -ltecio -pthread -lstdc++


# Source files
SRC = kind_defs.f90 \
	  grid_types.f90 \
	  flow_types.f90 \
	  namelist_definitions.f90 \
	  initialization.f90 \
	  gradients.f90 \
	  fluxes.f90 \
	  grid_properties.f90 \
	  dual_mesh.f90 \
	  write_output.f90 \
	  main.f90


# Object files
OBJ = $(SRC:.f90=.o)

# Executable name
EXEC= rans_solve

# Default target
all: $(EXEC)

# Rule to link object files to create the executable
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(TECIO_LIBS)

# Rule to compile .f90 files into .o files
%.o: %.f90
	$(FC) $(FFLAGS) -I$(TECIO_INC) -c $<

# Dependencies for object files
kind_defs.o: kind_defs.f90

grid_types.o: grid_types.f90 \
	          kind_defs.mod

flow_types.o: flow_types.f90 \
	          kind_defs.mod

namelist_definitions.o: namelist_definitions.f90 \
						kind_defs.mod

initialization.o: initialization.f90 \
	              kind_defs.mod \
				  grid_types.mod \
				  flow_types.mod \
				  namelist_definitions.mod

gradients.o: gradients.f90 \
	         kind_defs.mod \
			 grid_types.mod \
			 flow_types.mod

fluxes.o: fluxes.f90 \
	      kind_defs.mod \
		  grid_types.mod \
		  flow_types.mod \
		  namelist_definitions.mod \
		  initialization.mod

grid_properties.o: grid_properties.f90 \
				   kind_defs.mod \
				   grid_types.mod

dual_mesh.o: dual_mesh.f90 \
	         kind_defs.mod \
			 grid_types.mod

write_output.o: write_output.f90 \
	            kind_defs.mod \
				grid_types.mod \
				flow_types.mod \
				$(TECIO_INC)/tecio.f90

main.o: main.f90 \
	    grid_types.mod \
		flow_types.mod \
		namelist_definitions.mod \
		initialization.mod \
		gradients.mod \
		fluxes.mod \
	    grid_properties.mod \
		dual_mesh.mod \
		write_output.mod

# Clean target to remove object files and executable
clean:
	rm -f $(OBJ) $(EXEC)

# Phony targets
.PHONY: all clean

