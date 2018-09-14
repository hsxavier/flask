
############################
### Variable definitions ###
############################

# Option for debugging with valgrind:
DEBUG = -g -O0
# COMP = g++ $(DEBUG)

# Compiler and flask directory structure:
COMP = g++
BIN  = ../bin
OBJ  = ../obj
SRC  = ../src

# Healpix location:
HEALDIR  =  /home/skems/prog/Healpix_3.31
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib

# Location of other headers and libraries (if not in a default location):
# Uncomment and change path if required:
#CXXFITS  = -I/path/to/cfitsio/headers
#LDFITS   = -L/path/to/cfitsio/libraries
#CXXGSL   = -I/path/to/gsl/headers/folder
#LDGSL    = -L/path/to/gsl/libraries/folder

# Other compiler flags:
CXXOMP   = -fopenmp
LIBS     = -lgsl -lgslcblas -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lfftpack -lc_utils -lgomp

# Options for gina:
ifeq ($(HOSTNAME), gina.iag.usp.br)
HEALDIR  =  /sto/home/hsxavier/prog/Healpix_3.31
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib
CXXFITS  = -I/sto/home/hsxavier/prog/cfitsio
LDFITS   = -L/sto/home/hsxavier/prog/cfitsio
CXXGSL   = -I/sto/home/hsxavier/prog/gsl/include
LDGSL    = -L/sto/home/hsxavier/prog/gsl/lib
endif

# Options for GridUNESP computer:
ifeq ($(HOSTNAME), access.grid.unesp.br)
HEALDIR  = /home/rosenfel/FLASK/Healpix_3.31
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib
CXXFITS  = -I/home/rosenfel/FLASK/cfitsio/include
LDFITS   = -L/home/rosenfel/FLASK/cfitsio/lib
endif


# Options for Splinter computer:
ifeq ($(HOSTNAME), splinter-login.local)
#HEALDIR  = /home/hsxavier/Healpix_3.11
HEALDIR  = /share/splinter/cosmos/modules/may_2017/install_dir/Healpix_3.31
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib
CXXFITS  = -I/usr/include/cfitsio
endif

# Options for Alphacrucis computer:
ifeq ($(HOSTNAME), alphacrucis)
HEALDIR  = /sto/home/hsxavier/Healpix_3.20
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib
CXXFITS  = -I/sto/home/hsxavier/cfitsio
LDFITS   = -L/sto/home/hsxavier/cfitsio
endif

# Get Healpix version and check if it includes their own complex number class:
HEALVERSION := $(shell grep -o '[0-9]*\.[0-9]*' $(HEALDIR)/Version)
USEXCOMPLEX := $(shell echo "$(HEALVERSION) < 3.30" | bc) # Set to 1 if Healpix version < 3.30, 0 otherwise. 
# Get code version:
FLASKCOMMIT := $(shell git rev-parse --short HEAD)
ifeq ($(FLASKCOMMIT), )
FLASKCOMMIT = missing
endif

# Dynamic defines (#define that are set during compilation time, by passing instruction to compiler):
VCONTROL   = -D FLASKCOMMIT=\"$(FLASKCOMMIT)\"
DEFCOMPLEX = -D USEXCOMPLEX=$(USEXCOMPLEX)  


#############################
### Compiling intructions ###
#############################

# General instructions:

all: $(BIN)/flask $(BIN)/GenStarMask $(BIN)/XiCalculator $(BIN)/Dens2KappaCls

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*

test:


# Executables:

$(BIN)/Dens2KappaCls: $(OBJ)/Dens2KappaCls.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/Cosmology.o $(OBJ)/FieldsDatabase.o $(OBJ)/Integral.o \
			$(OBJ)/interpol.o $(OBJ)/Spline.o
	$(COMP) $(OBJ)/Dens2KappaCls.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
		$(OBJ)/Cosmology.o $(OBJ)/FieldsDatabase.o $(OBJ)/Integral.o \
		$(OBJ)/interpol.o $(OBJ)/Spline.o -o $@ -lgomp

$(BIN)/XiCalculator: $(OBJ)/XiCalculator.o $(OBJ)/Utilities.o 
	$(COMP) $(OBJ)/XiCalculator.o $(OBJ)/Utilities.o -o $@ -lgomp

$(BIN)/GenStarMask: $(OBJ)/GenStarMask.o
	$(COMP) $(LDHEAL) $(LDFITS) $(LDGSL) $(OBJ)/GenStarMask.o -o $@ $(LIBS)

$(BIN)/flask: $(OBJ)/flask.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/gsl_aux.o $(OBJ)/s2kit10_naive.o $(OBJ)/interpol.o \
			$(OBJ)/flask_aux.o $(OBJ)/Cosmology.o $(OBJ)/Integral.o \
			$(OBJ)/GeneralOutput.o $(OBJ)/SelectionFunc.o $(OBJ)/RegularizeCov.o \
			$(OBJ)/ClProcessing.o $(OBJ)/Maximize.o $(OBJ)/fitsfunctions.o \
			$(OBJ)/lognormal.o $(OBJ)/FieldsDatabase.o $(OBJ)/Spline.o
	$(COMP) $(LDHEAL) $(LDFITS) $(LDGSL) $(OBJ)/flask.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/gsl_aux.o $(OBJ)/s2kit10_naive.o $(OBJ)/interpol.o \
			$(OBJ)/flask_aux.o $(OBJ)/Cosmology.o $(OBJ)/Integral.o \
			$(OBJ)/GeneralOutput.o $(OBJ)/SelectionFunc.o $(OBJ)/RegularizeCov.o \
			$(OBJ)/ClProcessing.o $(OBJ)/Maximize.o $(OBJ)/fitsfunctions.o \
			$(OBJ)/lognormal.o $(OBJ)/FieldsDatabase.o $(OBJ)/Spline.o -o $@ $(LIBS)

# Objects:

$(OBJ)/Dens2KappaCls.o: $(SRC)/Dens2KappaCls.cpp $(SRC)/Utilities.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Cosmology.hpp $(SRC)/FieldsDatabase.hpp $(SRC)/definitions.hpp
	$(COMP) -c $(SRC)/Dens2KappaCls.cpp -o $@ $(CXXOMP)

$(OBJ)/XiCalculator.o: $(SRC)/XiCalculator.cpp Utilities.hpp
	$(COMP) $(CXXHEAL) -c $(SRC)/XiCalculator.cpp -o $@ $(CXXOMP)

$(OBJ)/GenStarMask.o: $(SRC)/GenStarMask.cpp
	$(COMP) $(CXXHEAL) $(CXXGSL) -c $(SRC)/GenStarMask.cpp -o $@ $(CXXOMP)

$(OBJ)/flask.o: $(SRC)/flask.cpp $(SRC)/flask_aux.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Utilities.hpp $(SRC)/gsl_aux.hpp $(SRC)/s2kit10_naive.hpp $(SRC)/Cosmology.hpp \
			$(SRC)/GeneralOutput.hpp $(SRC)/fitsfunctions.hpp $(SRC)/definitions.hpp \
			$(SRC)/lognormal.hpp $(SRC)/FieldsDatabase.hpp $(SRC)/SelectionFunc.hpp
	$(COMP) $(VCONTROL) $(DEFCOMPLEX) $(CXXHEAL) $(CXXGSL) -c $(SRC)/flask.cpp -o $@ $(CXXOMP)

$(OBJ)/Spline.o:	$(SRC)/Spline.cpp $(SRC)/Spline.hpp $(SRC)/Utilities.hpp $(SRC)/FieldsDatabase.hpp
	$(COMP) -c $(SRC)/Spline.cpp -o $@

$(OBJ)/FieldsDatabase.o: $(SRC)/FieldsDatabase.cpp $(SRC)/FieldsDatabase.hpp $(SRC)/Utilities.hpp
	 $(COMP) -c $(SRC)/FieldsDatabase.cpp -o $@ $(CXXOMP)

$(OBJ)/ClProcessing.o: $(SRC)/ClProcessing.cpp $(SRC)/ClProcessing.hpp $(SRC)/Utilities.hpp \
			$(SRC)/s2kit10_naive.hpp $(SRC)/GeneralOutput.hpp $(SRC)/RegularizeCov.hpp \
			$(SRC)/flask_aux.hpp $(SRC)/FieldsDatabase.hpp $(SRC)/Spline.hpp \
			$(SRC)/fitsfunctions.hpp
	$(COMP)  $(CXXHEAL) $(CXXGSL) -c $(SRC)/ClProcessing.cpp -o $@ $(CXXOMP)

$(OBJ)/flask_aux.o: $(SRC)/flask_aux.cpp $(SRC)/flask_aux.hpp \
			$(SRC)/gsl_aux.hpp $(SRC)/Cosmology.hpp $(SRC)/definitions.hpp \
			$(SRC)/fitsfunctions.hpp $(SRC)/GeneralOutput.hpp $(SRC)/Utilities.hpp \
			$(SRC)/lognormal.hpp
	$(COMP) $(DEFCOMPLEX) $(CXXHEAL) $(CXXGSL) -c $(SRC)/flask_aux.cpp -o $@ $(CXXOMP)

$(OBJ)/fitsfunctions.o: $(SRC)/fitsfunctions.cpp $(SRC)/fitsfunctions.hpp $(SRC)/definitions.hpp
	sed 's|#define HEALPIX_DATA .*|#define HEALPIX_DATA \"$(HEALDATA)\"|g' $(SRC)/fitsfunctions.cpp > $(OBJ)/fitsfunctions-made.cpp
	$(COMP) $(CXXHEAL) $(CXXFITS) $(CXXGSL) -c $(OBJ)/fitsfunctions-made.cpp -o $@

$(OBJ)/lognormal.o: $(SRC)/lognormal.cpp $(SRC)/lognormal.hpp
	$(COMP) -c $(SRC)/lognormal.cpp -o $@

$(OBJ)/RegularizeCov.o: $(SRC)/RegularizeCov.cpp $(SRC)/RegularizeCov.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Utilities.hpp
	$(COMP) $(CXXGSL) -c $(SRC)/RegularizeCov.cpp -o $@

$(OBJ)/SelectionFunc.o: $(SRC)/SelectionFunc.cpp $(SRC)/SelectionFunc.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Utilities.hpp $(SRC)/flask_aux.hpp $(SRC)/definitions.hpp
	$(COMP) $(CXXHEAL) $(CXXGSL) -c $(SRC)/SelectionFunc.cpp -o $@

$(OBJ)/GeneralOutput.o: $(SRC)/GeneralOutput.cpp $(SRC)/GeneralOutput.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/flask_aux.hpp $(SRC)/definitions.hpp $(SRC)/FieldsDatabase.hpp
	$(COMP) $(DEFCOMPLEX) $(CXXHEAL) $(CXXGSL) -c $(SRC)/GeneralOutput.cpp -o $@

$(OBJ)/Cosmology.o: $(SRC)/Cosmology.cpp $(SRC)/Cosmology.hpp $(SRC)/Utilities.hpp $(SRC)/Integral.hpp \
	$(SRC)/ParameterList.hpp $(SRC)/interpol.h
	$(COMP) -c $(SRC)/Cosmology.cpp -o $@

$(OBJ)/Integral.o: $(SRC)/Integral.cpp $(SRC)/Integral.hpp $(SRC)/Utilities.hpp $(SRC)/Cosmology.hpp
	$(COMP) -c $(SRC)/Integral.cpp -o $@

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@

$(OBJ)/gsl_aux.o: $(SRC)/gsl_aux.cpp $(SRC)/gsl_aux.hpp
	$(COMP) $(CXXGSL) -c $(SRC)/gsl_aux.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

$(OBJ)/interpol.o: $(SRC)/interpol.c $(SRC)/interpol.h
	$(COMP) -c $(SRC)/interpol.c -o $@

$(OBJ)/Maximize.o: $(SRC)/Maximize.c $(SRC)/Maximize.h $(SRC)/interpol.h
	$(COMP) -c $(SRC)/Maximize.c -o $@

$(OBJ)/s2kit10_naive.o: $(SRC)/s2kit10_naive.cpp $(SRC)/s2kit10_naive.hpp $(SRC)/interpol.h $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/s2kit10_naive.cpp -o $@ $(CXXOMP)
