# Variable definitions:
PROGRAM_NAME = corrlnfields
COMP = g++
BIN = ../bin
OBJ = ../obj
SRC = ../src

CXXOMP = -fopenmp
CXXHEAL = -I/home/skems/prog/Healpix_3.11/src/cxx/generic_gcc/include
LDHEAL = -L/home/skems/prog/Healpix_3.11/src/cxx/generic_gcc/lib
LIBS = -lgsl -lgslcblas -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lfftpack -lc_utils -lgomp

# General instructions:
all: $(BIN)/$(PROGRAM_NAME)

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*

$(BIN)/$(PROGRAM_NAME): $(OBJ)/corrlnfields.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/gsl_aux.o $(OBJ)/s2kit10_naive.o $(OBJ)/interpol.o \
			$(OBJ)/corrlnfields_aux.o $(OBJ)/Cosmology.o $(OBJ)/Integral.o \
			$(OBJ)/GeneralOutput.o $(OBJ)/SelectionFunc.o $(OBJ)/RegularizeCov.o \
			$(OBJ)/ClProcessing.o
	$(COMP) $(LDHEAL) $(OBJ)/corrlnfields.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/gsl_aux.o $(OBJ)/s2kit10_naive.o $(OBJ)/interpol.o \
			$(OBJ)/corrlnfields_aux.o $(OBJ)/Cosmology.o $(OBJ)/Integral.o \
			$(OBJ)/GeneralOutput.o $(OBJ)/SelectionFunc.o $(OBJ)/RegularizeCov.o \
			$(OBJ)/ClProcessing.o -o $@ $(LIBS)


$(OBJ)/corrlnfields.o: $(SRC)/corrlnfields.cpp $(SRC)/corrlnfields_aux.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Utilities.hpp $(SRC)/gsl_aux.hpp $(SRC)/s2kit10_naive.hpp $(SRC)/Cosmology.hpp \
			$(SRC)/GeneralOutput.hpp
	$(COMP) $(CXXHEAL) -c $(SRC)/corrlnfields.cpp -o $@

$(OBJ)/ClProcessing.o: $(SRC)/ClProcessing.cpp $(SRC)/ClProcessing.hpp $(SRC)/Utilities.hpp \
			$(SRC)/s2kit10_naive.hpp $(SRC)/GeneralOutput.hpp $(SRC)/RegularizeCov.hpp \
			$(SRC)/corrlnfields_aux.hpp
	$(COMP)  $(CXXHEAL) -c $(SRC)/ClProcessing.cpp -o $@ $(CXXOMP)

$(OBJ)/corrlnfields_aux.o: $(SRC)/corrlnfields_aux.cpp $(SRC)/corrlnfields_aux.hpp \
			$(SRC)/gsl_aux.hpp $(SRC)/Cosmology.hpp
	$(COMP) $(CXXHEAL) -c $(SRC)/corrlnfields_aux.cpp -o $@

$(OBJ)/RegularizeCov.o: $(SRC)/RegularizeCov.cpp $(SRC)/RegularizeCov.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/RegularizeCov.cpp -o $@

$(OBJ)/SelectionFunc.o: $(SRC)/SelectionFunc.cpp $(SRC)/SelectionFunc.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/Utilities.hpp $(SRC)/corrlnfields_aux.hpp
	$(COMP) $(CXXHEAL) -c $(SRC)/SelectionFunc.cpp -o $@

$(OBJ)/GeneralOutput.o: $(SRC)/GeneralOutput.cpp $(SRC)/GeneralOutput.hpp $(SRC)/ParameterList.hpp \
			$(SRC)/corrlnfields_aux.hpp
	$(COMP) $(CXXHEAL) -c $(SRC)/GeneralOutput.cpp -o $@

$(OBJ)/Cosmology.o: $(SRC)/Cosmology.cpp $(SRC)/Cosmology.hpp $(SRC)/Utilities.hpp $(SRC)/Integral.hpp \
	$(SRC)/ParameterList.hpp $(SRC)/interpol.h
	$(COMP) -c $(SRC)/Cosmology.cpp -o $@

$(OBJ)/Integral.o: $(SRC)/Integral.cpp $(SRC)/Integral.hpp $(SRC)/Utilities.hpp $(SRC)/Cosmology.hpp
	$(COMP) -c $(SRC)/Integral.cpp -o $@

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@

$(OBJ)/gsl_aux.o: $(SRC)/gsl_aux.cpp $(SRC)/gsl_aux.hpp
	$(COMP) -c $(SRC)/gsl_aux.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

$(OBJ)/interpol.o: $(SRC)/interpol.c $(SRC)/interpol.h
	$(COMP) -c $(SRC)/interpol.c -o $@

$(OBJ)/s2kit10_naive.o: $(SRC)/s2kit10_naive.cpp $(SRC)/s2kit10_naive.hpp $(SRC)/interpol.h $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/s2kit10_naive.cpp -o $@ $(CXXOMP)

