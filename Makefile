EXE  = experiment
MPI_EXE = mpi_experiment
SRC_DIR = src
OBJ_DIR = obj
RESULT_DIR=results
INC_DIR = include
LDFLAGS += -lboost_system -lboost_filesystem -lm -lIOH

SRC:= $(shell find src/ ! -name "experiment.cc" ! -name "mpi_experiment.cc" -name "*.cc")
OBJ = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
INC = -I $(INC_DIR)

CC      = g++
CFLAGS  = -Wall -std=c++20 -g# -O2

.PHONY: all
all: $(OBJ_DIR) $(RESULT_DIR) $(EXE)

.PHONY: mpi
mpi: $(OBJ_DIR) $(RESULT_DIR) $(MPI_EXE)

$(EXE): $(OBJ) $(OBJ_DIR)/experiment.o
	${CC} ${CFLAGS} -o $(EXE) $(OBJ) $(OBJ_DIR)/experiment.o ${LDFLAGS}

$(MPI_EXE): $(OBJ) $(OBJ_DIR)/mpi_experiment.o
	mpiCC ${CFLAGS} -o $(EXE) $(OBJ) $(OBJ_DIR)/mpi_experiment.o ${LDFLAGS}

$(OBJ_DIR)/mpi_experiment.o: $(SRC_DIR)/mpi_experiment.cc
	mpiCC -c $(CFLAGS) $(INC) -o $(OBJ_DIR)/mpi_experiment.o $(SRC_DIR)/mpi_experiment.cc

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/*
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(OBJ_DIR):
	mkdir $(OBJ_DIR)
$(RESULT_DIR):
	mkdir $(RESULT_DIR)

.PHONY:  clean
clean:
	rm -f $(OBJ_DIR)/*.o $(EXE) $(MPI_EXE)

.PHONY: cleanall
cleanall:
	rm -rf $(OBJ_DIR) $(EXE) $(MPI_EXE) results
