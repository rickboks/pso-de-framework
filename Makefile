EXE  = experiment
MPI_EXE = mpi_experiment
SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include
LDFLAGS += -lm

SRC = $(wildcard $(SRC_DIR)/*.cc)
OBJ = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
INC = -I $(INC_DIR)/

CC      = g++
CFLAGS  = -Wall -std=c++14 -g

COCOFLAGS = -g -ggdb -std=c89 -pedantic -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion

.PHONY: all
all: $(EXE)

.PHONY: mpi
mpi: $(MPI_EXE)

$(EXE): $(OBJ) $(OBJ_DIR)/coco.o $(OBJ_DIR)/experiment.o
	${CC} ${CFLAGS} -o $(EXE) $(OBJ) $(OBJ_DIR)/coco.o $(OBJ_DIR)/experiment.o ${LDFLAGS}

$(MPI_EXE): $(OBJ) $(OBJ_DIR)/coco.o $(OBJ_DIR)/mpi.o
	mpiCC ${CFLAGS} -o $(MPI_EXE) $(OBJ) $(OBJ_DIR)/coco.o $(OBJ_DIR)/mpi_experiment.o $(LDFLAGS)
	
$(OBJ_DIR)/mpi.o: $(INC_DIR)/coco.h $(SRC_DIR)/coco.c $(SRC_DIR)/mpi_experiment.c | $(OBJ_DIR)
	mpiCC -c $(CFLAGS) $(INC) -o $(OBJ_DIR)/mpi_experiment.o $(SRC_DIR)/mpi_experiment.c

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/* | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(OBJ_DIR)/coco.o: $(INC_DIR)/coco.h $(SRC_DIR)/coco.c | $(OBJ_DIR)
	gcc -c $(COCOFLAGS) $(INC) -o $(OBJ_DIR)/coco.o $(SRC_DIR)/coco.c

$(OBJ_DIR)/experiment.o: $(INC_DIR)/coco.h $(SRC_DIR)/coco.c $(SRC_DIR)/experiment.c | $(OBJ_DIR)
	${CC} -c $(CFLAGS) $(INC) -o $(OBJ_DIR)/experiment.o $(SRC_DIR)/experiment.c

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

.PHONY:  clean
clean:
	rm -f $(OBJ_DIR)/*.o $(EXE) $(MPI_EXE)

