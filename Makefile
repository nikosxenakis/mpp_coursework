#COMPILER OPTIONS
CFLAGS := -g -Wall -o3

CC := mpicc $(CFLAGS)

#DIRECTORIES
BASE_DIR := .

SRC_DIR := $(BASE_DIR)/src

HEADER_DIR := $(BASE_DIR)/include

BUILD_DIR := $(BASE_DIR)/build

BIN_DIR := $(BASE_DIR)/bin

#FILES
BIN := $(BIN_DIR)/mpp_coursework

SRC_FILES := $(wildcard $(SRC_DIR)/*.c) $(BASE_DIR)/main.c

OBJ_FILES := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES))

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADER_DIR)/%.h
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $<

$(BIN): $(OBJ_FILES)
	@echo " Linking..."
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@

all: $(BIN)
	@echo " $(BIN) ready."

run: $(BIN)
	mpirun -n 4 ./$(BIN)

clean: $(BIN)
	@echo " Cleaning..."
	rm -rf $(BIN_DIR) $(BUILD_DIR)

.PHONY: clean
