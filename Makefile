#COMPILER OPTIONS
CFLAGS := -g -Wall -O3
CC := mpicc $(CFLAGS)

# If the first argument is "run"...
ifeq (run,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

#DIRECTORIES
BASE_DIR := .
SRC_DIR := $(BASE_DIR)/src
HEADER_DIR := $(BASE_DIR)/include
BUILD_DIR := $(BASE_DIR)/build
BIN_DIR := $(BASE_DIR)/bin

#FILES
BIN := $(BIN_DIR)/imagenew
TEMPLATE := $(BIN_DIR)/image
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES))

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADER_DIR)/%.h
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $<

$(BIN): $(OBJ_FILES)
	@echo " Linking..."
	@mkdir -p $(BIN_DIR)
	$(CC) $^ $(BASE_DIR)/main.c -o $@

all: $(BIN)
	@echo " $(BIN) ready."

run: $(BIN)
	mpirun -n 4 ./$(BIN) $(RUN_ARGS)

test: test.py
	python test.py

clean:
	@echo " Cleaning..."
	rm -rf $(BIN_DIR) $(BUILD_DIR)

.PHONY: clean all run
