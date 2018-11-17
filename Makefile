#COMPILER OPTIONS
CFLAGS := -g -Wall -O3 -cc=icc 
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

#FILES
BIN := imagenew
SERIAL := serial
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES))

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADER_DIR)/%.h
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $<

$(BIN): $(OBJ_FILES)
	@echo " Linking..."
	$(CC) $^ $(BASE_DIR)/main.c -o $@

$(SERIAL): $(BUILD_DIR)/pgmio.o
	@echo " Linking..."
	$(CC) $^ $(BASE_DIR)/serial.c -o $@

all: $(BIN) $(SERIAL)
	@echo " executables ready."

run: $(BIN)
	mpirun -n 4 ./$(BIN) $(RUN_ARGS)

test: test.py
	python test.py

clean:
	@echo " Cleaning..."
	rm -rf $(BIN) $(SERIAL) $(BUILD_DIR)

.PHONY: clean all run
