#COMPILER OPTIONS
CFLAGS := -g -Wall -O3 -cc=icc 
CC := mpicc $(CFLAGS)

#DIRECTORIES
BASE_DIR := .
SRC_DIR := $(BASE_DIR)/src
HEADER_DIR := $(BASE_DIR)/include
BUILD_DIR := $(BASE_DIR)/build
SERIAL_DIR := $(BASE_DIR)/serial

#FILES
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES))

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADER_DIR)/%.h
	@mkdir -p $(BUILD_DIR)
	$(CC) -D $(mode) -c -o $@ $<

clean:
	@echo " Cleaning..."
	@rm -rf serial imagenew imagenew_timing imagenew_average_pixel imagenew_timing_intervals $(BUILD_DIR)

clean_build:
	@echo " Cleaning build files..."
	@rm -rf $(BUILD_DIR)

serial:
	@echo " Creating serial executable..."
	make $(BUILD_DIR)/pgmio.o
	@$(CC) $^ $(SERIAL_DIR)/serial.c -o $@

imagenew:
	@echo " Creating imagenew executable..."
	make $(OBJ_FILES) mode=PROD
	@$(CC) $(OBJ_FILES) $(BASE_DIR)/main.c -o $@

imagenew_timing:
	@echo " Creating imagenew_timing executable..."
	make $(OBJ_FILES) mode=TIMING_TEST
	@$(CC) $(OBJ_FILES) $(BASE_DIR)/main.c -o $@

imagenew_average_pixel:
	@echo " Creating imagenew_average_pixel executable..."
	make $(OBJ_FILES) mode=AVERAGE_PIXEL_TEST
	@$(CC) $(OBJ_FILES) $(BASE_DIR)/main.c -o $@

imagenew_timing_intervals:
	@echo " Creating imagenew_timing_intervals executable..."
	make $(OBJ_FILES) mode=TIMING_WITH_INTERVALS_TEST
	@$(CC) $(OBJ_FILES) $(BASE_DIR)/main.c -o $@

.PHONY: clean clean_build all test serial imagenew imagenew_timing imagenew_average_pixel imagenew_timing_intervals
