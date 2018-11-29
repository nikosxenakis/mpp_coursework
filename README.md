# mpp_coursework B136013

## Prerequisites

* Compiler: [GNU g++](https://www.gnu.org/software/gcc/)
* Build tool: [GNU Make](https://www.gnu.org/software/make/)

## Usage

### Building

To build this project run the following

```
make all
```

### Cleaning

To clean this project run the following

```
make clean
```

### Running

To run this project run the following

```
mpirun -n $PROC_NO ./bin/imagenew $INPUT_PATH
```

Where $PROC_NO is the number of processes and $INPUT_PATH is the path of the input edge image.
The output new image is stored in the output/ folder.
The record of the average pixel value is stored in the data/ folder.

### Testing

To test the output run the following

```
python test.py
```
