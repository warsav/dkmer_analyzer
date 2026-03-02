CUDA_PATH ?= /usr/local/cuda

CC = gcc
NVCC = nvcc

CFLAGS = -std=c11 -O3 -Wall -Wextra -pthread
NVCCFLAGS = -O3 -arch=sm_86 
LDFLAGS = -L$(CUDA_PATH)/lib64 -Wl,-rpath,$(CUDA_PATH)/lib64 -lm -lcudart -lstdc++

TARGET = dkmer_analyzer

SRC_C = main.c genome.c bitpacker.c analysis.c checkpoint.c
SRC_CU = gpu_bridge.cu
OBJ = $(SRC_C:.c=.o) $(SRC_CU:.cu=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
