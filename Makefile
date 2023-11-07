CC=g++  
INCLUDE=-I/share/home/liangwy/biotool/zlib
CFLAGS=-std=c++20 -g -O3 -Wall -Wextra -fopenmp -march=native
TARGETS=saluscbmap
OBJS=sequence_batch.o edlib.o align.o seed.o index.o main.o
LDFLAGS=-L/share/home/liangwy/biotool/zlib/lib -lm -lz 

all: $(TARGETS)

$(TARGETS): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJS) $(TARGETS)

format:
	clang-format -i *.cpp *.h