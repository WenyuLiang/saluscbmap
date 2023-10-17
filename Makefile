CC=g++  
INCLUDE=-I/share/home/liangwy/biotool/zlib
CFLAGS=-std=c++17 -g -O3 -Wall -Wextra -march=native
TARGETS=saluscbmap
OBJS=sequence_batch.o edlib.o align.o seed.o index.o main.o
LDFLAGS=-lm -lz 

all: $(TARGETS)

$(TARGETS): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJS) $(TARGETS)