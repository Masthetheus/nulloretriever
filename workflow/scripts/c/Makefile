.POSIX:
.SUFFIXES:
CC = gcc
CFLAGS = -O3 -Wall -Wextra
SRC = main.c kmer_encoding.c kmer_extraction.c kmerio.c
OBJ = $(patsubst %.c, build/%.o, $(SRC))

all: kmer_extractor

kmer_extractor: $(OBJ)
	$(CC) -o $@ $^
build/%.o: src/%.c | build
	$(CC) -c $(CFLAGS) $< -o $@
build:
	mkdir -p build/
clean:
	rm -f kmer_extractor build/*.o
