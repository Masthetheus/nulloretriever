CC = gcc
CFLAGS = -Wall -Wextra -O3 -Iinclude
LDFLAGS =

# Makes sure bin exists
$(shell mkdir -p bin)

# --------------------------------------------------------------------
# EXPLICIT RULES (made to allow future c scripts in the same makefile)
# --------------------------------------------------------------------

bin/null_order: src/nullomer_order.c | bin
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

# Target to compile a single file following the rule above
null_order: bin/null_order

clean:
	rm -f bin/*

help:
	@echo "Useful commands:"
	@echo "  make null_order        - Compiles nullomer order script."
	@echo "  make clean             - Cleans the bin folder from its files."
