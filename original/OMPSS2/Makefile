include ../config.mk
include flags.mk

all:
	$(CC) $(CFLAGS) $(COMMON_FLAGS) -c ompss2_driver.c
	$(CC) $(CFLAGS) $(COMMON_FLAGS) -c ompss2_propagate.c
	$(CC) $(CFLAGS) $(COMMON_FLAGS) -c ompss2_insertsource.c

clean:
	rm -f *.o *.a
