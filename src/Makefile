OFILES = main-direct.o tools.o xmalloc.o
TARGET = dynamics

AOFILES = pqr2flat.o xmalloc.o
ATARGET = pqr2flat 

BOFILES = txt2bin.o xmalloc.o
BTARGET = txt2bin 

CC = mpicc
CFLAGS = -Wall -Wextra  -O2

$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

$(ATARGET): $(AOFILES)
	$(CC) $(AOFILES) -o $@

$(BTARGET): $(BOFILES)
	$(CC) $(BOFILES) -o $@

clean:
	rm -f $(OFILES) $(TARGET) 
	      $(AOFILES) $(ATARGET)
	      $(BOFILES) $(BTARGET)


tools.o: tools.c array.h xmalloc.h
xmalloc.o: xmalloc.c xmalloc.h
main-direct.o: main-direct.c array.h xmalloc.h tools.h
pqr2flat.o: pqr2flat.c array.h xmalloc.h
txt2bin.o: txt2bin.c array.h xmalloc.h
