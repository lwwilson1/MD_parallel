OFILES = main.o routines.o tools.o xmalloc.o \
	 tree-pc.o tree-pc-yuk.o partition.o treecomm.o
TARGET = dynamics

AOFILES = pqr2flat.o xmalloc.o
ATARGET = pqr2flat 

BOFILES = txt2bin.o xmalloc.o
BTARGET = txt2bin 

CC = mpicc
CFLAGS = -Wall -Wextra -O3

$(TARGET): $(OFILES)
	$(CC) $(CFLAGS) $(OFILES) -o $@

$(ATARGET): $(AOFILES)
	$(CC) $(CFLAGS) $(AOFILES) -o $@

$(BTARGET): $(BOFILES)
	$(CC) $(CFLAGS) $(BOFILES) -o $@

clean:
	rm -f $(OFILES) $(TARGET) $(AOFILES) $(ATARGET) $(BOFILES) $(BTARGET)


tools.o: tools.c array.h xmalloc.h structs.h
xmalloc.o: xmalloc.c xmalloc.h
main.o: main.c array.h xmalloc.h treecode.h treecomm.h tools.h structs.h globtreevars.h
routines.o: routines.c array.h xmalloc.h tools.h structs.h

partition.o: partition.c partition.h tools.h structs.h
tree-pc.o: tree-pc.c treecode.h xmalloc.h array.h partition.h tools.h \
	       structs.h globtreevars.h 
tree-pc-yuk.o: tree-pc-yuk.c treecode.h xmalloc.h array.h partition.h tools.h \
	       structs.h globtreevars.h
treecomm.o: treecomm.c treecomm.h xmalloc.h array.h treecode.h \
           structs.h globtreevars.h

pqr2flat.o: pqr2flat.c array.h xmalloc.h
txt2bin.o: txt2bin.c array.h xmalloc.h
