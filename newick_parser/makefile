CC = gcc
OBJS = seqMain.o seqUtil.o Newickform.o

newicktree: ${OBJS}
	${CC} -o newicktree ${OBJS}
	
libnewick.a: ${OBJS}
	ar -cvq libnewick.a ${OBJS}
	
seqMain.o: seqMain.c
seqUtil.o: seqUtil.c seqUtil.h
Newickform.o: Newickform.c

clean:
	-rm *.o newicktree
