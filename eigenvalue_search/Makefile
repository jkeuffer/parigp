
DBGFLAGS   = -g -Wall
CFLAGS     = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer
#CFLAGS    = $(DBGFLAGS)

CC         = /usr/bin/gcc
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -fno-strict-aliasing -Wl,--export-dynamic  -fomit-frame-pointer 

RUNPTH     = -Wl,-rpath "/usr/local/lib"
LIBS       = -lm -L/usr/local/lib -lpari

RM = rm -f

ALL = essai-dyn

essai-dyn: essai.o bsgs.o exhaustive_search.o vp_mms.o
	$(LD) -o $@ essai.o bsgs.o exhaustive_search.o vp_mms.o $(LDFLAGS) $(RUNPTH) $(LIBS)

essai.o: essai.c bsgs.h exhaustive_search.h vp_mms.h
	$(CC) -c $(CFLAGS) $<
bsgs.o: bsgs.c bsgs.h
	$(CC) -c $(CFLAGS) $<
exhaustive_search.o: exhaustive_search.c exhaustive_search.h
	$(CC) -c $(CFLAGS) $<
vp_mms.o: vp_mms.c vp_mms.h
	$(CC) -c $(CFLAGS) $<

clean:
	-$(RM) *.o $(ALL) *~
