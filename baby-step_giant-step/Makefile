# Generic Makefile for PARI programs -- amd64 running linux (x86-64 kernel) 64-bit version
#
#  This file was created by Configure. Any change made to it will be
#  lost when Configure is run.
#
#
#

# change this TARGET to compile your own programs
TARGET = essai
SHELL  = /bin/sh

DBGFLAGS   = -g -Wall
CFLAGS     = -O3 -Wall -fno-strict-aliasing  -fomit-frame-pointer
#CFLAGS    = $(DBGFLAGS)

CC         = /usr/bin/gcc
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -fno-strict-aliasing -Wl,--export-dynamic  -fomit-frame-pointer 

RUNPTH     = -Wl,-rpath "/usr/local/lib"
LIBS       = -lm -L/usr/local/lib -lpari

RM = rm -f


OBJS = $(TARGET).o
ALL = $(TARGET)-dyn

all: $(ALL)

dyn: $(TARGET)-dyn

$(TARGET)-dyn: $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $<  $(RUNPTH) $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<
clean:
	-$(RM) *.o $(ALL)
