# Makefile for Plateau Border Program 'plat'

CC = gcc
RM = rm
MV = mv
OBJ = xplat.o plinit.o plequil.o pldiffus.o plminor.o plio.o mgl.o
SRC = xplat.c plinit.c plequil.c pldiffus.c plminor.c plio.c mgl.c
P_OBJ = plat.o plinit.o plequil.o pldiffus.o plminor.o plio.o mgl.o
P_SRC = plat.c plinit.c plequil.c pldiffus.c plminor.c plio.c mgl.c
CFLAGS = -DDEBUG -g #-DFIXEDSEED

# for X11
GRAPHLIB= -L/usr/lib -lXm -lXt -lX11

plat: $(P_OBJ) include.h machine.h
	$(CC) $(CFLAGS) $(P_OBJ) -lm -o plat

# xplat: $(OBJ) include.h machine.h
# 	$(CC) $(CFLAGS) $(OBJ) -L/usr/lib/X11R5 -lXm -lXt -lX11 -lm -o xplat
xplat: $(OBJ) include.h machine.h
	$(CC) $(CFLAGS) $(OBJ) $(GRAPHLIB) -lm -o xplat

debug: $(OBJ) include.h machine.h
	$(CC) -g $(CFLAGS) $(OBJ) $(GRAPHLIB) -lm -o plat

.c.o:
	$(CC) -c $(CFLAGS) $<

$(SRC): include.h

include.h: machine.h

clean:
	rm $(OBJ) $(P_OBJ) xplat plat
