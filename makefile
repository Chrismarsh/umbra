#Makefile for libmaw
#Chris Marsh 2011

#MATLAB setup
TMW_ROOT=/usr/local/Matlab2011a
Arch=glnxa64

CC=icpc 
CFLAGS=-O3 -openmp -std=c++0x
SRC=src

#brittle, fix
LDFLAGS=-L$(TMW_ROOT)/bin/$(Arch) -leng -lmx -L/usr/lib -llapack -Wl,-rpath-link,$(TMW_ROOT)/bin/$(Arch) -L../libmaw -lmaw
INCLUDES=-I$(TMW_ROOT)/extern/include -I../armadillo-2.0.1/include -I../libmaw

all: main

main: main.o bounding_rect.o triangle.o triangulation.o
	$(CC) $(CFLAGS) main.o bounding_rect.o triangle.o triangulation.o $(LDFLAGS)  -o umbra

main.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRC)/main.cpp
	
bounding_rect.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRC)/bounding_rect.cpp
	
triangle.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRC)/triangle.cpp
	
triangulation.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SRC)/triangulation.cpp	
	
	
	
clean:
	rm -rf *o umbra
