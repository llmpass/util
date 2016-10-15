
CC = g++
COPT = -ggdb -I.
LIBS =

all : test1_complex test2_complex test3_complex

%.o: %.c *.h Makefile
	$(CC) $(COPT) -c -o $@ $< 

test1_complex : test1_complex.o complex.o  mesh.o mesh_base.o Makefile
	$(CC) $(COPT) -o test1_complex test1_complex.o complex.o mesh.o mesh_base.o $(LIBS)

test2_complex : test2_complex.o complex.o  mesh.o mesh_base.o Makefile
	$(CC) $(COPT) -o test2_complex test2_complex.o complex.o mesh.o mesh_base.o $(LIBS)

test3_complex : test3_complex.o complex.o  mesh.o mesh_base.o Makefile
	$(CC) $(COPT) -o test3_complex test3_complex.o complex.o mesh.o mesh_base.o $(LIBS)

clean : 
	rm *~ *.o test1_complex test2_complex test3_complex 

