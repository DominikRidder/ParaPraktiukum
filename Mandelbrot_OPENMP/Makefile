
#--------------------------------
# insert name of your executable
#--------------------------------
PROGRAM = mandelseq2

#-------------------------------------------------
# insert compiler: 
#        icc          (serial)
#        mpicc        (mpi)
#        icc -qopenmp (openmp)
#-------------------------------------------------
CC  = icc


OBJ     = ${PROGRAM:=.o}

ARCH = _LINUX
INCLUDE_DIR = -I/homea/hpclab/train120/kurs/Tools
LIBS =  -L/homea/hpclab/train120/kurs/Tools -lppmwrite -lm
USER_CFLAGS =  $(INCLUDE_DIR) -D${ARCH} -g 

USER_FFLAGS =  $(INCLUDE_DIR) -r8  -g
srcdir=.

VPATH=.:$(srcdir)

$(PROGRAM): $(OBJ) 
	$(CC) $(USER_CFLAGS) -o $(PROGRAM) $(OBJ) $(LIBS) -qopenmp

.c.o:
	$(CC) $(USER_CFLAGS) -c $<

clean: 
	rm -f *.o $(PROGRAM)


