########################################################################
# Compiler and external dependences
########################################################################
CC        = mpic++
#CC 	   = CC #Wrapper for NERSC

########################################################################
# Compiling and linking options
########################################################################

#GCC
COPTS     = -g -Wno-unused-result -O3 -march=native -fopenmp -Drestrict=__restrict__

#INTEL Stellar
#COPTS     = -g -O3 -xhost -Wno-unknown-pragmas -restrict -qopenmp

#NVIDIA Traverse - CPU
#COPTS     = -g -w -fast -mp

#NVIDIA Traverse - GPU
#COPTS     = -g -w -fast -mp -acc=gpu #-Minfo=all

#CRAY Perlmutter - GPU
#COPTS     = -g -w -fast -mp -acc=gpu -target-accel=nvidia80 #-Minfo=all

CINCLUDES =
CDEFS     =
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)


LINKOPTS  = $(COPTS)
LIBS      = -lm 
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c

.c.o:
	$(CC) $(CFLAGS) -c $<

########################################################################
# List of all programs to be compiled
########################################################################
ALLPROGS = pic

all: $(ALLPROGS)

default: pic


########################################################################
# Example 5
########################################################################
pic: main.o particles.o collisions.o numeric.o timing.o input.o des.o desprng.o
	$(CC) -o $@ $^ $(LFLAGS)	


########################################################################
# Clean up
########################################################################
clean:
	rm -f *.o
distclean: clean
	rm -f $(ALLPROGS) $(ALLPROGS:=*~)
