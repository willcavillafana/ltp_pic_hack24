# LTP-PIC - Collision Kernel

This source code is a "mini-app" from the LTP-PIC code which encompasses the particle phase space update (Boris algorithm) and Monte-Carlo Collision (MCC) module of the particle-in-cell algorithm.

## Compiling

### Local system with GNU

The Makefile should by default include appropriate flags to compile with GNU. Simply run,
```
make clean
make
```


### Traverse with V100 GPUs

To compile the code, comment out line 12 `COPTS` and uncomment line 21 `COPTS` in the `Makefile`.

The standard modules being used are:
```
module load nvhpc/22.5
module load openmpi/nvhpc-22.5/4.1.3/64
module load cudatoolkit/11.7
```

Then run,
```
make clean
make
```

## Running

To run the code locally simple use the MPI execution command:
```
mpiexec -n <n> ./pic <input_file_name>.dat
```

Where `<n>` is the number of MPI tasks. Note that this code is embarassingly parallel, each task performs exactly the same operation on different randomly generated particles. There is no inter-task communication, except to sum and report the total number of particles being simulated.

`<input_file_name>.dat` is the input file which must have the `.dat` extension.

An example input file, called `example.dat` is in the top directory. This file includes self explanatory labels.

### Running on Traverse with slurm

A batch script titled `btrav` is included in the top directory and should submit a job to the Traverse cluster requesting 10 minutes on 1 GPU.


## Proposed Hackathon Tasks
1. Optimize MCC loop (`CollideParticlesNew` function in `collisions.c`)
   1. Investigate improved algorithms to create the random list of particles indice which will be collided in the main collision kernel.
   2. Consider storing collision particles into a buffer for continguous memory access.
   3. Consider moving collision particles to CPU for improved kernel performance (will need Grace-Hopper)
   4. Reduce kernel size to avoid register spilling
   5. Ivestigate collision performance with different phase space data precision (i.e. double vs single floating point precision)
   6. Implement sub-functions to clean up the code and improve readability
1. Modify the Boris algorithm (`PushParticles` in `particles.c`) to be able to utilize TF32 tensor cores.
2. Ideas for implementing binary collision algorithms?
3. Ideas fo mized fluid/kinetic particle methods?

