#include "main.h"

/*
FUNCTION DESCRIPTION: The main (top level) function for executing the code. Contains the initialization, time loop, and cleanup routines.
LAST EDIT: Re-arranged position of diagnostics to make sure everything is printed at the same time level. Re-arranged position of creation and collisions.
EDIT DATE: 1/11/22
NAME OF EDITOR: Tasman Powis
*/
int main(int argc, char **argv)
{

	//Counters
	int k;
	int my_rank;

	//Time markers
	double t0, t1;

	//Initializing allspecies, allcreation, grid, timer, output and control structures
	struct AllSpecies allspecies;
	struct AllCreation allcreation;
	struct AllCollisions allcollisions;
	struct Timing timer;
	struct Output output;
	struct Control control;


	//Initializing MPI and getting process rank
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	//Getting the input file name
	GetInputFileName(&argc,&argv,&control);

	//Reading the input file for details of setup for the timing, control, grid, BCs and output. Counts number of species.
	ProcessingInput(&control, &allspecies, &allcreation, &allcollisions, &timer, &output);

	//Initializing the Pseudo-Random-Number-Generator seed from /dev/urandom.
	InitializeGlobalSeed(&control);

	//Initializing the timer structure.
	InitializeTimer(&timer);

	//Allocating memory for the different species structures
	AllocateAllSpecies(&allspecies);

	//Allocating memory for the different species creation routines
	AllocateAllCreation(&allcreation);

	//Allocating memory for the species collision modules
	AllocateAllCollisions(&allcollisions);

	//Reading the input file, reading individual species information and storing details for initialization
	ProcessingSpeciesInput(&control, &allspecies);

	//Reading the input file, reading individual species creation information and storing details for initialization
	ProcessingCreationInput(&control, &allcreation, &allspecies, &timer);

	//Reading the input file, reading individual species collision module information and storing details for initialization
	ProcessingCollisionInput(&control, &allcollisions, &allspecies);

	//Initializing all of the species from the stored data
	InitializeAllSpecies(&control, &allspecies, &allcreation, &timer);

	//Initializing all of the particles from the "plasma" creation routine and preparing data for source routines
	InitializeAllCreation(&control, &allspecies, &allcreation, &timer);

	//Initializing all of the collision modules
	InitializeAllCollisions(&allspecies, &allcollisions);

	//Mapping the collision indices to the species structures
	MapCollisionsToSpecies(&allspecies, &allcollisions);

    //Computing the maximum probabilities of collision for the null-collision algorithm
	ComputeMaxProbabilities(&allspecies, &allcollisions, timer.dt);

	//Printing simulation information to the log file
	timer.flog = fopen("pic.log","w");

    //Synchronizing tasks prior to the time loop
    MPI_Barrier(MPI_COMM_WORLD);


	//MPI_Finalize(); exit(0);


	//Copying initialized data to the device
	MemoryToDeviceAllParticles(&allspecies);
	MemoryToDeviceAllCollisions(&allcollisions);


	//Time loop
	for (k = timer.start_step; k < timer.Nsteps + 1; k++)
	{

		//Beginning timer for total runtime
		t0 = MPI_Wtime();

		//Passing the iteration number to the time structure
		timer.step = k;

		//Printing the step to screen and log file
		if (my_rank == 0) {PrintStep(&timer, k);}


		//Counting and printing total number of particles
		{
			int i;
			long long np = 0, npt;

			for (i = 0; i < allspecies.num_species; i++)
			{
				np += (long long) allspecies.species[i].Npart_local;
			}

			MPI_Reduce(&np,&npt, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

			if (my_rank == 0) {printf(", Number of Particles = %lld",npt);}

		}



		// //Pushing the particles
		PushAllParticles(&allspecies, &timer);


		//Handling all collision algorithms - currently a memory increase bug
		CollideAllParticles(&allspecies,&allcollisions,&timer);


		//Blocking until each process has finished this iteratation
		Barrier(&timer);

		//Total time
		timer.time += MPI_Wtime() - t0;
		//exit(0);
    }

	//Returning final data to host following completion of the run
	MemoryDeleteDeviceAllCollisions(&allcollisions);
	MemoryDeleteDeviceAllParticles(&allspecies);



	//Communicating and printing timers
	CommunicateTiming(&timer);
	PrintTiming(&timer);


	//Clearing memory
	FinalizeControl(&control);
	FinalizeAllCollisions(&allcollisions);
	FinalizeAllCreation(&allcreation, &control);
	FinalizeAllSpecies(&allspecies);
	fclose(timer.flog);

    MPI_Finalize();

    return 0;

}
