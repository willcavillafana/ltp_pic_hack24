#include "main.h"

/*
FUNCTION DESCRIPTION: Allocate the memory for all species structures. The allspecies structure contains an array of species structures. These structures are agnostic of the creation routines.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void AllocateAllSpecies(struct AllSpecies *allspecies)
{
	int nspecies = allspecies->num_species;

	allspecies->species = (struct Species *) malloc(nspecies * sizeof(struct Species));
}



/*
FUNCTION DESCRIPTION: Free memory for all species structures. This is a multi-level process.
LAST EDIT: Created
EDIT DATE: 29/10/22
NAME OF EDITOR: Tasman Powis
*/
void FinalizeAllSpecies(struct AllSpecies *allspecies)
{
	int i;
	int nspecies = allspecies->num_species;

	//Clearing memory from each species struct
	for (i = 0; i < nspecies; i++)
	{
		FinalizeSpecies(&allspecies->species[i]);
	}

	free(allspecies->species);
}



/*
FUNCTION DESCRIPTION: Allocate the memory for all creation structures. The creation structure contains arrays of plasma and source structures. Distinguish between plasma and source creation routines.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void AllocateAllCreation(struct AllCreation *creation)
{
	int nplasma = creation->num_plasma;

	creation->plasma = (struct CreatePlasma *) malloc(nplasma * sizeof(struct CreatePlasma));
}



/*
FUNCTION DESCRIPTION: Free memory for all creation structures.
LAST EDIT: Created
EDIT DATE: 29/10/22
NAME OF EDITOR: Tasman Powis
*/
void FinalizeAllCreation(struct AllCreation *creation, struct Control *control)
{
	int i;
	int nplasma = creation->num_plasma;

	//Clearing out each struct
	if (control->restart == 0) {
		for (i = 0; i < nplasma; i++)
		{
			free(creation->plasma[i].rncount);
			free(creation->plasma[i].process_data);
			free(creation->plasma[i].thread_data);
		}
	}


	//Freeing the array of structs
	free(creation->plasma);
}



/*
FUNCTION DESCRIPTION: Initialize all species structures for a 2D simulation.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void InitializeAllSpecies(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *creation, struct Timing *timer)
{
	int i;
	int nspecies = allspecies->num_species;

	for (i = 0; i < nspecies; i++)
	{
		allspecies->species[i].seed = global_seed + i*10;

		InitializeSpecies(&allspecies->species[i], creation, timer);
	}

}



/*
FUNCTION DESCRIPTION: Initialize all of the creation routines for a 2D simulation.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void InitializeAllCreation(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *creation, struct Timing *timer)
{
	int i, j;
	int nspecies = allspecies->num_species;
	int nplasma = creation->num_plasma;
	int snum, sindex;


	//Creating the initial plasma particles.
	for (i = 0; i < nplasma; i++)
	{
		creation->plasma[i].seed = global_seed + (nspecies + i)*10;

		sindex = creation->plasma[i].sindex;

		InitializePlasma(&allspecies->species[sindex], &creation->plasma[i], timer);
	}



	//Cycling through the species to count the global number of initial particles
	for (i = 0; i < nspecies; i++)
	{
		long long int npl, np;
		struct Species *species;
		species = &allspecies->species[i];

		npl = (long long int) species->Npart_local;

		MPI_Reduce(&npl, &np, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		species->Npart = np;

	}

}











/*
FUNCTION DESCRIPTION: Copy all collision structures from the host to an accelerator.
LAST EDIT: Updated
EDIT DATE: 1/8/2022
NAME OF EDITOR: Tasman Powis
*/
void MemoryToDeviceAllCollisions(struct AllCollisions *allcollisions)
{
	int i, nl;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;

	#pragma acc enter data copyin(allcollisions[:1])
	#pragma acc enter data copyin(allcollisions->elastic[:nelastic])
	#pragma acc enter data copyin(allcollisions->excite[:nexcite])
	#pragma acc enter data copyin(allcollisions->ionize[:nionize])
	#pragma acc enter data copyin(allcollisions->cxchange[:ncxchange])

	//Cycling through the collisions
	for (i = 0; i < nelastic; i++)
	{
	//struct ElasticCollision *elastic = &allcollisions->elastic[i];
	nl = allcollisions->elastic[i].nlines;

	#pragma acc enter data copyin(allcollisions->elastic[i].energy_vals[:nl])
	#pragma acc enter data copyin(allcollisions->elastic[i].xsection_vals[:nl])

	#pragma acc enter data copyin(allcollisions->elastic[i].rncount[0:1])
	#pragma acc enter data copyin(allcollisions->elastic[i].process_data[0:1])
	#pragma acc enter data copyin(allcollisions->elastic[i].thread_data[0:1])
	}

	for (i = 0; i < nexcite; i++)
	{
	//struct ExcitationCollision *excite = &allcollisions->excite[i];
	nl = allcollisions->excite[i].nlines;

	#pragma acc enter data copyin(allcollisions->excite[i].energy_vals[:nl])
	#pragma acc enter data copyin(allcollisions->excite[i].xsection_vals[:nl])

	#pragma acc enter data copyin(allcollisions->excite[i].rncount[0:1])
	#pragma acc enter data copyin(allcollisions->excite[i].process_data[0:1])
	#pragma acc enter data copyin(allcollisions->excite[i].thread_data[0:1])
	}

	for (i = 0; i < nionize; i++)
	{
	//struct IonizationCollision *ionize = &allcollisions->ionize[i];
	nl = allcollisions->ionize[i].nlines;

	#pragma acc enter data copyin(allcollisions->ionize[i].energy_vals[:nl])
	#pragma acc enter data copyin(allcollisions->ionize[i].xsection_vals[:nl])

	#pragma acc enter data copyin(allcollisions->ionize[i].rncount[0:1])
	#pragma acc enter data copyin(allcollisions->ionize[i].process_data[0:1])
	#pragma acc enter data copyin(allcollisions->ionize[i].thread_data[0:1])
	}

	for (i = 0; i < ncxchange; i++)
	{
	//struct ChargeExchangeCollision *cxchange = &allcollisions->cxchange[i];
	nl = allcollisions->cxchange[i].nlines;

	#pragma acc enter data copyin(allcollisions->cxchange[i].energy_vals[:nl])
	#pragma acc enter data copyin(allcollisions->cxchange[i].xsection_vals[:nl])

	#pragma acc enter data copyin(allcollisions->cxchange[i].rncount[0:1])
	#pragma acc enter data copyin(allcollisions->cxchange[i].process_data[0:1])
	#pragma acc enter data copyin(allcollisions->cxchange[i].thread_data[0:1])
	}
}




/*
FUNCTION DESCRIPTION: Copy all collision structures from the host to an accelerator.
LAST EDIT: Updated
EDIT DATE: 1/8/2022
NAME OF EDITOR: Tasman Powis
*/
void MemoryDeleteDeviceAllCollisions(struct AllCollisions *allcollisions)
{
	int i, nl;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;

	//Cycling through the collisions
	for (i = 0; i < nelastic; i++)
	{
	struct ElasticCollision *elastic = &allcollisions->elastic[i];
	nl = elastic->nlines;

	#pragma acc exit data delete(elastic->energy_vals[:nl])
	#pragma acc exit data delete(elastic->xsection_vals[:nl])

	#pragma acc exit data delete(elastic->rncount[0:1])
	#pragma acc exit data delete(elastic->process_data[0:1])
	#pragma acc exit data delete(elastic->thread_data[0:1])
	}

	for (i = 0; i < nexcite; i++)
	{
	struct ExcitationCollision *excite = &allcollisions->excite[i];
	nl = excite->nlines;

	#pragma acc exit data delete(excite->energy_vals[:nl])
	#pragma acc exit data delete(excite->xsection_vals[:nl])

	#pragma acc exit data delete(excite->rncount[0:1])
	#pragma acc exit data delete(excite->process_data[0:1])
	#pragma acc exit data delete(excite->thread_data[0:1])
	}

	for (i = 0; i < nionize; i++)
	{
	struct IonizationCollision *ionize = &allcollisions->ionize[i];
	nl = ionize->nlines;

	#pragma acc exit data delete(ionize->energy_vals[:nl])
	#pragma acc exit data delete(ionize->xsection_vals[:nl])

	#pragma acc exit data delete(ionize->rncount[0:1])
	#pragma acc exit data delete(ionize->process_data[0:1])
	#pragma acc exit data delete(ionize->thread_data[0:1])
	}

	for (i = 0; i < ncxchange; i++)
	{
	struct ChargeExchangeCollision *cxchange = &allcollisions->cxchange[i];
	nl = cxchange->nlines;

	#pragma acc exit data delete(cxchange->energy_vals[:nl])
	#pragma acc exit data delete(cxchange->xsection_vals[:nl])

	#pragma acc exit data delete(cxchange->rncount[0:1])
	#pragma acc exit data delete(cxchange->process_data[0:1])
	#pragma acc exit data delete(cxchange->thread_data[0:1])
	}

	#pragma acc exit data delete(allcollisions->elastic[:nelastic])
	#pragma acc exit data delete(allcollisions->excite[:nexcite])
	#pragma acc exit data delete(allcollisions->ionize[:nionize])
	#pragma acc exit data delete(allcollisions->cxchange[:ncxchange])
	#pragma acc exit data delete(allcollisions[:1])

}






/*
FUNCTION DESCRIPTION: Reallocate particle double array memory on host, and on device if used.
LAST EDIT: Created
EDIT DATE: 6/15/2020
NAME OF EDITOR: Tasman Powis
*/
double* ReallocDouble(double *x, int mold, int mnew)
{
	#ifdef _OPENACC
	double *y;

	//Creating the new memory chunk on the GPU for the double array
	y = (double*) acc_malloc(mnew*sizeof(double));

	//Copying values from the old memory chunk into the new memory chunk
	acc_memcpy_device(y,acc_deviceptr(x),mold*sizeof(double));

	//Deleting the old memory chunk on the GPU
	acc_delete(x,mold*sizeof(double));
	#endif
	//Reallocating the new memory chunk on the CPU
	x = (double *) realloc(x,mnew*sizeof(double));
	#ifdef _OPENACC
	//Mapping the new GPU memory chunk to the pointer on the CPU
	acc_map_data(x,y,mnew*sizeof(double));
	#endif
	//Returning the new pointer
	return x;

}




/*
FUNCTION DESCRIPTION: Reallocate particle custom ptype array memory on host, and on device if used.
LAST EDIT: Created
EDIT DATE: 11/24/2020
NAME OF EDITOR: Tasman Powis
*/
ptype* ReallocPtype(ptype *x, int mold, int mnew)
{
	#ifdef _OPENACC
	ptype *y;

	//Creating the new memory chunk on the GPU for the double array
	y = (ptype*) acc_malloc(mnew*sizeof(ptype));

	//Copying values from the old memory chunk into the new memory chunk
	acc_memcpy_device(y,acc_deviceptr(x),mold*sizeof(ptype));

	//Deleting the old memory chunk on the GPU
	acc_delete(x,mold*sizeof(ptype));
	#endif
	//Reallocating the new memory chunk on the CPU
	x = (ptype *) realloc(x,mnew*sizeof(ptype));
	#ifdef _OPENACC
	//Mapping the new GPU memory chunk to the pointer on the CPU
	acc_map_data(x,y,mnew*sizeof(ptype));
	#endif
	//Returning the new pointer
	return x;

}





/*
FUNCTION DESCRIPTION: Reallocate particle int array memory on host, and on device if used.
LAST EDIT: Created
EDIT DATE: 6/15/2020
NAME OF EDITOR: Tasman Powis
*/
int* ReallocInt(int *x, int mold, int mnew)
{
	#ifdef _OPENACC
	int *y;

	//Creating the new memory chunk on the GPU for the double array
	y = (int*) acc_malloc(mnew*sizeof(int));

	//Copying values from the old memory chunk into the new memory chunk
	acc_memcpy_device(y,acc_deviceptr(x),mold*sizeof(int));

	//Deleting the old memory chunk on the GPU
	acc_delete(x,mold*sizeof(int));
	#endif
	//Reallocating the new memory chunk on the CPU
	x = (int *) realloc(x,mnew*sizeof(int));
	#ifdef _OPENACC
	//Mapping the new GPU memory chunk to the pointer on the CPU
	acc_map_data(x,y,mnew*sizeof(int));
	#endif
	//Returning the new pointer
	return x;

}




/*
FUNCTION DESCRIPTION: Reallocate particle short int array memory on host, and on device if used.
LAST EDIT: Created
EDIT DATE: 6/15/2020
NAME OF EDITOR: Tasman Powis
*/
short int* ReallocShortInt(short int *x, int mold, int mnew)
{
	#ifdef _OPENACC
	short int *y;

	//Creating the new memory chunk on the GPU for the double array
	y = (short int*) acc_malloc(mnew*sizeof(short int));

	//Copying values from the old memory chunk into the new memory chunk
	acc_memcpy_device(y,acc_deviceptr(x),mold*sizeof(short int));

	//Deleting the old memory chunk on the GPU
	acc_delete(x,mold*sizeof(short int));
	#endif
	//Reallocating the new memory chunk on the CPU
	x = (short int *) realloc(x,mnew*sizeof(short int));
	#ifdef _OPENACC
	//Mapping the new GPU memory chunk to the pointer on the CPU
	acc_map_data(x,y,mnew*sizeof(short int));
	#endif
	//Returning the new pointer
	return x;

}





/*
FUNCTION DESCRIPTION: Reallocate particle unsigned short int array memory on host, and on device if used.
LAST EDIT: Created
EDIT DATE: 6/15/2020
NAME OF EDITOR: Tasman Powis
*/
unsigned short int* ReallocUnsignedShortInt(unsigned short int *x, int mold, int mnew)
{
	#ifdef _OPENACC
	unsigned short int *y;

	//Creating the new memory chunk on the GPU for the double array
	y = (unsigned short int*) acc_malloc(mnew*sizeof(unsigned short int));

	//Copying values from the old memory chunk into the new memory chunk
	acc_memcpy_device(y,acc_deviceptr(x),mold*sizeof(unsigned short int));

	//Deleting the old memory chunk on the GPU
	acc_delete(x,mold*sizeof(unsigned short int));
	#endif
	//Reallocating the new memory chunk on the CPU
	x = (unsigned short int *) realloc(x,mnew*sizeof(unsigned short int));
	#ifdef _OPENACC
	//Mapping the new GPU memory chunk to the pointer on the CPU
	acc_map_data(x,y,mnew*sizeof(unsigned short int));
	#endif
	//Returning the new pointer
	return x;

}




/*
FUNCTION DESCRIPTION: Pushing all particles from all species.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void PushAllParticles(struct AllSpecies *allspecies, struct Timing *timer)
{
	int i, myid;
	int nspecies = allspecies->num_species;
	int subc;
	double t0 = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Cycling through the species
	for (i = 0; i < nspecies; i++)
	{
		subc = allspecies->species[i].subc;

		//Only pushing at the species sub-cycle interval
		if (timer->step % subc == 0)
		{
			double dt;
			struct Species *species;

			dt = timer->dt*subc;

			species = &allspecies->species[i];

			PushParticles(species, dt);

		}
	}

	timer->timePush += MPI_Wtime() - t0;

}







/*
FUNCTION DESCRIPTION: Streamlined MPI Barrier function which incorporates timing.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void Barrier(struct Timing *timer)
{
	double t0 = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);

	timer->timeBarrier += MPI_Wtime() - t0;
}






/*
FUNCTION DESCRIPTION: Initializing the species structure for a 2D simulation. This routine also looks into the plasma and source creation routines to make sure that there is sufficient initial memory allocated.
LAST EDIT: Updated for new beam creation routine
EDIT DATE: 12/05/2022
NAME OF EDITOR: Tasman Powis
*/
void InitializeSpecies(struct Species *species, struct AllCreation *creation, struct Timing *timer)
{
    int i;
    int myid, myid_region, nproc, nproc_region;
	int nplasma = creation->num_plasma;
	int psnum, isnum, mem_source = 0, mem_beam = 0;
	int Npart_local_mem = 0;
	unsigned long iprn;

	double Npart_plasma_float;
	long long Npart_plasma_int = 0;


	//Setting up the PRNG
	species->process_data = (desprng_common_t *) malloc(sizeof(desprng_common_t));
	species->thread_data = (desprng_individual_t *) malloc(sizeof(desprng_individual_t));

	species->rncount = 0;

	initialize_common(species->process_data);
	create_identifier(&species->seed);
	initialize_individual(species->process_data, species->thread_data, species->seed);


	//Determining the total number and Region specific number of initial plasma particles to create. The method is based on area weighting and therefore requires grabbing and setting the Region limits on plasma creation for this task.
	for (i = 0; i < nplasma; i++)
	{

		struct CreatePlasma *plasma;

		plasma = &creation->plasma[i];

		if (plasma->snum == species->snum)
		{

			Npart_plasma_float = (plasma->xmax_glob - plasma->xmin_glob)*(plasma->ymax_glob - plasma->ymin_glob)*plasma->n0_plasma/species->N;
			Npart_plasma_int += (long long) Npart_plasma_float;

			plasma->Npart_plasma_local = (int) Npart_plasma_float;

		}
	}



	//printf("\nNpart_plasma_int = %lld\n",Npart_plasma_int);

	//Setting the initial particle counters to zero
	species->Npart = 0;
	species->Npart_local = 0;

	//Initialize memory sufficient for all plasma creation and Npart_init times all sourceed creation for this species
	Npart_local_mem = 2 * Npart_plasma_int;
	species->Npart_local_mem = Npart_local_mem;

	//printf("\nNpart_local_mem = %d\n",Npart_local_mem);

    //Allocating memory 2D particle phase space coordinates.
    species->x = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->y = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->vx = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->vy = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->vz = (ptype *) calloc(Npart_local_mem, sizeof(ptype));

    //Allocating memory 2D particle electric field values.
    species->Ex = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->Ey = (ptype *) calloc(Npart_local_mem, sizeof(ptype));

	//Allocating memory 3D particle magnetic field values.
    species->Bx = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->By = (ptype *) calloc(Npart_local_mem, sizeof(ptype));
    species->Bz = (ptype *) calloc(Npart_local_mem, sizeof(ptype));

    //Allocating memory 2D particle cell number (for sorting algorithm)
    species->cell = (int *) calloc(Npart_local_mem, sizeof(int));
    species->cell_stat = (int *) calloc(Npart_local_mem, sizeof(int));

    //Allocating memory 2D particle BC tags
    species->xbc = (short int *) calloc(Npart_local_mem, sizeof(short int));
    species->ybc = (short int *) calloc(Npart_local_mem, sizeof(short int));

    //Allocating memory 2D particle tag and printing logic
    species->tag = (unsigned short int *) calloc(Npart_local_mem, sizeof(unsigned short int));
    species->print_tag = (unsigned short int *) calloc(Npart_local_mem, sizeof(unsigned short int));

	//Allocating memory to array which stores particle indices which will NOT be deleted when applying BCs
	species->pbc_keep = (int *) calloc(Npart_local_mem, sizeof(int));
	//species->pbc_test = (int *) calloc(Npart_local_mem, sizeof(int));
	species->pbc_test = (int *) calloc(4*Npart_buff, sizeof(int));
	species->pbc_dtf = (ptype *) calloc(4*Npart_buff, sizeof(ptype));
	species->pbc_rem = (int *) calloc(4*Npart_buff, sizeof(int));
	species->rem_left = (int *) calloc(4*Npart_buff, sizeof(int));
	species->keep_right = (int *) calloc(4*Npart_buff, sizeof(int));

	species->pindex = (int *) calloc(Npart_local_mem, sizeof(int));

	//Allocating memory for particle source buffer, necessary since randomly created particles must be produced on host before being transmitted to device
	species->buff_source = (struct Species_buffer *) malloc(sizeof(struct Species_buffer));
	AllocateParticleBuffer(species->buff_source, Npart_buff);

	//Allocating memory to particle exit locations (in memory) and sorted location
	species->buff = (struct Species_buffer *) malloc(sizeof(struct Species_buffer));
	AllocateParticleBuffer(species->buff, 4*Npart_buff);

}



/*
FUNCTION DESCRIPTION: Clear all memory for a 2D species struct.
LAST EDIT: Created
EDIT DATE: 29/10/22
NAME OF EDITOR: Tasman Powis
*/
void FinalizeSpecies(struct Species *species)
{
	free(species->x);
	free(species->y);
	free(species->vx);
	free(species->vy);
	free(species->vz);

	free(species->Ex);
	free(species->Ey);

	free(species->Bx);
	free(species->By);
	free(species->Bz);

	free(species->cell);
	free(species->cell_stat);

	free(species->xbc);
	free(species->ybc);

	free(species->tag);
	free(species->print_tag);

	free(species->pbc_keep);
	free(species->pbc_test);
	free(species->pbc_dtf);
	free(species->pbc_rem);
	free(species->rem_left);
	free(species->keep_right);
	free(species->pindex);

	free(species->cindex);
	free(species->cdefs);

	free(species->process_data);
	free(species->thread_data);

	FinalizeParticleBuffer(species->buff_source);
	FinalizeParticleBuffer(species->buff);
}




/*
FUNCTION DESCRIPTION: Sub-routine for allocating memory for particle send and receive buffers for inter-Region communication.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void AllocateParticleBuffer(struct Species_buffer *buff, int mem)
{

	buff->npart_send = 0;
	buff->npart_recv = 0;
	buff->npart_keep = 0;
	buff->npart_recv_new = 0;

	buff->dtf_send = (ptype *) calloc(mem, sizeof(ptype));
    buff->x_send = (ptype *) calloc(mem, sizeof(ptype));
    buff->y_send = (ptype *) calloc(mem, sizeof(ptype));
    buff->vx_send = (ptype *) calloc(mem, sizeof(ptype));
    buff->vy_send = (ptype *) calloc(mem, sizeof(ptype));
    buff->vz_send = (ptype *) calloc(mem, sizeof(ptype));
    buff->xbc_send = (short int *) calloc(mem, sizeof(short int));
    buff->ybc_send = (short int *) calloc(mem, sizeof(short int));
    buff->cell_send = (int *) calloc(mem, sizeof(int));
    buff->tag_send = (unsigned short int *) calloc(mem, sizeof(unsigned short int));
    buff->print_tag_send = (unsigned short int *) calloc(mem, sizeof(unsigned short int));

	buff->dtf_recv = (ptype *) calloc(mem, sizeof(ptype));
    buff->x_recv = (ptype *) calloc(mem, sizeof(ptype));
    buff->y_recv = (ptype *) calloc(mem, sizeof(ptype));
    buff->vx_recv = (ptype *) calloc(mem, sizeof(ptype));
    buff->vy_recv = (ptype *) calloc(mem, sizeof(ptype));
    buff->vz_recv = (ptype *) calloc(mem, sizeof(ptype));
    buff->xbc_recv = (short int *) calloc(mem, sizeof(short int));
    buff->ybc_recv = (short int *) calloc(mem, sizeof(short int));
    buff->cell_recv = (int *) calloc(mem, sizeof(int));
    buff->tag_recv = (unsigned short int *) calloc(mem, sizeof(unsigned short int));
    buff->print_tag_recv = (unsigned short int *) calloc(mem, sizeof(unsigned short int));

}



/*
FUNCTION DESCRIPTION: Free memory for particle buffer.
LAST EDIT: Created
EDIT DATE: 29/10/22
NAME OF EDITOR: Tasman Powis
*/
void FinalizeParticleBuffer(struct Species_buffer *buff)
{
	free(buff->dtf_send);
	free(buff->x_send);
	free(buff->y_send);
	free(buff->vx_send);
	free(buff->vy_send);
	free(buff->vz_send);
	free(buff->xbc_send);
	free(buff->ybc_send);
	free(buff->cell_send);
	free(buff->tag_send);
	free(buff->print_tag_send);

	free(buff->dtf_recv);
	free(buff->x_recv);
	free(buff->y_recv);
	free(buff->vx_recv);
	free(buff->vy_recv);
	free(buff->vz_recv);
	free(buff->xbc_recv);
	free(buff->ybc_recv);
	free(buff->cell_recv);
	free(buff->tag_recv);
	free(buff->print_tag_recv);

	free(buff);
}



/*
FUNCTION DESCRIPTION: Initialize particles for a plasma creation routine for a 2D simulation.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void InitializePlasma(struct Species *species, struct CreatePlasma *plasma, struct Timing *timer)
{
    int i, j;
    int istart, iend;
    ptype xcenter, ycenter;
    ptype xwidth, ywidth;
    double print_rand;
    double V0;
    //unsigned long iprn;
    //double xprn;

	//Total number of particles to try and create
	istart = species->Npart_local;
	iend = istart + plasma->Npart_plasma_local;

	//Setting up the PRNG
	plasma->rncount = (unsigned long *) malloc(sizeof(unsigned long));
	plasma->process_data = (desprng_common_t *) malloc(sizeof(desprng_common_t));
	plasma->thread_data = (desprng_individual_t *) malloc(sizeof(desprng_individual_t));

	plasma->rncount[0] = 0;

	initialize_common(plasma->process_data);
	create_identifier(&plasma->seed);
	initialize_individual(plasma->process_data, plasma->thread_data, plasma->seed);

	//Setting up the creation domain for the particles
	plasma->xwidth = plasma->xmax - plasma->xmin;
	plasma->ywidth = plasma->ymax - plasma->ymin;
	xwidth = (ptype) plasma->xwidth;
	ywidth = (ptype) plasma->ywidth;

    plasma->xcenter = 0.5*(plasma->xmin + plasma->xmax);
    plasma->ycenter = 0.5*(plasma->ymin + plasma->ymax);
	xcenter = (ptype) plasma->xcenter;
	ycenter = (ptype) plasma->ycenter;

    //Computing the initial thermal velocity
    V0 = sqrt(e_div_me*plasma->T0/species->M);

	//Grabbing the inital drift velocities
	plasma->Vx = (ptype) sgn((double) plasma->Kx)*sqrt(e_div_me*fabs(plasma->Kx)/species->M);
	plasma->Vy = (ptype) sgn((double) plasma->Ky)*sqrt(e_div_me*fabs(plasma->Ky)/species->M);
	plasma->Vz = (ptype) sgn((double) plasma->Kz)*sqrt(e_div_me*fabs(plasma->Kz)/species->M);

	//Setting the initial count for placing particles
	j = istart;

	//#pragma omp parallel
	{
	unsigned long icount, iprn;
	int ii, jj, n;
	ptype x, y;

	//#pragma omp for //simd - throws up a warning that loop is not vectorized
    for (i = istart; i < iend; i++)
    {
    	//Need to provide a uniqie int for each PRNG generated. We allow 10 posibilities for each normally distributed number.
    	icount = plasma->rncount[0] + 33*i;

        //Initializing particle position in phase space
        x = (ptype) (get_uniform_prn(plasma->process_data, plasma->thread_data, icount, &iprn) - 0.5)*xwidth + xcenter;
        y = (ptype) (get_uniform_prn(plasma->process_data, plasma->thread_data, icount+1, &iprn) - 0.5)*ywidth + ycenter;

		//If it's in a vacuum cell, continue setting up the particle
		species->x[j] = x;
		species->y[j] = y;

        species->vx[j] = (ptype) RanGaussianDesprng(plasma->process_data, plasma->thread_data, icount + 2 , V0) + plasma->Vx;
        species->vy[j] = (ptype) RanGaussianDesprng(plasma->process_data, plasma->thread_data, icount + 12, V0) + plasma->Vy;
        species->vz[j] = (ptype) RanGaussianDesprng(plasma->process_data, plasma->thread_data, icount + 22, V0) + plasma->Vz;

		//Setting the particle boundary status (0 means within the region)
		species->xbc[j] = 0;
		species->ybc[j] = 0;

		//Setting the particle tag (0 means initial plasma species)
		species->tag[j] = 0;

		//Setting print tag for each particle (1 means print this particle in phase space output)
        print_rand = get_uniform_prn(plasma->process_data, plasma->thread_data, icount + 32, &iprn);
        if (print_rand < species->print_frac) {
                species->print_tag[j] = 1;
        }
        else {
                species->print_tag[j] = 0;
        }

		//Iterate counter
		j++;

	}

	}

	species->Npart_local = j;

}



/*
FUNCTION DESCRIPTION: Increasing the memory for particles storage of type "species" to new size "new_mem_size".
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void IncreaseSpeciesMemory(struct Species *species, int new_mem_size)
{

	int mold, mnew;
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	mold = species->Npart_local_mem;
	mnew = new_mem_size;

	//Resetting pointers
	species->x = ReallocPtype(species->x,mold,mnew);
	species->y = ReallocPtype(species->y,mold,mnew);
	species->vx = ReallocPtype(species->vx,mold,mnew);
	species->vy = ReallocPtype(species->vy,mold,mnew);
	species->vz = ReallocPtype(species->vz,mold,mnew);

	species->Ex = ReallocPtype(species->Ex,mold,mnew);
	species->Ey = ReallocPtype(species->Ey,mold,mnew);

	species->Bx = ReallocPtype(species->Bx,mold,mnew);
	species->By = ReallocPtype(species->By,mold,mnew);
	species->Bz = ReallocPtype(species->Bz,mold,mnew);

	species->cell = ReallocInt(species->cell,mold,mnew);
	species->cell_stat = ReallocInt(species->cell_stat,mold,mnew);

	species->xbc = ReallocShortInt(species->xbc,mold,mnew);
	species->ybc = ReallocShortInt(species->ybc,mold,mnew);

	species->tag = ReallocUnsignedShortInt(species->tag,mold,mnew);
	species->print_tag = ReallocUnsignedShortInt(species->print_tag,mold,mnew);

	species->pbc_keep = ReallocInt(species->pbc_keep,mold,mnew);
	//species->pbc_test = ReallocInt(species->pbc_test,mold,mnew);
	species->pindex = ReallocInt(species->pindex,mold,mnew);

}





/*
FUNCTION DESCRIPTION: Copy the particle structure from the host to an accelerator.
LAST EDIT: Updated
EDIT DATE: 1/8/2022
NAME OF EDITOR: Tasman Powis
*/
void MemoryToDeviceAllParticles(struct AllSpecies *allspecies)
{
	int i;
	int nspecies = allspecies->num_species;
	int npbm = 4*Npart_buff;
	int npim = Npart_buff;

	#pragma acc enter data copyin(allspecies[0:1])
	#pragma acc enter data copyin(allspecies->species[0:nspecies])

	//Cycling through the species
	for (i = 0; i < nspecies; i++)
	{

		int npm = allspecies->species[i].Npart_local_mem;
		int cnum = allspecies->species[i].cnum;
		int cinum = allspecies->species[i].coll_index_len;

		//Copying in the phase space data
		#pragma acc enter data copyin(allspecies->species[i].x[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].y[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].vx[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].vy[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].vz[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].Ex[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].Ey[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].Bx[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].By[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].Bz[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].cell[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].cell_stat[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].xbc[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].ybc[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].tag[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].print_tag[0:npm])

		#pragma acc enter data copyin(allspecies->species[i].pbc_keep[0:npm])
		#pragma acc enter data copyin(allspecies->species[i].pindex[0:npm])

		#pragma acc enter data copyin(allspecies->species[i].pbc_rem[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].pbc_test[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].pbc_dtf[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].rem_left[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].keep_right[0:npbm])

		//Copying in the source buffer
		#pragma acc enter data copyin(allspecies->species[i].buff_source[0:1])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->x_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->y_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->vx_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->vy_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->vz_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->xbc_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->ybc_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->cell_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->tag_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->print_tag_send[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->x_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->y_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->vx_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->vy_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->vz_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->xbc_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->ybc_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->cell_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->tag_recv[0:npim])
		#pragma acc enter data copyin(allspecies->species[i].buff_source->print_tag_recv[0:npim])

		//Copying in the communication buffer
		#pragma acc enter data copyin(allspecies->species[i].buff[0:1])
		#pragma acc enter data copyin(allspecies->species[i].buff->dtf_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->x_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->y_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->vx_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->vy_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->vz_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->xbc_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->ybc_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->cell_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->tag_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->print_tag_send[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->x_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->y_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->vx_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->vy_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->vz_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->xbc_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->ybc_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->cell_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->tag_recv[0:npbm])
		#pragma acc enter data copyin(allspecies->species[i].buff->print_tag_recv[0:npbm])

		//Copying in PRNG data
		#pragma acc enter data copyin(allspecies->species[i].rncount[0:1])
		#pragma acc enter data copyin(allspecies->species[i].process_data[:1])
		#pragma acc enter data copyin(allspecies->species[i].thread_data[:1])

		//Copying in species collision data
		#pragma acc enter data copyin(allspecies->species[i].cindex[0:cnum])
		#pragma acc enter data copyin(allspecies->species[i].cdefs[0:cnum])
		#pragma acc enter data copyin(allspecies->species[i].coll_index[0:cinum])
	}
}







/*
FUNCTION DESCRIPTION: Copy the particle structure from the an accelerator back to the host.
LAST EDIT: Updated
EDIT DATE: 1/8/2022
NAME OF EDITOR: Tasman Powis
*/
void MemoryDeleteDeviceAllParticles(struct AllSpecies *allspecies)
{
	int i;
	int nspecies = allspecies->num_species;
	int npbm = 4*Npart_buff;
	int npim = Npart_buff;

	//Cycling through the species
	for (i = 0; i < nspecies; i++)
	{

		#pragma acc update host(allspecies->species[i].Npart_local_mem)
		#pragma acc update host(allspecies->species[i].cnum)
		#pragma acc update host(allspecies->species[i].coll_index_len)

		int npm = allspecies->species[i].Npart_local_mem;
		int cnum = allspecies->species[i].cnum;
		int cinum = allspecies->species[i].coll_index_len;

		//Deleting the phase space data
		#pragma acc exit data delete(allspecies->species[i].x[0:npm])
		#pragma acc exit data delete(allspecies->species[i].y[0:npm])
		#pragma acc exit data delete(allspecies->species[i].vx[0:npm])
		#pragma acc exit data delete(allspecies->species[i].vy[0:npm])
		#pragma acc exit data delete(allspecies->species[i].vz[0:npm])
		#pragma acc exit data delete(allspecies->species[i].Ex[0:npm])
		#pragma acc exit data delete(allspecies->species[i].Ey[0:npm])
		#pragma acc exit data delete(allspecies->species[i].Bx[0:npm])
		#pragma acc exit data delete(allspecies->species[i].By[0:npm])
		#pragma acc exit data delete(allspecies->species[i].Bz[0:npm])
		#pragma acc exit data delete(allspecies->species[i].cell[0:npm])
		#pragma acc exit data delete(allspecies->species[i].cell_stat[0:npm])
		#pragma acc exit data delete(allspecies->species[i].xbc[0:npm])
		#pragma acc exit data delete(allspecies->species[i].ybc[0:npm])
		#pragma acc exit data delete(allspecies->species[i].tag[0:npm])
		#pragma acc exit data delete(allspecies->species[i].print_tag[0:npm])

		#pragma acc exit data delete(allspecies->species[i].pbc_keep[0:npm])
		#pragma acc exit data delete(allspecies->species[i].pindex[0:npm])

		#pragma acc exit data delete(allspecies->species[i].pbc_rem[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].pbc_test[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].pbc_dtf[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].rem_left[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].keep_right[0:npbm])

		//Deleting the source buffer
		#pragma acc exit data delete(allspecies->species[i].buff_source->x_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->y_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->vx_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->vy_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->vz_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->xbc_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->ybc_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->cell_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->tag_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->print_tag_send[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->x_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->y_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->vx_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->vy_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->vz_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->xbc_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->ybc_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->cell_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->tag_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source->print_tag_recv[0:npim])
		#pragma acc exit data delete(allspecies->species[i].buff_source[0:1])

		//Deleting the communication buffer
		#pragma acc exit data delete(allspecies->species[i].buff->dtf_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->x_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->y_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->vx_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->vy_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->vz_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->xbc_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->ybc_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->cell_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->tag_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->print_tag_send[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->x_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->y_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->vx_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->vy_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->vz_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->xbc_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->ybc_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->cell_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->tag_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff->print_tag_recv[0:npbm])
		#pragma acc exit data delete(allspecies->species[i].buff[0:1])

		//Deleting PRNG data
		#pragma acc exit data delete(allspecies->species[i].rncount[0:1])
		#pragma acc exit data delete(allspecies->species[i].process_data[:1])
		#pragma acc exit data delete(allspecies->species[i].thread_data[:1])

		//Deleting species collision data
		#pragma acc exit data delete(allspecies->species[i].cindex[0:cnum])
		#pragma acc exit data delete(allspecies->species[i].cdefs[0:cnum])
		#pragma acc exit data delete(allspecies->species[i].coll_index[0:cinum])
	}

	#pragma acc exit data delete(allspecies->species[0:nspecies])
	#pragma acc exit data delete(allspecies[0:1])

}





/*
FUNCTION DESCRIPTION: Copy the particle structure from the an accelerator back to the host.
LAST EDIT: Updated
EDIT DATE: 1/8/2022
NAME OF EDITOR: Tasman Powis
*/
void UpdateHostAllParticles(struct AllSpecies *allspecies, struct Timing *timer)
{
	int i;
	int nspecies = allspecies->num_species;
	double t0 = MPI_Wtime();

	//Cycling through the species
	for (i = 0; i < nspecies; i++)
	{
		//struct Species *species = &allspecies->species[i];
		int npm;

		#pragma acc update host(allspecies->species[i].Npart_local_mem)
		npm = allspecies->species[i].Npart_local_mem;

		#pragma acc update host(allspecies->species[i].x[0:npm])
		#pragma acc update host(allspecies->species[i].y[0:npm])
		#pragma acc update host(allspecies->species[i].vx[0:npm])
		#pragma acc update host(allspecies->species[i].vy[0:npm])
		#pragma acc update host(allspecies->species[i].vz[0:npm])
		#pragma acc update host(allspecies->species[i].Ex[0:npm])
		#pragma acc update host(allspecies->species[i].Ey[0:npm])
		#pragma acc update host(allspecies->species[i].Bx[0:npm])
		#pragma acc update host(allspecies->species[i].By[0:npm])
		#pragma acc update host(allspecies->species[i].Bz[0:npm])
		#pragma acc update host(allspecies->species[i].xbc[0:npm])
		#pragma acc update host(allspecies->species[i].ybc[0:npm])
		#pragma acc update host(allspecies->species[i].cell[0:npm])
		#pragma acc update host(allspecies->species[i].tag[0:npm])
		#pragma acc update host(allspecies->species[i].print_tag[0:npm])
	}

	timer->timeMemCopy += MPI_Wtime() - t0;

}






/*
FUNCTION DESCRIPTION: Pushing all particles in "species". Uses the Boris algorithm, analytic magnetic field and interpolated electric fields.
LAST EDIT: Replaced OpenMP threading with OpenACC threading for multicore and GPUs
EDIT DATE: 5/29/2020
NAME OF EDITOR: Tasman Powis
*/
void PushParticles(struct Species *species, double dt)
{
	//Particle properties
	const ptype q_div_m = (ptype) e_div_me*species->Q/species->M;
	const ptype q_dth_div_m = (ptype) 0.5*dt*q_div_m;

	//Loop properties
	int i;
	const int Npart_local = species->Npart_local;
	const int npl = species->Npart_local;

	//Restricted pointers for species array quantities
	ptype *restrict X, *restrict Y, *restrict VX, *restrict VY, *restrict VZ, *restrict EX, *restrict EY, *restrict BXp, *restrict BYp, *restrict BZp;

	//Assigning restricted pointers
	X = species->x;
	Y = species->y;
	VX = species->vx;
	VY = species->vy;
	VZ = species->vz;
	EX = species->Ex;
	EY = species->Ey;
	BXp = species->Bx;
	BYp = species->By;
	BZp = species->Bz;


	#ifdef _OPENACC
	#pragma acc parallel present(species[0:1],X[0:npl],Y[0:npl],VX[0:npl],VY[0:npl],VZ[0:npl],EX[0:npl],EY[0:npl],BXp[0:npl],BYp[0:npl],BZp[0:npl])
	#else
	#pragma omp parallel
	#endif
	{

	ptype x, y, z;
	ptype vxp, vyp, vzp, vxm, vym, vzm, vxprime, vyprime, vzprime;
	ptype Ex, Ey;
	ptype Bx, By, Bz;
	ptype tx, ty, tz;
	ptype sx, sy, sz;
	ptype divtabs2;

	#ifdef _OPENACC
	#pragma acc loop
	#else
	#pragma omp for simd
	#endif
	for (i = 0; i < Npart_local; i++)
	{
		//Particle position
		x = X[i];
		y = Y[i];

		//Local magnetic field from analytic expression - found to be faster to do this in-situ, i.e. not pre-calculated
		// Bx = BX(x,y);
		// By = BY(x,y);
		// Bz = BZ(x,y);
		Bx = BXp[i];
		By = BYp[i];
		Bz = BZp[i];

		//Computing relevant parameters - NEED TO OPTIMIZE THIS
		tx = q_dth_div_m*Bx;
		ty = q_dth_div_m*By;
		tz = q_dth_div_m*Bz;

		divtabs2 = 2.0/(1 + tx*tx + ty*ty + tz*tz);

		sx = tx*divtabs2;
		sy = ty*divtabs2;
		sz = tz*divtabs2;

		//Getting the local electric field.
		Ex = EX[i];
		Ey = EY[i];

		//Velocity space push
		vxm = VX[i] + q_dth_div_m*Ex;
		vym = VY[i] + q_dth_div_m*Ey;
		vzm = VZ[i];

		vxprime = vxm + (vym*tz - vzm*ty);
		vyprime = vym + (vzm*tx - vxm*tz);
		vzprime = vzm + (vxm*ty - vym*tx);

		vxp = vxm + (vyprime*sz - vzprime*sy);
		vyp = vym + (vzprime*sx - vxprime*sz);
		vzp = vzm + (vxprime*sy - vyprime*sx);

		vxm = vxp + q_dth_div_m*Ex;
		vym = vyp + q_dth_div_m*Ey;
		vzm = vzp;

		VX[i] = vxm;
		VY[i] = vym;
		VZ[i] = vzm;

		//Position space push
		x += dt*vxm;
		y += dt*vym;

		//Saving positions
		X[i] = x;
		Y[i] = y;

    }
    }


}



/*
FUNCTION DESCRIPTION: Sub-routine for storing all particles from the newly created collision particles the main species structure.
LAST EDIT: Created
EDIT DATE: 2/8/2021
NAME OF EDITOR: Tasman Powis
*/
void StoreNewCollisionParticlesBuffer(struct Species *species, struct Species_buffer *buff, int npbm)
{

	//Number of particles to receive
	int npart = buff->npart_recv;

	//Restricted pointers for species array quantities
	ptype *restrict X;
	ptype *restrict Y;
	ptype *restrict VX;
	ptype *restrict VY;
	ptype *restrict VZ;
	short int *restrict XBC;
	short int *restrict YBC;
	unsigned short int *restrict T;
	unsigned short int *restrict TP;

	//Initializing and assigning pointers for buffer
	ptype *restrict Xb = species->buff_source->x_recv;
	ptype *restrict Yb = species->buff_source->y_recv;
	ptype *restrict VXb = species->buff_source->vx_recv;
	ptype *restrict VYb = species->buff_source->vy_recv;
	ptype *restrict VZb = species->buff_source->vz_recv;
	unsigned short int *restrict Tb = species->buff_source->tag_recv;
	unsigned short int *restrict TPb = species->buff_source->print_tag_recv;

	//Size of memory chunks on GPU
	int npm = species->Npart_local_mem;


	//Storing particles
	if (npart > 0) {
		int i, j = 0;
		int istart, iend;
		int Npart_local_mem;
		int new_mem;


		istart = species->Npart_local;
		iend = istart + npart;

		Npart_local_mem = species->Npart_local_mem;

		while (iend >= Npart_local_mem) {

			if(Npart_local_mem < 10) {Npart_local_mem += 10;}

			Npart_local_mem *= 2;

			IncreaseSpeciesMemory(species,Npart_local_mem);

		}

		species->Npart_local_mem = Npart_local_mem;
		npm = Npart_local_mem;

		//#pragma acc update device(species[0:1])
		#pragma acc update device(species->Npart_local_mem)

		//Assigning pointers post possible memory increase
		X = species->x;
		Y = species->y;
		VX = species->vx;
		VY = species->vy;
		VZ = species->vz;
		XBC = species->xbc;
		YBC = species->ybc;
		T = species->tag;
		TP = species->print_tag;


		#ifdef _OPENACC
		#pragma acc parallel loop default(present)
		#else
		#pragma omp parallel for
		#endif
		for (i = istart; i < iend; i++) {

			j = i - istart;

			X[i] = Xb[j];
			Y[i] = Yb[j];
			VX[i] = VXb[j];
			VY[i] = VYb[j];
			VZ[i] = VZb[j];
			XBC[i] = 0;
			YBC[i] = 0;
			T[i] = Tb[j];
			TP[i] = TPb[j];


		}

		//Storing the new particle count
		species->Npart_local = iend;

	}

	//#pragma acc update device(species[0:1])
	#pragma acc update device(species->Npart_local)

}