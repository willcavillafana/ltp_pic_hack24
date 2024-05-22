#include "main.h"



/*
FUNCTION DESCRIPTION: Allocate memory for each of the collision modules before they are read in from the input file.
LAST EDIT: 2/2/2021
EDIT DATE: Creation
NAME OF EDITOR: Tasman Powis
*/
void AllocateAllCollisions(struct AllCollisions *allcollisions)
{

	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;


	allcollisions->elastic = (struct ElasticCollision *) malloc(nelastic * sizeof(struct ElasticCollision));
	allcollisions->excite = (struct ExcitationCollision *) malloc(nexcite * sizeof(struct ExcitationCollision));
	allcollisions->ionize = (struct IonizationCollision *) malloc(nionize * sizeof(struct IonizationCollision));
	allcollisions->cxchange = (struct ChargeExchangeCollision *) malloc(ncxchange * sizeof(struct ChargeExchangeCollision));
}



/*
FUNCTION DESCRIPTION: Initialize each of the collision modules using data read from the input file.
LAST EDIT: 2/2/2021
EDIT DATE: Creation
NAME OF EDITOR: Tasman Powis
*/
//void InitializeAllCollisions(struct Control *control, struct AllSpecies *allspecies, struct AllCollisions *creation, struct Timing *timer);
void InitializeAllCollisions(struct AllSpecies *allspecies, struct AllCollisions *allcollisions)
{
	int i;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;
	int sindex;


	for (i = 0; i < nelastic; i++)
	{
		allcollisions->elastic[i].seed = global_seed + 0; //NEED TO UPDATE BASED ON PREVIOUSLY ASSIGNED SEEDS

		sindex = allcollisions->elastic[i].species_in_index;

		InitializeElasticCollision(&allspecies->species[sindex], &allcollisions->elastic[i]);
	}

	for (i = 0; i < nexcite; i++)
	{
		allcollisions->excite[i].seed = global_seed + 0; //NEED TO UPDATE BASED ON PREVIOUSLY ASSIGNED SEEDS

		sindex = allcollisions->excite[i].species_in_index;

		InitializeExcitationCollision(&allspecies->species[sindex], &allcollisions->excite[i]);
	}

	for (i = 0; i < nionize; i++)
	{
		allcollisions->ionize[i].seed = global_seed + 0; //NEED TO UPDATE BASED ON PREVIOUSLY ASSIGNED SEEDS

		sindex = allcollisions->ionize[i].electron_in_index;

		InitializeIonizationCollision(allspecies, &allcollisions->ionize[i]);
	}

	for (i = 0; i < ncxchange; i++)
	{
		allcollisions->cxchange[i].seed = global_seed + 0; //NEED TO UPDATE BASED ON PREVIOUSLY ASSIGNED SEEDS

		sindex = allcollisions->cxchange[i].species_in_index;

		InitializeChargeExchangeCollision(&allspecies->species[sindex], &allcollisions->cxchange[i]);
	}


}









/*
FUNCTION DESCRIPTION: Initialize an elastic collision module
LAST EDIT: 2/2/2021
EDIT DATE: Creation
NAME OF EDITOR: Tasman Powis
*/
void InitializeElasticCollision(struct Species *species, struct ElasticCollision *elastic)
{

	int i;

	//Filling in the collision cross section data
	elastic->nlines = GetXSectionFileNumLines(elastic->xsection_file);

	elastic->energy_vals = (double *) calloc(elastic->nlines, sizeof(double));
	elastic->xsection_vals = (double *) calloc(elastic->nlines, sizeof(double));

	GetXSectionFileData(elastic->xsection_file, elastic->energy_vals, elastic->xsection_vals);

	// printf("\n\nenergy[%d] = %lf\n",0,elastic->energy_vals[0]);
	// printf("energy[%d] = %lf\n",elastic->nlines-1,elastic->energy_vals[elastic->nlines-1]);
	// printf("xsection[%d] = %.6e\n\n",elastic->nlines-1,elastic->xsection_vals[elastic->nlines-1]);


	//Setting up the PRNG
	elastic->process_data = (desprng_common_t *) malloc(sizeof(desprng_common_t));
	elastic->thread_data = (desprng_individual_t *) malloc(sizeof(desprng_individual_t));

	elastic->rncount = 0;

	initialize_common(elastic->process_data);
	create_identifier(&elastic->seed);
	initialize_individual(elastic->process_data, elastic->thread_data, elastic->seed);


}




/*
FUNCTION DESCRIPTION: Initialize an excitation collision module
LAST EDIT: 2/2/2021
EDIT DATE: Creation
NAME OF EDITOR: Tasman Powis
*/
void InitializeExcitationCollision(struct Species *species, struct ExcitationCollision *excite)
{

	//Informing the species which collisions they are associated with - TO DO LATER DURING OPTIMIZATION


	//Filling in the collision cross section data
	excite->nlines = GetXSectionFileNumLines(excite->xsection_file);

	excite->energy_vals = (double *) calloc(excite->nlines, sizeof(double));
	excite->xsection_vals = (double *) calloc(excite->nlines, sizeof(double));

	GetXSectionFileData(excite->xsection_file, excite->energy_vals, excite->xsection_vals);

	//printf("\n\nenergy[1] = %lf\n\n",elastic->energy_vals[1]);


	//Setting up the PRNG
	excite->process_data = (desprng_common_t *) malloc(sizeof(desprng_common_t));
	excite->thread_data = (desprng_individual_t *) malloc(sizeof(desprng_individual_t));

	excite->rncount = 0;

	initialize_common(excite->process_data);
	create_identifier(&excite->seed);
	initialize_individual(excite->process_data, excite->thread_data, excite->seed);


}






/*
FUNCTION DESCRIPTION: Initialize a charge exchange collision module
LAST EDIT: 2/13/2021
EDIT DATE: Created
NAME OF EDITOR: Tasman Powis
*/
void InitializeChargeExchangeCollision(struct Species *species, struct ChargeExchangeCollision *cxchange)
{

	//Informing the species which collisions they are associated with - TO DO LATER DURING OPTIMIZATION


	//Filling in the collision cross section data
	cxchange->nlines = GetXSectionFileNumLines(cxchange->xsection_file);

	cxchange->energy_vals = (double *) calloc(cxchange->nlines, sizeof(double));
	cxchange->xsection_vals = (double *) calloc(cxchange->nlines, sizeof(double));

	GetXSectionFileData(cxchange->xsection_file, cxchange->energy_vals, cxchange->xsection_vals);

	//printf("\n\nenergy[1] = %lf\n\n",elastic->energy_vals[1]);


	//Setting up the PRNG
	cxchange->process_data = (desprng_common_t *) malloc(sizeof(desprng_common_t));
	cxchange->thread_data = (desprng_individual_t *) malloc(sizeof(desprng_individual_t));

	cxchange->rncount = 0;

	initialize_common(cxchange->process_data);
	create_identifier(&cxchange->seed);
	initialize_individual(cxchange->process_data, cxchange->thread_data, cxchange->seed);


}





/*
FUNCTION DESCRIPTION: Initialize an ionization collision module
LAST EDIT: 2/2/2021
EDIT DATE: Creation
NAME OF EDITOR: Tasman Powis
*/
void InitializeIonizationCollision(struct AllSpecies *allspecies, struct IonizationCollision *ionize)
{

	int sindex, eindex, iindex;


	//Filling in the collision cross section data
	ionize->nlines = GetXSectionFileNumLines(ionize->xsection_file);

	ionize->energy_vals = (double *) calloc(ionize->nlines, sizeof(double));
	ionize->xsection_vals = (double *) calloc(ionize->nlines, sizeof(double));

	GetXSectionFileData(ionize->xsection_file, ionize->energy_vals, ionize->xsection_vals);

	//printf("\n\nbion = %lf\n\n",ionize->bion);

	//Pre-computing the weight ratio for the species to be created
	sindex = ionize->electron_in_index;
	eindex = ionize->electron_out_index;
	iindex = ionize->ion_out_index;

	//ionize->enew_int = (int) allspecies->species[sindex].N/allspecies->species[eindex].N;
	//ionize->enew_frac = allspecies->species[sindex].N/allspecies->species[eindex].N - (double) ionize->enew_int;

	ionize->inew_int = (int) allspecies->species[sindex].N/allspecies->species[iindex].N;
	ionize->inew_frac = allspecies->species[sindex].N/allspecies->species[iindex].N - (double) ionize->inew_int;

	//printf("\nenew_int = %d, enew_frac = %f\n",ionize->enew_int,ionize->enew_frac);
	//printf("\ninew_int = %d, inew_frac = %f\n",ionize->inew_int,ionize->inew_frac);

	//Setting up the PRNG
	ionize->process_data = (desprng_common_t *) malloc(sizeof(desprng_common_t));
	ionize->thread_data = (desprng_individual_t *) malloc(sizeof(desprng_individual_t));

	ionize->rncount = 0;

	initialize_common(ionize->process_data);
	create_identifier(&ionize->seed);
	initialize_individual(ionize->process_data, ionize->thread_data, ionize->seed);


}










/*
FUNCTION DESCRIPTION: Free memory from each of the collision modules.
LAST EDIT: Created
EDIT DATE: 29/10/22
NAME OF EDITOR: Tasman Powis
*/
void FinalizeAllCollisions(struct AllCollisions *allcollisions)
{
	int i;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;
	int sindex;


	//Freeing data from each collision struct
	for (i = 0; i < nelastic; i++)
	{
		free(allcollisions->elastic[i].energy_vals);
		free(allcollisions->elastic[i].xsection_vals);
		free(allcollisions->elastic[i].process_data);
		free(allcollisions->elastic[i].thread_data);
	}

	for (i = 0; i < nexcite; i++)
	{
		free(allcollisions->excite[i].energy_vals);
		free(allcollisions->excite[i].xsection_vals);
		free(allcollisions->excite[i].process_data);
		free(allcollisions->excite[i].thread_data);
	}

	for (i = 0; i < nionize; i++)
	{
		free(allcollisions->ionize[i].energy_vals);
		free(allcollisions->ionize[i].xsection_vals);
		free(allcollisions->ionize[i].process_data);
		free(allcollisions->ionize[i].thread_data);
	}

	for (i = 0; i < ncxchange; i++)
	{
		free(allcollisions->cxchange[i].energy_vals);
		free(allcollisions->cxchange[i].xsection_vals);
		free(allcollisions->cxchange[i].process_data);
		free(allcollisions->cxchange[i].thread_data);
	}


	//Freeing the struct arrays
	free(allcollisions->elastic);
	free(allcollisions->excite);
	free(allcollisions->ionize);
	free(allcollisions->cxchange);


}










/*
FUNCTION DESCRIPTION: Map the collision indices in to the species structure.
LAST EDIT: 4/2/2021
EDIT DATE: Created
NAME OF EDITOR: Tasman Powis
*/
void MapCollisionsToSpecies(struct AllSpecies *allspecies, struct AllCollisions *allcollisions)
{
	int i, s;
	int nspecies = allspecies->num_species;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;
	int sindex, sindex_eout, sindex_iout;

	struct Species *species;


	for (s = 0; s < nspecies; s++)
	{

		//int nel = 0, nex = 0, nio = 0;
		int nc;

		//First we count the number of each type of collision
		species = &allspecies->species[s];

		//species->cnum_elastic = 0;
		//species->cnum_excite = 0;
		//species->cnum_ionize = 0;

		species->cnum = 0;


		for (i = 0; i < nelastic; i++)
		{
			sindex = allcollisions->elastic[i].species_in_index;

			if (s == sindex) {
				//species->cnum_elastic++;
				species->cnum++;
			}

		}

		for (i = 0; i < nexcite; i++)
		{

			sindex = allcollisions->excite[i].species_in_index;

			if (s == sindex) {
				//species->cnum_excite++;
				species->cnum++;
			}

		}

		for (i = 0; i < nionize; i++)
		{

			sindex = allcollisions->ionize[i].electron_in_index;
			//sindex_eout = allcollisions->ionize[i].electron_out_index;
			//sindex_iout = allcollisions->ionize[i].ion_out_index;

			//if (s == sindex || s == sindex_eout || s == sindex_iout) {
			if (s == sindex) {
				//species->cnum_ionize++;
				species->cnum++;
			}

		}

		for (i = 0; i < ncxchange; i++)
		{

			sindex = allcollisions->cxchange[i].species_in_index;

			if (s == sindex) {
				//species->cnum_excite++;
				species->cnum++;
			}

		}


		//Then allocate space to store the collision module indices
		//species->cindex_elastic = (int *) malloc(species->cnum_elastic*sizeof(int));
		//species->cindex_excite = (int *) malloc(species->cnum_excite*sizeof(int));
		//species->cindex_ionize = (int *) malloc(species->cnum_ionize*sizeof(int));

		species->cindex = (int *) malloc(species->cnum*sizeof(int));
		species->cdefs = (ctype *) malloc(species->cnum*sizeof(ctype));


		//Then we store the collision module indices
		nc = 0;
		for (i = 0; i < nelastic; i++)
		{
			sindex = allcollisions->elastic[i].species_in_index;

			if (s == sindex) {
				//species->cindex_elastic[nel] = i;
				//nel++;
				species->cindex[nc] = i;
				species->cdefs[nc] = ELASTIC;
				nc++;
			}

		}

		for (i = 0; i < nexcite; i++)
		{

			sindex = allcollisions->excite[i].species_in_index;

			if (s == sindex) {
				//species->cindex_excite[nex] = i;
				//nex++;
				species->cindex[nc] = i;
				species->cdefs[nc] = EXCITE;
				nc++;
			}

		}

		for (i = 0; i < nionize; i++)
		{

			sindex = allcollisions->ionize[i].electron_in_index;
			//sindex_eout = allcollisions->ionize[i].electron_out_index;
			//sindex_iout = allcollisions->ionize[i].ion_out_index;

			//if (s == sindex || s == sindex_eout || s == sindex_iout) {
			if (s == sindex) {
				//species->cindex_ionize[nio] = i;
				//nio++;
				species->cindex[nc] = i;
				species->cdefs[nc] = IONIZE;
				nc++;
			}

		}

		for (i = 0; i < ncxchange; i++)
		{

			sindex = allcollisions->cxchange[i].species_in_index;

			if (s == sindex) {
				//species->cindex_excite[nex] = i;
				//nex++;
				species->cindex[nc] = i;
				species->cdefs[nc] = CXCHANGE;
				nc++;
			}

		}


		// printf("\nSpecies %d, collision indices and types:\n",s);
		// for (i = 0; i < species->cnum; i++) {
		// 	printf("type = %d, num = %d\n",species->cdefs[i],species->cindex[i]);
		// }


/*
		printf("\nSpecies %d, elastic collision indices: ",s);
		for (i = 0; i < species->cnum_elastic; i++) {
			printf("%d ",species->cindex_elastic[i]);
		}


		printf("\nSpecies %d, excitation collision indices: ",s);
		for (i = 0; i < species->cnum_excite; i++) {
			printf("%d ",species->cindex_excite[i]);
		}

		printf("\nSpecies %d, ionization collision indices: ",s);
		for (i = 0; i < species->cnum_ionize; i++) {
			printf("%d ",species->cindex_ionize[i]);
		}
		printf("\n");
*/

	}
	//printf("\n\n");



}










/*
FUNCTION DESCRIPTION: Map the collision indices in to the species structure.
LAST EDIT: 4/2/2021
EDIT DATE: Created
NAME OF EDITOR: Tasman Powis
*/
void ComputeMaxProbabilities(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, double dt)
{
	int s, myid;
	int nspecies = allspecies->num_species;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;
	int sindex, sindex_eout, sindex_iout;

	struct Species *species;

	//Getting the collision type structs
	struct ElasticCollision *elastic = allcollisions->elastic;
	struct ExcitationCollision *excite = allcollisions->excite;
	struct IonizationCollision *ionize = allcollisions->ionize;
	struct ChargeExchangeCollision *cxchange = allcollisions->cxchange;

	int nsamples = 100000;
	double ns0 = -4.0;
	double ns1 = 4.0;
	double dke;


	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	dke = (ns1 - ns0)/nsamples;

	//Running through the species loop
	for (s = 0; s < nspecies; s++)
	{

		//int nel = 0, nex = 0, nio = 0;
		int i, I, k;
		double numax_tot = 0.0, pmax_tot = 0.0, ng_max = 0.0;
		ctype colltype;
		double *KE = (double *) calloc(nsamples+1, sizeof(double));
		double *VS = (double *) calloc(nsamples+1, sizeof(double));
		double *XST = (double *) calloc(nsamples+1, sizeof(double));

		//First we count the number of each type of collision
		species = &allspecies->species[s];


		//Setting up vectors to track data
		for (k = 0; k < nsamples+1; k++)
		{
			KE[k] = pow(10.0,(ns0 + k*dke));
			VS[k] = pow((2.0*e_div_me/species->M*KE[k]),0.5);

			//if (myid == 0)
			//printf("K[%d] = %f, VS[%d] = %f\n",k,KE[k],k,VS[k]);

		}

		//Running through each collision associated with the particle
		for (i = 0; i < species->cnum; i++)
		{
			I = species->cindex[i];
        	colltype = species->cdefs[i];

			//Elastic collisions
			if (colltype == ELASTIC)
			{
				double ke, xs;

				// if (myid == 0) {
				// 	printf("ELASTIC\n");
				// }

				if (elastic[I].neutral_n0 > ng_max) {
					ng_max = elastic[I].neutral_n0;
				}

				for (k = 0; k < nsamples+1; k++)
				{
					int j;
					int nl = elastic[I].nlines;
					double e0, e1, xs0, xs1;

					ke = KE[k];

					if (ke < elastic[I].energy_vals[0]) {
						xs = elastic[I].xsection_vals[0];
					}

					else if (ke >= elastic[I].energy_vals[nl-1]){
						xs = elastic[I].xsection_vals[nl-1];
					}

					else{
						for (j = 0; j < nl-1; j++)
						{
							if (ke >= elastic[I].energy_vals[j] && ke < elastic[I].energy_vals[j+1]) {
								e0 = elastic[I].energy_vals[j];
								e1 = elastic[I].energy_vals[j+1];
								xs0 = elastic[I].xsection_vals[j];
								xs1 = elastic[I].xsection_vals[j+1];

								xs = xs0 + (ke - e0)*(xs1 - xs0)/(e1 - e0);
								break;
							}
						}
					}

					XST[k] += xs;

				}

			}


			//Elastic collisions
			if (colltype == EXCITE)
			{
				double ke, xs;

				// if (myid == 0) {
				// 	printf("EXCITE\n");
				// }

				if (excite[I].neutral_n0 > ng_max) {
					ng_max = excite[I].neutral_n0;
				}

				for (k = 0; k < nsamples+1; k++)
				{
					int j;
					int nl = excite[I].nlines;
					double e0, e1, xs0, xs1;

					ke = KE[k];

					if (ke < excite[I].energy_vals[0]) {
						xs = excite[I].xsection_vals[0];
					}

					else if (ke >= excite[I].energy_vals[nl-1]){
						xs = excite[I].xsection_vals[nl-1];
					}

					else{
						for (j = 0; j < nl-1; j++)
						{
							if (ke >= excite[I].energy_vals[j] && ke < excite[I].energy_vals[j+1]) {
								e0 = excite[I].energy_vals[j];
								e1 = excite[I].energy_vals[j+1];
								xs0 = excite[I].xsection_vals[j];
								xs1 = excite[I].xsection_vals[j+1];

								xs = xs0 + (ke - e0)*(xs1 - xs0)/(e1 - e0);
								break;
							}
						}
					}

					XST[k] += xs;

				}

			}


			//Elastic collisions
			if (colltype == IONIZE)
			{
				double ke, xs;

				// if (myid == 0) {
				// 	printf("IONIZE\n");
				// }

				if (ionize[I].neutral_n0 > ng_max) {
					ng_max = ionize[I].neutral_n0;
				}

				for (k = 0; k < nsamples+1; k++)
				{
					int j;
					int nl = ionize[I].nlines;
					double e0, e1, xs0, xs1;

					ke = KE[k];

					if (ke < ionize[I].energy_vals[0]) {
						xs = ionize[I].xsection_vals[0];
					}

					else if (ke >= ionize[I].energy_vals[nl-1]){
						xs = ionize[I].xsection_vals[nl-1];
					}

					else{
						for (j = 0; j < nl-1; j++)
						{
							if (ke >= ionize[I].energy_vals[j] && ke < ionize[I].energy_vals[j+1]) {
								e0 = ionize[I].energy_vals[j];
								e1 = ionize[I].energy_vals[j+1];
								xs0 = ionize[I].xsection_vals[j];
								xs1 = ionize[I].xsection_vals[j+1];

								xs = xs0 + (ke - e0)*(xs1 - xs0)/(e1 - e0);
								break;
							}
						}
					}

					XST[k] += xs;

				}

			}


			//Elastic collisions
			if (colltype == CXCHANGE)
			{
				double ke, xs;

				// if (myid == 0) {
				// 	printf("CXCHANGE\n");
				// }

				if (cxchange[I].neutral_n0 > ng_max) {
					ng_max = cxchange[I].neutral_n0;
				}

				for (k = 0; k < nsamples+1; k++)
				{
					int j;
					int nl = cxchange[I].nlines;
					double e0, e1, xs0, xs1;

					ke = KE[k];

					if (ke < cxchange[I].energy_vals[0]) {
						xs = cxchange[I].xsection_vals[0];
					}

					else if (ke >= cxchange[I].energy_vals[nl-1]){
						xs = cxchange[I].xsection_vals[nl-1];
					}

					else{
						for (j = 0; j < nl-1; j++)
						{
							if (ke >= cxchange[I].energy_vals[j] && ke < cxchange[I].energy_vals[j+1]) {
								e0 = cxchange[I].energy_vals[j];
								e1 = cxchange[I].energy_vals[j+1];
								xs0 = cxchange[I].xsection_vals[j];
								xs1 = cxchange[I].xsection_vals[j+1];

								xs = xs0 + (ke - e0)*(xs1 - xs0)/(e1 - e0);
								break;
							}
						}
					}

					XST[k] += xs;

				}

			}

		}


		//Computing the maximum collision frequency
		for (k = 0; k < nsamples+1; k++)
		{
			if (ng_max*VS[k]*XST[k] > numax_tot) {
				numax_tot = ng_max*VS[k]*XST[k];
			}

			//if (myid == 0)
			//printf("KE[%d] = %f, XS[k] = %.8e, NU[k] = %.8e\n",k,KE[k],k,XST[k],k,ng_max*VS[k]*XST[k]);

		}


		species->numax_tot = numax_tot;
		species->pmax_tot = 1.0 - exp(-dt*species->numax_tot);

		//Allocating the initial collision index array
		species->coll_index_len = 1000;
		species->coll_index = (int *) calloc(species->coll_index_len,sizeof(int));



		// if (myid == 0) {
		//   	printf("numax_tot[%d] = %lf, pmax_tot[%d] = %lf\n",s,species->numax_tot,s,species->pmax_tot);
		// }

		free(KE);
		free(VS);
		free(XST);


	}
	//printf("\n\n");

}





/*
FUNCTION DESCRIPTION: Execute all collision algorithms at the given time step
LAST EDIT: 8/8/2022
EDIT DATE: Updated to reproduce Turner benchmark
NAME OF EDITOR: Tasman Powis
*/
void CollideAllParticles(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, struct Timing *timer)
{

	int i, myid;
	int nspecies = allspecies->num_species;

	double t0 = MPI_Wtime();


	//Cycling through the species
	for (i = 0; i < nspecies; i++)
	{

		//Collisions of species
		CollideParticlesNew(allspecies,allcollisions,i,timer->dt);

	}


	//Storing created particles
	for (i = 0; i < nspecies; i++)
	{


		#pragma acc update host(allspecies->species[i].buff_source->npart_recv)
		//printf("New species %d: %d\n",i,allspecies->species[i].buff_source->npart_recv);
		StoreNewCollisionParticlesBuffer(&allspecies->species[i], allspecies->species[i].buff_source, Npart_buff);

		allspecies->species[i].buff_source->npart_recv = 0;
		#pragma acc update device(allspecies->species[i].buff_source->npart_recv)


	}

	timer->timeCollide += MPI_Wtime() - t0;

}













/*
FUNCTION DESCRIPTION: Execute all collision algorithms for a specific species
LAST EDIT: 30/1/2023
EDIT DATE: Updated to precompute the colliding particles
NAME OF EDITOR: Tasman Powis
*/
void CollideParticlesNew(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, int sindex, double dt)
{

	//Getting the species struct
	int snum = allspecies->num_species;
	struct Species *species = &allspecies->species[sindex];
	int cnum = species->cnum;
	int ncoll, ncoll_est, new_len;
	//double ncoll_frac;
	//int istart, iend;

	//Getting the collision type structs
	struct ElasticCollision *elastic = allcollisions->elastic;
	struct ExcitationCollision *excite = allcollisions->excite;
	struct IonizationCollision *ionize = allcollisions->ionize;
	struct ChargeExchangeCollision *cxchange = allcollisions->cxchange;

	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;

	//Collision counts
	int ccoll = 0;
	int celastic = 0;
	int cexcite = 0;
	int cionize = 0;
	int ccxchange = 0;
	int cnull = 0;

	//Species properties
	double M = species->M;
	double N = species->N;
	double numax = species->numax_tot;
	double pmax = species->pmax_tot;

	//Loop properties
	//unsigned long i;
	int i, jj = 0;
        int Npart_local = species->Npart_local;

	//Restricted pointers for species array quantities
	ptype *restrict X, *restrict Y, *restrict Z, *restrict VX, *restrict VY, *restrict VZ;

	//Assigning restricted pointers
	X = species->x;
	Y = species->y;
	VX = species->vx;
	VY = species->vy;
	VZ = species->vz;


	//Calculate the estimated number of collision to occur
	ncoll_est = (int) (((double) Npart_local)*pmax);

	//Increasing the size of the collision index array if needed - with added buffer of 25%
	new_len = (int) (1.5*ncoll_est);
	if (new_len > species->coll_index_len) {
		species->coll_index = ReallocInt(species->coll_index, species->coll_index_len, new_len);
		species->coll_index_len = new_len;
	}
	#pragma acc update device(species->coll_index_len)

	//printf("\nspecies->coll_index_len = %d",species->coll_index_len);

	//Creating the collision index list
	#pragma acc parallel default(present) copy(jj) present(species[0:1],species->coll_index[:species->coll_index_len])
	{
	unsigned long icount, iprn;
	int kk;

	#pragma acc loop
	for (i = 0; i < Npart_local; i++)
    	{
		icount = species->rncount + i;
		if (get_uniform_prn(species->process_data, species->thread_data, icount+0, &iprn) >= pmax) {
			continue;
		}
		else {
			#pragma acc atomic capture
			{kk = jj; jj++;}
			species->coll_index[kk] = i;
			//printf("\njj = %lu",jj);
		}
	}
	}

	//Actual number of collisions
	ncoll = jj;

	//printf("\nncoll = %d",jj);

	/*
	#pragma acc parallel present(species->coll_index[:species->coll_index_len])
	for (i = 0; i < species->coll_index_len; i++)
	{
		printf("\ncoll_index[%d] = %d",i,species->coll_index[i]);
	}
	*/

	//return;

	//Updating PRNG counter
	species->rncount += Npart_local;
	#pragma acc update device(species->rncount)


	//Collision algorithm
	#ifdef _OPENACC
	#pragma acc parallel default(present) copy(ccoll,celastic,cexcite,cionize,ccxchange,cnull) present(elastic[:nelastic],excite[:nexcite],ionize[:nionize],cxchange[:ncxchange],species->coll_index[:species->coll_index_len])
	//#pragma acc parallel default(present) copy(ccoll,celastic,cexcite,cionize,ccxchange,cnull)
	//#pragma acc parallel default(present) copy(ccoll,celastic,cexcite,cionize,ccxchange,cnull) present(elastic[:nelastic],excite[:nexcite],ionize[:nionize],cxchange[:ncxchange]) vector_length(1024)
	#else
	#pragma omp parallel
	#endif
    {

	int j, I, J, k, nl, m;
	double x, vx, vy, vz, v, v2, kev;
	double xs, e0, e1, xs0, xs1;
	double nu, ng, pcoll0, pcoll1;
	ctype colltype;
	unsigned long icount, jcount, iprn;
	double srand, rand;

	double MG, TG, V0;
	double gx, gy, gz, g, vscale;
	double gx1, gy1, gz1, g1;
	double gx2, gy2, gz2, g2;
	double theta, phi, chi, eta;
	double chi1, chi2, eta1, eta2;
	double coschi, sinchi, cosphi, sinphi;
	double st, ct, sp, cp;
	double sc, cc, se, ce;
	double sc1, cc1, se1, ce1;
	double sc2, cc2, se2, ce2;
	double r11, r12, r13, r21, r22, r23, r31, r32, r33;
	double up1, up2, up3, mag;

	double vxn, vyn, vzn;
	double vxion, vyion, vzion;
	double vxrel, vyrel, vzrel, v2rel, vrel, kerel;
	double wex, wion, bion, ke, ke_ev, ke1, ke2, ke1_ev, ke2_ev, mu;
	double vscale1, vscale2;

	int le, li, seo, sio, inew;
	struct Species_buffer *buff_eout;
	struct Species_buffer *buff_iout;


	#ifdef _OPENACC
	#pragma acc loop
	#else
	#pragma omp for
	#endif
	for (i = 0; i < ncoll; i++)
    {

		//Getting the particle index
		I = species->coll_index[i];

		//Get particle data
		vx = VX[I];
		vy = VY[I];
		vz = VZ[I];

		v2 = vx*vx + vy*vy + vz*vz;
		v = pow(v2,0.5);
		kev = 0.5*me_div_e*M*v2;

		//Set initial collision probability to zero and draw a new number
		pcoll0 = 0.0;
		pcoll1 = 0.0;

		//Drawing a random number for testing which collision occurs
		icount = species->rncount + i;
		srand = get_uniform_prn(species->process_data, species->thread_data, icount+0, &iprn);

		//Cycling through the collisions of this species to determine which collision occurs
		#pragma acc loop seq
		for (j = 0; j < species->cnum; j++)
		{
			J = species->cindex[j];
        	colltype = species->cdefs[j];

			//Elastic collisions - electrons
			if ((colltype == ELASTIC) && (M == 1.0))
			{

				nl = elastic[J].nlines;
				ng = elastic[J].neutral_n0;

				if (kev < elastic[J].energy_vals[0]) {
					xs = elastic[J].xsection_vals[0];
				}

				else if (kev >= elastic[J].energy_vals[nl-1]){
					xs = elastic[J].xsection_vals[nl-1];
				}

				else{
					#pragma acc loop seq
					for (k = 0; k < nl-1; k++)
					{
						if (kev >= elastic[J].energy_vals[k] && kev < elastic[J].energy_vals[k+1]) {
							e0 = elastic[J].energy_vals[k];
							e1 = elastic[J].energy_vals[k+1];
							xs0 = elastic[J].xsection_vals[k];
							xs1 = elastic[J].xsection_vals[k+1];

							xs = xs0 + (kev - e0)*(xs1 - xs0)/(e1 - e0);
							//break; //Can't have break in parallelisation
						}
					}
				}

				//Updating the propability
				pcoll1 += ng*v*xs/numax;

				if ((srand >= pcoll0) && (srand < pcoll1)) {

					MG = elastic[J].neutral_M0;

					//Adding to collision count
					#ifdef _OPENACC
					#pragma acc atomic
					#else
					#pragma omp atomic
					#endif
					celastic++;



					//COLLISION ALGORITHM
					jcount = elastic[J].rncount + i*2;

					gx = vx;
					gy = vy;
					gz = vz;

					g = sqrt(gx*gx + gy*gy + gz*gz);

					theta = atan2(sqrt(gy*gy + gz*gz),gx);
					phi = atan2(gz,gy);

					st = sin(theta);
					ct = cos(theta);
					sp = sin(phi);
					cp = cos(phi);

					//Scattering angles
					rand = get_uniform_prn(elastic[J].process_data, elastic[J].thread_data, jcount+0, &iprn);
					#ifdef BENCH_TURNER
					chi = acos(1.0 - 2.0*rand);
					#else
					if (kev != 0.0) {
						chi = acos((2.0 + kev - 2.0*pow((1.0 + kev),rand))/kev);
					}
					else {
						chi = acos(1.0 - 2.0*rand);
					}
					#endif

					rand = get_uniform_prn(elastic[J].process_data, elastic[J].thread_data, jcount+1, &iprn);
					eta = 2.0*PI*rand;

					sc  = sin(chi);
					cc  = cos(chi);
					se  = sin(eta);
					ce  = cos(eta);

					gx = g*(ct*cc - st*sc*ce);
					gy = g*(st*cp*cc + ct*cp*sc*ce - sp*sc*se);
					gz = g*(st*sp*cc + ct*sp*sc*ce + cp*sc*se);

					vscale = sqrt(1.0-2.0*M/MG*(1.0 - cc));

					VX[I] = vscale*gx;
					VY[I] = vscale*gy;
					VZ[I] = vscale*gz;

				}

				//Copying over the probability for the next test
				pcoll0 = pcoll1;

			}




			//Elastic collisions - ions
			else if ((colltype == ELASTIC) && (M > 1.0))
			{
				MG = M; //elastic[J].neutral_M0; - Only allowing collisions with neutrals of the same weight
				TG = elastic[J].neutral_T0;

				//Sampling a neutral particle
				jcount = elastic[J].rncount + i*40;
				V0 = sqrt(e_div_me*TG/MG);

				vxn = (ptype) RanGaussianDesprng(elastic[J].process_data, elastic[J].thread_data, jcount + 0 , V0);
				vyn = (ptype) RanGaussianDesprng(elastic[J].process_data, elastic[J].thread_data, jcount + 10 , V0);
				vzn = (ptype) RanGaussianDesprng(elastic[J].process_data, elastic[J].thread_data, jcount + 20 , V0);

				//Relative velocity of ion to sampled neutral
				vxrel = vx - vxn;
				vyrel = vy - vyn;
				vzrel = vz - vzn;
				v2rel = vxrel*vxrel + vyrel*vyrel + vzrel*vzrel;
				vrel = sqrt(v2rel);
				mu = 0.5*M;
				kerel = 0.5*me_div_e*mu*v2rel;

				nl = elastic[J].nlines;
				ng = elastic[J].neutral_n0;

				if (kerel < elastic[J].energy_vals[0]) {
					xs = elastic[J].xsection_vals[0];
				}

				else if (kerel >= elastic[J].energy_vals[nl-1]){
					xs = elastic[J].xsection_vals[nl-1];
				}

				else{
					#pragma acc loop seq
					for (k = 0; k < nl-1; k++)
					{
						if (kerel >= elastic[J].energy_vals[k] && kerel < elastic[J].energy_vals[k+1]) {
							e0 = elastic[J].energy_vals[k];
							e1 = elastic[J].energy_vals[k+1];
							xs0 = elastic[J].xsection_vals[k];
							xs1 = elastic[J].xsection_vals[k+1];

							xs = xs0 + (kerel - e0)*(xs1 - xs0)/(e1 - e0);
							//break; //Can't have break in parallelisation
						}
					}
				}

				//Updating the propability
				pcoll1 += ng*vrel*xs/numax;

				if ((srand >= pcoll0) && (srand < pcoll1)) {

					//Adding to collision count
					#ifdef _OPENACC
					#pragma acc atomic
					#else
					#pragma omp atomic
					#endif
					celastic++;



					//COLLISION ALGORITHM
					gx = vx - vxn;
					gy = vy - vyn;
					gz = vz - vzn;

					g = sqrt(gx*gx + gy*gy + gz*gz);

					//Scattering angles
					rand = get_uniform_prn(elastic[J].process_data, elastic[J].thread_data, jcount+30, &iprn);
					coschi = sqrt(rand);
					if (coschi*coschi > 1.0) {
						coschi = 1.0;
					}
					sinchi = sqrt(1.0 - coschi*coschi);

					rand = get_uniform_prn(elastic[J].process_data, elastic[J].thread_data, jcount+31, &iprn);
					phi = 2.0*PI*rand;
					cosphi = cos(phi);
					sinphi = sin(phi);

					if (g > 0.0) {
						r13 = gx/g;
						r23 = gy/g;
						r33 = gz/g;
					}
					else {
						r13 = 0.0;
						r23 = 0.0;
						r33 = 0.0;
					}

					up1 = 0.0;
					up2 = 0.0;
					up3 = 0.0;

					if (r33 == 1.0) {
						up2 = 1.0;
					}
					else {
						up3 = 1.0;
					}

					r12 = r23*up3 - r33*up2;
					r22 = r33*up1 - r13*up3;
					r32 = r13*up2 - r23*up1;
					mag = sqrt(r12*r12 + r22*r22 + r32*r32);

					r12 = r12/mag;
					r22 = r22/mag;
					r32 = r32/mag;

					r11 = r22*r33 - r32*r23;
					r21 = r32*r13 - r12*r33;
					r31 = r12*r23 - r22*r13;

					g = g*coschi;

					gx = g*(r11*sinchi*cosphi + r12*sinchi*sinphi + r13*coschi);
					gy = g*(r21*sinchi*cosphi + r22*sinchi*sinphi + r23*coschi);
					gz = g*(r31*sinchi*cosphi + r32*sinchi*sinphi + r33*coschi);

					VX[I] = gx + vxn;
					VY[I] = gy + vyn;
					VZ[I] = gz + vzn;

				}

				//Copying over the probability for the next test
				pcoll0 = pcoll1;

			}


			//Excitation collisions
			else if (colltype == EXCITE)
			{
				nl = excite[J].nlines;
				ng = excite[J].neutral_n0;

				if (kev < excite[J].energy_vals[0]) {
					xs = excite[J].xsection_vals[0];
				}

				else if (kev >= excite[J].energy_vals[nl-1]){
					xs = excite[J].xsection_vals[nl-1];
				}

				else{
					#pragma acc loop seq
					for (k = 0; k < nl-1; k++)
					{
						if (kev >= excite[J].energy_vals[k] && kev < excite[J].energy_vals[k+1]) {
							e0 = excite[J].energy_vals[k];
							e1 = excite[J].energy_vals[k+1];
							xs0 = excite[J].xsection_vals[k];
							xs1 = excite[J].xsection_vals[k+1];

							xs = xs0 + (kev - e0)*(xs1 - xs0)/(e1 - e0);
							//break; //Can't have break in parallelisation
						}
					}
				}

				pcoll1 += ng*v*xs/numax;

				if ((srand >= pcoll0) && (srand < pcoll1)) {

					MG = excite[J].neutral_M0;
					wex = excite[J].neutral_wex;

					//Adding to collision count
					#ifdef _OPENACC
					#pragma acc atomic
					#else
					#pragma omp atomic
					#endif
					cexcite++;

					//COLLISION ALGORITHM
					jcount = excite[J].rncount + i*2;

					gx = vx;
					gy = vy;
					gz = vz;

					g = sqrt(gx*gx + gy*gy + gz*gz);

					theta = atan2(sqrt(gy*gy + gz*gz),gx);
					phi = atan2(gz,gy);

					st = sin(theta);
					ct = cos(theta);
					sp = sin(phi);
					cp = cos(phi);

					//Adjusting the kinetic energy for excitation
					ke = 0.5*M*mele*g*g;
					ke = fabs(ke - wex*qele);
                    ke_ev = ke/qele;
					g = sqrt(2.0*ke/mele/M);

					//Scattering angles
					rand = get_uniform_prn(excite[J].process_data, excite[J].thread_data, jcount+0, &iprn);

					#ifdef BENCH_TURNER
					chi = acos(1.0 - 2.0*rand);
					#else
					if (ke_ev > 1.0e-8) {
						chi = acos((2.0 + ke_ev - 2.0*pow((1.0 + ke_ev),rand))/ke_ev);
                        //chi = 0.0;
					}
					else {
						chi = acos(1.0 - 2.0*rand);
					}
					#endif

					rand = get_uniform_prn(excite[J].process_data, excite[J].thread_data, jcount+1, &iprn);
					eta = 2.0*PI*rand;

					sc  = sin(chi);
					cc  = cos(chi);
					se  = sin(eta);
					ce  = cos(eta);

					gx = g*(ct*cc - st*sc*ce);
					gy = g*(st*cp*cc + ct*cp*sc*ce - sp*sc*se);
					gz = g*(st*sp*cc + ct*sp*sc*ce + cp*sc*se);

					#ifdef BENCH_TURNER
					vscale = 1.0;
					#else
					vscale = sqrt(1.0-2.0*M/MG*(1.0 - cc));
					#endif

					VX[I] = vscale*gx;
					VY[I] = vscale*gy;
					VZ[I] = vscale*gz;
				}

				//Copying over the probability for the next test
				pcoll0 = pcoll1;

			}



			//Ionization collisions
			else if (colltype == IONIZE)
			{
				nl = ionize[J].nlines;
				ng = ionize[J].neutral_n0;

				if (kev < ionize[J].energy_vals[0]) {
					xs = ionize[J].xsection_vals[0];
				}

				else if (kev >= ionize[J].energy_vals[nl-1]){
					xs = ionize[J].xsection_vals[nl-1];
				}

				else{
					#pragma acc loop seq
					for (k = 0; k < nl-1; k++)
					{
						if (kev >= ionize[J].energy_vals[k] && kev < ionize[J].energy_vals[k+1]) {
							e0 = ionize[J].energy_vals[k];
							e1 = ionize[J].energy_vals[k+1];
							xs0 = ionize[J].xsection_vals[k];
							xs1 = ionize[J].xsection_vals[k+1];

							xs = xs0 + (kev - e0)*(xs1 - xs0)/(e1 - e0);
							//break; //Can't have break in parallelisation
						}
					}
				}

				pcoll1 += ng*v*xs/numax;

				if ((srand >= pcoll0) && (srand < pcoll1)) {

					//double MG = ionize[J].neutral_M0;
					TG = ionize[J].neutral_T0;
					wion = ionize[J].neutral_wion;
					bion = ionize[J].bion;

					seo = ionize[J].electron_out_index;
					sio = ionize[J].ion_out_index;
					MG = allspecies->species[sio].M; //Gas atom mass is equal to that of the ion being created.
					buff_eout = allspecies->species[seo].buff_source;
					buff_iout = allspecies->species[sio].buff_source;

					//Adding to collision count
					#ifdef _OPENACC
					#pragma acc atomic
					#else
					#pragma omp atomic
					#endif
					cionize++;

					//COLLISION ALGORITHM
					//jcount = ionize[J].rncount + (i-istart)*40;
					jcount = ionize[J].rncount + i*(10+(ionize->inew_int+1)*40);

					gx = vx;
					gy = vy;
					gz = vz;

					g = sqrt(gx*gx + gy*gy + gz*gz);

					theta = atan2(sqrt(gy*gy + gz*gz),gx);
					phi = atan2(gz,gy);

					st = sin(theta);
					ct = cos(theta);
					sp = sin(phi);
					cp = cos(phi);

					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+0, &iprn);
					ke = 0.5*M*mele*g*g;
					ke1 = fabs(ke - wion*qele);
					#ifdef BENCH_TURNER
					ke1 = 0.5*ke1;
					ke2 = ke1;
					#else
					ke1 = bion*tan(rand*atan(0.5*ke1/bion));
					ke2 = ke - ke1 - wion*qele;
					#endif
                    ke1_ev = ke1/qele;
                    ke2_ev = ke2/qele;

					g1 = sqrt(2.0*ke1/mele/M);
					g2 = sqrt(2.0*ke2/mele/M);

					//Computing scattering angles
					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+1, &iprn);

					#ifdef BENCH_TURNER
					chi1 = acos(1.0 - 2.0*rand);
					#else
					if (ke1_ev > 1.0e-8) {
						chi1 = acos((2.0 + ke1_ev - 2.0*pow((1.0 + ke1_ev),rand))/ke1_ev);
					}
					else {
						chi1 = acos(1.0 - 2.0*rand);
					}
					#endif

					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+2, &iprn);
					#ifdef BENCH_TURNER
					chi2 = acos(1.0 - 2.0*rand);
					#else
					if (ke2_ev > 1.0e-8) {
						chi2 = acos((2.0 + ke2_ev - 2.0*pow((1.0 + ke2_ev),rand))/ke2_ev);
					}
					else {
						chi2 = acos(1.0 - 2.0*rand);
					}
					#endif

					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+3, &iprn);
					eta1 = 2.0*PI*rand;

					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+4, &iprn);
					eta2 = 2.0*PI*rand;

					sc1  = sin(chi1);
					cc1  = cos(chi1);
					se1  = sin(eta1);
					ce1  = cos(eta1);

					sc2  = sin(chi2);
					cc2  = cos(chi2);
					se2  = sin(eta2);
					ce2  = cos(eta2);

					gx1 = g1*(ct*cc1 - st*sc1*ce1);
					gy1 = g1*(st*cp*cc1 + ct*cp*sc1*ce1 - sp*sc1*se1);
					gz1 = g1*(st*sp*cc1 + ct*sp*sc1*ce1 + cp*sc1*se1);

					gx2 = g2*(ct*cc2 - st*sc2*ce2);
					gy2 = g2*(st*cp*cc2 + ct*cp*sc2*ce2 - sp*sc2*se2);
					gz2 = g2*(st*sp*cc2 + ct*sp*sc2*ce2 + cp*sc2*se2);

					#ifdef BENCH_TURNER
					vscale1 = 1.0;
					vscale2 = 1.0;
					#else
					vscale1 = sqrt(1.0-2.0*M/MG*(1.0 - cc1));
					vscale2= sqrt(1.0-2.0*M/MG*(1.0 - cc2));
					#endif

					//Updated velocity of colliding electron
					VX[I] = vscale1*gx1;
					VY[I] = vscale1*gy1;
					VZ[I] = vscale1*gz1;

					//printf("before le: %d\n",le);
					//printf("before npart_recv: %d\n",buff_eout->npart_recv);

					//Velocity of new electron
					#ifdef _OPENACC
					#pragma acc atomic capture
					#else
					#pragma omp atomic capture
					#endif
					{le = buff_eout->npart_recv; buff_eout->npart_recv++;}
					//printf("nionize = %d\n",le);
					//{le = allspecies->species[seo].buff_source->npart_recv; allspecies->species[seo].buff_source->npart_recv++;}

					//Storing the created electrons in to a buffer
					buff_eout->x_recv[le] = X[I];
					buff_eout->y_recv[le] = Y[I];
					buff_eout->vx_recv[le] = vscale2*gx2;
					buff_eout->vy_recv[le] = vscale2*gy2;
					buff_eout->vz_recv[le] = vscale2*gz2;

					buff_eout->tag_recv[le] = 2; //Ionization tag

					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+5, &iprn);

					if (rand < allspecies->species[seo].print_frac) {
					  buff_eout->print_tag_recv[le] = 1;
					}
					else {
					  buff_eout->print_tag_recv[le] = 0;
					}


					//Determining the number of new ions to create
					inew = ionize->inew_int;

					rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+6, &iprn);
					if (rand < ionize->inew_frac) {
						inew++;
					}

					//Ion thermal velocity
					V0 = sqrt(e_div_me*TG/MG);
					jcount += 7;

					//Creating new ions
					#pragma acc loop seq
					for (m = 0; m < inew; m++) {

						//Velocity components of new ion
						vxion = (ptype) RanGaussianDesprng(ionize[J].process_data, ionize[J].thread_data, jcount , V0);
						vyion = (ptype) RanGaussianDesprng(ionize[J].process_data, ionize[J].thread_data, jcount + 10, V0);
						vzion = (ptype) RanGaussianDesprng(ionize[J].process_data, ionize[J].thread_data, jcount + 20, V0);

						#ifdef _OPENACC
						#pragma acc atomic capture
						#else
						#pragma omp atomic capture
						#endif
						{li = buff_iout->npart_recv; buff_iout->npart_recv++;}

						//Storing the created ions in to a buffer
						buff_iout->x_recv[li] = X[I];
						buff_iout->y_recv[li] = Y[I];
						buff_iout->vx_recv[li] = vxion;
						buff_iout->vy_recv[li] = vyion;
						buff_iout->vz_recv[li] = vzion;

						buff_iout->tag_recv[li] = 2;

		      			rand = get_uniform_prn(ionize[J].process_data, ionize[J].thread_data, jcount+30, &iprn);

						if (rand < allspecies->species[sio].print_frac) {
						  buff_iout->print_tag_recv[li] = 1;
						}
						else {
						  buff_iout->print_tag_recv[li] = 0;
						}

						jcount += 31;

					}

				}

				//Copying over the probability for the next test
				pcoll0 = pcoll1;

			}


			//Charge-exchange collisions
			else if (colltype == CXCHANGE)
			{
				MG = M; //cxchange[J].neutral_M0; M; //Only allowing collisions with neutrals of the same weight
				TG = cxchange[J].neutral_T0;

				//Sampling a neutral particle
				jcount = cxchange[J].rncount + i*50;
				V0 = sqrt(e_div_me*TG/MG);

				vxn = (ptype) RanGaussianDesprng(cxchange[J].process_data, cxchange[J].thread_data, jcount + 0 , V0);
				vyn = (ptype) RanGaussianDesprng(cxchange[J].process_data, cxchange[J].thread_data, jcount + 10 , V0);
				vzn = (ptype) RanGaussianDesprng(cxchange[J].process_data, cxchange[J].thread_data, jcount + 20 , V0);

				//Relative velocity of ion to sampled neutral
				vxrel = vx - vxn;
				vyrel = vy - vyn;
				vzrel = vz - vzn;
				v2rel = vxrel*vxrel + vyrel*vyrel + vzrel*vzrel;
				vrel = sqrt(v2rel);
                #ifdef BENCH_TURNER
				mu = 0.5*M;
				kerel = 0.5*me_div_e*mu*v2rel;
                #else
                kerel = 0.5*me_div_e*M*v2rel;
                #endif

				nl = cxchange[J].nlines;
				ng = cxchange[J].neutral_n0;

				if (kerel < cxchange[J].energy_vals[0]) {
					xs = cxchange[J].xsection_vals[0];
				}

				else if (kerel >= cxchange[J].energy_vals[nl-1]){
					xs = cxchange[J].xsection_vals[nl-1];
				}

				else{
					#pragma acc loop seq
					for (k = 0; k < nl-1; k++)
					{
						if (kerel >= cxchange[J].energy_vals[k] && kerel < cxchange[J].energy_vals[k+1]) {
							e0 = cxchange[J].energy_vals[k];
							e1 = cxchange[J].energy_vals[k+1];
							xs0 = cxchange[J].xsection_vals[k];
							xs1 = cxchange[J].xsection_vals[k+1];

							xs = xs0 + (kerel - e0)*(xs1 - xs0)/(e1 - e0);
							//break; //Can't have break in parallelisation
						}
					}
				}

				pcoll1 += ng*vrel*xs/numax;

				if ((srand >= pcoll0) && (srand < pcoll1)) {

					//Adding to collision count
					#ifdef _OPENACC
					#pragma acc atomic
					#else
					#pragma omp atomic
					#endif
					ccxchange++;

					//COLLISION ALGORITHM
					VX[I] = vxn;
					VY[I] = vyn;
					VZ[I] = vzn;

				}

				//Copying over the probability for the next test
				pcoll0 = pcoll1;

			}

		} //End loop over collision types

		//If we did not have a collision, then we had a null collision
		if (srand > pcoll0) {
			//Adding to collision count
			#ifdef _OPENACC
			#pragma acc atomic
			#else
			#pragma omp atomic
			#endif
			cnull++;
		}

	} //End loop over particles


	//Update the coll-pnrg indexes
	#pragma acc loop seq
	for (i = 0; i < species->cnum; i++)
	{

	    J = species->cindex[i];
	    colltype = species->cdefs[i];


	    //Elastic collisions
	    if (colltype == ELASTIC)
	    {
			elastic[J].rncount += 40*ncoll;
	    }

	    //Excitation collisions
	    if (colltype == EXCITE)
	    {
			excite[J].rncount += 3*ncoll;
	    }

	    //Ionization collisions
	    if (colltype == IONIZE)
	    {
			ionize[J].rncount += 500*ncoll;
	    }

		//Exchange collisions
	    if (colltype == CXCHANGE)
	    {
		    cxchange[J].rncount += 30*ncoll;
	    }

	}

	} //End parallel region

	species->rncount += 2*ncoll;
	//#pragma acc update device(species[0:1])
	#pragma acc update device(species->rncount)

}