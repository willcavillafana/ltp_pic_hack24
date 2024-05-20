#include "main.h"



/*
FUNCTION DESCRIPTION: Reads in a string and finds the extension of the name
LAST EDIT: Modified to be c++ compliant
EDIT DATE: 20/11/2022
NAME OF EDITOR: Tasman Powis
*/
const char *get_filename_ext(char *filename)
{
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}




/*
FUNCTION DESCRIPTION: Reads in the command line arguments, finds the file name with *.dat and sets this as the input file. Default is "input.dat"
LAST EDIT: Modified to be c++ compliant
EDIT DATE: 20/11/2022
NAME OF EDITOR: Tasman Powis
*/
void GetInputFileName(int *argc, char ***argv, struct Control *control)
{
	int i, myid, ndat = 0;
	char *str;
	const char *val;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	for (i = 0; i < *argc; i++)
	{
		str = argv[0][i];

		val = get_filename_ext(str);

		if (strcmp(val,"dat") == 0) {
			strcpy(control->input,str);
			ndat++;
		}
	}

	if (ndat == 0) {
		if (myid == 0) {printf("No input file listed. Must be a single file of format *.dat\nEXITING\n\n");}
		MPI_Finalize();
		exit(0);
	}

	if (ndat > 1) {
		if (myid == 0) {printf("Too many input files listed. Must be a single file of format *.dat\nEXITING\n\n");}
		MPI_Finalize();
		exit(0);
	}

}




/*
FUNCTION DESCRIPTION: Checkes whether a string is empty or not
LAST EDIT: Added
EDIT DATE: 1/12/22
NAME OF EDITOR: Copied from: https://codeforwin.org/2018/02/c-program-remove-empty-lines-from-file.html
*/
int IsEmpty(const char *str)
{
    char ch;
    do
    {
        ch = *(str++);

        // Check non whitespace character
        if(ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' && ch != '\0')
            return 0;

    } while (ch != '\0');

    return 1;
}









/*
FUNCTION DESCRIPTION: Reading the input file for details of setup for the timing, control, grid, BCs and output. Counts number of species, creations, blocks, collisions, phase space diagnostics and probes
LAST EDIT: Updated with a more robust way to handle input
EDIT DATE: 2/12/22
NAME OF EDITOR: Tasman Powis
*/
void ProcessingInput(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *allcreation, struct AllCollisions *allcollisions, struct Timing *timer, struct Output *output)
{

	int i, j, k, cc, nlines_all, nlines, slength = 256;
	int *input_loc;
	char **input_all, **input;
	char **parameters, **args;
	int myid;
	int nspecies = 0, nplasma = 0, nsource = 0, nbeam = 0, nblock = 0, nelastic = 0, nexcite = 0, nionize = 0, ncxchange = 0, nphase = 0, nprobe = 0, nsee = 0;
	//int *snum_list, snum;
	FILE *finput;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Opening input.dat to read
	finput = fopen(control->input,"r");

	if (finput == NULL) {
		if (myid == 0) {printf("Input file does not exist\nEXITING\n");}
		MPI_Finalize();
		exit(0);
	}

	//Counting the number of lines
	nlines_all = CountLines(finput);

	//Pointing back to the start of the file
	rewind(finput);

	//Allocating space for the array of lines from input.dat
	input_all = (char **) malloc(nlines_all * sizeof(char *));
	for (i = 0; i < nlines_all; i++) {
		input_all[i] = (char *) calloc(slength,sizeof(char));
	}

	//Loading in the lines
	for (i = 0; i < nlines_all; i++) {
		fgets(input_all[i], slength, finput);
	}

	//Closing input.dat
	fclose(finput);


	//Removing comments (anything including and following a #)
	for (j = 0; j < nlines_all; j++) {

		//Checking if there is a comment in a given line and finding location of '#'
		cc = -1;
		for (i = 0; i < slength; i++) {
			if (input_all[j][i] == '#') {
				cc = i;
				break;
			}
		}

		//If a comment is found, we remove all characters including the '#'
		if (cc >= 0) {
			input_all[j][cc] = '\n';

			for (i = cc+1; i < slength-1; i++) {
				input_all[j][i] = '\0';
			}
		}
	}


	//Counting number of lines which are not empty
	nlines = nlines_all;
	for (i = 0; i < nlines_all; i++) {
		nlines -= IsEmpty(input_all[i]);
	}


	//Creating a new array to store only the lines which are not empty
	input = (char **) malloc(nlines * sizeof(char *));
	for (i = 0; i < nlines; i++) {
		input[i] = (char *) malloc(slength * sizeof(char));
	}

	//Creating an array which points to the line number of the input file
	input_loc = (int *) malloc(nlines * sizeof(int));


	//Copying lines and populating line number array
	k = 0;
	for (i = 0; i < nlines_all; i++) {
		if (!IsEmpty(input_all[i])) {
			for (j = 0; j < slength; j++) {
				input[k][j] = input_all[i][j];
			}
			input_loc[k] = i;
			k++;
		}
	}


	//Initializing arrays for parameters and args
	parameters = (char **) malloc(nlines * sizeof(char *));
	for (i = 0; i < nlines; i++) {
		parameters[i] = (char *) malloc(slength * sizeof(char));
	}

	args = (char **) malloc(nlines * sizeof(char *));
	for (i = 0; i < nlines; i++) {
		args[i] = (char *) malloc(slength * sizeof(char));
	}


	//Separating the lines into parameters and args
	for (i = 0; i < nlines; i++) {

		if (sscanf(input[i],"%s %s",parameters[i],args[i]) != 2) {
			if (myid == 0) {printf("Input line %d is not in the correct format. Ensure the line has one parameter and one argument\nEXITING\n",input_loc[i]+1);}
			MPI_Finalize();
			exit(0);
		}

	}

	//Freeing stored arrays
	for (i = 0; i < nlines_all; i++) {
		free(input_all[i]);
	}
	free(input_all);

	for (i = 0; i < nlines; i++) {
		free(input[i]);
	}
	free(input);

	free(input_loc);


	//Storing the parameters and args into the control struct
	control->input_nlines = nlines;
	control->input_params = parameters;
	control->input_args = args;

	//Printing lines for debug
	// for (i = 0; i < nlines; i++) {
	// 	printf("%s %s\n",control->input_params[i],control->input_args[i]);
	// }




	//Scanning lines for time stepping parameters
	timer->dt = GetInputValueDouble(control, 0, control->input_nlines, "time_step_size", 0.0, "MISSING: time_step_size\nEXITING\n", 1);

	//printf("dt = %.8e\n",timer->dt);


	timer->Nsteps = GetInputValueInt(control, 0, control->input_nlines, "number_steps", 0, "MISSING: number_steps\nEXITING\n", 1);

	timer->print_interval = GetInputValueInt(control, 0, control->input_nlines, "print_interval", 0, "MISSING: print_interval\nSet to default of 0 (none)\n", 0);

	timer->sort_interval = 0;

	timer->checkpoint_interval = GetInputValueInt(control, 0, control->input_nlines, "checkpoint_interval", 0, "MISSING: checkpoint_interval\nSet to default of 0 (none)\n", 0);


	//Scanning lines for models
	control->epsr = GetInputValueDouble(control, 0, control->input_nlines, "relative_permittivity", 1.0, "MISSING: relative_permittivity\nSetting to default 1.0\n", 0);
	//control->fsubc = GetInputValueInt(&str, 0, nlines, "field_subcycle", 1, "MISSING: Field sub-cycle setting\nSetting to default 1 (solve very time step)\n", 0);
	control->fsubc = 1;
	control->sid = GetInputValueInt(control, 0, control->input_nlines, "solver_id", 31, "MISSING: Field solver method (solver_id)\nSetting to default 31 - Hypre PFMG preconditioner with GMRES\n", 0);
	//grid->hypre_gpu = GetInputValueInt(control, 0, control->input_nlines, "hypre_gpu", 0, "MISSING: Field Hypre GPU setting (hypre_gpu)\nSetting to default 0 (use CPU solver)\n", 0);
	control->seed_input = GetInputValueInt(control, 0, control->input_nlines, "random_seed", -1, "Setting random seed from dev/urandom\n\n", 0);


	//Scanning the lines for number of blocks, species, plasmas, sources, collisions and diagnostics
	for (i = 0; i < nlines; i++) {

		if (strcmp(control->input_params[i],"species") == 0) {
			nspecies++;
		}
		if (strcmp(control->input_params[i],"plasma") == 0) {
			nplasma++;
		}
		if (strcmp(control->input_params[i],"source") == 0) {
			nsource++;
		}
		if (strcmp(control->input_params[i],"beam") == 0) {
			nbeam++;
		}
		if (strcmp(control->input_params[i],"elastic") == 0) {
			nelastic++;
		}
		if (strcmp(control->input_params[i],"excitation") == 0) {
			nexcite++;
		}
		if (strcmp(control->input_params[i],"ionization") == 0) {
			nionize++;
		}
		if (strcmp(control->input_params[i],"cxchange") == 0) {
			ncxchange++;
		}
		if (strcmp(control->input_params[i],"phase") == 0) {
			nphase++;
		}
		if (strcmp(control->input_params[i],"probe") == 0) {
			nprobe++;
		}
		if (strcmp(control->input_params[i],"secondary") == 0) {
			nsee++;
		}		
	}

	//Storing species and creation counts to relevant structures
	allspecies->num_species = nspecies;
	allcreation->num_plasma = nplasma;
	allcollisions->num_elastic = nelastic;
	allcollisions->num_excite = nexcite;
	allcollisions->num_ionize = nionize;
	allcollisions->num_cxchange = ncxchange;
	//allwallcollisions->num_see = nsee;
	output->num_phase = nphase;
	output->num_probe = nprobe;


	//Scanning the lines for output settings
	output->TotalParticlesFlag = GetInputValueInt(control, 0, control->input_nlines, "print_total_particles", 0, "MISSING: print_total_particles flag\nSet to default of 0 (none)\n", 0);
	output->TotalMomentumFlag = GetInputValueInt(control, 0, control->input_nlines, "print_total_momentum", 0, "MISSING: print_total_momentum flag\nSet to default of 0 (none)\n", 0);
	output->TotalEnergyFlag = GetInputValueInt(control, 0, control->input_nlines, "print_total_energy", 0, "MISSING: print_total_energy flag\nSet to default of 0 (none)\n", 0);


	//Allocating memory for the phase space diagnostics (these will be read in later)
	output->phase = (struct PhaseDiagnostic *) malloc(output->num_phase*sizeof(struct PhaseDiagnostic));

}




/*
FUNCTION DESCRIPTION: Removing memory allocated in the control struct
LAST EDIT: Created
EDIT DATE: 4/12/22
NAME OF EDITOR: Tasman Powis
*/
void FinalizeControl(struct Control *control)
{
	int i;
	int nlines = control->input_nlines;

	for (i = 0; i < nlines; i++)
	{
		free(control->input_params[i]);
		free(control->input_args[i]);
	}

	free(control->input_params);
	free(control->input_args);
}





/*
FUNCTION DESCRIPTION: Function which processes the species defined in the input file
LAST EDIT: Updated for new input format
EDIT DATE: 2/12/22
NAME OF EDITOR: Tasman Powis
*/
void ProcessingSpeciesInput(struct Control *control, struct AllSpecies *allspecies)
{
	int i, j, cc, nlines, slength = 128;
	char **str;
	int val, myid;
	int nspecies = allspecies->num_species, snum, *slines;
	FILE *finput;

	nlines = control->input_nlines;
	slines = (int *) malloc((nspecies + 1) * sizeof(int));

	//Scanning the lines for the location of "species" identifier. We then only search between each identifiers to determine properties.
	j = 0;
	for (i = 0; i < nlines; i++) {
		if (strcmp(control->input_params[i],"species") == 0) {
			slines[j] = i;
			j++;
		}
	}

	slines[nspecies] = nlines;



	//Scanning the lines for information on each species
	for (i = 0; i < nspecies; i++) {
		allspecies->species[i].snum = GetInputValueInt(control, slines[i], slines[i+1], "species", 0, "MISSING: Species number (species)\nEXITING\n", 1);
		allspecies->species[i].N = (ptype) GetInputValueDouble(control, slines[i], slines[i+1], "clumping", 1.0, "MISSING: Species clumping factor (clumping)\nEXITING\n", 1);
		allspecies->species[i].Q = (ptype) GetInputValueDouble(control, slines[i], slines[i+1], "charge", 1.0, "MISSING: Species charge (charge)\nEXITING\n", 1);
		allspecies->species[i].M = (ptype) GetInputValueDouble(control, slines[i], slines[i+1], "mass", 1.0, "MISSING: Species mass (mass)\nEXITING\n", 1);
		allspecies->species[i].mag = GetInputValueInt(control, slines[i], slines[i+1], "magnetized", 1, "MISSING: Species magnetization setting (magnetized)\nSetting to default 1 (magnetized)\n", 0);
		//allspecies->species[i].subc = GetInputValueInt(control, slines[i], slines[i+1], "subcycle", 1, "MISSING: Species subcycle setting (subcycle)\nSetting to default 1 (push every time step)\n", 0);
		allspecies->species[i].subc = GetInputValueInt(control, slines[i], slines[i+1], "subcycle", 1, "", 0);
		//allspecies->species[i].print_frac = GetInputValueDouble(control, slines[i], slines[i+1], "print_frac", 0.0, "MISSING: Species print_frac\nSetting to default of 0.0 (none)\n", 0); //Not needed now we have a phase space diagnostic. But might be needed in the future for particle tracing.
		allspecies->species[i].print_frac = GetInputValueDouble(control, slines[i], slines[i+1], "print_frac", 1.0, "", 0);
	}

	//Checking for repeated numbers. Exiting if found.
	for (i = 0; i < nspecies; i++)
	{
		int n1, n2;

		n1 = allspecies->species[i].snum;

		for (j = 0; j < nspecies; j++)
		{
			if (j == i) {continue;}

			n2 = allspecies->species[j].snum;

			if (n1 == n2) {
				if (myid == 0) {printf("Species number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}

	free(slines);

}



/*
FUNCTION DESCRIPTION: Function which processes the creation routines defined in the input file
LAST EDIT: Updated for new input format
EDIT DATE: 2/12/22
NAME OF EDITOR: Tasman Powis
*/
void ProcessingCreationInput(struct Control *control, struct AllCreation *allcreation, struct AllSpecies *allspecies, struct Timing *timer)
{
	int i, j, np, ni, nb, sf, nlines;
	int val, myid;
	int nplasma = allcreation->num_plasma;
	int ncreation = nplasma;
	int cnum, *clines, *cdefs;
	int nspecies = allspecies->num_species;
	FILE *finput;

	nlines = control->input_nlines;
	clines = (int *) malloc((ncreation + 1) * sizeof(int));
	cdefs = (int *) malloc((ncreation + 1) * sizeof(int));


	//Scanning the lines for the location of creation identifiers. We then only search between each identifier to determine properties.
	j = 0;
	for (i = 0; i < nlines; i++) {
		if (strcmp(control->input_params[i],"plasma") == 0) {
			clines[j] = i;
			cdefs[j] = 0;
			j++;
		}
		if (strcmp(control->input_params[i],"source") == 0) {
			clines[j] = i;
			cdefs[j] = 1;
			j++;
		}
		if (strcmp(control->input_params[i],"beam") == 0) {
			clines[j] = i;
			cdefs[j] = 2;
			j++;
		}
	}

	clines[ncreation] = nlines;


	//Scanning the lines for information on each creation identifier.
	np = 0;
	ni = 0;
	nb = 0;
	for (i = 0; i < ncreation; i++) {
		//printf("cdefs = %d\n",cdefs[i]);
		switch(cdefs[i]) {
			//Information for a plasma creation identifier
			case 0:
				//printf("Loading Plasma %d\n",np);
				allcreation->plasma[np].plasma_num = GetInputValueInt(control, clines[i], clines[i+1], "plasma", 0, "MISSING: plasma number\nEXITING\n", 1);
				allcreation->plasma[np].snum = GetInputValueInt(control, clines[i], clines[i+1], "plasma_species", 0, "MISSING: plasma_species\nEXITING\n", 1);
				allcreation->plasma[np].n0_plasma = GetInputValueDouble(control, clines[i], clines[i+1], "plasma_density", 0, "MISSING: plasma_density\nEXITING\n", 1);
				allcreation->plasma[np].xmin_glob = GetInputValueDouble(control, clines[i], clines[i+1], "plasma_xmin", 0.0, "MISSING: plasma_xmin\nEXITING\n", 1);
				allcreation->plasma[np].xmax_glob = GetInputValueDouble(control, clines[i], clines[i+1], "plasma_xmax", 0.0, "MISSING: plasma_xmax\nEXITING\n", 1);
				allcreation->plasma[np].ymin_glob = GetInputValueDouble(control, clines[i], clines[i+1], "plasma_ymin", 0.0, "MISSING: plasma_ymin\nEXITING\n", 1);
				allcreation->plasma[np].ymax_glob = GetInputValueDouble(control, clines[i], clines[i+1], "plasma_ymax", 0.0, "MISSING: plasma_ymax\nEXITING\n", 1);
				allcreation->plasma[np].T0 = GetInputValueDouble(control, clines[i], clines[i+1], "plasma_temperature", 0.0, "MISSING: plasma_temperature\nSetting to 0.0\n", 0);
				allcreation->plasma[np].Kx = (ptype) GetInputValueDouble(control, clines[i], clines[i+1], "plasma_kinetic_x", 0.0, "MISSING: plasma_kinetic_x\nSetting to 0.0\n", 0);
				allcreation->plasma[np].Ky = (ptype) GetInputValueDouble(control, clines[i], clines[i+1], "plasma_kinetic_y", 0.0, "MISSING: plasma_kinetic_y\nSetting to 0.0\n", 0);
				allcreation->plasma[np].Kz = (ptype) GetInputValueDouble(control, clines[i], clines[i+1], "plasma_kinetic_z", 0.0, "MISSING: plasma_kinetic_z\nSetting to 0.0\n", 0);

				sf = 0;
				for (j = 0; j < nspecies; j++) {
					if (allcreation->plasma[np].snum == allspecies->species[j].snum) {
						allcreation->plasma[np].sindex = j;
						sf = 1;
					}
				}

				if (sf == 0) {
					if (myid == 0) {printf("Plasma %d species number could not be found\nEXITING\n",allcreation->plasma[np].plasma_num);}
					MPI_Finalize();
					exit(0);
				}

				np++;
				break;
		}
	}


	//printf("nplasma = %d, nsource = %d\n",np,ni);
	free(clines);
	free(cdefs);


	//Checking for repeated numbers. Exiting if found. Plasmas
	for (i = 0; i < nplasma; i++)
	{
		int n1, n2;

		n1 = allcreation->plasma[i].plasma_num;

		for (j = 0; j < nplasma; j++)
		{
			if (j == i) {continue;}

			n2 = allcreation->plasma[j].plasma_num;

			if (n1 == n2) {
				if (myid == 0) {printf("Plasma number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}

}







/*
FUNCTION DESCRIPTION: Function which processes the collision modules defined in the input file
LAST EDIT: Created
EDIT DATE: 2/3/2020
NAME OF EDITOR: Tasman Powis
*/
void ProcessingCollisionInput(struct Control *control, struct AllCollisions *allcollisions, struct AllSpecies *allspecies)
{
	int i, j, cc, sf, sfe, sfi, nel, nex, nio, ncx, nlines;
	int val, myid;
	int nelastic = allcollisions->num_elastic;
	int nexcite = allcollisions->num_excite;
	int nionize = allcollisions->num_ionize;
	int ncxchange = allcollisions->num_cxchange;
	int ncoll = nelastic + nexcite + nionize + ncxchange;
	int cnum, *clines;
	int snum;
	double smass;
	ctype *cdefs;
	int nspecies = allspecies->num_species;
	FILE *finput;

	nlines = control->input_nlines;
	clines = (int *) malloc((ncoll + 1) * sizeof(int));
	cdefs = (ctype *) malloc((ncoll + 1) * sizeof(ctype));

	//Scanning the lines for the location of collision identifiers. We then only search between each identifier to determine properties.
	j = 0;
	for (i = 0; i < nlines; i++) {
		if (strcmp(control->input_params[i],"elastic") == 0) {
			clines[j] = i;
			cdefs[j] = ELASTIC;
			j++;
		}
		if (strcmp(control->input_params[i],"excitation") == 0) {
			clines[j] = i;
			cdefs[j] = EXCITE;
			j++;
		}
		if (strcmp(control->input_params[i],"ionization") == 0) {
			clines[j] = i;
			cdefs[j] = IONIZE;
			j++;
		}
		if (strcmp(control->input_params[i],"cxchange") == 0) {
			clines[j] = i;
			cdefs[j] = CXCHANGE;
			j++;
		}
	}

	clines[ncoll] = nlines;


	//Scanning the lines for information on each creation identifier.
	nel = 0;
	nex = 0;
	nio = 0;
	ncx = 0;

	for (i = 0; i < ncoll; i++) {
		//printf("cdefs = %d\n",cdefs[i]);
		switch(cdefs[i]) {
			//Information for an elastic collision identifier
			case ELASTIC:
				//printf("Loading Elastic %d\n",np);
				allcollisions->elastic[nel].col_num = GetInputValueInt(control, clines[i], clines[i+1], "elastic", 0, "MISSING: elastic number\nEXITING\n", 1);
				allcollisions->elastic[nel].species_in_num = GetInputValueInt(control, clines[i], clines[i+1], "species_in", 0, "MISSING: elastic species_in\nEXITING\n", 1);
				allcollisions->elastic[nel].neutral_n0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_density", 0.0, "MISSING: elastic neutral_density\nEXITING\n", 1);
				allcollisions->elastic[nel].neutral_T0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_temperature", 0.0, "MISSING: elastic neutral_temperature\nEXITING\n", 1);
				//allcollisions->elastic[nel].neutral_M0 = GetInputValueDouble(&str, clines[i], clines[i+1], "neutral_mass", 0.0, "MISSING: elastic neutral_mass\nEXITING\n", 1);
				allcollisions->elastic[nel].xsection_file = GetInputValueStr(control, clines[i], clines[i+1], "xsection_file", "MISSING: elastic xsection_file name\nEXITING\n", 1);

				sf = 0;
				for (j = 0; j < nspecies; j++) {
					if (allcollisions->elastic[nel].species_in_num == allspecies->species[j].snum) {
						allcollisions->elastic[nel].species_in_index = j;
						sf = 1;
					}
				}

				if (sf == 0) {
					if (myid == 0) {printf("Elastic %d species_in number could not be found\nEXITING\n",allcollisions->elastic[nel].col_num);}
					MPI_Finalize();
					exit(0);
				}

				snum = allcollisions->elastic[nel].species_in_index;
				smass = allspecies->species[snum].M;

				if (smass == 1.0) {
					allcollisions->elastic[nel].neutral_M0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_mass", 0.0, "MISSING: electron-neutral elastic collision neutral_mass\nEXITING\n", 1);
				}

				nel++;
				break;

			//Information for an excitation collision identifier
			case EXCITE:
				//printf("Loading Excitation %d\n",np);
				allcollisions->excite[nex].col_num = GetInputValueInt(control, clines[i], clines[i+1], "excitation", 0, "MISSING: excitation number\nEXITING\n", 1);
				allcollisions->excite[nex].species_in_num = GetInputValueInt(control, clines[i], clines[i+1], "species_in", 0, "MISSING: excitation species_in\nEXITING\n", 1);
				allcollisions->excite[nex].neutral_n0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_density", 0.0, "MISSING: excitation neutral_density\nEXITING\n", 1);
				allcollisions->excite[nex].neutral_T0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_temperature", 0.0, "MISSING: excitation neutral_temperature\nEXITING\n", 1);
				allcollisions->excite[nex].neutral_M0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_mass", 0.0, "MISSING: excitation neutral_mass\nEXITING\n", 1);
				allcollisions->excite[nex].neutral_wex = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_excitation_energy", 0.0, "MISSING: excitation neutral_excitation_energy\nEXITING\n", 1);
				allcollisions->excite[nex].xsection_file = GetInputValueStr(control, clines[i], clines[i+1], "xsection_file", "MISSING: excitation xsection_file name\nEXITING\n", 1);

				sf = 0;
				for (j = 0; j < nspecies; j++) {
					if (allcollisions->excite[nex].species_in_num == allspecies->species[j].snum) {
						allcollisions->excite[nex].species_in_index = j;
						sf = 1;
					}
				}

				if (sf == 0) {
					if (myid == 0) {printf("Excitation %d species_in number could not be found\nEXITING\n",allcollisions->excite[nex].col_num);}
					MPI_Finalize();
					exit(0);
				}

				nex++;
				break;

			//Information for an ionization collision identifier
			case IONIZE:
				//printf("Loading Ionization %d\n",np);
				allcollisions->ionize[nio].col_num = GetInputValueInt(control, clines[i], clines[i+1], "ionization", 0, "MISSING: ionization number\nEXITING\n", 1);
				allcollisions->ionize[nio].electron_in_num = GetInputValueInt(control, clines[i], clines[i+1], "electron_in", 0, "MISSING: ionization electron_in\nEXITING\n", 1);
				//allcollisions->ionize[nio].electron_out_num = GetInputValueInt(&str, clines[i], clines[i+1], "electron_out", 0, "MISSING: ionization electron_out\nEXITING\n", 1);
				allcollisions->ionize[nio].ion_out_num = GetInputValueInt(control, clines[i], clines[i+1], "ion_out", 0, "MISSING: ionization ion_out\nEXITING\n", 1);
				allcollisions->ionize[nio].neutral_n0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_density", 0.0, "MISSING: ionization neutral_density\nEXITING\n", 1);
				allcollisions->ionize[nio].neutral_T0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_temperature", 0.0, "MISSING: ionization neutral_temperature\nEXITING\n", 1);
				//allcollisions->ionize[nio].neutral_M0 = GetInputValueDouble(&str, clines[i], clines[i+1], "neutral_mass", 0.0, "MISSING: ionization neutral_mass\nEXITING\n", 1);
				allcollisions->ionize[nio].neutral_wion = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_ionization_energy", 0.0, "MISSING: ionization neutral_ionization_energy\nEXITING\n", 1);
				allcollisions->ionize[nio].bion = GetInputValueDouble(control, clines[i], clines[i+1], "ionization_b_value", 10.0, "MISSING: ionization_b_value. Setting to 10.0, the value for Argon ionization.\n", 0);
				allcollisions->ionize[nio].xsection_file = GetInputValueStr(control, clines[i], clines[i+1], "xsection_file", "MISSING: ionization xsection_file name\nEXITING\n", 1);

				sfe = 0;
				sfi = 0;
				for (j = 0; j < nspecies; j++) {
					if (allcollisions->ionize[nio].electron_in_num == allspecies->species[j].snum) {
						allcollisions->ionize[nio].electron_in_index = j;
						sfe = 1;
					}
					// if (allcollisions->ionize[nio].electron_out_num == allspecies->species[j].snum) {
					// 	allcollisions->ionize[nio].electron_out_index = j;
					// }
					if (allcollisions->ionize[nio].ion_out_num == allspecies->species[j].snum) {
						allcollisions->ionize[nio].ion_out_index = j;
						sfi = 1;
					}
				}

				if (sfe == 0) {
					if (myid == 0) {printf("Ionization %d electron_in number could not be found\nEXITING\n",allcollisions->ionize[nio].col_num);}
					MPI_Finalize();
					exit(0);
				}

				if (sfi == 0) {
					if (myid == 0) {printf("Ionization %d ion_out number could not be found\nEXITING\n",allcollisions->ionize[nio].col_num);}
					MPI_Finalize();
					exit(0);
				}

				//Electron out must be the same as electron in. We do not yet have an algorithm implemented to split particles in the event of unequal electron weights.
				allcollisions->ionize[nio].electron_out_num = allcollisions->ionize[nio].electron_in_num;
				allcollisions->ionize[nio].electron_out_index = allcollisions->ionize[nio].electron_in_index;

				nio++;
				break;

			//Information for an excitation collision identifier
			case CXCHANGE:
				//printf("Loading Excitation %d\n",np);
				allcollisions->cxchange[ncx].col_num = GetInputValueInt(control, clines[i], clines[i+1], "cxchange", 0, "MISSING: cxchange number\nEXITING\n", 1);
				allcollisions->cxchange[ncx].species_in_num = GetInputValueInt(control, clines[i], clines[i+1], "ion_in", 0, "MISSING: charge exchange ion_in\nEXITING\n", 1);
				allcollisions->cxchange[ncx].neutral_n0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_density", 0.0, "MISSING: charge exchange neutral_density\nEXITING\n", 1);
				allcollisions->cxchange[ncx].neutral_T0 = GetInputValueDouble(control, clines[i], clines[i+1], "neutral_temperature", 0.0, "MISSING: charge exchange neutral_temperature\nEXITING\n", 1);
				//allcollisions->cxchange[ncx].neutral_M0 = GetInputValueDouble(&str, clines[i], clines[i+1], "neutral_mass", 0.0, "MISSING: charge exchange neutral_mass\nEXITING\n", 1);
				allcollisions->cxchange[ncx].xsection_file = GetInputValueStr(control, clines[i], clines[i+1], "xsection_file", "MISSING: charge exchange xsection_file name\nEXITING\n", 1);

				sf = 0;
				for (j = 0; j < nspecies; j++) {
					if (allcollisions->cxchange[ncx].species_in_num == allspecies->species[j].snum) {
						allcollisions->cxchange[ncx].species_in_index = j;
						sf = 1;
					}
				}

				if (sf == 0) {
					if (myid == 0) {printf("CXchange %d ion_in number could not be found\nEXITING\n",allcollisions->cxchange[ncx].col_num);}
					MPI_Finalize();
					exit(0);
				}

				ncx++;
				break;

		}
	}

	//printf("nelastic = %d, nexcite = %d, nionize = %d\n",nel,nex,nio);
	free(clines);
	free(cdefs);


	//Checking for repeated numbers. Exiting if found. Elastic
	for (i = 0; i < nelastic; i++)
	{
		int n1, n2;

		n1 = allcollisions->elastic[i].col_num;

		for (j = 0; j < nelastic; j++)
		{
			if (j == i) {continue;}

			n2 = allcollisions->elastic[j].col_num;

			if (n1 == n2) {
				if (myid == 0) {printf("Elastic number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}


	//Checking for repeated numbers. Exiting if found. Excitation
	for (i = 0; i < nexcite; i++)
	{
		int n1, n2;

		n1 = allcollisions->excite[i].col_num;

		for (j = 0; j < nexcite; j++)
		{
			if (j == i) {continue;}

			n2 = allcollisions->excite[j].col_num;

			if (n1 == n2) {
				if (myid == 0) {printf("Excitation number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}


	//Checking for repeated numbers. Exiting if found. Ionization
	for (i = 0; i < nionize; i++)
	{
		int n1, n2;

		n1 = allcollisions->ionize[i].col_num;

		for (j = 0; j < nionize; j++)
		{
			if (j == i) {continue;}

			n2 = allcollisions->ionize[j].col_num;

			if (n1 == n2) {
				if (myid == 0) {printf("Ionization number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}


	//Checking for repeated numbers. Exiting if found. Cxchange
	for (i = 0; i < ncxchange; i++)
	{
		int n1, n2;

		n1 = allcollisions->cxchange[i].col_num;

		for (j = 0; j < ncxchange; j++)
		{
			if (j == i) {continue;}

			n2 = allcollisions->cxchange[j].col_num;

			if (n1 == n2) {
				if (myid == 0) {printf("CXchange number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}


}








/*
FUNCTION DESCRIPTION: Function which processes the phases space diagnostics defined in the input file
LAST EDIT: Created
EDIT DATE: 16/9/2022
NAME OF EDITOR: Tasman Powis
*/
void ProcessingPhaseDiagnosticInput(struct Control *control, struct Output *output, struct AllSpecies *allspecies)
{
	int i, j, cc, sf, nlines;
	int val, myid;
	int nphase = output->num_phase, snum, *slines;
	int nspecies = allspecies->num_species;
	FILE *finput;

	nlines = control->input_nlines;
	slines = (int *) malloc((nphase + 1) * sizeof(int));

	//Scanning the lines for the location of "phase" identifier. We then only search between each identifiers to determine properties.
	j = 0;
	for (i = 0; i < nlines; i++) {
		if (strcmp(control->input_params[i],"phase") == 0) {
			slines[j] = i;
			j++;
		}
	}

	slines[nphase] = nlines;


	//Scanning the lines for information on each species
	for (i = 0; i < nphase; i++) {
		output->phase[i].pnum = GetInputValueInt(control, slines[i], slines[i+1], "phase", 0, "MISSING: Phase diagnostic number (phase)\nEXITING\n", 1);
		output->phase[i].snum = GetInputValueInt(control,slines[i], slines[i+1], "phase_species", 0, "MISSING: phase_species\nEXITING\n", 1);
		output->phase[i].print_interval = GetInputValueInt(control, slines[i], slines[i+1], "phase_interval", 0, "MISSING: phase_interval\nEXITING\n", 1);
		output->phase[i].xmin = GetInputValueDouble(control, slines[i], slines[i+1], "phase_xmin", 0.0, "MISSING: phase_xmin\nEXITING\n", 1);
		output->phase[i].xmax = GetInputValueDouble(control,slines[i], slines[i+1], "phase_xmax", 0.0, "MISSING: phase_xmax\nEXITING\n", 1);
		output->phase[i].ymin = GetInputValueDouble(control, slines[i], slines[i+1], "phase_ymin", 0.0, "MISSING: phase_ymin\nEXITING\n", 1);
		output->phase[i].ymax = GetInputValueDouble(control,slines[i], slines[i+1], "phase_ymax", 0.0, "MISSING: phase_ymax\nEXITING\n", 1);

		sf = 0;
		for (j = 0; j < nspecies; j++) {
			if (output->phase[i].snum == allspecies->species[j].snum) {
				output->phase[i].sindex = j;
				sf = 1;
			}
		}

		if (sf == 0) {
			if (myid == 0) {printf("Phase diagnostic %d species number could not be found\nEXITING\n",output->phase[i].pnum);}
			MPI_Finalize();
			exit(0);
		}

	}

	free(slines);


	//Checking for repeated numbers. Exiting if found.
	for (i = 0; i < nphase; i++)
	{
		int n1, n2;

		n1 = output->phase[i].pnum;

		for (j = 0; j < nphase; j++)
		{
			if (j == i) {continue;}

			n2 = output->phase[j].pnum;

			if (n1 == n2) {
				if (myid == 0) {printf("Phase number %d is repeated\nEXITING\n",n1);}
				MPI_Finalize();
				exit(0);
			}
		}
	}


}



/*
FUNCTION DESCRIPTION: Checks that the file exists and counts the number of un-commented lines.
LAST EDIT: Created
EDIT DATE: 2/3/2021
NAME OF EDITOR: Tasman Powis
*/
int GetXSectionFileNumLines(char *fname)
{
	int i, flen, sl = 256;
	int nl;
	char str[sl];
	FILE *finput;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Loading collision cross-section file.
	finput = fopen(fname,"r");

	if (finput == NULL) {
		printf("\n");
		if (myid == 0) {printf("Collision input file \"%s\" could not be found\nEXITING\n",fname);}
		MPI_Finalize();
		exit(0);
	}

	//Get the number of lines in the file
	flen = CountLines(finput);
	rewind(finput);


	//Counting the number of lines which are not comments
	nl = 0;
	for (i = 0; i < flen; i++)
	{
		fgets(str,sl,finput);

		if (str[0] == '#')
			continue;

		nl++;
	}

	fclose(finput);

	return nl+1;
}





/*
FUNCTION DESCRIPTION: Get the cross-section data from file, ignoring comments. Exits if the data is in the wrong format.
LAST EDIT: Created
EDIT DATE: 2/3/2021
NAME OF EDITOR: Tasman Powis
*/
void GetXSectionFileData(char *fname, double *E, double *XS)
{
	int i, flen, sl = 256;
	int myid;
	int nc;
	double a, b;
	char str[sl];
	FILE *finput;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Load the file
	finput = fopen(fname,"r");

	//Get the number of lines in the file
	flen = CountLines(finput);
	rewind(finput);

	//Reading each line in turn (skipping those which are comments only)
	nc = 0;
	for (i = 0; i < flen+1; i++)
	{
		fgets(str,sl,finput);

		if (str[0] == '#')
			continue;



		if(sscanf(str,"%lf %lf", &E[nc], &XS[nc]) != 2) {
			printf("\n");
			if (myid == 0) {printf("Collision input file \"%s\" is not formatted correctly\nEXITING\n",fname);}
			MPI_Finalize();
			exit(0);
		}

		nc++;

	}

	fclose(finput);

}







/*
FUNCTION DESCRIPTION: Reads an the input between lines "istart" and "iend" and searches for a matching parameter string, then the corresponding argument string is converted to the appropriate format. If the functions fails to identify the line it will spit out the "fail_message" and then execute the "fail response" which is either to continue or stop the program. If the program continues the value will be that of "def" or default.
LAST EDIT: Updated for new input read in format
EDIT DATE: 2/11/2022
NAME OF EDITOR: Tasman Powis
*/
double GetInputValueDouble(struct Control *control, int istart, int iend, char const *line, double def, char const *fail_message, int fail_response)
{
	int i, loc = -1;
	bool fail;
	double val;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Locating the string in the parameters
	for (i = istart; i < iend; i++)
	{
		if(strcmp(control->input_params[i],line) == 0) {
			loc = i;
			break;
		}
	}

	//Exiting if the value parameter is not found
	if(loc < 0) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
		return val;
	}

	//Converting the argument string to the correct format and checking for failures
	if(sscanf(control->input_args[loc],"%lf",&val) != 1) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
	}

	return val;
}





/*
FUNCTION DESCRIPTION: Reads an the input between lines "istart" and "iend" and searches for a matching parameter string, then the corresponding argument string is converted to the appropriate format. If the functions fails to identify the line it will spit out the "fail_message" and then execute the "fail response" which is either to continue or stop the program. If the program continues the value will be that of "def" or default.
LAST EDIT: Updated for new input read in format
EDIT DATE: 2/11/2022
NAME OF EDITOR: Tasman Powis
*/
float GetInputValueFloat(struct Control *control, int istart, int iend, char const *line, float def, char const *fail_message, int fail_response)
{
	int i, loc = -1;
	bool fail;
	float val;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Locating the string in the parameters
	for (i = istart; i < iend; i++)
	{
		if(strcmp(control->input_params[i],line) == 0) {
			loc = i;
			break;
		}
	}

	//Exiting if the value parameter is not found
	if(loc < 0) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
		return val;
	}

	//Converting the argument string to the correct format and checking for failures
	if(sscanf(control->input_args[loc],"%f",&val) != 1) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
	}

	return val;
}





/*
FUNCTION DESCRIPTION: Reads an the input between lines "istart" and "iend" and searches for a matching parameter string, then the corresponding argument string is converted to the appropriate format. If the functions fails to identify the line it will spit out the "fail_message" and then execute the "fail response" which is either to continue or stop the program. If the program continues the value will be that of "def" or default.
LAST EDIT: Updated for new input read in format
EDIT DATE: 2/11/2022
NAME OF EDITOR: Tasman Powis
*/
int GetInputValueInt(struct Control *control, int istart, int iend, char const *line, int def, char const *fail_message, int fail_response)
{
	int i, loc = -1;
	bool fail;
	int val;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Locating the string in the parameters
	for (i = istart; i < iend; i++)
	{
		if(strcmp(control->input_params[i],line) == 0) {
			loc = i;
			break;
		}
	}

	//Exiting if the value parameter is not found
	if(loc < 0) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
		return val;
	}

	//Converting the argument string to the correct format and checking for failures
	if(sscanf(control->input_args[loc],"%d",&val) != 1) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
	}

	return val;
}






/*
FUNCTION DESCRIPTION: Reads an the input between lines "istart" and "iend" and searches for a matching parameter string, then the corresponding argument string is converted to the appropriate format. If the functions fails to identify the line it will spit out the "fail_message" and then execute the "fail response" which is either to continue or stop the program. If the program continues the value will be that of "def" or default.
LAST EDIT: Updated for new input read in format
EDIT DATE: 2/11/2022
NAME OF EDITOR: Tasman Powis
*/
long long int GetInputValueLongInt(struct Control *control, int istart, int iend, char const *line, long long int def, char const *fail_message, int fail_response)
{
	int i, loc = -1;
	bool fail;
	long long int val;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//Locating the string in the parameters
	for (i = istart; i < iend; i++)
	{
		if(strcmp(control->input_params[i],line) == 0) {
			loc = i;
			break;
		}
	}

	//Exiting if the value parameter is not found
	if(loc < 0) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
		return val;
	}

	//Converting the argument string to the correct format and checking for failures
	if(sscanf(control->input_args[loc],"%lld",&val) != 1) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
	}

	return val;
}





/*
FUNCTION DESCRIPTION: Reads an the input between lines "istart" and "iend" and searches for a matching parameter string, then the corresponding argument string is converted to the appropriate format. If the functions fails to identify the line it will spit out the "fail_message" and then execute the "fail response" which is either to continue or stop the program. If the program continues the value will be that of "def" or default.
LAST EDIT: Updated for new input read in format
EDIT DATE: 2/11/2022
NAME OF EDITOR: Tasman Powis
*/
char* GetInputValueStr(struct Control *control, int istart, int iend, char const *line, char const *fail_message, int fail_response)
{
	int i, loc = -1;
	bool fail;
	char *val;
	char *def;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	def = (char *) calloc(1,sizeof(char));

	//Locating the string in the parameters
	for (i = istart; i < iend; i++)
	{
		if(strcmp(control->input_params[i],line) == 0) {
			loc = i;
			break;
		}
	}

	//Exiting if the value parameter is not found
	if(loc < 0) {
		if (myid == 0) {printf("%s",fail_message);}
		if (fail_response) {MPI_Finalize(); exit(0);}
		if (!fail_response) {val = def;}
		return val;
	}

	val = control->input_args[loc];

	return val;
}





/*
FUNCTION DESCRIPTION: Counts the number of lines in "finput"
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
int CountLines(FILE *finput)
{
	int nl = 0;
	int c;

	for (c = getc(finput); c != EOF; c = getc(finput)) {
		if (c == '\n') {
			nl++;
		}
	}

	return nl;

}
