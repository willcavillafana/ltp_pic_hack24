/* Header file for PIC code */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "mpi.h"
#include "desprng.h"
#ifdef _OPENACC
#include "openacc.h"
#endif

//Checking for CUDA-aware-MPI
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  #define CUDA_AWARE_MPI
#endif /* MPIX_CUDA_AWARE_SUPPORT */

//Defining pi
#ifdef M_PI
  #define PI M_PI
#else
  #define PI 3.14159265358979
#endif

//Setting physical ratios
#define e_div_me 1.75889402239543308815e11
#define me_div_e 5.68563010356572208005e-12
#define e_div_eps0 1.80937165798621600464e-8
#define mele 9.10938356e-31
#define qele 1.60217662e-19


//Setting the precision for the particles - double precision is default
#ifdef SINGLE_PRECISION
typedef float ptype;
#else
typedef double ptype;
#endif


//Setting the number of dimensions global variable
#define Npart_init 100000
#define Npart_buff 100000
#define Npart_buff_corner 2500



//Macros for Quicksort algorithm
#define pivot_index() (begin+(end-begin)/2)
#define swap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))


//Global PRNG seed - created from "true" random numbers
extern unsigned long global_seed;

//Type definitions for the collisions.
typedef enum {ELASTIC, EXCITE, IONIZE, CXCHANGE} ctype;


//Control structure
struct Control {
	double epsr;
	int sid;
	//Restart logic: 0 - fresh sart, 1 - restarting from checkpoint
	int restart;
	int fsubc; //Field solve every this number of steps
	char input[100];
	char **input_params, **input_args;
	int input_nlines;
	int seed_input;
};



//Structure for all species
struct AllSpecies {
	int num_species; //Total number of species
	struct Species *species;

};

//Structure for all creation methods
struct AllCreation {
	int num_plasma;  //Total number of plasma creation routines
	struct CreatePlasma *plasma;
};



//Structure for a species buffer
struct Species_buffer {

	//Send buffers for particles
	ptype *dtf_send;
	ptype *vx_send, *vy_send, *vz_send;
	ptype *x_send;
	short int *xbc_send;
	ptype *y_send;
	short int *ybc_send;
	int *cell_send;
	unsigned short int *tag_send, *print_tag_send;


	//Receive buffers for particles
	ptype *dtf_recv;
	ptype *vx_recv, *vy_recv, *vz_recv;
	ptype *x_recv;
	short int *xbc_recv;
	ptype *y_recv;
	short int *ybc_recv;
	int *cell_recv;
	unsigned short int *tag_recv, *print_tag_recv;


	//Send and receive particle counts
	int npart_send, npart_recv, npart_keep;
	int npart_recv_new, npart_keep_new;
};




//Species structure in 2D
struct Species {
	int myid, myid_region;
	int nproc, nproc_region;
	//MPI_Comm comm_part;
	//int rid;

	char name; //Name of species
	int snum; //Species number
	ptype N; //Clumping parameter.
	ptype Q; //Charge in units of electron charge.
	ptype M; //Mass in units of electron mass.
	int mag; //0 if species is unmagnetized, 1 if it is magnetized
	int subc; //Push every this number of steps
	long long Npart; //Total global number of particles in this species structure
	long long Npart_mem; //Total global number of particles which can fit in to allocated memory within this species structure
	int Npart_region; //Total region number of particle in this species structure
	int Npart_region_mem; //Total region number of particle which can fit in to allocated memory within this species structure
	int Npart_local; //Number of active particles on the local process
	int Npart_local_mem; //Number of memory blocks allocated for particles on this process
	double print_frac;

	ptype *vx, *vy, *vz;
	//ptype *vxh, *vyh, *vzh;
	double *Bx, *By, *Bz;
	ptype *x;
	ptype *Ex;
	short int *xbc;
	ptype *y;
	ptype *Ey;
	short int *ybc;
	int *cell, *cell_stat;
	unsigned short int *tag, *print_tag;
	int *pbc_keep;
	int *pbc_test;
	ptype *pbc_dtf;
	int *pbc_rem;
	int npart_keep;
	int npart_test;
	int npart_rem;
	int *rem_left;
	int *keep_right;
	int *pindex;

	int *coll_index; //Array for storing which particles will collide
	int coll_index_len;


	//Lists of collisions which pertain to this species
	int cnum;
	int *cindex;
	ctype *cdefs;

	double numax_tot;
	double pmax_tot;


	//DESPRNG for each species
	unsigned long seed;
	unsigned long rncount;
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;


	//Particle source buffers
	struct Species_buffer *buff_source;

	//Particle boundary condition buffers
	struct Species_buffer *buff;
};




//Structure storing information for initial plasma creation of given species
struct CreatePlasma {
	int plasma_num;
	int snum;
	int sindex;
	double n0_plasma;
	long long int Npart_plasma;
	long long int Npart_plasma_region;
	int Npart_plasma_local;
	double T0;
	double xmin_glob, xmax_glob;
	double xmin, xmax, xwidth, xcenter;
	double ymin_glob, ymax_glob;
	double ymin, ymax, ywidth, ycenter;
	ptype Kx, Ky, Kz; //Initial uniform plasma energies
	ptype Vx, Vy, Vz; //Initial uniform plasma velocities

	//DESPRNG for each plasma creation routine
	unsigned long seed;
	unsigned long *rncount;
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;
};



//Structure for all collision modules
struct AllCollisions {
	int num_elastic;
	int num_excite;
	int num_ionize;
	int num_cxchange;

	struct ElasticCollision *elastic;
	struct ExcitationCollision *excite;
	struct IonizationCollision *ionize;
	struct ChargeExchangeCollision *cxchange;
};




//Structure for elastic collisions
struct ElasticCollision {
	int col_num;

	int species_in_num;
	int species_in_index;

	double neutral_n0;
	double neutral_T0;
	double neutral_M0;

	char *xsection_file;
	double *energy_vals;
	double *xsection_vals;
	int nlines;

	//DESPRNG for each collision module
	unsigned long seed;
	unsigned long rncount;
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;
};



struct ExcitationCollision {
	int col_num;

	int species_in_num;
	int species_in_index;

	double neutral_n0;
	double neutral_T0;
    double neutral_M0;
    double neutral_wex;


	char *xsection_file;
	double *energy_vals;
	double *xsection_vals;
	int nlines;

	//DESPRNG for each collision module
	unsigned long seed;
	unsigned long rncount;
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;
};




struct ChargeExchangeCollision {
	int col_num;

	int species_in_num;
	int species_in_index;

	double neutral_n0;
	double neutral_T0;
    double neutral_M0;

    char *xsection_file;
	double *energy_vals;
	double *xsection_vals;
	int nlines;

	//DESPRNG for each collision module
	unsigned long seed;
	unsigned long rncount;
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;
};



//Structure for elastic collisions
struct IonizationCollision {
	int col_num;

	int electron_in_num;
	int electron_in_index;

	int electron_out_num;
	int electron_out_index;
	int ion_out_num;
	int ion_out_index;

	int enew_int, inew_int;
	double enew_frac, inew_frac;

	double neutral_n0;
	double neutral_T0;
    double neutral_M0;
    double neutral_wion;
	double bion;

	char *xsection_file;
	double *energy_vals;
	double *xsection_vals;
	int nlines;

	//DESPRNG for each collision module
	unsigned long seed;
	unsigned long rncount;
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;
};



//Timing structure
struct Timing {
	//Simulation timing variables
	int step, start_step, Nsteps, print_interval, sort_interval, checkpoint_interval, print_step, diagnostic_step;
	double dt;

	//Accumulated times
	double time, timeGrid, timeGridCommReduce, timeGridCommEdges, timeInterpToGrid, timeInterpToPart, timePrepareRHS, timeFieldSolve, timeGetSolution, timePush, timeSource, timePartBC, timePartComm, timePartSort, timeDiagnose, timeBarrier, timeMemCopy, timeCollide;
	double time_avg, timeGrid_avg, timeGridCommReduce_avg, timeGridCommEdges_avg, timeInterpToGrid_avg, timeInterpToPart_avg, timePrepareRHS_avg, timeFieldSolve_avg, timeGetSolution_avg, timePush_avg, timeSource_avg, timePartBC_avg, timePartComm_avg, timePartSort_avg, timeDiagnose_avg, timeBarrier_avg, timeMemCopy_avg, timeCollide_avg;

	double timePartBC_Remove, timePartBC_Replace, timePartBC_Sort, timePartBC_Communicate, timePartBC_Combine, timePartBC_Store;
	double timePartBC_Remove_avg, timePartBC_Replace_avg, timePartBC_Sort_avg, timePartBC_Communicate_avg, timePartBC_Combine_avg, timePartBC_Store_avg;

	double timePartBC_Test, timePartBC_SortRecv, timePartBC_first, timePartBC_second;
	double timePartBC_Test_avg, timePartBC_SortRecv_avg, timePartBC_first_avg, timePartBC_second_avg;

	FILE *flog;

};



//Phase space diagnostic struct
struct PhaseDiagnostic {
	int pnum;
	int snum;
	int sindex;
	int print_interval;
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;
	//double print_frac;
};



//Output structure
struct Output {
	//Output data flags
	//int ParticlePhaseFlag;
	int DensityFlag;
	int ChargeDensityFlag;
	int PotentialFlag;
	int ElectricFieldFlag;
	int CurrentDensityFlag;
	int TemperatureFlag;
	int VelocityMomentsFlag;
	int TotalParticlesFlag;
	int TotalMomentumFlag;
	int TotalEnergyFlag;

	//Arrays for storing simulation totals
	double *NP, *PX, *PY, *PZ, *KE, PE;

	//Array of phases space diagnostics
	int num_phase;
	struct PhaseDiagnostic *phase;

	//Array of probe diagnostics
	int num_probe;
	struct ProbeDiagnostic *probe;

	//Size and extents for MPI parallel IO
	MPI_Datatype region_array;
	MPI_Datatype print_array;
	//MPI_Datatype region_array_cells;
};











/********** FUNCTION REFERENCES **********/


//Input functions
int CountLines(FILE *finput);
const char *get_filename_ext(char *filename);
//char *get_filename_ext(char *filename);
void GetInputFileName(int *argc, char ***argv, struct Control *control);
void ProcessingInput(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *allcreation, struct AllCollisions *allcollisions, struct Timing *timer, struct Output *output);
void ProcessingSpeciesInput(struct Control *control, struct AllSpecies *allspecies);
void ProcessingCreationInput(struct Control *control, struct AllCreation *allcreation, struct AllSpecies *allspecies, struct Timing *timer);
void ProcessingCollisionInput(struct Control *control, struct AllCollisions *allcollisions, struct AllSpecies *allspecies);
void ProcessingPhaseDiagnosticInput(struct Control *control, struct Output *output, struct AllSpecies *allspecies);
double GetInputValueDouble(char ***input, int start, int end, char const *line, double def, char const *fail_message, int fail_response);
ptype GetInputValuePtype(char ***input, int start, int end, char const *line, ptype def, char const *fail_message, int fail_response);
int GetInputValueInt(char ***input, int start, int end, char const *line, int def, char const *fail_message, int fail_response);
long long int GetInputValueLongInt(char ***input, int start, int end, char const *line, long long int def, char const *fail_message, int fail_response);
char* GetInputValueStr(char ***input, int start, int end, char const *line, char const *fail_message, int fail_response);


//Preparation functions
void PrepareOutputDirectories(struct Output *output, struct Control *control, struct Grid *grid, struct AllSpecies *allspecies, struct Timing *timer);
void PrintInitializeInformation(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *allcreation, struct AllCollisions *allcollisions, struct AllWallCollisions *allwallcollisions, struct Grid *grid, struct Timing *timer, struct Output *output);

//Timing functions
void InitializeTimer(struct Timing *timer);
void PrintStep(struct Timing *timer, int k);
void CommunicateTiming(struct Timing *timer);
void PrintTiming(struct Timing *timer);
void PrintTimingToFile(struct Timing *timer, char *fname);


//Empirical magnetic field functions
void GetMagneticFieldGridData(struct Grid *grid);
void GetMagneticFieldRegionData(struct Grid *grid);



//Species initialization functions
void AllocateAllSpecies(struct AllSpecies *allspecies);
void InitializeAllSpecies(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *creation, struct Timing *timer);
void InitializeSpecies(struct Species *species, struct AllCreation *creation, struct Timing *timer);
void AllocateParticleBuffer(struct Species_buffer *buff, int mem);


//Creation initialization functions
void AllocateAllCreation(struct AllCreation *creation);
void InitializeAllCreation(struct Control *control, struct AllSpecies *allspecies, struct AllCreation *creation, struct Timing *timer);
void InitializePlasma(struct Species *species, struct CreatePlasma *plasma, struct Timing *timer);


//Collision initialization functions
void AllocateAllCollisions(struct AllCollisions *allcollisions);
void AllocateAllWallCollisions(struct AllWallCollisions *allwallcollisions);
void InitializeAllCollisions(struct AllSpecies *allspecies, struct AllCollisions *creation);
void InitializeElasticCollision(struct Species *species, struct ElasticCollision *elastic);
void InitializeExcitationCollision(struct Species *species, struct ExcitationCollision *excite);
void InitializeIonizationCollision(struct AllSpecies *allspecies, struct IonizationCollision *ionize);
void InitializeChargeExchangeCollision(struct Species *species, struct ChargeExchangeCollision *cxchange);
void MapCollisionsToSpecies(struct AllSpecies *allspecies, struct AllCollisions *allcollisions);
void ComputeMaxProbabilities(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, double dt);
int GetXSectionFileNumLines(char *fname);
void GetXSectionFileData(char *fname, double *E, double *XS);



//Source of species at each time step
void AllSourceCreation(struct AllSpecies *allspecies, struct AllCreation *creation, struct Region *region, struct Timing *timer);
void SourceCreation(struct Species *species, struct CreateSource *source, struct Region *region);


//Memory increase functions
void IncreaseSpeciesMemory(struct Species *species, int new_mem_size);
double* ReallocDouble(double *x, int mold, int mnew);
ptype* ReallocPtype(ptype *x, int mold, int mnew);
int* ReallocInt(int *x, int mold, int mnew);
short int* ReallocShortInt(short int *x, int mold, int mnew);
unsigned short int* ReallocUnsignedShortInt(unsigned short int *x, int mold, int mnew);


//Memory copy functions from Host to Device and visa-versa
void MemoryToDeviceAllParticles(struct AllSpecies *allspecies);
void MemoryDeleteDeviceAllParticles(struct AllSpecies *allspecies);
void UpdateHostAllParticles(struct AllSpecies *allspecies, struct Timing *timer);

void MemoryToDeviceAllSource(struct AllCreation *allspecies);
void MemoryDeleteDeviceAllSource(struct AllCreation *allspecies);

void MemoryToDeviceAllCollisions(struct AllCollisions *allcollisions);
void MemoryDeleteDeviceAllCollisions(struct AllCollisions *allcollisions);




//Evolution of Particles
void PushAllParticles(struct AllSpecies *allspecies, struct Timing *timer);
void PushParticles(struct Species *species, double dt);


//Handling of particle collisions
void CollideAllParticles(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, struct Timing *timer);
void CollideParticles(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, struct Region *region, int sindex, double dt);
void CollideParticlesNew(struct AllSpecies *allspecies, struct AllCollisions *allcollisions, int sindex, double dt);
void StoreNewCollisionParticlesBuffer(struct Species *species, struct Species_buffer *buff, int npbm);



//Magnetic field functions
#pragma acc routine seq
ptype BX(ptype x, ptype y);
#pragma acc routine seq
ptype BY(ptype x, ptype y);
#pragma acc routine seq
ptype BZ(ptype x, ptype y);

//Timed barrier function
void Barrier(struct Timing *timer);




//Numeric
void InitializeGlobalSeed(struct Control *control);
#pragma acc routine(RanGaussianDesprng) seq
double RanGaussianDesprng(desprng_common_t *process_data, desprng_individual_t *thread_data, unsigned long rcount, double sigma);
#pragma acc routine(HalfRanGaussianDesprng) seq
double HalfRanGaussianDesprng(desprng_common_t *process_data, desprng_individual_t *thread_data, unsigned long rcount, double sigma);

void WeightedHighestTwoFactors(int N, int nx, int ny, int *f1, int *f2);
int PowerOfTwo(int n);
void HighestTwoFactors(int N, int *f1, int *f2);
double dmin(double a, double b);
double dmax(double a, double b);
int sgn(double val);



//Functions for printing diagnostics
void PrintAllDiagnostics(struct AllSpecies *allspecies, struct Grid *grid, struct Timing *timer, struct Output *output);
void DiagnosePlasma(struct AllSpecies *allspecies, struct Grid *grid, struct Timing *timer, struct Output *output);
void PrintAllPhaseDiagnostics(struct AllSpecies *allspecies, struct Timing *timer, struct Output *output);
void PrintPhaseDiagnostic(struct PhaseDiagnostic *phase, struct Species *species, int step);
void ComputeTotals(struct AllSpecies *allspecies, struct Grid *grid, struct Output *output);
void PrintTotals(struct AllSpecies *allspecies, struct Timing *timer, struct Output *output);


//Prototype function for parallel I/O with MPI
void InitializeIO(struct Output *output, struct Grid *grid);


//Memory clearing functions
void FinalizeControl(struct Control *control);
void FinalizeOutput(struct Output *output);
void FinalizeAllCollisions(struct AllCollisions *allcollisions);
void FinalizeAllCreation(struct AllCreation *creation, struct Control *control);
void FinalizeAllSpecies(struct AllSpecies *allspecies);
void FinalizeSpecies(struct Species *species);
void FinalizeParticleBuffer(struct Species_buffer *buff);
void FinalizeGrid(struct Grid *grid);



//Other functions
int IsEmpty(const char *str);
double GetInputValueDouble(struct Control *control, int istart, int iend, char const *line, double def, char const *fail_message, int fail_response);
float GetInputValueFloat(struct Control *control, int istart, int iend, char const *line, float def, char const *fail_message, int fail_response);
int GetInputValueInt(struct Control *control, int istart, int iend, char const *line, int def, char const *fail_message, int fail_response);
long long int GetInputValueLongInt(struct Control *control, int istart, int iend, char const *line, long long int def, char const *fail_message, int fail_response);
char* GetInputValueStr(struct Control *control, int istart, int iend, char const *line, char const *fail_message, int fail_response);
