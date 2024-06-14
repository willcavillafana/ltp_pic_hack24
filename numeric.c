#include "main.h"
#include <stdint.h>

//Declaration of global_seed variable
unsigned long global_seed;


/*
FUNCTION DESCRIPTION: Value for normally distributed PRNG algorithm.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
/* position of right-most step for normal distribution*/
#define PARAM_R 3.44428647676



/*
FUNCTION DESCRIPTION: Tabulated values for normally distributed PRNG algorithm.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
/* tabulated values for the heigt of the Ziggurat levels */
static const double ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};
#pragma acc declare copyin(ytab[0:128])


/*
FUNCTION DESCRIPTION: Tabulated values for normally distributed PRNG algorithm.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
/* tabulated values for 2^24 times x[i]/x[i+1],
 *  * used to accept for U*x[i+1]<=x[i] without any floating point operations */
static const unsigned long ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};
#pragma acc declare copyin(ktab[0:128])



/*
FUNCTION DESCRIPTION: Tabulated values for normally distributed PRNG algorithm.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
/* tabulated values of 2^{-24}*x[i] */
static const double wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};
#pragma acc declare copyin(wtab[0:128])




/*
FUNCTION DESCRIPTION: Initializes the PRNG seed from a "true" random number sampled from /dev/urandom.
LAST EDIT: Created
EDIT DATE: 1/20/2021
NAME OF EDITOR: Tasman Powis
*/
void InitializeGlobalSeed(struct Control *control)
{
	FILE *frand;

	//Obtaining a true random number for the seed
	frand = fopen("/dev/urandom","r");
	fread(&global_seed, sizeof(global_seed), 1, frand);



	if (control->seed_input >= 0) {
		int myid;

		MPI_Comm_rank(MPI_COMM_WORLD, &myid);

		global_seed = (unsigned long) control->seed_input + myid;
	}

	//printf("global_seed = %lu\n\n", global_seed);

}



#pragma acc routine(philox4_32) seq
inline void philox4_32(uint32_t k0, uint32_t k1, uint32_t count[4], uint32_t output[4]) {

    uint32_t R0 = count[0];
    uint32_t L0 = count[1];
    uint32_t R1 = count[2];
    uint32_t L1 = count[3];

    // Can replace __umulhi(a, b) with uint32_t(uint64_t(a)*b >> 32) if intrinsic is unavailable

    #pragma unroll
    for (int i = 0; i < 10; ++i) {
        uint32_t L1_new = R1 * 0xD2511F53;
        uint32_t R1_new = __umulhi(R0, 0xCD9E8D57) ^ k0 ^ L0;
        uint32_t L0_new = R0 * 0xCD9E8D57;
        uint32_t R0_new = __umulhi(R1, 0xD2511F53) ^ k1 ^ L1;

        k0 += 0xBB67AE85;
        k1 += 0x9E3779B9;

        R0 = R0_new;
        R1 = R1_new;
        L0 = L0_new;
        L1 = L1_new;
    }

    output[0] = R0;
    output[1] = L0;
    output[2] = R1;
    output[3] = L1;
}

#pragma acc routine(get_uniform_prn_new) seq
inline double get_uniform_prn_new(uint64_t seed, uint64_t countin) {
    uint32_t output[4];
    uint32_t count[] = {uint32_t(countin), uint32_t(countin >> 32), 0, 0};
    philox4_32(uint32_t(seed), uint32_t(seed >> 32), count, output);

    return ((uint64_t(output[1]) << 32) | output[0]) / (1.0 + ULONG_MAX);
}


/*
FUNCTION DESCRIPTION: Gaussian random number generated by DESPRNG.
LAST EDIT: Created
EDIT DATE: 1/20/2021
NAME OF EDITOR: Tasman Powis
*/
double RanGaussianDesprng(unsigned long seed, unsigned long rcount, double sigma)
//double RanGaussianDesprng(desprng_common_t *process_data, desprng_individual_t *thread_data, unsigned long rcount, double sigma)
{
	unsigned long  U, sign, i, j;//;, iprn;
	double  x, y;

	//printf("wtab[0] = %f\n",wtab[0]);

	while (1) {
		U = (unsigned long) (get_uniform_prn_new(seed , rcount)*4294967295); rcount++;
		//U = (unsigned long) (get_uniform_prn(process_data, thread_data, rcount, &iprn)*4294967295); rcount++;
		i = U & 0x0000007F;		/* 7 bit to choose the step */
		sign = U & 0x00000080;	/* 1 bit for the sign */
		j = U>>8;			/* 24 bit for the x-value */

		x = j*wtab[i];
		if (j < ktab[i])  break;

		if (i<127) {
			double  y0, y1;
			y0 = ytab[i];
			y1 = ytab[i+1];
			y = y1+(y0-y1)*get_uniform_prn_new(seed, rcount); rcount++;
			//y = y1+(y0-y1)*get_uniform_prn(process_data, thread_data, rcount, &iprn); rcount++;
		} else {
			x = PARAM_R - log(1.0-get_uniform_prn_new(seed, rcount))/PARAM_R; rcount++;
			y = exp(-PARAM_R*(x-0.5*PARAM_R))*get_uniform_prn_new(seed, rcount); rcount++;
			//x = PARAM_R - log(1.0-get_uniform_prn(process_data, thread_data, rcount, &iprn))/PARAM_R; rcount++;
			//y = exp(-PARAM_R*(x-0.5*PARAM_R))*get_uniform_prn(process_data, thread_data, rcount, &iprn); rcount++;
		}
		if (y < exp(-0.5*x*x))  break;
	}

	//printf("rcount = %d\n",rcount);

	return  sign ? sigma*x : -sigma*x;
}





/*
FUNCTION DESCRIPTION: Half Gaussian random number generated by DESPRNG.
LAST EDIT: Created
EDIT DATE: 1/20/2021
NAME OF EDITOR: Tasman Powis
*/
double HalfRanGaussianDesprng(desprng_common_t *process_data, desprng_individual_t *thread_data, unsigned long rcount, double sigma)
{
	unsigned long  U, sign, i, j, iprn;
	double  x, y;
	double out = -1.0;

	//printf("wtab[0] = %f\n",wtab[0]);

	while (out < 0.0) {
		while (1) {
			U = (unsigned long) (get_uniform_prn(process_data, thread_data, rcount, &iprn)*4294967295); rcount++;
			i = U & 0x0000007F;		/* 7 bit to choose the step */
			sign = U & 0x00000080;	/* 1 bit for the sign */
			j = U>>8;			/* 24 bit for the x-value */

			x = j*wtab[i];
			if (j < ktab[i])  break;

			if (i<127) {
				double  y0, y1;
				y0 = ytab[i];
				y1 = ytab[i+1];
				y = y1+(y0-y1)*get_uniform_prn(process_data, thread_data, rcount, &iprn); rcount++;
			} else {
				x = PARAM_R - log(1.0-get_uniform_prn(process_data, thread_data, rcount, &iprn))/PARAM_R; rcount++;
				y = exp(-PARAM_R*(x-0.5*PARAM_R))*get_uniform_prn(process_data, thread_data, rcount, &iprn); rcount++;
			}
			if (y < exp(-0.5*x*x))  break;
		}

		//printf("rcount = %d\n",rcount);

		out = sign ? sigma*x : -sigma*x;
		//printf("out = %.6e\n",out);
	}

	return  out;
}





/*
FUNCTION DESCRIPTION: Computes the factors for "N" which are both of the largest value (i.e. which are closest in value). This is used to compute the optimum processor layout for a square grid.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void HighestTwoFactors(int N, int *f1, int *f2)
{
	int i, lim;
	int a, b;


	*f1 = 1;

	lim = (int) sqrt((double) N);

	for(i = 1; i <= lim; i++) {
		if (N % i == 0) {
			a = i;
		}
	}
	b = N/a;

	*f1 = b;
	*f2 = a;

}



/*
FUNCTION DESCRIPTION: Computes the factors for "N" which are both of the largest value (i.e. which are closest in value), weighted by the factors "nx" and "ny". This is used to compute the optimum processor layout for a rectangular grid of dimensions nx x ny
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void WeightedHighestTwoFactors(int N, int nx, int ny, int *f1, int *f2)
{
	int i, lim;
	int a, b, at, bt;
	double ratio, ratiot, diff, difft;

	ratio = (double) nx / (double) ny;
	diff = (double) nx;

	for (i = 1; i <= nx; i++) {
		if (N % i == 0) {
			at = i;
			bt = N/i;
			ratiot = (double) at / (double) bt;
			difft = fabs(ratio - ratiot);
			if (difft < diff) {
				a = at;
				b = bt;
				diff = difft;
			}
		}
	}


	*f1 = a;
	*f2 = b;

}


/*
FUNCTION DESCRIPTION: Computes the closest power of two to the input integer "n".
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
int PowerOfTwo(int n)
{
	if (n == 0) {return 0;}

	while (n != 1) {
		if (n % 2 != 0) {return 0;}
		n = n/2;
	}
	return 1;
}



/*
FUNCTION DESCRIPTION: Returns the sign (or 0) of an integer.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
int sgn(double val) {
	if (val > 0) return 1;
	if (val < 0) return -1;
	return 0;
}


/*
FUNCTION DESCRIPTION: Returns the minimum value of two doubles (only used for initialization so doesn't need to be efficient).
LAST EDIT: Created
EDIT DATE: 5/21/2020
NAME OF EDITOR: Tasman Powis
*/
double dmin(double a, double b)
{
	if (a >= b) return b;
	if (a < b) return a;
	return 0;
}



/*
FUNCTION DESCRIPTION: Returns the maximum value of two doubles (only used for initialization so doesn't need to be efficient).
LAST EDIT: Created
EDIT DATE: 5/21/2020
NAME OF EDITOR: Tasman Powis
*/
double dmax(double a, double b)
{
	if (a >= b) return a;
	if (a < b) return b;
	return 0;
}
