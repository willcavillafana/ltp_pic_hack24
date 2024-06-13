#include "main.h"
/*
FUNCTION DESCRIPTION: Zeroing initial timer values. Also creating log file.
LAST EDIT: Added creating of log file
EDIT DATE: 31/10/2022
NAME OF EDITOR: Tasman Powis
*/
void InitializeTimer(struct Timing *timer)
{

	timer->start_step = 0;
    timer->time = 0.0;

	timer->timePush = 0.0;
	timer->timeBarrier = 0.0;
	timer->timeCollide = 0.0;
}




/*
FUNCTION DESCRIPTION: Rule for printing the current step. Reducing the print interval as the number of steps increases.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void PrintStep(struct Timing *timer, int k)
{

	fprintf(timer->flog,"\nStep %d, time %.2e ns", k, k*timer->dt*1.0e9);

	if (k <= 10) {
		printf("\nStep %d, time %.2e ns", k, k*timer->dt*1.0e9);
		timer->print_step = 1;
	}
	else if (k > 10 && k <= 100 && k % 10 == 0) {
		printf("\nStep %d, time %.2e ns", k, k*timer->dt*1.0e9);
		timer->print_step = 1;
	}
	else if (k > 100 && k <= 1000 && k % 100 == 0) {
		printf("\nStep %d, time %.2e ns", k, k*timer->dt*1.0e9);
		timer->print_step = 1;
	}
	else if (k > 1000 && k % 1000 == 0) {
		printf("\nStep %d, time %.2e ns", k, k*timer->dt*1.0e9);
		timer->print_step = 1;
	}
	else {
		timer->print_step = 0;
	}

	if (timer->print_interval > 0 && timer->step % timer->print_interval == 0) {
		timer->diagnostic_step = 1;
	}
	else {
		timer->diagnostic_step = 0;
	}

}




/*
FUNCTION DESCRIPTION: Communicating all final timing values for then computing and printings stats.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void CommunicateTiming(struct Timing *timer)
{


		MPI_Reduce(&timer->time,&timer->time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&timer->timePush,&timer->timePush_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&timer->timeBarrier,&timer->timeBarrier_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Reduce(&timer->timeCollide,&timer->timeCollide_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

}




/*
FUNCTION DESCRIPTION: Printing average times for each sub-routine.
LAST EDIT: Updated and clarified labeling
EDIT DATE: 4/1/2020
NAME OF EDITOR: Tasman Powis
*/
void PrintTiming(struct Timing *timer)
{

	int my_rank, nproc;


        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        MPI_Comm_size(MPI_COMM_WORLD, &nproc);


        if (my_rank == 0) {

		timer->time_avg *= 1.0/nproc;
		timer->timePush_avg *= 1.0/nproc;
		timer->timeCollide_avg *= 1.0/nproc;
		timer->timeBarrier_avg *= 1.0/nproc;


		//Printing timers
		printf("\n\n");
		printf("timePush %f s\n", timer->timePush_avg);
		printf("timeCollide %f s\n", timer->timeCollide_avg);
		printf("timeBarrier %f s\n", timer->timeBarrier_avg);
		printf("time %f s\n\n", timer->time_avg);

        }

}
