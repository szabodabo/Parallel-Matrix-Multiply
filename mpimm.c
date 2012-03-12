#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <rdtsc.h>
#include <mtrand.h>
#include <mpi.h>

#define MATRIX_SIZE 16
#define SEND 0
#define RECV 1

double CLOCK_RATE = 700000000.0; // change to 700000000.0 for Blue Gene/L
int MY_RANK;

int index_translate( int r, int c ) {
//	printf("[%d] Translating r=%d, c=%d to %d\n", MY_RANK, r, c, (r * MATRIX_SIZE) + c);
	return (r * MATRIX_SIZE) + c;
}

double matrix_mult( double *A, double *B, int cRow, int cCol ) {
	int idx;
	double sum = 0;
	for (idx = 0; idx < MATRIX_SIZE; idx++) {
		sum += A[index_translate(cRow, idx)] * B[index_translate(idx, cCol)];
	}
	return sum;
}

int main(int argc, char **argv)
{
	unsigned int matrix_size=MATRIX_SIZE;
	unsigned long rng_init_seeds[6]={0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
	unsigned long rng_init_length=6;

  int i, j;
	int myRank, commSize;
	int partitionSize, startIdx;
	int col_offset;
	int partitionsRemaining;
	int Bsource;
	//TODO: Find out what these values really represent
	double my_irecv_time = 0.0;
	double my_isend_time = 0.0;
	double my_compute_time = 0.0;
	int irecv_bytes = 0;
	int isend_bytes = 0;
	unsigned long long startSend, startRecv, startComp, finishSend, finishRecv, finishComp;

	double *A_ROWS;
	double *B_COLS;
	double *C_ROWS;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MY_RANK = myRank;
	Bsource = myRank;

	//TODO: Rewrite this to use tabs, not spaces
	rng_init_seeds[0] = myRank;
  init_by_array(rng_init_seeds, rng_init_length);	

	partitionSize = matrix_size / commSize;
	startIdx = partitionSize * myRank;
	col_offset = partitionSize * myRank;
	partitionsRemaining = commSize;

	MPI_Request *requests = malloc(2 * sizeof(MPI_Request));

	double *recvBuffer = malloc( matrix_size * partitionSize * sizeof(double) );

	MPI_Datatype MPI_Partition;
	MPI_Type_contiguous(partitionSize * matrix_size, MPI_DOUBLE, &MPI_Partition);
	MPI_Type_commit(&MPI_Partition);

	A_ROWS = (double *) malloc( matrix_size * partitionSize * sizeof(double) );
	B_COLS = (double *) malloc( partitionSize * matrix_size * sizeof(double) );
	C_ROWS = (double *) malloc( matrix_size * partitionSize * sizeof(double) );

	printf("I am MPI Rank %d and I have rows %d - %d\n", 
			myRank, startIdx, startIdx+partitionSize-1);

	for( i = 0; i < partitionSize; i++ ) {
  	for( j = 0; j < matrix_size; j++ ) {
			A_ROWS[index_translate(i, j)] = genrand_res53();
	    B_COLS[index_translate(j, i)] = genrand_res53();
	  }
	}

	//I will use my N/P full rows to compute a N/P * N/P segment of C
	//By passing around columns of B, I can compute a row of C in
	//  N operations... 
	//Because of that, this rank will only store one row of C.

	printf("[%d] Matrix A:\n", myRank);
	for (i = 0; i < partitionSize; i++) {
		for (j = 0; j < matrix_size; j++) {
			printf("%5.5lf ", A_ROWS[index_translate(i, j)]);
		}
		printf("\n");
	}

	printf("[%d] Matrix B (%d - %d):\n", myRank, col_offset, col_offset+partitionSize-1);
	for (i = 0; i < matrix_size; i++) {
		for (j = 0; j < partitionSize; j++) {
			printf("%5.2lf ", B_COLS[index_translate(i, j)]);
		}
		printf("\n");
	}


	while (partitionsRemaining > 0) {
		printf("[%2d] %d / %d partitions remaining.\n", myRank, partitionsRemaining, commSize);
		if (partitionsRemaining != commSize) { //This isn't the first go	
			//Get the next B offset and partition
			startRecv = rdtsc();
			int sourceRank = myRank -1;
			if (sourceRank == -1) { sourceRank = commSize-1; }
			MPI_Irecv(recvBuffer, 1, MPI_Partition, sourceRank, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[RECV]);
			Bsource--;
			if (Bsource == -1) { Bsource = commSize-1; }
			col_offset = partitionSize * Bsource;
			printf("[%d] Posted recv for BSrc %d\n", myRank, Bsource);

			//Post a send to get the current B partition out of here
			startSend = rdtsc();
			int destRank = (myRank+1) % commSize;
			MPI_Isend(B_COLS, 1, MPI_Partition, destRank, 1, MPI_COMM_WORLD, &requests[SEND]);
			printf("[%d] Posted send of current data\n", myRank);
			int completedIdx;
			for (i = 0; i < 2; i++) {
				printf("[%d] Beginning wait...\n", myRank);
				MPI_Waitany(2, requests, &completedIdx, MPI_STATUS_IGNORE);
				if (completedIdx == RECV) {
					finishRecv = rdtsc();
					printf("[%d] Completed recv\n", myRank);
				} else if (completedIdx == SEND) {
					finishSend = rdtsc();
					printf("[%d] Completed send\n", myRank);
				} else { printf("[%d] PROBLEM!\n", myRank); }
				requests[completedIdx] = MPI_REQUEST_NULL;
			}
			printf("[%d] Done with comm\n", myRank);
			
			irecv_bytes += (matrix_size * partitionSize * sizeof(double));
			isend_bytes += (matrix_size * partitionSize * sizeof(double));
			my_isend_time += (finishSend - startSend) / CLOCK_RATE;
			my_irecv_time += (finishRecv - startRecv) / CLOCK_RATE;

			printf("[%d] Moving buffer...\n", myRank);

			memcpy(B_COLS, recvBuffer, matrix_size * partitionSize * sizeof(double));
			printf("[%d] Done moving buffer\n", myRank);
			//Move the next B from buffer to current
		/*	for (i = 0; i < matrix_size; i++) {
				for (j = 0; j < partitionSize; j++) {
					B_COLS[index_translate(i, j)] = recvBuffer[index_translate(i, j)];
				}
			}*/
		} //End if not first go

/*		if (myRank == 0) {
		printf("[%d] Matrix B (%d - %d):\n", myRank, col_offset, col_offset+partitionSize-1);
		for (i = 0; i < matrix_size; i++) {
			for (j = 0; j < partitionSize; j++) {
				printf("%5.5lf ", B_COLS[i][j]);
			}
			printf("\n");
		}
		}*/

		//Now compute a block of C for the current B partition
		printf("[%d] Starting compute\n", myRank);
		startComp = rdtsc();
		for (i = 0; i < partitionSize; i++) {
			for (j = 0; j < partitionSize; j++) {
				C_ROWS[index_translate(i, j+col_offset)] = matrix_mult( A_ROWS, B_COLS, i, j );
				
			}
		}
		finishComp = rdtsc();
		printf("[%d] Finished compute\n", myRank);
		my_compute_time += (finishComp - startComp) / CLOCK_RATE;
		partitionsRemaining--;
	}
/*
	MPI_Cancel(&nextBuffer.offsetReq);
	MPI_Cancel(&sendStatus.offset);
	for (i = 0; i < matrix_size; i++) {
		MPI_Cancel(&nextBuffer.rowReq[i]);
		MPI_Cancel(&sendStatus.bufs[i]);
	}
*/
/*	printf("[%d] Matrix C:\n", myRank);
	for (i = 0; i < partitionSize; i++) {
		for (j = 0; j < matrix_size; j++) {
			printf("%5.5lf ", C_ROWS[i][j]);
		}
		printf("\n");
	}
*/

	printf("[%d] Beginning stat analysis\n", myRank);

	int PRINT_CORRECTNESS_CHECK = 1;

	MPI_Barrier(MPI_COMM_WORLD);
	if (myRank == 0 && PRINT_CORRECTNESS_CHECK) {
		//Output some of C for correctness check
		printf("Rank 0's part of Matrix C (first %d rows):\n", partitionSize);
		for (i = 0; i < partitionSize; i++) {
			for (j = 0; j < matrix_size; j++) {
				printf("%5.5lf ", C_ROWS[index_translate(i, j)]);
			}
			printf("\n");
		}
	}

	double my_total_time = my_isend_time + my_irecv_time + my_compute_time;

	double my_isend_bw = (double) isend_bytes / my_isend_time;
	double my_irecv_bw = (double) irecv_bytes / my_irecv_time;
	double total_exec_time = 0.0;
	double max_isend_bw, min_isend_bw, sum_isend_bw, avg_isend_bw;
	double max_irecv_bw, min_irecv_bw, sum_irecv_bw, avg_irecv_bw;

	printf("[%d] Variables initialized.\n", myRank);
	printf("[%d] MyTotalTime: %8.7lf sec\n", myRank, my_total_time);

	//TODO: Not sure if max should be used here
	MPI_Allreduce(&my_total_time, &total_exec_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	MPI_Allreduce(&my_isend_bw, &max_isend_bw, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&my_isend_bw, &min_isend_bw, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&my_isend_bw, &sum_isend_bw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(&my_irecv_bw, &max_irecv_bw, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&my_irecv_bw, &min_irecv_bw, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&my_irecv_bw, &sum_irecv_bw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	printf("[%d] Allreduces done!\n", myRank);

	if (myRank == 0) { //Yay me
		printf("Rank 0 printing stats!\n");
		avg_isend_bw = sum_isend_bw / (double) commSize;
		avg_irecv_bw = sum_irecv_bw / (double) commSize;

		printf("Total exec time: %8.7lf seconds\n", total_exec_time);
		printf("Max Isend BW:    %8.7lf MB/sec\n", max_isend_bw / 1000 / 1000);
		printf("Min Isend BW:    %8.7lf MB/sec\n", min_isend_bw / 1000 / 1000);
		printf("Avg Isend BW:    %8.7lf MB/sec\n", avg_isend_bw / 1000 / 1000);
		printf("Max Irecv BW:    %8.7lf MB/sec\n", max_irecv_bw / 1000 / 1000);
		printf("Min Irecv BW:    %8.7lf MB/sec\n", min_irecv_bw / 1000 / 1000);
		printf("Avg Irecv BW:    %8.7lf MB/sec\n", avg_irecv_bw / 1000 / 1000);
		
		FILE *dataFile = fopen("out.csv", "a");
		//num_nodes  exec_time  min_isend  min_irecv  max_isend  max_irecv  avg_isend  avg_irecv
		fprintf(dataFile, "%d %8.7lf %8.7lf %8.7lf %8.7lf %8.7lf %8.7lf %8.7lf\n",
				commSize, total_exec_time, 
				min_isend_bw / 1000000.0, min_irecv_bw / 1000000.0, 
				max_isend_bw / 1000000.0, max_irecv_bw / 1000000.0,
				avg_isend_bw / 1000000.0, avg_irecv_bw / 1000000.0);
		printf("[0] Closing data file...\n");
		fclose(dataFile);
	}
	printf("[%d] Breakpoint!\n", myRank);

	MPI_Barrier(MPI_COMM_WORLD);

	printf("[%d] Passed barrier\n", myRank); fflush(stdout);

//	printf("[%d] SendWaits: %d; RecvWaits: %d\n", myRank, sendWait, recvWait);

	printf("[%d] Freeing matrix buffers\n", myRank); fflush(stdout);
//	free(A_ROWS);
//	free(B_COLS);
//	free(C_ROWS);
	printf("[%d] Freeing receive buffer\n", myRank); fflush(stdout);
//	free(recvBuffer);
//	free(requests);
	printf("[%d] Done freeing.\n", myRank); fflush(stdout);

	int rc = MPI_Finalize();
	printf("[%d] Return code of finalize: %d\n", myRank, rc); fflush(stdout);
	return 0;
}
