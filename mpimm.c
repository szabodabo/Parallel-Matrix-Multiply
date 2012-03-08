#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <rdtsc.h>
#include <mtrand.h>
#include <mpi.h>

double CLOCK_RATE = 700000000.0; // change to 700000000.0 for Blue Gene/L

double matrix_mult( double **A, double **B, int size, int cRow, int cCol ) {
	int idx;
	double sum = 0;
	for (idx = 0; idx < size; idx++) {
		sum += A[cRow][idx] * B[idx][cCol];
	}
	return sum;
}

int main(int argc, char **argv)
{
	unsigned int matrix_size=32;
	unsigned long rng_init_seeds[6]={0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
	unsigned long rng_init_length=6;

  int i, j;
	int myRank, commSize;
	int partitionSize, startIdx;
	int col_offset;
	int partitionsRemaining;
	int OFFSET_TAG = matrix_size+1;
	//TODO: Find out what these values really represent
	double my_irecv_time = 0.0;
	double my_isend_time = 0.0;
	double my_compute_time = 0.0;
	int irecv_bytes = 0;
	int isend_bytes = 0;
	unsigned long long startSend, startRecv, startComp, finishSend, finishRecv, finishComp;
	//int recvWait = 0;
	//int sendWait = 0;

	double **A_ROWS;
	double **B_COLS;
	double **C_ROWS;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	//TODO: Rewrite this to use tabs, not spaces
	rng_init_seeds[0] = myRank;
  init_by_array(rng_init_seeds, rng_init_length);	
	
	partitionSize = matrix_size / commSize;
	startIdx = partitionSize * myRank;
	col_offset = partitionSize * myRank;
	partitionsRemaining = commSize;

	typedef struct {
		MPI_Request offset;
		MPI_Request *bufs; //matrix_size
	} send_req_struct;
	send_req_struct sendStatus;
	sendStatus.bufs = calloc(matrix_size, sizeof(MPI_Request));

	typedef struct {
		int offset;
		double **rows;
		MPI_Request offsetReq;
		MPI_Request *rowReq;
	} buffer_struct;
	buffer_struct nextBuffer;
	nextBuffer.rowReq = calloc(matrix_size, sizeof(MPI_Request));
	nextBuffer.rows = calloc(matrix_size, sizeof(double *));
	for (i = 0; i < matrix_size; i++) {
		nextBuffer.rows[i] = calloc(partitionSize, sizeof(double));
	}
/*
	MPI_Datatype MPI_BRow;
	MPI_Type_contiguous(partitionSize, MPI_DOUBLE, &MPI_BRow);
	MPI_Type_commit(&MPI_BRow);
*/

	A_ROWS = (double **) calloc( partitionSize, sizeof(double *) );
	for (i = 0; i < partitionSize; i++) {
		A_ROWS[i] = (double *) calloc( matrix_size, sizeof(double *) );
	}

	B_COLS = (double **) calloc( matrix_size, sizeof(double *) );
	for (i = 0; i < matrix_size; i++) {
		B_COLS[i] = (double *) calloc( partitionSize, sizeof(double *) );
	}

	C_ROWS = (double **) calloc( partitionSize, sizeof(double *) );
	for (i = 0; i < partitionSize; i++) {
		C_ROWS[i] = (double *) calloc( matrix_size, sizeof(double *) );
	}

	printf("I am MPI Rank %d and I have rows %d - %d\n", 
			myRank, startIdx, startIdx+partitionSize-1);

	for( i = 0; i < partitionSize; i++ ) {
  	for( j = 0; j < matrix_size; j++ ) {
    	A_ROWS[i][j] = genrand_res53(); // (double)i*j;
	    B_COLS[j][i] = genrand_res53(); //A[i][j] * A[i][j];
	  }
	}

	//I will use my N/P full rows to compute a N/P * N/P segment of C
	//By passing around columns of B, I can compute a row of C in
	//  N operations... 
	//Because of that, this rank will only store one row of C.

/*	printf("[%d] Matrix A:\n", myRank);
	for (i = 0; i < partitionSize; i++) {
		for (j = 0; j < matrix_size; j++) {
			printf("%5.5lf ", A_ROWS[i][j]);
		}
		printf("\n");
	}*/
/*
	printf("[%d] Matrix B (%d - %d):\n", myRank, col_offset, col_offset+partitionSize-1);
	for (i = 0; i < matrix_size; i++) {
		for (j = 0; j < partitionSize; j++) {
			printf("%5.2lf ", B_COLS[i][j]);
		}
		printf("\n");
	}
*/

	while (partitionsRemaining > 0) {
		if (partitionsRemaining != commSize) { //This isn't the first go	

			//Get the next B offset and partition
			startRecv = rdtsc();
			int sourceRank = myRank -1;
			if (sourceRank == -1) { sourceRank = commSize-1; }
			MPI_Irecv(&nextBuffer.offset, 1, MPI_INT, sourceRank, OFFSET_TAG, MPI_COMM_WORLD, &nextBuffer.offsetReq);
			for (i = 0; i < matrix_size; i++) {
				MPI_Irecv(&nextBuffer.rows[i][0], partitionSize, MPI_DOUBLE, sourceRank, i, MPI_COMM_WORLD, &nextBuffer.rowReq[i]);
			}


			//Post a send to get the current B partition out of here
			startSend = rdtsc();
			int destRank = (myRank+1) % commSize;
			MPI_Isend(&col_offset, 1, MPI_INT, destRank, OFFSET_TAG, MPI_COMM_WORLD, &sendStatus.offset);
			for (i = 0; i < matrix_size; i++) {
				MPI_Isend(&B_COLS[i][0], partitionSize, MPI_DOUBLE, destRank, i, MPI_COMM_WORLD, &sendStatus.bufs[i]);
			}

			int sentOffset = 0;
			int recvOffset = 0;
			int sentData = 0;
			int recvData = 0;
			int recvJustCompleted = 0;
			int sendJustCompleted = 0;

			while (!sentOffset || !recvOffset || !sentData || !recvData) {
				if (!sentOffset) {
					MPI_Test(&sendStatus.offset, &sentOffset, MPI_STATUS_IGNORE);
					if (sentOffset) {
						sendJustCompleted = 1;
					}
				}
				if (!sentData) {
					MPI_Testall(matrix_size, sendStatus.bufs, &sentData, MPI_STATUSES_IGNORE);
					if (sentData) {
						sendJustCompleted = 1;
					} else {
		//				printf("[%d] Waiting for send...\n", myRank);
//						sendWait++;
					}
				}
				if (sendJustCompleted && sentOffset && sentData) {
					finishSend = rdtsc();
					isend_bytes += sizeof(int);
					isend_bytes += (matrix_size * partitionSize * sizeof(double));
					my_isend_time += (finishSend - startSend) / CLOCK_RATE;
				}
				if (!recvOffset) {
					MPI_Test(&nextBuffer.offsetReq, &recvOffset, MPI_STATUS_IGNORE);
					if (recvOffset) {
						recvJustCompleted = 1;
					}
				}
				if (!recvData) {
					MPI_Testall(matrix_size, nextBuffer.rowReq, &recvData, MPI_STATUSES_IGNORE);
					if (recvData) {
						recvJustCompleted = 1;
					}
					else {
	//					recvWait++;
			//			printf("[%d] Waiting for recv...\n", myRank);
					}
				}
				if (recvJustCompleted && recvOffset && recvData) {
					finishRecv = rdtsc();
					irecv_bytes += sizeof(int);
					irecv_bytes += (matrix_size * partitionSize * sizeof(double));
					my_irecv_time += (finishRecv - startRecv) / CLOCK_RATE;
				}
				sendJustCompleted = 0;
				recvJustCompleted = 0;
			}

			//TODO: Time these waits
			//Make sure we have the next B info before copying it from buffer
//			MPI_Wait(&nextBuffer.offsetReq, MPI_STATUS_IGNORE);
//			MPI_Waitall(matrix_size, nextBuffer.rowReq, MPI_STATUSES_IGNORE);

			//Make sure our send has completed before we remove its data
//			MPI_Wait(&sendStatus.offset, MPI_STATUS_IGNORE);
//			MPI_Waitall(matrix_size, sendStatus.bufs, MPI_STATUSES_IGNORE);



			//Move the next B from buffer to current
			col_offset = nextBuffer.offset;
			for (i = 0; i < matrix_size; i++) {
				for (j = 0; j < partitionSize; j++) {
					B_COLS[i][j] = nextBuffer.rows[i][j];
				}
			}
		}
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
		startComp = rdtsc();
		for (i = 0; i < partitionSize; i++) {
			printf("[%d] Computing row %d of Matrix C...\n", myRank, startIdx+i);
			for (j = 0; j < partitionSize; j++) {
				printf("[%d] Computing C[%d][%d]...\n", myRank, startIdx+i, j+col_offset);
				C_ROWS[i][j+col_offset] = matrix_mult( A_ROWS, B_COLS, matrix_size, i, j );
				
			}
			printf("[%d] Done computing row %d of Matrix C.\n", myRank, startIdx+i);
		}
		finishComp = rdtsc();
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

	int PRINT_CORRECTNESS_CHECK = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	if (myRank == 0 && PRINT_CORRECTNESS_CHECK) {
		//Output some of C for correctness check
		printf("Rank 0's part of Matrix C (first %d rows):\n", partitionSize);
		for (i = 0; i < partitionSize; i++) {
			for (j = 0; j < matrix_size; j++) {
				printf("%5.5lf ", C_ROWS[i][j]);
			}
			printf("\n");
		}
	}

	double my_total_time = my_isend_time + my_irecv_time + my_compute_time;

	double my_isend_bw = (double) isend_bytes / my_isend_time;
	double my_irecv_bw = (double) irecv_bytes / my_irecv_time;
	double total_exec_time;
	double max_isend_bw, min_isend_bw, sum_isend_bw, avg_isend_bw;
	double max_irecv_bw, min_irecv_bw, sum_irecv_bw, avg_irecv_bw;

	//TODO: Not sure if max should be used here
	MPI_Reduce(&my_total_time, &total_exec_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	MPI_Reduce(&my_isend_bw, &max_isend_bw, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&my_isend_bw, &min_isend_bw, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&my_isend_bw, &sum_isend_bw, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&my_irecv_bw, &max_irecv_bw, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&my_irecv_bw, &min_irecv_bw, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&my_irecv_bw, &sum_irecv_bw, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myRank == 0) { //Yay me
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
		fclose(dataFile);
	}

//	MPI_Barrier(MPI_COMM_WORLD);

//	printf("[%d] SendWaits: %d; RecvWaits: %d\n", myRank, sendWait, recvWait);

	MPI_Finalize();
	return 0;
}
