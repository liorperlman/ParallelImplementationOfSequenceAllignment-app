#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include "myProto.h"
#include "Sequence.h"
#include "myMacro.h"

/*
MPI+OpenMP+CUDA Integration project 
Initially the array is known to process 0.
It sends half of the positions to process 1.
Both processes start to make the mutant array for half of the positions each.
Half of the mutants and signs arrays generated using OpenMP and half generated using CUDA.
The best mutant details found by process 1 is sended from process 1 to process 0, which perform a test to check if the mutant is better (gets better score) than the best mutant found by process 0.
*/
 

int main(int argc, char *argv[]) {
   	MPI_Status  status;
   	int size = 0, rank, i, numOfProc, position, sizeToSend = 0, print = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

	char *fileName;
	char mod[7]; char line1[10000]; char line2[5000];
	double weights[4], otherScore = -10000, myScore = -10000;
	int repositions[5001];// max seq1 = 10000, min seq2 = 1 -> max repositions = 10000/2 = 5000+1 (flag)
	int current, numOfRepositions, flag = -1, myBestPosition, otherBestPosition, modValue;
	Sequence *seq1, *seq2;
	SignType minimum;
   	
	if (numOfProc != 2) {
		printf("Run the example with two processes only\n");
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Divide the tasks between both processes
	if (rank == 0) {
		
		// Allocate memory for the whole array and send a half of the array to other process
		if ((seq1 = (Sequence *)malloc(sizeof(Sequence))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		if ((seq2 = (Sequence *)malloc(sizeof(Sequence))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		
		// Assign weights, mod and sequences values from input.txt and update variables according to mod
		readFromFile("input.txt", weights, line1, line2, mod);
		minimum = updateMin(weights);
		if (strcmp(mod, "maximum") == 0){
			modValue = 1;
			myScore = -10000;
			otherScore = -10000;
		}	
		else{
			modValue = 0;
			myScore = 10000;
			otherScore = 10000;
		}			
		seq1->count = strlen(line1);
		seq2->count = strlen(line2);		
		if ((seq1->data = (char *)malloc(seq1->count)) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		if ((seq2->data = (char *)malloc(seq2->count)) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		strcpy(seq1->data, line1);
		strcpy(seq2->data, line2);
		
		numOfRepositions = seq1->count - seq2->count + 1;
		sizeToSend = seq2->count/2;

		// Broadcast both sequences and weights array
		MPI_Bcast(&seq1->count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&seq2->count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&line1, seq1->count, MPI_CHAR, 0, MPI_COMM_WORLD);
		MPI_Bcast(&line2, seq2->count, MPI_CHAR, 0, MPI_COMM_WORLD);
		MPI_Bcast(&weights, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&sizeToSend, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numOfRepositions, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&modValue, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// Divide work by positions
		int index = 0;
		for (position = 0; position < numOfRepositions; position++){
			if (position % 2 == 1)			
				MPI_Send(&position, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			else
				repositions[index++] = position;		
		}
		repositions[index] = flag; // master's stop flag
		MPI_Send(&flag, 1, MPI_INT, 1, 0, MPI_COMM_WORLD); // slave's stop flag
	} 
	else {
		// Allocate memory and recieve both sequences from other process
		if ((seq1 = (Sequence *)malloc(sizeof(Sequence))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		if ((seq2 = (Sequence *)malloc(sizeof(Sequence))) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		MPI_Bcast(&seq1->count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&seq2->count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if ((seq1->data = (char *)malloc(seq1->count)) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		if ((seq2->data = (char *)malloc(seq2->count)) == NULL)
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		
		MPI_Bcast(seq1->data, seq1->count, MPI_CHAR, 0, MPI_COMM_WORLD);
		MPI_Bcast(seq2->data, seq2->count, MPI_CHAR, 0, MPI_COMM_WORLD);
		MPI_Bcast(&weights, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&sizeToSend, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numOfRepositions, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&modValue, 1, MPI_INT, 0, MPI_COMM_WORLD);
		minimum = updateMin(weights);
		if (modValue == 1)	// maximum
			myScore = -10000;
			
		else 	// modValue = 0 => minimum
			myScore = 10000;
		
		MPI_Recv(&position, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		int i = 0;
		while (position > 0){
			
			repositions[i++] = position;	
			MPI_Recv(&position, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
		repositions[i] = position; // slave's stop flag
		
		if (print == 1){
			printf("in slave -> seq1->count = %d\n", seq1->count);
			printf("in slave -> seq2->count = %d\n", seq2->count);
			printf("in slave -> seq1->data = %s\n", seq1->data);
			printf("in slave -> seq2->data = %s\n", seq2->data);
			printf("in slave -> weights = w1=%lf, w2=%lf, w3=%lf, w4=%lf\n", weights[0], weights[1],
																				weights[2], weights[3]);
			printf("in slave -> sizeToSend = %d\n", sizeToSend);
			printf("in slave -> position = %d\n", position);
			printf("in slave -> modValue = %d, myScore = %lf\n", modValue, myScore);
		}
	}

	char myBestMutant[seq2->count]; char otherBestMutant[seq2->count];

	// On each process - perform a first half of its task with OpenMP   
#pragma omp parallel num_threads(numOfRepositions)
{ 
	#pragma omp for
	for (int i = 0; i < numOfRepositions; i++){
		
		if (repositions[i] != -1 && !(i!= 0 && repositions[i] == 0)){
			
			
			SignType signs[seq2->count];
			char mutant[seq2->count];
			// On each thread - perform half of its task with OpenMP (creating half of the mutant and updated signs array)
			createMutant(signs, seq1->data + repositions[i], seq2->data, sizeToSend, mutant, modValue, minimum);
			createSigns(signs, seq1->data + repositions[i], mutant, sizeToSend);
			// On each thread - perform a second half of its task with CUDA
			if (computeOnGPU(signs + sizeToSend, seq1->data + repositions[i] + sizeToSend, seq2->data + sizeToSend,
								seq2->count - sizeToSend, mutant + sizeToSend, modValue, minimum) != 0)
		  		MPI_Abort(MPI_COMM_WORLD, __LINE__);
			
			// Calculating the new mutant score
			double newScore = getCount(signs, seq2->count, weights);
			
	       	//if we found a better score, update values
	       	if ((newScore > myScore && modValue == 1) || (newScore < myScore && modValue == 0)){
	       		myScore = newScore;
	       		myBestPosition = repositions[i];
	       		strcpy(myBestMutant, mutant);	
	       	}
	       	// Print for debugging
	       	if (print == 1){
	       		printf("repositions[%d] = %d with thread #%d\n", i, repositions[i], omp_get_thread_num());
				printf("mutant for position %d is: %s with the score of %lf\n", repositions[i], mutant, newScore);
			}
		}
	}
}
	// Waiting for all to finish
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank!=0){	//sending slave's results
		MPI_Send(&myScore, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&myBestPosition, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&myBestMutant, seq2->count, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	// Collect the result on one of processes
	if (rank == 0) { //receiving slave's results
		MPI_Recv(&otherScore, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&otherBestPosition, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&otherBestMutant, seq2->count, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
		
		// Final compare between the best mutant each process made
		if ((otherScore > myScore && modValue == 1) || (otherScore < myScore && modValue == 0)){
			myScore = otherScore;
      		myBestPosition = otherBestPosition;
      		strcpy(myBestMutant, otherBestMutant);
      		
		}
		// Print for debugging
		if (print == 1){
			printf("bestmutant recieved from slave to rank %d best position is %d with a score of %lf\n", rank, otherBestPosition, otherScore);
			printf("the overall best mutant is:\n%s\nwith the position of %d and a score of %lf\n", myBestMutant, myBestPosition, myScore);
		}
		// Writing results to output file
		writeOutputToFile(myBestMutant, myBestPosition, myScore);
		     	
	}

	freeAll(seq1, seq2);
	MPI_Finalize();
	
	return 0;
}
