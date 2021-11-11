#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"
#include "myMacro.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


__device__ SignType defineSignCuda(char c1, char c2);
__device__ char *myStrchr(const char *s, int c);
__device__ char getCharReplacementCuda(char c1, char c2, int modValue, SignType minimum);


__global__ void createMutantCuda(SignType *d_signs, char *d_seq1, char *d_seq2, int numElements, char *d_mutant, int modValue, SignType minimum){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < numElements){
		//if (d_seq2[i] == letterToReplace){
		if (modValue == 1)
			d_mutant[i] = getCharReplacementCuda(d_seq1[i], d_seq2[i], modValue, minimum);
		else{ //minimum
			if (defineSignCuda(d_seq1[i], d_seq2[i]) != minimum)	
				d_mutant[i] = getCharReplacementCuda(d_seq1[i], d_seq2[i], modValue, minimum);
			else
				d_mutant[i] = d_seq2[i];
		}
		
		d_signs[i] = defineSignCuda(d_seq1[i], d_mutant[i]);
	}		
}

__device__ char getCharReplacementCuda(char c1, char c2, int modValue, SignType minimum){
	int i;
	char *ptr1, *ptr2;
	SignType sign, sign1, sign2;
	const char *similar[SIM_SIZE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	int similarSize = SIM_SIZE;
	
	if (modValue == 1) {	// maximum
		// Check if both c1 and c2 are in one of Similar groups
		for (i = 0; i < similarSize; i++) {
			ptr1 = myStrchr(similar[i], c1);
			ptr2 = myStrchr(similar[i], c2);
			if (ptr1 != NULL && ptr2 != NULL)
				return c2;
		}
		// If not we can replace it with the char from seq1 to increase the score
		return c1;
	}
	else {	// modValue = 0 => minimum
		// Find a char who's giving the worst score in comparison with c1 and not similar in comparison with c2 
		sign = defineSignCuda(c1, c2);
		if (sign != minimum){
			for (i = 0; i < 26; i++){
				sign1 = defineSignCuda(c1, i+65);
				sign2 = defineSignCuda(c2, i+65);
				if (sign1 == minimum && sign2 != Similar)
					return i+65;	
			}
			return c2;
		}
		else{ // sign = minimum
			return c2;
		}
	}
}

__device__ SignType defineSignCuda(char c1, char c2) {
	int i;
	char *ptr1, *ptr2;
	const char *similar[SIM_SIZE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	int similarSize = SIM_SIZE;
	const char *almostSimilar[ALMOST_SIM_SIZE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK",
							 			"NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };
	int almostSimSize = ALMOST_SIM_SIZE;
	
	// Check if Equal
	if (c1 == c2)
		return Equal;

	// Check if both c1 and c2 are in one of Similar groups
	for (i = 0; i < similarSize; i++) {
		ptr1 = myStrchr(similar[i], c1);
		ptr2 = myStrchr(similar[i], c2);
		if (ptr1 != NULL && ptr2 != NULL)
			return Similar;
	}

	// Check if both c1 and c2 are in one of Almost Similar groups
	for (i = 0; i < almostSimSize; i++) {
		ptr1 = myStrchr(almostSimilar[i], c1);
		ptr2 = myStrchr(almostSimilar[i], c2);
		if (ptr1 != NULL && ptr2 != NULL)
			return AlmostSimilar;
	}
	// Not Equal and Not found in Similar or AlmostSimilar groups 
	return NotEqual;
}

__device__ char *myStrchr(const char *s, int c)
{
    while (*s != (char)c)
        if (!*s++)
            return NULL;
    return (char *)s;
}

int computeOnGPU(SignType *h_signs, char *h_seq1, char *h_seq2, int numElements, char *h_mutant, int modValue, SignType minimum) {
    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;

    size_t charSize = numElements * sizeof(char);
    size_t signSize = numElements * sizeof(SignType);

    // Allocate memory on GPU to copy the data from the host
    char *d_seq1;
    char *d_seq2;
    char *d_mutant;
    SignType *d_signs;
    err = cudaMalloc((void **)&d_seq1, charSize);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory for seq1 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    } 
    err = cudaMalloc((void **)&d_seq2, charSize);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory for seq2 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&d_mutant, charSize);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory for seq2 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&d_signs, signSize);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory for signs - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_seq1, h_seq1, charSize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device for seq1 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(d_seq2, h_seq2, charSize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device for seq2 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Launch the Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
    createMutantCuda<<<blocksPerGrid, threadsPerBlock>>>(d_signs, d_seq1, d_seq2, numElements, d_mutant, modValue, minimum);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the result from GPU to the host memory.
    err = cudaMemcpy(h_mutant, d_mutant, charSize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(h_signs, d_signs, signSize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free allocated memory on GPU
    if (cudaFree(d_seq1) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    if (cudaFree(d_seq2) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    if (cudaFree(d_mutant) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
	if (cudaFree(d_signs) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    return 0;
}
