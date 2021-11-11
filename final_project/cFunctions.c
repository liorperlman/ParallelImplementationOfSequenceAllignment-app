#include "myProto.h"
#include <stdio.h>
#include <stdlib.h>
#include "myMacro.h"
#include <string.h>


char *similar[SIM_SIZE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
int similarSize = SIM_SIZE;
char *almostSimilar[ALMOST_SIM_SIZE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK",
							 			"NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };
int almostSimSize = ALMOST_SIM_SIZE;
	
void readFromFile(char const *fileName, double *weights, char *line1, char *line2, char *mod){
	
	FILE *fp = fopen(fileName,"r");
	if(fp == NULL){
		fprintf(stderr,"faild to open input.txt\n");
		exit(EXIT_FAILURE);
	}
	fscanf(fp,"%lf %lf %lf %lf", &weights[0], &weights[1], &weights[2], &weights[3]);
	fscanf(fp, "%s %s %s", line1, line2, mod);

	fclose(fp);
}

SignType defineSign(char c1, char c2) {
	int i;
	char *ptr1, *ptr2;
	
	// Check if Equal
	if (c1 == c2)
		return Equal;

	// Check if both c1 and c2 are in one of Similar groups
	for (i = 0; i < similarSize; i++) {
		ptr1 = strchr(similar[i], c1);
		ptr2 = strchr(similar[i], c2);
		if (ptr1 != NULL && ptr2 != NULL)
			return Similar;
	}

	// Check if both c1 and c2 are in one of Almost Similar groups
	for (i = 0; i < almostSimSize; i++) {
		ptr1 = strchr(almostSimilar[i], c1);
		ptr2 = strchr(almostSimilar[i], c2);
		if (ptr1 != NULL && ptr2 != NULL)
			return AlmostSimilar;
	}
	// Not Equal and Not found in Similar or AlmostSimilar groups 
	return NotEqual;
}
void createSigns(SignType *signs, char *seq1, char *seq2, int size){

	// Assign SignType values to signs array 
	for (int j = 0; j< size; j++){ 
  		signs[j] = defineSign(seq1[j], seq2[j]);
  	}
}

double getCount(SignType *signs, int size, double *weights) {
	int i, numOfStars = 0, numOfColons = 0, numOfDots = 0, numOfSpaces = 0;

	for (i = 0; i < size; i++) {
		// Analyze signs according to the criteria
			if (signs[i] == Equal)
				numOfStars++;
			else if (signs[i] == Similar)
				numOfColons++;
			else if (signs[i] == AlmostSimilar)
				numOfDots++;
			else
				numOfSpaces++;
	}
	return numOfStars*weights[0] - numOfColons*weights[1] - numOfDots*weights[2] - numOfSpaces*weights[3];
}

char getCharReplacement(char c1, char c2, int modValue, SignType minimum){
	int i;
	char *ptr1, *ptr2;
	SignType sign, sign1, sign2;
	if (modValue == 1) {	// maximum
		// Check if both c1 and c2 are in one of Similar groups
		for (i = 0; i < similarSize; i++) {
			ptr1 = strchr(similar[i], c1);
			ptr2 = strchr(similar[i], c2);
			if (ptr1 != NULL && ptr2 != NULL)
				return c2;
		}
		// If not we can replace it with the char from seq1 to increase the score
		return c1;
	}
	else {	// modValue = 0 => minimum
		// Find a char who's giving the worst score in comparison with c1 and not Conservative in comparison with c2
		sign = defineSign(c1, c2);
		if (sign != minimum){
			for (i = 0; i < 26; i++){
				sign1 = defineSign(c1, i+'A');
				sign2 = defineSign(c2, i+'A');
				if (sign1 == minimum && sign2 != Similar)
					return i+'A';	
			}
			return c2;
		}
		else {//sign = minimum
			return c2;
		}
	}
}

// For single letter replacement version
char getCharToReplace(SignType *signs, char *seq2, int size, int mod, double *weights){
	double letters[27] = {0};
	double maxScore = 0;
	char mostCommonEqualOrNotEqualLetter;
		
	if (mod == 1) {	// maximum
		for (int i = 0; i < size; i++){
			// Calculating only NotEqual and AlmostSimilar since Similar can't be replaced
			if (seq2[i] != '-'){
				if (signs[i] == NotEqual)
					letters[seq2[i]-65] += (weights[3]+weights[0]);
					
				else if (signs[i] == AlmostSimilar)
					letters[seq2[i]-65] += (weights[2]+weights[0]);
				
				if (letters[seq2[i]-65] > maxScore){
						maxScore = letters[seq2[i]-65];
						mostCommonEqualOrNotEqualLetter = seq2[i];
				}
				
			}
			else { // seq[i] = '-'
				if (signs[i] == NotEqual)
					letters[26] += (weights[3]+weights[0]);
				else if (signs[i] == AlmostSimilar)
					letters[26] += (weights[2]+weights[0]);
				if (letters[26] > maxScore){
					maxScore = letters[26];
					mostCommonEqualOrNotEqualLetter = seq2[i];
				}
			}
		}  	
	}
	else {	// mod = 0 => minimum
		for (int i = 0; i < size; i++){
			if (seq2[i] != '-'){
				if (signs[i] == Equal)
					letters[seq2[i]-65]+=(weights[3]+weights[0]);
				else if (signs[i] == Similar)
					letters[seq2[i]-65]+=(weights[3]-weights[1]);
				else if (signs[i] == AlmostSimilar)
					letters[seq2[i]-65]+=(weights[3]-weights[2]);
				
				if (letters[seq2[i]-65] > maxScore){
					maxScore = letters[seq2[i]-65];
					mostCommonEqualOrNotEqualLetter = seq2[i];
				}
			}
			else { // seq[i] = '-'
				if (signs[i] == Equal)
					letters[26]+=(weights[3]+weights[0]);
				else if (signs[i] == Similar)
					letters[26]+=(weights[3]-weights[1]);
				else if (signs[i] == AlmostSimilar)
					letters[26]+=(weights[3]-weights[2]);
					
				if (letters[26] > maxScore){
					maxScore = letters[26];
					mostCommonEqualOrNotEqualLetter = seq2[i];
				}
			}
		}  	
	}
	
	return mostCommonEqualOrNotEqualLetter;
}

void createMutant(SignType *signs, char *seq1, char *seq2, int size, char *mutant, int modValue, SignType minimum){
	
	for(int i = 0; i < size; i++)
		mutant[i] = getCharReplacement(seq1[i], seq2[i], modValue, minimum);
	
}

void writeOutputToFile(char *mutant, int position, double score){
	FILE *out = fopen("Output.txt", "w");
	if (out == NULL){
		printf("Error opening output file!\n");
		exit(1);
	}
	fprintf(out, "Mutant: %s\n", mutant);
	fprintf(out, "Offset: %d, Alignment Score: %lf\n", position, score);
	fclose(out);
}

// This function finds the worst weights for the minimum mod
SignType updateMin(double *weights){
	int min = 1;
	for (int i = 2; i < 4; i++){
		if (weights[i] >= weights[i-1])
			min = i;
	}
	if (min == 1)
		return Similar;
	else if (min == 2)
		return AlmostSimilar;
	else	// min == 3
		return NotEqual;
}

// Free allocated memory
void freeAll(Sequence *seq1,Sequence *seq2){
	free(seq1->data);
	free(seq2->data);
	free(seq1);
	free(seq2);
}
