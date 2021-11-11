#pragma once
#include "myMacro.h"
#include "Sequence.h"

int computeOnGPU(SignType *h_signs, char *h_seq1, char *h_seq2, int numElements, char *mutant, int modValue, SignType minimum);
void readFromFile(char const *fileName, double *weights, char *line1, char *line2, char *mod);
double cpu_Max(double *tmpArr, int iterationSize);
void calcSumOfVecB(double *h_B, int numElements);
SignType defineSign(char c1, char c2);
void createSigns(SignType *signs, char *seq1, char *seq2, int size);
double getCount(SignType *signs, int size, double *weights);
char getCharReplacement(char c1, char c2, int modValue, SignType minimum);
char getCharToReplace(SignType *signs, char *seq2, int size, int modValue, double *weights);
void createMutant(SignType *signs, char *seq1, char *seq2, int size, char *mutant, int modValue, SignType minimum);
void writeOutputToFile(char *mutant, int position, double score);
SignType updateMin(double *weights);
void freeAll(Sequence *seq1,Sequence *seq2);
