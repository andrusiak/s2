using namespace std;
/*
 ============================================================================
 Name        : s2calc.c
 Author      : Yurii Andrusyak <yuriy.andrusyak@gmail.com>
 Version     : 1.0
 Copyright   : Your copyright notice
 Description : Enhanced algorithm of generalized order parameter S2
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

#define NMAX 10000
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct AutocorrFunction{
	float* t;
	float* C1;
	long frames_number;
};

struct SimpleModelState{
	float s2;
	float tau_s;
};

struct ExtendedModelState{
	float s2;
	float sf2;
	float tau_f;
	float tau_s;
};

void printACF(float * t, float * C1, int n){
    for(int i = 0; i<n; i++){
     printf("time = %5.1f , C1 = %f\n", *(t+i),*(C1+i));
    }
}

SimpleModelState doIteration(AutocorrFunction* acf, SimpleModelState initialState) {
	/**
	 * Use Newton's method for optimization
	 *
	 * X[t+1]=X[t] - H^-1 * g
	 * 		where H(s2,ts) - hessian of f(x,y)= (1 - s2) Exp[-ti/ts] + s2 in point(s2, ts)
	 */

	float s2 = initialState.s2, tau_s = initialState.tau_s;

	Matrix hessian(2, 2, "H");
	Matrix gradient(2, 1, "grad");

	for (long i = 0; i < acf->frames_number; i++) {

	float ti=*(acf)->t, Ci = *(acf)->C1;
	hessian(0,0) += 0;
	hessian(0,1) += -(ti*exp(-ti/tau_s))/(tau_s*tau_s);
	hessian(1,0) += -(ti* exp(-ti/tau_s))/(tau_s*tau_s);
	hessian(1,1) += -(ti*(s2-1)*(ti-2*tau_s)*exp(-ti/tau_s))/(tau_s*tau_s*tau_s*tau_s);

	gradient(0,0) += 2.0f*(exp(-ti/tau_s)-1)*(Ci + exp(-ti/tau_s)*(s2 - 1) - s2); // 1 - exp(-ti / tau_s);
	gradient(1,0) += 2.0f*exp(-2*ti/tau_s)*(s2 - 1)*ti*(-1 + exp(ti/tau_s)*(Ci - s2) + s2)/(tau_s*tau_s); //-(ti*(s2-1)*exp(-ti/tau_s))/(tau_s*tau_s);
	}
	return initialState;
}

SimpleModelState calcSimpleApproad(AutocorrFunction* acf, SimpleModelState* initialValue) {
	float eps = 0.0001;
//	float s2i = 0.5, s2next = s2i + 10 * eps;
//	float tau_s = 50.0;
	SimpleModelState currentValue = *initialValue;
	do {
		*initialValue = currentValue;

		currentValue = doIteration(acf, *initialValue);
	} while (fabs(currentValue.s2 - initialValue->s2) > eps);

	return currentValue;
}

void calcExtendedApproad(float* t, float* C1, float& s2, float& sf2, float& tau_f, float& tau_s){

}

float processACF(float * t, float * C1, int n){
    int t_e, it_e= -1;
    float A0;
    int kpr = 4;
    for(int i = 0; i<n; i++){
        if(C1[i]<1/exp(1)){
            t_e =(int)t[i]; it_e=i;
   //         printf("t_e =%d ", t_e);
            break;
        }
    }
    if(it_e== -1){ it_e = n-1; t_e=(int)t[it_e];}

    int itmax = MIN(kpr * it_e, n-1);
    A0 = C1[itmax];
    float s2 = A0 + (1-A0)*exp(-t[itmax]/t_e);
   // printf("A0=C1[%d] =%f \n",itmax, A0);
//    printf("s2=%f \n", s2);

    return s2;
}

int main(void) {
    static const char acfFile[] = "rotacf1.xvg";
    static const char s2File[] = "s2kan.xvg";
    char * header ="@    title \"Rotational Correlation Function\"\n@    xaxis  label \"Time (ps)\" \n@    yaxis  label \"C(t)\" @TYPE xy\n";

    FILE *file = fopen( acfFile, "r");
    if (file != NULL) {
        printf("File %s has been successfully opened\n", acfFile);
        char line[100];
        float time_, C1_;
        int numValuesRead;
        float t [NMAX];
        float C1[NMAX];
        float s2;
        long i = 0, nFrames, nResidue = 0;
        FILE *output = fopen(s2File, "w");
        fprintf(output,"%s", header);

        while(fgets(line, sizeof line, file)!= NULL) /* read a line */
        {
            if(strchr(line, '#')== NULL && strchr(line, '@')== NULL){
                numValuesRead = sscanf(line,  "%f %f", &time_, &C1_);
                if(numValuesRead == 0 && strchr(line, '&')!= NULL){
                  //  puts("Switched to the next residue.");
                    nFrames = i;
                    i = 0;
                    // printACF(t, C1, nFrames);
                    s2 = processACF(t, C1, nFrames);

                    AutocorrFunction acf;
                    acf.t=t;
                    acf.C1 = C1;
                    acf.frames_number = nFrames;

                    SimpleModelState* state;
					state->s2 = 0.5;
					state->tau_s = 50;
                    calcSimpleApproad(&acf, state);

                    ++nResidue;
                    fprintf(output,"%d \t%f \t%f \n", nResidue, state->s2, state->tau_s);
                  //  exit(1);
                }else{
                    t[i] = time_; C1[i] = C1_; i++;
            }}
        }
        fclose(file);
        fclose(output);
    } else {
        perror(acfFile); /* why didn't the file open? */
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


