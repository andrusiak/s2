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

#define NMAX 10000
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct AutocorrFunction{
	float* t;
	float* C1;
	long frames_number;
};

void printACF(float * t, float * C1, int n){
    for(int i = 0; i<n; i++){
     printf("time = %5.1f , C1 = %f\n", *(t+i),*(C1+i));
    }
}

void calcSimpleApproad(AutocorrFunction* acf, float& s2, float& tau_s){

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
                    AutocorrFunction acf = {t1, c2, nFrames};
                    float tau_s;
                    calcSimpleApproad(&acf, s2, tau_s);
                    ++nResidue;
                    fprintf(output,"%d \t%f \n", nResidue, s2);
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


