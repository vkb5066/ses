#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "def.h"

int main(int argc, char** argv){
    uint seed;
    job runparams; uchar err;

    //we'll need a random seed - get it and print it if not supplied by cla
    seed = (argc == 2)? (uint)atoi(argv[1]) : time(NULL);
    srand(seed);

    writeheader(seed);

    //Read main input file, figure out what to do
    if((err = readjob(&runparams))) goto quitjob;

    switch(runparams.runmode){
        //---Normal run modes---------------------------------------------------
        ///??   
        case (uchar)0:
            break;
        ///static bandstructure (no relaxation, output eigenvals to file)
        case (uchar)1:
            runstaticbs(&runparams);
            break;
        ///atomic relaxation (no band structure, output final positions to file)
        case (uchar)2:
            runmd(&runparams);
            break;

        //---Testing modes------------------------------------------------------
#if(COMPILE_TEST_ROUTINES)
        ///read input files
        case (uchar)50:
            testio(runparams);
            break;
        ///test hermitian qr eigensolver
        case (uchar)60:
            testqrh();
            break;
        ///test recip-space hamil gen
        case (uchar)70:
            testexphamilgen(runparams);
            break;
#endif
        //---Other--------------------------------------------------------------
#if(COMPILE_PSO_FILE)
        ///pso fitting
        case (uchar)100:
            psomin(&runparams);
            break;
#endif
        ///Unknown code error
        default:
            goto quitbadmode;
    }


    
    puts("");
    return 0;

    quitjob: 
    puts("main(): cannot read "INFILE_RUNPARAMS);
    return 1;
    quitbadmode:
    puts("main(): unrecognized runmode in "INFILE_RUNPARAMS);
    return 2;
}