#include <stdio.h>
#include "def.h"

int main(int argc, char** argv){
    job runparams; uchar err;
 
    writeheader();

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
        ///Unknown code error
        default:
            goto quitbadmode;
    }


    

    return 0;

    quitjob: 
    puts("main(): cannot read "INFILE_RUNPARAMS);
    return 1;
    quitbadmode:
    puts("main(): unrecognized runmode in "INFILE_RUNPARAMS);
    return 2;
}