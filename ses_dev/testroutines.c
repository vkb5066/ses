//Contains main function routines that are called from main
#include "def.h"


//---Testing routines-----------------------------------------------------------
#if(COMPILE_TEST_ROUTINES)

#include <stdlib.h>
#include <stdio.h>

static forceinline void _tstcmul(const cfp a, const cfp b, cfp* res){
    res->re = a.re*b.re - a.im*b.im;
    res->im = a.re*b.im + a.im*b.re;
}

//routine id 50: test basic io functions
void testio(job runparams){
    uchar i; ushort j, acc;
    lat lattice; pot* pots; ushort nkpts; kpt* kpts;
    uchar err;

    printf("parsing run param file (after setting defaults) ...\n");
    INIT_RUNPARAMS(runparams);
    err = readjob(&runparams);
    printf("readjob() returned with exit code %hhu\n", err);
    if(!err){ 
        reportjob(&runparams);   
        ///nothing to free here
    }

    printf("reading kpoints file ...\n");
    err = readkpoints(&nkpts, &kpts);
    printf("readkpoints() returned with exit code %hhu\n", err);
    if(!err){
        reportkpoints(nkpts, kpts);
        free(kpts);
    }

    printf("reading lattice file ...\n");
    err = readlattice(&lattice);    
    printf("readlattice() returned with exit code %hhu\n", err);
    if(!err){
        reportlattice(&lattice);
        for(i = (uchar)0, acc=(ushort)0; i < lattice.nspecs; ++i){
            for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc){
                free(lattice.sites[acc].spec);
            }
        }
        free(lattice.sites);
        free(lattice.speccounts);
    }

    printf("reading potential file ...\n");
    err = readpotential(lattice.nspecs, &pots);    
    printf("readpotential() returned with exit code %hhu\n", err);
    if(!err){
        reportpotential(lattice.nspecs, pots);
        for(i = (uchar)0; i < lattice.nspecs; ++i){
            free(pots[i].samples);
        }
        free(pots);
    }
}


//routine id 60: test single-thread hermitian qr eigenvalue algo
void testqrh(){
#define SIZE 6 ///6x6
    int i, j, k;
    fp KNOWNVALS[SIZE];

    fp*restrict vals = malloc(SIZE*sizeof(fp));             // all of this 
    cfp*restrict*restrict vecs = malloc(SIZE*sizeof(cfp*)); // must be malloced
    for(i = 0; i < SIZE; ++i){                              // to match the
        vecs[i] = malloc(SIZE*sizeof(cfp));                 // normal use case
    }                                                       // of my function
    cfp* M = malloc(SIZE*SIZE*sizeof(cfp));
    cfp* O = malloc(SIZE*SIZE*sizeof(cfp));

    cfp Mv[SIZE];
    cfp tmp;
    cfp sum;


    //Set hermitian matrix with known eigenvalues
M[0].re = (fp)+0.32411008980085400; M[0].im = (fp)+0.000000000000000000;
M[1].re = (fp)-5.33266261725691000; M[1].im = (fp)-9.772591442025960000;
M[2].re = (fp)-9.03176921351199100; M[2].im = (fp)+6.702850351076897000;
M[3].re = (fp)+1.89533626208986840; M[3].im = (fp)-0.994990514117287000;
M[4].re = (fp)-8.70554688311609800; M[4].im = (fp)-2.148895546957707000;
M[5].re = (fp)+9.11241246255341900; M[5].im = (fp)-9.456381740239134000;
M[6].re = (fp)-5.33266261725691000; M[6].im = (fp)+9.772591442025960000;
M[7].re = (fp)-6.13695200069311200; M[7].im = (fp)+0.000000000000000000;
M[8].re = (fp)-9.25364347589143700; M[8].im = (fp)-0.035430525813229250;
M[9].re = (fp)-1.21284983509976390; M[9].im = (fp)-3.495627163999025000;
M[10].re = (fp)+8.1988465577028930; M[10].im = (fp)+7.25894099959833000;
M[11].re = (fp)-0.5799012334070497; M[11].im = (fp)+4.85096421383668600;
M[12].re = (fp)-9.0317692135119910; M[12].im = (fp)-6.70285035107689700;
M[13].re = (fp)-9.2536434758914370; M[13].im = (fp)+0.03543052581322925;
M[14].re = (fp)-2.1776102305040830; M[14].im = (fp)+0.00000000000000000;
M[15].re = (fp)-6.9795423967650830; M[15].im = (fp)+0.23312622066087485;
M[16].re = (fp)-7.7407052887145330; M[16].im = (fp)+0.06488893521261296;
M[17].re = (fp)+7.1455233286559190; M[17].im = (fp)-5.44665833516633200;
M[18].re = (fp)+1.8953362620898684; M[18].im = (fp)+0.99499051411728700;
M[19].re = (fp)-1.2128498350997639; M[19].im = (fp)+3.49562716399902500;
M[20].re = (fp)-6.9795423967650830; M[20].im = (fp)-0.23312622066087485;
M[21].re = (fp)+9.9961304390338130; M[21].im = (fp)+0.00000000000000000;
M[22].re = (fp)-8.7876822080343490; M[22].im = (fp)-6.14201669814857500;
M[23].re = (fp)+3.7768454368288220; M[23].im = (fp)-1.77043715050100130;
M[24].re = (fp)-8.7055468831160980; M[24].im = (fp)+2.14889554695770700;
M[25].re = (fp)+8.1988465577028930; M[25].im = (fp)-7.25894099959833000;
M[26].re = (fp)-7.7407052887145330; M[26].im = (fp)-0.06488893521261296;
M[27].re = (fp)-8.7876822080343490; M[27].im = (fp)+6.14201669814857500;
M[28].re = (fp)+1.8314309230066605; M[28].im = (fp)+0.00000000000000000;
M[29].re = (fp)+4.2014919783389630; M[29].im = (fp)+3.65335513916172160;
M[30].re = (fp)+9.1124124625534190; M[30].im = (fp)+9.45638174023913400;
M[31].re = (fp)-0.5799012334070497; M[31].im = (fp)-4.85096421383668600;
M[32].re = (fp)+7.1455233286559190; M[32].im = (fp)+5.44665833516633200;
M[33].re = (fp)+3.7768454368288220; M[33].im = (fp)+1.77043715050100130;
M[34].re = (fp)+4.2014919783389630; M[34].im = (fp)-3.65335513916172160;
M[35].re = (fp)-2.1015271892900690; M[35].im = (fp)+0.00000000000000000;

    for(i = 0; i < SIZE*SIZE; ++i){
        O[i].re = M[i].re; O[i].im = M[i].im; ///matrix gets destroyed in qrh 
    }                                         /// so need a deepcopy

    KNOWNVALS[0] = (fp)-30.2449695770542560;
    KNOWNVALS[1] = (fp)-15.8871382822438400;
    KNOWNVALS[2] = (fp)-03.2038564181751448;
    KNOWNVALS[3] = (fp)+04.4926944551792385;
    KNOWNVALS[4] = (fp)+17.1532735244540700;
    KNOWNVALS[5] = (fp)+29.4255783291940420;


    //Actually use the eigensolver
    qrh((uint)SIZE, &O, &vals, &vecs, 25u, 1.0e-5F, (uchar)1);

    //Print out differences in eigenvalues between known and test
    printf("diff in eigs:\n");
    for(i = 0; i < SIZE; ++i){
        printf(FPFORMAT"\n", ABSF(vals[i] - KNOWNVALS[i]));
    }
    //Print out the error in Mv - lv = 0 ... for (l, v) = an eigenpair
    for(i = 0; i < SIZE; ++i){ ///eigenval loop
        printf("err in Mv - lv for eig %i: ", i);
        for(j = 0; j < SIZE; ++j){ ///rows of M
            Mv[j].re = ZERO; Mv[j].im = ZERO;
            for(k = 0; k < SIZE; ++k){ ///"cols" of vecs (exc. vecs is trnspsd)
                _tstcmul(M[ACC2(SIZE, j, k)], vecs[i][k], &tmp);
                Mv[j].re += tmp.re; 
                Mv[j].im += tmp.im;
            }
        }
        sum.re = ZERO; sum.im = ZERO;
        for(j = 0; j < SIZE; ++j){
            sum.re += ABSF(Mv[j].re - vals[i]*vecs[i][j].re);
            sum.im += ABSF(Mv[j].im - vals[i]*vecs[i][j].im);
        }
        printf(FPFORMAT" + "FPFORMAT"i\n", sum.re, sum.im);
    }



    //Clean up
    free(vals);
    for(i = 0; i < SIZE; ++i){
        free(vecs[i]);
    }
    free((cfp*)vecs);
    free(M); free(O);

#undef SIZE
}





//routine id 70: test the eigenvalues of explicit hamiltonian generation
void testexphamilgen(job runparams){
    uint i;
    uchar ii; ushort acc, j;
    lat lattice; pot* pots; ushort nkpts; kpt* kpts;
    hamil ham;
    ushort maxmill[3]; uint maxbasis;    

    fp gcut; 
    cfp*restrict H; fp*restrict va; cfp*restrict*restrict ve;

    readlattice(&lattice);
    readpotential(lattice.nspecs, &pots);
    readkpoints(&nkpts, &kpts);

    tocrdsc(&lattice);
    torecipbasis(lattice.A);

    gcut = SQRT(runparams.encut/HBARSQD_OVER_TWOM);
    setmaxdims(gcut, runparams.fftgmul, lattice.A, nkpts, kpts,
               maxmill, ham.dims, ham.l2dims, &maxbasis);

    //setup the hamiltonian
    ham.vloc = calloc(ham.dims[0]*ham.dims[1]*ham.dims[2], sizeof(cfp));
    sethamvloc(lattice, pots, &ham);

    ham.kp = kpts;   
    ham.mills = malloc(3u*maxbasis*sizeof(short));
    ham.ke = malloc(maxbasis*sizeof(fp));
    sethamkin(gcut*gcut, maxmill, lattice.A, runparams.meff, &ham);
 


    H = malloc(ham.npw*ham.npw*sizeof(cfp));
    va = malloc(ham.npw*sizeof(fp));
    ve = malloc(ham.npw*sizeof(cfp*)); ///really, these should be set to kpts
    for(i = 0u; i < ham.npw; ++i){     ///but im lazy
        ve[i] = malloc(ham.npw*sizeof(cfp));                 
    }                 
    buildexphamil(ham, &H);
    uint err = qrh(ham.npw, &H, &va, &ve, 25u, 1.0e-6F, (uchar)0);
    (err++);


    free(va);
    for(i = 0u; i < ham.npw; ++i){
        free(ve[i]);
    }
    free((cfp*)ve);
    for(ii = (uchar)0, acc=(ushort)0; ii < lattice.nspecs; ++ii){
        for(j = (ushort)0; j < lattice.speccounts[ii]; ++j, ++acc){
            free(lattice.sites[acc].spec);
        }
        free(pots[ii].samples);
    }
    free(lattice.sites);
    free(lattice.speccounts);
    free(pots);
    free(kpts);
    free(H);
    free(ham.vloc);
    free(ham.ke);
    free(ham.mills);
}



#endif