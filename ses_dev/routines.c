#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "def.h"

//Sets sensible defaults for job control, given some other info
static forceinline void initjobdefs(const uchar nelems, 
                                    const ushort*restrict elemcounts,
                                    const pot*restrict pots, 
                                    const uchar diagmode, 
                                    job*restrict runparams){
    uchar i;    
    fp g2max;

    //if no bands are set, a reasonable default is the number of val electrons
    if(runparams->nbands == (ushort)0){
        for(i = (uchar)0; i < nelems; ++i){
            runparams->nbands += elemcounts[i] * (ushort)pots[i].nvel;
        }
    }
    //if no energy cutoff is set, a reasonable default is the energy 
    //corresponding to the maximum sampled G value in the pot files
    if(runparams->encut < ZERO){
        g2max = pots[0].q2hi;
        for(i = (uchar)1; i < nelems; ++i){
            if(pots[i].q2hi > g2max) g2max = pots[i].q2hi;
        }
        runparams->encut = HBARSQD_OVER_TWOM * g2max;
    }
    //if no energy tol is set, a reasonable value depends on the diag mode
    if(runparams->entol < ZERO){
        switch(diagmode){
            case (uchar)0: ///qr
                runparams->entol = DEF_QR_EPS;
                break;
            case (uchar)1: ///jd
                runparams->entol = DEF_JD_EPS;
                break;
            case (uchar)3: ///lopcg
                runparams->entol = DEF_CG_EPS;
                break;
            default: ///??? just set to the smallest value we have
                runparams->entol = DEF_QR_EPS;
        }
    }

    //similarly, if no itr limit is set, reasonable values depend on diag mode
    if(runparams->itrlim == (ushort)0){
        switch(diagmode){
            case (uchar)0: ///qr
                runparams->itrlim = DEF_QR_ITRLIM;
                break;
            case (uchar)1: ///jd
                runparams->itrlim = DEF_JD_ITRLIM;
                break;
            case (uchar)3: ///lopcg
                runparams->itrlim = DEF_CG_ITRLIM;
                break;
            default: ///??? just set to the smallest value we have
                runparams->entol = DEF_QR_EPS;
        }
    }
}



//routine 1: static bandstructure
void runstaticbs(job* runparams){
    uchar sc; uint i; ushort j, acc;

    lat latt; fp B[3][3]; pot*restrict pots; ushort nkpts; kpt*restrict kpts;
    hamil ham; fp gcut; ushort maxmill[3]; uint maxbasis;    
    uchar err;


    puts("static band structure");
    if(runparams->nfstargets) puts("(w/ folded spectrum method)");
    

    //io
    if((err = readlattice(&latt))) goto quitlat;
    if((err = readpotential(latt.nspecs, &pots))) goto quitpot;
    memcpy(B, latt.A, 9*sizeof(fp));
    torecipbasis(B); ///latt.A -> B   
    if((err = readkpoints(B, &nkpts, &kpts))) goto quitkpt;

    //finish initializing everything needed for main routine
    initjobdefs(latt.nspecs, latt.speccounts, pots, runparams->diagmode, 
                runparams);
    tocrdsc(latt);
    gcut = SQRT(runparams->encut/HBARSQD_OVER_TWOM);
    setmaxdims(gcut, runparams->fftgmul, B, nkpts, kpts,
               maxmill, ham.dims, ham.l2dims, &maxbasis);

    //allocate, set V(G)
    ham.vloc = calloc(ham.dims[0]*ham.dims[1]*ham.dims[2], sizeof(cfp));
    sethamvloc(latt, B, pots, &ham);
    
    //we no longer need a lot of this stuff - get rid of it
    for(sc = (uchar)0, acc=(ushort)0; sc < latt.nspecs; ++sc){
        for(j = (ushort)0; j < latt.speccounts[sc]; ++j, ++acc)
            free(latt.sites[acc].spec);
        free(pots[sc].samples);
    }
    free(latt.sites);
    free(latt.speccounts);
    free(pots);

    //here is the important part: actually calculating the eigenpairs
    //exactly how this is done depends strongly on the 'diagmode'
    switch(runparams->diagmode){
    ///hermitian qr algo, builds explicit hamiltonians for each kpt
    case (uchar)0:{//-----------------------------------------------------------
        cfp*restrict A;

        ///allocate all the reusable space we'll need
        ham.mills = malloc(3u*maxbasis*sizeof(short));
        ham.ke = malloc(maxbasis*sizeof(fp));
        A = malloc(maxbasis*maxbasis*sizeof(cfp)); ///will hold explicit ham
    
        ///solve for the eigenvals of each k-point
        fputs("working on k-point               ", stdout);
        for(j = (ushort)0; j < nkpts; ++j){
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05hu of %05hu", 
                   j + (ushort)1, nkpts);
            ///setup hamiltonian, kpoints
            ham.kp = kpts + j;
            sethamkin(gcut*gcut, maxmill, B, runparams->meff, &ham);
            kpts[j].evals = malloc(ham.npw*sizeof(fp));
            ////direct diagonalization
            buildexphamil(ham, A);
            qrh(ham.npw, A, kpts[j].evals, kpts[j].evecs, 
                runparams->itrlim, runparams->entol, (uchar)0);
        }
        puts("");

        writebands('0', B, MIN(runparams->nbands, maxbasis), nkpts, kpts);

        ///clean up the stuff we don't need anymore (basically all but kpts)
        free(ham.vloc);
        free(ham.mills);
        free(ham.ke);
        free(A);

        break;
    }

    case(uchar)1:{//------------------------------------------------------------
        ///TODO: this will be jacobi-davidson method
        break;

    }

    //LOPCG iterative routine (now with fancy preconditioning!)
    case (uchar)3:{//-----------------------------------------------------------
        uchar emul;
        cfp*restrict*restrict vecbuf1, *restrict*restrict vecbuf2;
        uchar*restrict IWS;
        fp*restrict RWS; cfp*restrict CWS; cfp*restrict*restrict CCWS;        
        uint prevbasissize;

        ///if we're not doing fsm, need to set the target to a dummy value
        fp dum = ZERO;
        if(!(runparams->nfstargets)) runparams->fstargets = &dum;

        ///allocate all the reusable space we'll need
        ham.mills = malloc(3u*maxbasis*sizeof(short));
        ham.ke = malloc(maxbasis*sizeof(fp));
#if USE_BAD_PRECONDITIONER
        if(runparams->nfstargets) ham.la = malloc(maxbasis*sizeof(fp));
#endif
        for(j = (ushort)0; j < nkpts; ++j){
            kpts[j].evals = malloc(maxbasis*sizeof(fp));
        }
        IWS = malloc(runparams->nbands*sizeof(uchar));
        RWS = malloc((5*runparams->nbands + maxbasis)*sizeof(fp));
        CWS = malloc((18*runparams->nbands*runparams->nbands +
                      6*runparams->nbands*maxbasis +
                      maxbasis + 
                      2*ham.dims[0]*ham.dims[1]*ham.dims[2])*sizeof(cfp));
        CCWS = malloc(3*runparams->nbands*sizeof(cfp*));
        for(j = (ushort)0; j < 3*runparams->nbands; ++j){
            CCWS[j] = malloc(3*runparams->nbands*sizeof(cfp));
        }
        vecbuf1 = malloc(runparams->nbands*sizeof(cfp*));
        for(j = (ushort)0; j < runparams->nbands; ++j) 
            vecbuf1[j] = malloc(maxbasis*sizeof(cfp));
        vecbuf2 = malloc(runparams->nbands*sizeof(cfp*));
        for(j = (ushort)0; j < runparams->nbands; ++j) 
            vecbuf2[j] = malloc(maxbasis*sizeof(cfp));

        //need to get V(r) from V(G)
        fftconv3dif(ham.dims, ham.l2dims, ham.vloc, CWS); 

        //fsm loop: even if fsm isn't being used, this will run at least once
        sc = (uchar)0;
        do{
            printf("loop %hhu of %hhu\n", sc + (uchar)1, runparams->nfstargets);

            ///to begin, do a cold-start on the first k-point with random init
            ///guess vectors
            printf("working on k-point 00001 of %05hu\n", nkpts);
            ham.kp = kpts;
            sethamkin(gcut*gcut, maxmill, B, runparams->meff, &ham);
#if USE_BAD_PRECONDITIONER
            if(runparams->nfstargets) sethamlap(B, runparams->meff, &ham);
#endif
            emul = lcg(runparams->nbands, ham.npw, ham, IWS, RWS, CWS, CCWS, 
                       0u, 0u, NULL, ham.kp->evals, vecbuf1, 
                       runparams->nfstargets, runparams->fstargets[sc], 
                       runparams->itrlim, runparams->entol, runparams->stab);

            ///now that we have an initial guess, do the rest of the k-points,
            ///exchanging initial guess buffers
            for(j = (ushort)1; j < nkpts; ++j){
                fputs("\33[2K\r", stdout); ///clear entire row (S.O. #1508490)
                fputs("\033[1A", stdout); ///move up a line (S.O. #26423537)
                printf("working on k-point %05hu of %05hu\n", 
                       j + (ushort)1, nkpts);
                prevbasissize = ham.npw;
                ham.kp = kpts + j;
                sethamkin(gcut*gcut, maxmill, B, runparams->meff, &ham);
#if USE_BAD_PRECONDITIONER
                if(runparams->nfstargets) sethamlap(B, runparams->meff, &ham);
#endif
                emul = emul? (uchar)0 : (uchar)1; ///no reuse vecs if cg failed
                //emul = (uchar)0; ///somethimes this works better???
                if(j%(ushort)2 == 1){
                    emul = lcg(runparams->nbands, ham.npw, ham, IWS, RWS, CWS, CCWS, 
                               emul*runparams->nbands, emul*prevbasissize, (const cfp*restrict*restrict)vecbuf1, ham.kp->evals, vecbuf2, 
                               runparams->nfstargets, runparams->fstargets[sc], 
                               runparams->itrlim, runparams->entol, runparams->stab);
                }
                else{
                    emul = lcg(runparams->nbands, ham.npw, ham, IWS, RWS, CWS, CCWS, 
                               emul*runparams->nbands, emul*prevbasissize, (const cfp*restrict*restrict)vecbuf2, ham.kp->evals, vecbuf1, 
                               runparams->nfstargets, runparams->fstargets[sc], 
                               runparams->itrlim, runparams->entol, runparams->stab);
                }
            }

            ///if we are doing FSM, we're about to overwrite the evals - print
            ///them out before that happens
            writebands('0' + sc, B, MIN(runparams->nbands, maxbasis), 
                       nkpts, kpts);

            sc++;
        } while(sc < runparams->nfstargets);


        ///clean up the stuff we don't need anymore (basically all but kpts)
        free(ham.vloc);
        free(ham.mills);
        free(ham.ke);
        if(runparams->nfstargets) free(runparams->fstargets);
        free(IWS); free(RWS); free(CWS);
        for(j = (ushort)0; j < runparams->nbands; ++j) free(vecbuf1[j]);
        free((cfp**)vecbuf1);
        for(j = (ushort)0; j < runparams->nbands; ++j) free(vecbuf2[j]);
        free((cfp**)vecbuf2);
        for(i = 0u; i < 3*runparams->nbands; ++i) free(CCWS[i]);
        free((cfp**)CCWS);

        break;
    }

    default: ///nothing needed here: diagmode is user-set or set-to-default 
        break;
    }
    
    //all thats left to be freed are the kpoint-related values
    for(j = (ushort)0; j < nkpts; ++j) free(kpts[j].evals);
    free(kpts);


    //end conditions
    ///TODO: write out some program end strings, maybe timing?
    return; ///everything is okay

    quitlat: 
    puts("runstaticbs(): cannot open "INFILE_LATTICE);
    exit(EXIT_FAILURE);
    quitpot: 
    puts("runstaticbs(): cannot open "INFILE_POTENTIAL);
    exit(EXIT_FAILURE);
    quitkpt:
    if(err == (uchar)1)      puts("runstaticbs(): cannot open "INFILE_KPOINTS);
    else if(err == (uchar)2) puts("runstaticbs(): unrecognized control "
                                  "parameter in "INFILE_KPOINTS);
    exit(EXIT_FAILURE);
}




void runmd(job* runparams){
    uchar i; ushort j, nsites;

    lat latt; pot*restrict potentials; keat*restrict keatings;
    ushort nxsites; site* xsites;
    uchar err;
    fp*restrict RWS;
    fp Ai[3][3];

    puts("keating molecular dynamics");

    //io
    if((err = readlattice(&latt))) goto quitlat;
    if((err = readpotential(latt.nspecs, &potentials))) goto quitpot;
    if((err = readkeating(latt.nspecs, &keatings))) goto quitkeat;

    //get pairs, triplets necessary for keating potential
    tocrdsc(latt); ///dont forget to do this or your program will explode :)
    setnbrs(latt, &nxsites, &xsites, 
            runparams->rcut, DEG2RAD(runparams->acut), DEF_MD_EPS);

    //the work space necessary for the verlet integration
    for(i = (uchar)0, nsites = (ushort)0; i < latt.nspecs; ++i)
    for(j = (ushort)0; j < latt.speccounts[i]; ++j) nsites++;
    RWS = malloc(9u*nsites*sizeof(fp));

    relax(RWS, latt, potentials, keatings, 
          runparams->timestep, runparams->ftol, runparams->mditrlim);

    //print out the resulting positions
    memcpy(Ai, latt.A, 9u*sizeof(fp)); toinverse(Ai); tocrdsf(Ai, latt);
    backtocell(latt, (uchar)0);
    writelattice(latt, "Ga As");

    //end conditions
    ///TODO: write out some program end strings, maybe timing?
    return;


    quitlat: 
    puts("runmd(): cannot open "INFILE_LATTICE);
    exit(EXIT_FAILURE);
    quitpot: 
    puts("runmd(): cannot open "INFILE_POTENTIAL);
    exit(EXIT_FAILURE);
    quitkeat:
    puts("runmd(): cannot open "INFILE_KEATING);
    exit(EXIT_FAILURE);    
}

