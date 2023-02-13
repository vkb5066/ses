#include <stdlib.h>
#include <stdio.h>
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
    //(qr < jd < dav)
    if(runparams->entol < ZERO){
        switch(diagmode){
            case (uchar)0: ///qr
                runparams->entol = DEF_QR_EPS;
                break;
            case (uchar)1: ///jd
                runparams->entol = DEF_JD_EPS;
                break;
            case (uchar)2: ///dav
                runparams->entol = DEF_DA_EPS;
                break;
            default: ///??? just set to the smallest value we have
                runparams->entol = DEF_QR_EPS;
        }
    }
}





//routine 1: static bandstructure
void runstaticbs(job* runparams){
    uchar sc; uint i; ushort j, acc;

    lat latt; pot*restrict pots; ushort nkpts; kpt*restrict kpts;
    hamil ham; fp gcut; ushort maxmill[3]; uint maxbasis;    
    uchar err;


    puts("static band structure");

    //io
    if((err = readlattice(&latt))) goto quitlat;
    if((err = readpotential(latt.nspecs, &pots))) goto quitpot;   
    if((err = readkpoints(&nkpts, &kpts))) goto quitkpt;

    //finish initializing everything needed for main routine
    initjobdefs(latt.nspecs, latt.speccounts, pots, runparams->diagmode, 
                runparams);
    tocrdsc(&latt);
    torecipbasis(latt.A); ///latt.A -> latt.B
    gcut = SQRT(runparams->encut/HBARSQD_OVER_TWOM);
    setmaxdims(gcut, runparams->fftgmul, latt.A, nkpts, kpts,
               maxmill, ham.dims, ham.l2dims, &maxbasis);

    //allocate, set V(G)
    ham.vloc = calloc(ham.dims[0]*ham.dims[1]*ham.dims[2], sizeof(cfp));
    sethamvloc(latt, pots, &ham);
    
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
    case (uchar)0: //-----------------------------------------------------------
        cfp*restrict A;
#if(!PARALLEL)
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
            sethamkin(gcut*gcut, maxmill, latt.A, runparams->meff, &ham);
            kpts[j].evals = malloc(ham.npw*sizeof(fp));
            ////direct diagonalization
            buildexphamil(ham, A);
            qrh(ham.npw, A, kpts[j].evals, kpts[j].evecs, 
                QR_ITR_LIM, runparams->entol, (uchar)0);
        }
        puts("");

        ///clean up the stuff we don't need anymore (basically all but kpts)
        free(ham.vloc);
        free(ham.mills);
        free(ham.ke);
        free(A);
#elif(PARALLEL)

#endif
        break;


    case(uchar)1://-------------------------------------------------------------
        ///TODO: this will be jacobi-davidson method
        break;


    ///hermitian davidson routine, iterative solution of kpt n-1 used as a
    ///starting guess for kpt n
    case (uchar)2://------------------------------------------------------------
        uint maxdavbasis;
        uint prevbasissize;
        cfp*restrict V, *restrict Q, *restrict W;
        cfp*restrict*restrict L;
        cfp*restrict psig1, *restrict psig2;        
        cfp*restrict*restrict vecbuf1, *restrict*restrict vecbuf2;
        ///allocate all the reusable space we'll need
        ham.mills = malloc(3u*maxbasis*sizeof(short));
        ham.ke = malloc(maxbasis*sizeof(fp));
        maxdavbasis = runparams->mbsmul * runparams->nbands;
        V = malloc(maxdavbasis*maxbasis*sizeof(cfp));
        Q = malloc(maxdavbasis*maxbasis*sizeof(cfp));
        W = malloc(maxdavbasis*maxbasis*sizeof(cfp));
        L = malloc(maxdavbasis*sizeof(cfp*));
        for(i = 0u; i < maxdavbasis; ++i) L[i] = malloc(maxbasis*sizeof(cfp));
        psig1 = malloc(ham.dims[0]*ham.dims[1]*ham.dims[2]*sizeof(cfp));
        psig2 = malloc(ham.dims[0]*ham.dims[1]*ham.dims[2]*sizeof(cfp));
        vecbuf1 = malloc(runparams->nbands*sizeof(cfp*));
        for(j = (ushort)0; j < runparams->nbands; ++j) 
            vecbuf1[j] = malloc(maxbasis*sizeof(cfp));
        vecbuf2 = malloc(runparams->nbands*sizeof(cfp*));
        for(j = (ushort)0; j < runparams->nbands; ++j) 
            vecbuf2[j] = malloc(maxbasis*sizeof(cfp));

        //need to get V(r) from V(G)
        fftconv3dif(ham.dims, ham.l2dims, ham.vloc, psig1); 

        ///to begin, do a cold-start on the first k-point with random initial
        ///guess vectors
        fputs("working on k-point               ", stdout);
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05hu of %05hu", (ushort)1, nkpts);
        ham.kp = kpts;
        sethamkin(gcut*gcut, maxmill, latt.A, runparams->meff, &ham);
        kpts[0].evals = malloc(ham.npw*sizeof(fp));
        dav(ham.npw, ham, maxdavbasis, V, Q, W, L, psig1, psig2, 
            0u, 0u, NULL, 
            runparams->nbands, kpts[0].evals, vecbuf1, 
            0u, ZERO, DA_ITR_LIM, runparams->entol);
        prevbasissize = ham.npw;

        ///now that we have an initial guess, move through the rest of the
        ///k-points and interchange initial guess buffers
        for(j = (ushort)1; j < nkpts; ++j){
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05hu of %05hu", 
                   j + (ushort)1, nkpts);
            ///setup hamiltonian, kpoints
            ham.kp = kpts + j;
            sethamkin(gcut*gcut, maxmill, latt.A, runparams->meff, &ham);
            kpts[j].evals = malloc(ham.npw*sizeof(fp));
            ////iterative diagonalization
            if(j%(ushort)2 == (ushort)1){
                dav(ham.npw, ham, maxdavbasis, V, Q, W, L, psig1, psig2, 
                    runparams->nbands, prevbasissize, 
                    (const cfp*restrict*restrict)vecbuf1, 
                    runparams->nbands, kpts[j].evals, vecbuf2, 
                    0u, ZERO, DA_ITR_LIM, runparams->entol);
            }
            else{
                dav(ham.npw, ham, maxdavbasis, V, Q, W, L, psig1, psig2, 
                    runparams->nbands, prevbasissize, 
                    (const cfp*restrict*restrict)vecbuf2, 
                    runparams->nbands, kpts[j].evals, vecbuf1, 
                    0u, ZERO, DA_ITR_LIM, runparams->entol);
            }
            prevbasissize = ham.npw;
        }
        puts("");


        ///clean up the stuff we don't need anymore (basically all but kpts)
        free(ham.vloc);
        free(ham.mills);
        free(ham.ke);
        free(V); free(Q); free(W);
        for(i = 0u; i < maxdavbasis; ++i) free(L[i]);
        free((cfp**)L);
        free(psig1); free(psig2);
        for(j = (ushort)0; j < runparams->nbands; ++j){
            free(vecbuf1[j]);
            free(vecbuf2[j]);
        }
        free((cfp**)vecbuf1); free((cfp**)vecbuf2);


    default: ///nothing needed here: diagmode is user-set or set-to-default 
        break;
    }

    //ok, now we have all of the eigenvals.  print them out
    writebands(latt.A, MIN(runparams->nbands, maxbasis), nkpts, kpts);
    
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