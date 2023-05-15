//PSO for fitting force field params to data
//No fancy input/output - this is hard coded ... make sure you go through it
//and change stuff as needed
#include "def.h"
#if COMPILE_PSO_FILE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


//general constants-------------------------------------------------------------
#define STACK_ALLOC 0 ///1 = fast, dangerous.  0 = slow, safe (heap alloc)

#define N_SPECIES 2u
#define N_ATOMS 2u
#define N_MD_SAMPLES 10000u
#define INFILE_MD_SAMPLES "mdsamples" ///input file - see below for format

//simulation constants----------------------------------------------------------
#define MINIMIZE_AVERAGE 0     ///if 0, minimize the max error instead of <err>
#define N_DIM 2u               ///symmetry INEQUIV, IMPORTANT params to fit
                               ///search space bounds [lo, hi]
static const fp SS_BOUNDS[N_DIM][2] = { {(fp)+00.10, (fp)+5.00},   ///dim 0
                                        {(fp)+00.10, (fp)+5.00} }; ///dim 1

#define STRICT_BOUNDS 1u ///force all simulation particles to stay in bounds
#define N_ITER 1000u
#define N_PARTICLES 75u

#define PARAM_W (fp)0.7 ///inertia, must be 0 < w < 1
#define PARAM_P (fp)2.6 ///individual strength, usually 1 < p < 3
#define PARAM_G (fp)1.2 ///global strangth, usually 1 < g < 3

///kickout params---------------------------------------------------------------
#define IMPROV_EPS (fp)1.0e-4 ///min change in global best before kickout reset
#define N_B4_KICKOUT 25u      ///iterations w/o improvement (^) before kickout 
#define OPT_NBRHD_SQD (fp)5.5 ///square distance from optima for kickout
#define KICKOUT_PROB (fp)0.5  ///how likley particle near optima is reset

//map function------------------------------------------------------------------
//consider the harmonic bond strength interactions for covalent compound AB.
//in theory, there are 3 (4) unsymmetric (symmetric) interactions: 
//A-A, A-B, B-A, B-B.  In most situations, however, the only important 
//interaction is the 1st nearest neighbor one, which usually constrains the 
//possible interacations to only 1 (2) unsymmetric (symmetric).  
//We obviously only want to spend time optimizing the important interactions, so
//here we define a routine that maps a pso vector (only unsymmetric, important
//interactions) to a keating object (all interactions).
//In the example above, we'd need to map pso element 0 to keating element 1, 2

//the keating matrix is always laid out as keatings[i] = j is the interaction ij
//the pso vector should go through i, j alpha params (j the fast index) and
//then i, j sqrt(beta) params (j the fast index)

//below maps a pso vector p to a keating set k
//2-species system, only A-B interactions (e.x. GaAs) ... 
// ... 1 interaction each (Ga-As = As-Ga)
#if N_SPECIES == 2u
#define ALPHMAP(p, k) do{ k[0].al[1] =     p[0];\
                          k[1].al[0] =     p[0];\
                        }while(0)
#define BETAMAP(p, k) do{ k[0].sqrtbe[1] = p[1];\
                          k[1].sqrtbe[0] = p[1];\
                        }while(0)
//3-species system, only (A,B)-C interactions (e.x. (Cd,Zn)Te) ...
// ... 2 interactions each (Cd-Te = Te-Cd, Zn-Te = Te-Zn, no cation-cation)
//assumes cations are first, than anions
#elif N_SPECIES == 3u

#endif


//main object-------------------------------------------------------------------
typedef struct part{
    fp vc[N_DIM]; //velocity, current
    fp pc[N_DIM]; //position, current
    fp pb[N_DIM]; //position, best
    fp fb;        //f(pb)
} part;


//pso function definitions------------------------------------------------------
//n-dim distance between a and b, squared
static forceinline fp _dist2(const fp a[restrict N_DIM], 
                             const fp b[restrict N_DIM]){
    uint i; fp res;

    res = b[0] - a[0];
    res = res*res;
    for(i = 1u; i < N_DIM; ++i){
        res += (b[i]-a[i])*(b[i]-a[i]);
    }

    return res;
}

//randomize velocity and current position of particle p
static forceinline void _rset(part *restrict p){
    uint i;

    for(i = 0u; i < N_DIM; ++i)
        p->vc[i] = ABSF(SS_BOUNDS[i][1] - SS_BOUNDS[i][0]) * 
                   ((fp)2.0*((fp)rand()/(fp)RAND_MAX) - ONE);
    for(i = 0u; i < N_DIM; ++i)
        p->pc[i] = SS_BOUNDS[i][0] + ABSF(SS_BOUNDS[i][1] - SS_BOUNDS[i][0]) *
                   ((fp)rand()/(fp)RAND_MAX);
}


//io functions------------------------------------------------------------------
static void readmddata(fp data[restrict N_MD_SAMPLES*N_ATOMS*6u]){
    //input: x, y, z, fx, fy, fz for atom 0, md sample 0 (cartesian coords)
    //       x, y, z, fx, fy, fz for atom 1, md sample 0
    //       ...
    //       x, y, z, fx, fy, fz for atom N_ATOMS-1, md sample 0
    //       x, y, z, fx, fy, fz for atom 0, md sample 1
    //       ... etc until N_MD_SAMPLES - 1

    FILE* infile;
    uint i, acc;

    infile = fopen(INFILE_MD_SAMPLES, "r");
    if(!infile) printf("unable to open file\n");

    for(i = 0u, acc = 0u; i < N_MD_SAMPLES*N_ATOMS; ++i, acc += 6u){
        fscanf(infile, " "FPFORMAT" "FPFORMAT" "FPFORMAT
                       " "FPFORMAT" "FPFORMAT" "FPFORMAT,
               &data[acc+0u], &data[acc+1u], &data[acc+2u], &data[acc+3u], 
               &data[acc+4u], &data[acc+5u]);
    }

    fclose(infile);
    return;
}

//objective eval functions------------------------------------------------------
//calculate force on a site
//copied straight from md.c except this returns nada
//mind that 'res' must be zero before being supplied
#define DOT(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
static forceinline void _setforce(const uchar nspecs, 
                                  const site si, const keat ki, 
                                  fp res[restrict 3]){
    uchar j, k, jj, kk;
    fp rij[3], rik[3];
    fp cij, cijk, inv;

    //two-body term
    for(j = (uchar)0; j < si.npairs; ++j){
        rij[0] = si.pairs[j]->crdsc[0] - si.crdsc[0];
        rij[1] = si.pairs[j]->crdsc[1] - si.crdsc[1];
        rij[2] = si.pairs[j]->crdsc[2] - si.crdsc[2];

        ///4*a_ij * (|r_ij|^2 / r_0ij^2  -  1)
        jj = *(si.pairs[j]->spec);
        cij = (fp)4.0*ki.al[jj]*(DOT(rij, rij)*ki.r0inv[jj]*ki.r0inv[jj] - ONE);

        res[0] += cij*rij[0]; 
        res[1] += cij*rij[1]; 
        res[2] += cij*rij[2];

        ///three-body term
        inv = ONE / (fp)(si.ntrips[j]);
        for(k = (uchar)0; k < si.ntrips[j]; ++k){
            rik[0] = si.trips[j][k]->crdsc[0] - si.crdsc[0];
            rik[1] = si.trips[j][k]->crdsc[1] - si.crdsc[1];
            rik[2] = si.trips[j][k]->crdsc[2] - si.crdsc[2];

            ////2*sqrt(b_ij*b_ik) * (r_ij@r_ik / r_0ij*r_0ik - cos(theta_0jik))
            kk = *(si.trips[j][k]->spec);
            cijk = (fp)2.0*ki.sqrtbe[jj]*ki.sqrtbe[kk] *
                   (DOT(rij, rik)*ki.r0inv[jj]*ki.r0inv[kk] - 
                    ki.c0[ACC2(nspecs, jj, kk)]);
            cijk *= inv;

            res[0] += cijk*(rij[0] + rik[0]); ////x_j + x_k - 2*x_i
            res[1] += cijk*(rij[1] + rik[1]); ////y_j + y_k - 2*y_i
            res[2] += cijk*(rij[2] + rik[2]); ////z_j + z_k - 2*z_i
        }
    }
    inv = ONE / (fp)(si.npairs);    
    res[0] *= inv; res[1] *= inv; res[2] *= inv;

    return;
}
#undef DOT

//evaluate objective function: 
//sum_nmdsamples[ sum_natoms( |f_trial - f_ref|^2 ) / natoms ] / nmdsamples
static forceinline fp _feval(const fp psov[restrict N_DIM],
                             const fp*restrict mddata, site*restrict sites,
                             const keat*restrict keatings){
    uint i, j, accm, acca;
    ushort ii;
    fp del[3];
#if MINIMIZE_AVERAGE
    fp erratom, errsample;
#else
    fp errths, errmax;
#endif

    //map the symmetric pso vector into the unsymmetric keating data
    ALPHMAP(psov, keatings);
    BETAMAP(psov, keatings);

    //for each md sample...
#if MINIMIZE_AVERAGE
    errsample = ZERO;
#else
    errmax = ZERO;
#endif
    for(i = 0u, accm = 0u; i < N_MD_SAMPLES; ++i, accm += 6u*N_ATOMS){
        ///shift all atom coords to their sample points
        for(j = 0u, acca = accm; j < N_ATOMS; ++j, acca += 6u){
            if(sites[j].nimgs){
                del[0] = mddata[acca+0u] - sites[j].crdsc[0];
                del[1] = mddata[acca+1u] - sites[j].crdsc[1];
                del[2] = mddata[acca+2u] - sites[j].crdsc[2];
            }

            memcpy(sites[j].crdsc, mddata+acca, 3u*sizeof(fp));
            for(ii = (ushort)0; ii < sites[j].nimgs; ++ii){
                sites[j].imgs[ii]->crdsc[0] += del[0];
                sites[j].imgs[ii]->crdsc[1] += del[1];
                sites[j].imgs[ii]->crdsc[2] += del[2];
            }
        }

        ///eval < |trial force - ref force|^2 > over all atoms
#if MINIMIZE_AVERAGE
        erratom = ZERO;
#endif
        for(j = 0u, acca = accm; j < N_ATOMS; ++j, acca += 6u){
            memset(del, 0, 3u*sizeof(fp));
            _setforce(N_SPECIES, sites[j], keatings[*sites[j].spec], del);
#if MINIMIZE_AVERAGE
            erratom += (del[0] - mddata[acca+3u])*(del[0] - mddata[acca+3u]) +
                       (del[1] - mddata[acca+4u])*(del[1] - mddata[acca+4u]) +
                       (del[2] - mddata[acca+5u])*(del[2] - mddata[acca+5u]);
#else
            errths = (del[0] - mddata[acca+3u])*(del[0] - mddata[acca+3u]) +
                     (del[1] - mddata[acca+4u])*(del[1] - mddata[acca+4u]) +
                     (del[2] - mddata[acca+5u])*(del[2] - mddata[acca+5u]);
            errmax = MAX(errmax, errths);
#endif
        }
#if MINIMIZE_AVERAGE        
        errsample += (erratom / (fp)N_ATOMS);
#endif
    }
#if MINIMIZE_AVERAGE
    return errsample / (fp)N_MD_SAMPLES;
#else
    return errmax;
#endif
}


//main function-----------------------------------------------------------------
//base pso algorithm from 
//Clerc, "Standard Particle Swarm Optimization" (2012) - HAL Open Access Archive
//kickout/restart procedure by
//Barone, "A way to Hopefully make this Work Better" (2023) - My Brain Just Now
void psomin(job* runparams){
    uint i, j, nkicks, count, nslu;
    fp pbest[N_DIM], fbest, ftmp;
    lat lattice; keat*restrict keatings;
    fp*restrict mddata; mddata = malloc(N_MD_SAMPLES*N_ATOMS*6u*sizeof(fp));
    ushort nxsites; site* xsites;
    #if STACK_ALLOC
    part parts[N_PARTICLES];
    #else
    part*restrict parts; parts = malloc(N_PARTICLES*sizeof(part));
    #endif


    //initialization (for function eval, etc)
    printf("reading lattice file from %s ... ", INFILE_LATTICE);
    readlattice(&lattice); printf("ok\n");
    printf("reading md sample data from %s ... ", INFILE_MD_SAMPLES);
    readmddata(mddata); printf("ok\n");    
    printf("reading keating params from %s ... ", INFILE_KEATING);
    readkeating(N_SPECIES, &keatings); printf("ok\n");

    tocrdsc(lattice);
    setnbrs(lattice, &nxsites, &xsites, 
            runparams->rcut, DEG2RAD(runparams->acut), DEF_MD_EPS);

    //initialization (for pso)
    _rset(parts);
    memcpy(parts[0].pb, parts[0].pc, N_DIM*sizeof(fp));
    parts[0].fb = _feval(parts[0].pb, mddata, lattice.sites, keatings);

    memcpy(pbest, parts[0].pc, N_DIM*sizeof(fp));
    fbest = parts[0].fb;

    for(i = 1u; i < N_PARTICLES; ++i){
        _rset(parts + i);
        memcpy(parts[i].pb, parts[i].pc, N_DIM*sizeof(fp));
        parts[i].fb = _feval(parts[i].pb, mddata, lattice.sites, keatings);

        if(parts[i].fb < fbest){
            fbest = parts[i].fb;
            memcpy(pbest, parts[i].pb, N_DIM*sizeof(fp));
        }
    }


    //main loop
    for(count = 0u, nslu = 0u, nkicks = -1u; count < N_ITER; ++count, ++nslu){

        //convergence info
        printf("\33[2K\r %u "FPFORMAT" (last nkick=%u)", count, fbest, nkicks);
        fflush(stdout);

        //standard pso loop
        for(i = 0u; i < N_PARTICLES; ++i){
            ///update velocities, positions
            for(j = 0u; j < N_DIM; ++j){
                parts[i].vc[j] = PARAM_W*parts[i].vc[j] +
                                 PARAM_P*((fp)rand()/(fp)RAND_MAX)*
                                         (parts[i].pb[j] - parts[i].pc[j]) +
                                 PARAM_G*((fp)rand()/(fp)RAND_MAX)*
                                         (pbest[j] - parts[i].pc[j]);
#if STRICT_BOUNDS
                if(parts[i].pc[j] < SS_BOUNDS[j][0]){
                    parts[i].vc[j] = ZERO;
                    parts[i].pc[j] = SS_BOUNDS[j][0];
                }
                else if(parts[i].pc[j] > SS_BOUNDS[j][1]){
                    parts[i].vc[j] = ZERO;
                    parts[i].pc[j] = SS_BOUNDS[j][1];
                }
                else{
                    parts[i].pc[j] += parts[i].vc[j];
                }
#else
                parts[i].pc[j] += parts[i].vc[j];
#endif
            }

            ///update local, global optima
            ftmp = _feval(parts[i].pc, mddata, lattice.sites, keatings);
            if(ftmp < parts[i].fb){
                memcpy(parts[i].pb, parts[i].pc, N_DIM*sizeof(fp));
                parts[i].fb = ftmp;
                if(ftmp + IMPROV_EPS < fbest) nslu = -(1u); ////kickout update
                if(ftmp < fbest){
                    memcpy(pbest, parts[i].pb, N_DIM*sizeof(fp));
                    fbest = ftmp;
                }
            }
        }

        //restart procedure
        if(nslu == N_B4_KICKOUT){
            for(i = 0u, nkicks = 0u; i < N_PARTICLES; ++i){
                if((fp)rand()/(fp)RAND_MAX > KICKOUT_PROB) continue;
                if(_dist2(parts[i].pc, pbest) > OPT_NBRHD_SQD) continue;
                _rset(parts + i);
                nkicks++;
            }
            nslu = -(1u);
        }
    }

    printf("\ndone. f(opt)= "FPFORMAT"\nopt=", fbest);
    for(i = 0u; i < N_DIM; ++i){
        printf(" "FPFORMAT, pbest[i]);
    } 

#if(!STACK_ALLOC)
    free(parts);
#endif
}






#endif