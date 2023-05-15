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

//simulation constants----------------------------------------------------------
#define N_DIM 50u               ///symmetry INEQUIV params to fit
#define SS_BOUNDS_LO (fp)-10.0 ///search space minima, maxima in eV/A^2
#define SS_BOUNDS_HI (fp)+10.0 ///constant across all dimensions
#define SS_BOUNDS_DIF ABSF(SS_BOUNDS_HI - SS_BOUNDS_LO) ///hope eval @ cmpl time

#define N_ITER 100000u
#define N_PARTICLES 75u

#define PARAM_W (fp)0.5 ///inertia, must be 0 < w < 1
#define PARAM_P (fp)1.6 ///individual strength, usually 1 < p < 3
#define PARAM_G (fp)2.5 ///global strangth, usually 1 < g < 3

///kickout params---------------------------------------------------------------
#define IMPROV_EPS (fp)1.0e-4 ///min change in global best before kickout reset
#define N_B4_KICKOUT 25u      ///iterations w/o improvement (^) before kickout 
#define OPT_NBRHD_SQD (fp)0.5 ///square distance from optima for kickout
#define KICKOUT_PROB (fp)0.5  ///how likley particle near optima is reset


//main object-------------------------------------------------------------------
typedef struct part{
    fp vc[N_DIM]; //velocity, current
    fp pc[N_DIM]; //position, current
    fp pb[N_DIM]; //position, best
    fp fb;        //f(pb)
} part;


//pso function definitions------------------------------------------------------
//testing function to minimize - Rastringin function, true min = 0  @  x = 0
static forceinline fp _eval(const fp x[restrict N_DIM]){
    uint i; fp res;

    res = (fp)10.0*(fp)N_DIM;
    for(i = 0u; i < N_DIM; ++i){
        res += x[i]*x[i] - (fp)10.0*COS(TWOPI*x[i]);
    }

    return res;
}

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
        p->vc[i] = SS_BOUNDS_DIF * ((fp)2.0*((fp)rand()/(fp)RAND_MAX) - ONE);
    for(i = 0u; i < N_DIM; ++i)
        p->pc[i] = SS_BOUNDS_LO + ((fp)rand()/(fp)RAND_MAX)*SS_BOUNDS_DIF;
}


//main function-----------------------------------------------------------------
//base pso algorithm from 
//Clerc, "Standard Particle Swarm Optimization" (2012) - HAL Open Access Archive
//kickout/restart procedure by
//Barone, "A way to Hopefully make this Work Better" (2023) - My Brain Just Now
void psomin(){
    uint i, j, count, nslu;
    fp pbest[N_DIM], fbest, ftmp;
    #if STACK_ALLOC
    part parts[N_PARTICLES];
    #else
    part*restrict parts; parts = malloc(N_PARTICLES*sizeof(part));
    #endif


    //initialization
    _rset(parts);
    memcpy(parts[0].pb, parts[0].pc, N_DIM*sizeof(fp));
    parts[0].fb = _eval(parts[0].pb);

    memcpy(pbest, parts[0].pc, N_DIM*sizeof(fp));
    fbest = parts[0].fb;

    for(i = 1u; i < N_PARTICLES; ++i){
        _rset(parts + i);
        memcpy(parts[i].pb, parts[i].pc, N_DIM*sizeof(fp));
        parts[i].fb = _eval(parts[i].pb);

        if(parts[i].fb < fbest){
            fbest = parts[i].fb;
            memcpy(pbest, parts[i].pb, N_DIM*sizeof(fp));
        }
    }


    //main loop
    for(count = 0u, nslu = 0u; count < N_ITER; ++count, ++nslu){

        //convergence info
        printf("\33[2K\r %u "FPFORMAT, count, fbest);

        //standard pso loop
        for(i = 0u; i < N_PARTICLES; ++i){
            ///update velocities, positions
            for(j = 0u; j < N_DIM; ++j){
                parts[i].vc[j] = PARAM_W*parts[i].vc[j] +
                                 PARAM_P*((fp)rand()/(fp)RAND_MAX)*
                                         (parts[i].pb[j] - parts[i].pc[j]) +
                                 PARAM_G*((fp)rand()/(fp)RAND_MAX)*
                                         (pbest[j] - parts[i].pc[j]);
                parts[i].pc[j] += parts[i].vc[j];
            }

            ///update local, global optima
            ftmp = _eval(parts[i].pc);
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
            for(i = 0u; i < N_PARTICLES; ++i){
                if((fp)rand()/(fp)RAND_MAX > KICKOUT_PROB) continue;
                if(_dist2(parts[i].pc, pbest) > OPT_NBRHD_SQD) continue;
                _rset(parts + i);
            }
            nslu = -(1u);
        }
    }

    printf("\ndone. f(opt)="FPFORMAT, fbest);

#if(!STACK_ALLOC)
    free(parts);
#endif
}





#endif