#include <stdlib.h>
#include <stdio.h>
#include <alloca.h>
#include <string.h>

#include "def.h"
#include "kdtree-master/kdtree.h"


#define DOT(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

//change of basis x = Ay, stores result in x
//A should be stored with its mathamatical columns as rows
static forceinline void _cob(const fp A[restrict 3][3],
                             fp x[restrict 3], const fp y[restrict 3]){
    x[0] = A[0][0]*y[0] + A[1][0]*y[1] + A[2][0]*y[2];
    x[1] = A[0][1]*y[0] + A[1][1]*y[1] + A[2][1]*y[2];
    x[2] = A[0][2]*y[0] + A[1][2]*y[1] + A[2][2]*y[2];
}

//square distance between a 3d vector v and a point x, y, z
static forceinline fp _dist2vxyz(const fp v[restrict 3], 
                                const fp x, const fp y, const fp z){
    fp ba, res;
    
    ba = v[0] - x; res  = ba*ba;
    ba = v[1] - y; res += ba*ba;
    ba = v[2] - z;
    return res + ba*ba;
}



/*
   !!! CONTENT WARNING
   severe abuse of preprocessor macros below.  reader discretion is advised.
   CONTENT WARNING !!!
*/
//set pairs, trips, periodic images for a set of base atoms
//allocates extra atoms (the periodic images) -> xatoms
//rcut = cutoff dist between central atom and its pair (angstroms)
//acut = cutoff angle for a triplet, in radians [0, pi]
#define IS_IMAGE           (ushort)-1
#define IS_IMPORTANT_IMAGE (ushort)-2
void setnbrs(const lat lattice, 
             ushort*restrict nxatoms, site*restrict*restrict xatoms,
             const fp rcut, const fp acut, const fp eps){
    uint i, j, k, ii;
    int a, b, c;

    uint nimgstot; int nimgs[3];
    fp latlen;

#define atoms lattice.sites
    uint natoms, nallatoms, reset, caa, cta;
    site*restrict allatoms; ///!dyn mem
#if(!SAFE_NBR_ALLOC)
    site* reallocchecker;
#endif

    uint* indexmap; ///!dyn mem
    void* kdnode; ///!dyn mem
    int ressize;
    struct kdres* res; ///!dyn mem

    fp x, y, z, cosacut;
    fp ij[3], nij, ik[3];


    //find num of periodic images needed to satisfy rcut
    nimgstot = 1u;
    for(i = 0u; i < 3u; ++i){
        latlen = SQRT(DOT(lattice.A[i], lattice.A[i]));
        for(j = 1u;; ++j){
            if(rcut < (fp)j*latlen){
                nimgs[i] = (int)j;
                break;
            }
        }
        nimgstot *= (2u*(uint)nimgs[i]+1u);
    }

    //periodic image loop
    natoms = 0u; 
    for(i = 0u; i < (uint)lattice.nspecs; ++i) 
        natoms += (uint)lattice.speccounts[i];
    nallatoms = nimgstot*natoms;
    allatoms = malloc(nallatoms*sizeof(site));
    for(i = 0u, caa = 0u; i < natoms; ++i){
        ///ensure that this atom is in the cell
        reset = 0u;
        for(j = 0u; j < 3u; ++j){
            while(atoms[i].crdsf[j] < ZERO){
                atoms[i].crdsf[j] += ONE;
                reset = 1u;
            }
            while(atoms[i].crdsf[j] >= ONE){
                atoms[i].crdsf[j] -= ONE;
                reset = 1u;
            }
        }
        if(reset) _cob(lattice.A, atoms[i].crdsc, atoms[i].crdsf); 

        ///deal with a=b=c=0: this is an atom inside the cell and we must keep
        ///track of its periodic images
        ////deepcopy important info from 'atoms' ... ii indexes this atom
        ii = (nimgstot)*(i+1u) - 1u;
        allatoms[ii].spec = atoms[i].spec;
        memcpy(allatoms[ii].crdsf, atoms[i].crdsf, 3u*sizeof(fp));
        memcpy(allatoms[ii].crdsc, atoms[i].crdsc, 3u*sizeof(fp));
        allatoms[ii].nimgs = (ushort)(nimgstot - 1u); ////minus (0,0,0) "image"
        allatoms[ii].imgs = malloc((allatoms[ii].nimgs)*sizeof(site*));
        allatoms[ii].self = &atoms[i];

        cta = 0u;
        for(a = -nimgs[0]; a <= nimgs[0]; ++a)
        for(b = -nimgs[1]; b <= nimgs[1]; ++b)
        for(c = -nimgs[2]; c <= nimgs[2]; ++c){
            if(!a && !b && !c) continue;

            ////deepcopy important info from 'atoms' (+ coord shift) but keep
            ////a reference to the original species (we want to independently 
            ////change coords, but always keep species identical for atom 
            ////hopping later)
            allatoms[caa].spec = allatoms[ii].spec;
            allatoms[caa].crdsf[0] = atoms[i].crdsf[0] + (fp)a;
            allatoms[caa].crdsf[1] = atoms[i].crdsf[1] + (fp)b;
            allatoms[caa].crdsf[2] = atoms[i].crdsf[2] + (fp)c;
            _cob(lattice.A, allatoms[caa].crdsc, allatoms[caa].crdsf);
            allatoms[caa].nimgs = IS_IMAGE; ///use this entry as a flag

            ////add this shifted atom to the periodic images
            allatoms[ii].imgs[cta] = &allatoms[caa];
            cta++;
            caa++;
        }
    
        caa++;
    }

    //fill in a kd tree for quick nearest neighbor searches
    indexmap = malloc(nallatoms*sizeof(uint));
    for(i = 0u; i < nallatoms; ++i) indexmap[i] = i;
    kdnode = kd_create(3);
    for(i = 0u; i < nallatoms; ++i){
        kd_insert3(kdnode, allatoms[i].crdsc[0], allatoms[i].crdsc[1], 
                   allatoms[i].crdsc[2], (void*)(indexmap + i));
    }

    //neighbor setting loop ... this is going to get messy
    //looping indices guarantee that we only query atoms in the original cell
    //at the same time, we also reduce the num images of each base atom from
    //max possible -> minimum number
    *nxatoms = 0u; 
    *xatoms = malloc(nallatoms*sizeof(site));
    cosacut = COS(acut);
    for(i = nimgstot-1u, ii = 0u; i < nallatoms; i += nimgstot, ++ii){
        ///query the kd tree
        x = allatoms[i].crdsc[0];
        y = allatoms[i].crdsc[1];
        z = allatoms[i].crdsc[2];
        res = kd_nearest_range3(kdnode, x, y, z, rcut);
        ressize = kd_res_size(res);

        ///--- set pairs -------------------------------------------------------
        allatoms[i].npairs = (uint)ressize - 1u;
        atoms[ii].npairs = allatoms[i].npairs;
        allatoms[i].pairs = malloc(allatoms[i].npairs*sizeof(site*));
        atoms[ii].pairs = malloc(atoms[ii].npairs*sizeof(site*));

        for(j = 0u, k = 0u; j < (uint)ressize; ++j){
#if   PREC == 1
            cta = *(  (uint*)kd_res_item3f(res, &x, &y, &z)  );
#elif PREC == 2
            cta = *(  (uint*)kd_res_item3(res, &x, &y, &z)   );
#endif
            
            if(_dist2vxyz(allatoms[i].crdsc, x, y, z) < eps) goto skip;
            ///yet another flag - this one shows that this atom is either a base
            ///atom or the neighbor of one.  only set if it isnt already a 
            ///reasonable value
            if(allatoms[cta].nimgs == IS_IMAGE){
                ////ok, this is a unique periodic image that we want to keep
                ///this will go into the extra atoms array - deepcopy it there
                allatoms[cta].nimgs = IS_IMPORTANT_IMAGE;
                (*xatoms)[*nxatoms].spec = allatoms[cta].spec;
                memcpy((*xatoms)[*nxatoms].crdsf, allatoms[cta].crdsf, 
                        3u*sizeof(fp));
                memcpy((*xatoms)[*nxatoms].crdsc, allatoms[cta].crdsc, 
                        3u*sizeof(fp));
                ////connect kd tree's self refs to this extra atom's ref
                (*xatoms)[*nxatoms].self = &(*xatoms)[*nxatoms];
                allatoms[cta].self = &(*xatoms)[*nxatoms];
                (*nxatoms)++;
            }
            allatoms[i].pairs[k] = &(allatoms[cta]);
            atoms[ii].pairs[k] = allatoms[cta].self;

            k++;
            skip:
            kd_res_next(res);
        }
        kd_res_free(res);

        ///--- set trips -------------------------------------------------------
#define ai lattice.sites[ii]
#define aj lattice.sites[ii].pairs[j]
#define ak lattice.sites[ii].pairs[k]
        ///ntrips is an array w/ jth index = number of trips of pair ij
        ai.ntrips = calloc(ai.npairs, sizeof(uchar));
        ai.trips = malloc(ai.npairs*sizeof(site**));

        for(j = 0u; j < ai.npairs; ++j){
            ij[0] = aj->crdsc[0] - ai.crdsc[0];
            ij[1] = aj->crdsc[1] - ai.crdsc[1];
            ij[2] = aj->crdsc[2] - ai.crdsc[2];
            nij = SQRT(DOT(ij, ij));

            ////here, we begin filling trips of j
            ai.trips[j] = malloc(ai.npairs*sizeof(site*));
            for(k = 0u; k < ai.npairs; ++k){
                if(j == k) continue;
                ik[0] = ak->crdsc[0] - ai.crdsc[0];
                ik[1] = ak->crdsc[1] - ai.crdsc[1];
                ik[2] = ak->crdsc[2] - ai.crdsc[2];

                ////trip of j is made when j-i, k-i < rcut and angle jik < acut
                ////for j != k
                ////we've already dont the rcut check, now just do the acut 
                ////check by cos(theta) > cos(theta_cut) (not a typo!) as long 
                //// 0 <= theta <= pi
                //// cos(angle jik) = ( (j-i)@(k-i) )  /  ( |j-i||k-i| )
                if(cosacut < DOT(ij, ik)/(nij*SQRT(DOT(ik, ik)))){
                    ai.trips[j][ai.ntrips[j]] = ak; /////ak is already a ptr
                    ai.ntrips[j]++;
                }
            }

            ///we've allocated more memory than needed ... get some back
            ai.trips[j] = realloc(ai.trips[j], ai.ntrips[j]*sizeof(site*));
        }
#undef ai
#undef aj
#undef ak
    }

    //recuce, set nimgs
    for(i = nimgstot-1u, ii = 0u; i < nallatoms; i += nimgstot, ++ii){
        atoms[ii].nimgs = (ushort)0;
        atoms[ii].imgs = malloc((nimgstot-1u)*sizeof(site*));

        ///set imgs to keep by checking vs. the "throw me away" flag
        for(j = 0u; j < allatoms[i].nimgs; ++j){
            if(allatoms[i].imgs[j]->nimgs == IS_IMAGE) continue;
            ////if we're here, we want to keep this one
            atoms[ii].imgs[atoms[ii].nimgs] = allatoms[i].imgs[j]->self;
            atoms[ii].nimgs++;
        }
        ///yet again, we've allocated too much memory ... get some back
        atoms[ii].imgs = realloc(atoms[ii].imgs, atoms[ii].nimgs*sizeof(site*));
    }


    //finish up
#if(!SAFE_NBR_ALLOC)
    ///the call to realloc(), although the new size is always <= current size,
    ///technically may decide to "copy" *xatoms to a new memory address.
    ///this is an issue since the lattice sites 'img', 'pair', and 'trip' 
    ///members will then point to the old memory addresses, and altering them
    ///with the expectation of consistent changes will fuck. you. up.
    ///any self-respecting implementation of realloc() won't cause this problem
    ///if you're shrinking the array, but there is no guarentee.
    ///this is also useful to comment out because valgrind freaks out even
    ///if the memory doesn't end up getting moved.  
    reallocchecker = *xatoms;
    *xatoms = realloc(*xatoms, (*nxatoms)*sizeof(site)); 
    if(reallocchecker != *xatoms){
        puts("!!! WARNING !!!\n"
             "setnbrs(): realloc() moved *xatoms.  any routine that alters\n"
             "lattice site images, pairs, or triplets will result in errors.\n"
             "you may comment out the *xatoms = realloc(...) line to fix this\n"
             "error, at the expense of more memory use.  Don't blame me,\n"
             "blame whoever wrote your realloc() implementation.\n"
             "!!! WARNING !!!");
    }
#endif
    for(i = nimgstot-1u; i < nallatoms; i += nimgstot){  
        free(allatoms[i].imgs);
        free(allatoms[i].pairs);
    }
    free(allatoms);
    kd_free(kdnode);
    free(indexmap);
}
#undef atoms
#undef IS_IMAGE
#undef IS_IMPORTANT_IMAGE


//calculates net force vector on site i (si) -> res (eV/angstrm)
//also returns the square of the net force
//MIND: res must be set to zero before being sent to this function
static forceinline fp _setforce(const uchar nspecs, 
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

    return DOT(res, res);
}

//shifts an atom and its periodic images at time t to time t + dt
//via x(t+dt) = x(t) + v(t)dt + 1/2a(t)dt^2
static forceinline void _shift(const fp vtdt[restrict 3], 
                               const fp atdt2[restrict 3], 
                               site*restrict s){
    ushort i;
    fp dx[3];

    dx[0] = vtdt[0] + atdt2[0];
    dx[1] = vtdt[1] + atdt2[1];
    dx[2] = vtdt[2] + atdt2[2];

    //update main coords, then images
    s->crdsc[0] += dx[0];
    s->crdsc[1] += dx[1];
    s->crdsc[2] += dx[2];
    for(i = (ushort)0; i < s->nimgs; ++i){
        s->imgs[i]->crdsc[0] += dx[0];
        s->imgs[i]->crdsc[1] += dx[1];
        s->imgs[i]->crdsc[2] += dx[2];
    }
}


//prints out convergence progress, returns true if converged
static forceinline uchar _mdupdateconv(const uint step, const uint nsteps,
                                       const fp ctime, const fp maxf, 
                                       const fp bestf, const fp ftol){

    fputs("\33[2K\r", stdout); ///clears entire line (S.O. #1508490)
    printf(" (%u / %u)  t = "FPFORMAT" fs  max[f(t)] = "FPFORMAT" eV/A"
           "  min[f(t-dt)] = "FPFORMAT" eV/A",
           step, nsteps, ctime, maxf, bestf);
    fflush(stdout);    

    return (uchar)(maxf < ftol);
}

//replaces atom's crdsf with the crdsc that minimize the energy if
//current force < min force up to this point.  also sets min force
static forceinline void _mdupdatebests(const fp thsfrc, fp*restrict minfrc,
                                       const lat latt){
    uchar i; ushort j, k, acc;

    if(thsfrc < *minfrc) *minfrc = thsfrc;
    else return; ///no need to go through this if this force > minimum so far

    for(i = (uchar)0, acc = 0u; i < latt.nspecs; ++i){
    for(j = (ushort)0; j < latt.speccounts[i]; ++j, ++acc){
        memcpy(latt.sites[acc].crdsf, latt.sites[acc].crdsc, 3u*sizeof(fp));
        for(k = (ushort)0; k < latt.sites[acc].nimgs; ++k){
            memcpy(latt.sites[acc].imgs[k]->crdsf, 
                   latt.sites[acc].imgs[k]->crdsc, 3u*sizeof(fp));
        }
    }
    }
}

//quenches the velocity v(t) for structural minimization
//behavior of the function is changed in def.h
static forceinline void _mdquench(fp vtdt[restrict 3], 
                                  const fp atdt2[restrict 3]){
#if   MD_QUENCH_STYLE == 0 ///no quenching
    
#elif MD_QUENCH_STYLE == 1 ///vasp style
    //v(t)dt proj along a(t) unless v(t) antiparallel to a(t), where v(t) = 0
    fp fac, vtl, atli;
    
    vtl = SQRT(DOT(vtdt, vtdt)) + SMALL;
    atli = ONE / (SQRT(DOT(atdt2, atdt2)) + SMALL);

    //antiparallel check
    fac = DOT(vtdt, atdt2)/vtl*atli; /// = v(t)dt@a(t)dt^2 / |v(t)dt||a(t)dt^2|
    if(ABSF(fac + ONE) < DEF_MD_VASP_QUENCH_EPS){ ///not a mistake: cmp fac, -1
        memset(vtdt, 0, 3u*sizeof(fp));
        return;
    }
    //project v(t)dt along a(t)
    fac = fac*vtl*atli; /// = v(t)dt@a(t)dt^2 / |a(t)dt^2|^2
    vtdt[0] = fac*atdt2[0];
    vtdt[1] = fac*atdt2[1];
    vtdt[2] = fac*atdt2[2];

#elif MD_QUENCH_STYLE == 2 ///Mattoni style
    //v_i(t) = 0 if v_i points in the opposite direction of a_i (i = x,y,z)
    if(vtdt[0]*atdt2[0] < ZERO) vtdt[0] = ZERO;
    if(vtdt[1]*atdt2[1] < ZERO) vtdt[1] = ZERO;
    if(vtdt[2]*atdt2[2] < ZERO) vtdt[2] = ZERO;

#endif

    return;
}

//atomic relaxation via velocity verlet integration w/ keating potential
//the following must be allocated prior to use:
// ---workspace---
//  (9*natom) * sizeof(fp)    b for RWS
//note that the cartesian coordinates are updated at the end of this routine, 
//and the direct coordinates are filled with GARBAGE VALUES
//this looks more complicated than necessary since x(t), v(t) and a(t) have
//VERY different dimensionalities ... therefore, I instead work with
//x(t), v(t)dt, and a(t)dt^2 ... all in angstroms to avoid numerical issues
uchar relax(fp*restrict RWS, const lat lattice, const pot*restrict pots, 
            const keat*restrict keatings, const fp ts,
            const fp maxforcetol, const uint itrlim){

    uchar i; ushort j, k;
    uint cnt, acc, accm;
    fp ctime, tmpfrc, maxfrc, bstfrc;
    fp*restrict imasses = alloca(lattice.nspecs*sizeof(fp)); ///don't free 
    fp dt2;
    fp*restrict vtdt, *restrict acctdt2, *restrict acctpdtdt2;


    //workspace setup: all three arrays are natoms x 3 arrays  
    //v(t)*dt, a(t), and a(t+dt)
    for(i = (uchar)0, acc = 0u; i < lattice.nspecs; ++i)
    for(j = (ushort)0; j < lattice.speccounts[i]; ++j) acc++;    
    acc = 3u*acc;
    vtdt = RWS + 0; acctdt2 = vtdt + acc; acctpdtdt2 = acctdt2 + acc;

    //initialize: get inverse masses into useful units
    dt2 = ts*ts;
    for(i = (uchar)0; i < lattice.nspecs; ++i){
        imasses[i] = VERLET_ACC_CONV * ONE/pots[i].mass * dt2;
    }

    //initialize: find a(t=0) from initial positions, v(t=0) = 0
    maxfrc = ZERO;
    memset(acctdt2, 0, acc*sizeof(fp)); ///here, acc = 3*natoms
    for(i = (uchar)0, acc = 0u, accm = 0u; i < lattice.nspecs; ++i){
    for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc, accm += 3u){
        tmpfrc = _setforce(lattice.nspecs, lattice.sites[acc], keatings[i], 
                           acctdt2 + accm); ///F(t) -> acct[offset + 0, 1, 2]
        maxfrc = MAX(tmpfrc, maxfrc);
        acctdt2[accm + 0u] *= imasses[i]; ///F(t)/m = a(t)dt^2
        acctdt2[accm + 1u] *= imasses[i]; 
        acctdt2[accm + 2u] *= imasses[i];    
    }
    }

    if(_mdupdateconv(0u, itrlim, ZERO, maxfrc, maxfrc, maxforcetol)) 
        return (uchar)0;
    bstfrc = maxfrc + maxforcetol;
    _mdupdatebests(maxfrc, &bstfrc, lattice);
    memset(vtdt, 0, (accm)*sizeof(fp)); ///here, accm = 3 x # of atoms

    //main loop
    for(cnt = 1u, ctime = ts; cnt <= itrlim; ++cnt, ctime += ts){

        ///calculate x(t+dt)
        for(i = (uchar)0, acc = 0u, accm = 0u; i < lattice.nspecs; ++i){
        for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc, accm += 3u){
            _mdquench(vtdt + accm, acctdt2 + accm);
            _shift(vtdt + accm, acctdt2 + accm, &lattice.sites[acc]);
        }
        }   

        ///calculate a(t+dt)dt^2
        maxfrc = ZERO;
        memset(acctpdtdt2, 0, accm*sizeof(fp));
        for(i = (uchar)0, acc = 0u, accm = 0u; i < lattice.nspecs; ++i){
        for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc, accm += 3u){
            tmpfrc = _setforce(lattice.nspecs, lattice.sites[acc], 
                               keatings[i], acctpdtdt2 + accm); 
            maxfrc = MAX(tmpfrc, maxfrc);
            acctpdtdt2[accm + 0u] *= imasses[i];
            acctpdtdt2[accm + 1u] *= imasses[i];
            acctpdtdt2[accm + 2u] *= imasses[i];    
        }
        }

        ///check convergence, update best positions if needed
        if(_mdupdateconv(cnt, itrlim, ctime, maxfrc, bstfrc, maxforcetol)) 
            return (uchar)0;
        _mdupdatebests(maxfrc, &bstfrc, lattice);


        ///calculate v(t+dt)dt, store in v(t)dt 
        for(acc = 0u; acc < accm; ++acc){ ////accm still holds 3 x natoms
            vtdt[acc] += acctdt2[acc] + acctpdtdt2[acc];
        }

        ///replace a(t)dt^2 with a(t+dt)dt^2 to prep for next iteration
        memcpy(acctdt2, acctpdtdt2, accm*sizeof(fp)); ////accm = 3 x natoms
    }

    //if we're here, our md didnt reach the minimum force tol within the given
    //time :(
    //return the best coords back into crdsc
    for(i = (uchar)0, acc = 0u; i < lattice.nspecs; ++i){
    for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc){
        memcpy(lattice.sites[acc].crdsc, lattice.sites[acc].crdsf, 
               3u*sizeof(fp));
        for(k = (ushort)0; k < lattice.sites[acc].nimgs; ++k){
            memcpy(lattice.sites[acc].imgs[k]->crdsc, 
                   lattice.sites[acc].imgs[k]->crdsf, 3u*sizeof(fp));
        }
    }
    }
    return (uchar)1;
}






#undef DOT