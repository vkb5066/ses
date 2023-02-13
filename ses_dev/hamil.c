//Implementation of hamil operator init / manip
#include <stdlib.h>
#include "def.h"

//Give sample max index: for a given gcut2, recip lattice vector len, and k pt
//element, returns the max sample index needed to describe V(G)
//e.x. if bi = b3, then giving ki = k along b3 returns l_max 
static forceinline ushort _gsmi(const fp gcut2, const fp bi2, const fp ki){
    fp x;
    ushort arg1, arg2;
    
    x = SQRT(gcut2/bi2);
    arg1 = (ushort)abs((int)(+x - ki));
    arg2 = (ushort)abs((int)(-x - ki));

    return MAX(arg1, arg2);
}

//Give the power of 2 that is >= the input
static forceinline ushort _npo2(const ushort n){
    ushort p;

    p = (ushort)1;
    if(n && !(n&(n-(ushort)1))) return n;
    while(p < n) p <<= (ushort)1;

    return p;
}

//Gives the power that 2 is raised to to yield the input (the input MUST be
//a pow of 2)
static forceinline uchar _lo2(ushort n){
    ushort p;

    p = (ushort)0;
    while(n >>= (ushort)1) p++;

    return (uchar)p;
}


/*
*  Returns max possible sizes that will be needed for pw basis
*  Updates: hklmt   (max abs sampling indices for the k.e. part)
*           hklsv   (sizes of the V(G) grid)
*           hkll2sv (log_2(hklsv))
*           sizet   (upper limit on # of G vectors in any given basis)
*  kpoint coords should be fractional along b1, b2, b3
*/
void setmaxdims(const fp gcut, const fp gcutvmul, const fp B[restrict 3][3],
                const ushort nkpts, const kpt*restrict kpts,                  
                ushort hklmt[restrict 3], ushort hklsv[restrict 3],
                uchar hkll2sv[restrict 3], uint*restrict sizet){
    ushort i;
    ushort tst;
    fp b1l2, b2l2, b3l2; 
    fp gcut2, gcut2ext;


    gcut2 = gcut*gcut; gcut2ext = gcut2*gcutvmul*gcutvmul;   
    //square lengths of the B lattice vectors
    b1l2 = B[0][0]*B[0][0] + B[0][1]*B[0][1] + B[0][2]*B[0][2];
    b2l2 = B[1][0]*B[1][0] + B[1][1]*B[1][1] + B[1][2]*B[1][2];
    b3l2 = B[2][0]*B[2][0] + B[2][1]*B[2][1] + B[2][2]*B[2][2];
    for(i = (ushort)0; i < (ushort)3; ++i){
        hklmt[i] = (ushort)0;
        hklsv[i] = (ushort)0;
    }

    //scan over all kpoints and pick out the largest h, k, l
    for(i = (ushort)0; i < nkpts; ++i){
        ///kinetic energy: can work with normal gcut2.  test h->k->l
        tst = (short)_gsmi(gcut2, b1l2, kpts[i].crds[0]);
        hklmt[0] = MAX(tst, hklmt[0]);
        tst = (short)_gsmi(gcut2, b2l2, kpts[i].crds[1]);
        hklmt[1] = MAX(tst, hklmt[1]);
        tst = (short)_gsmi(gcut2, b3l2, kpts[i].crds[2]);
        hklmt[2] = MAX(tst, hklmt[2]);
        
        ///potential energy: needs the (possibly modified) gcut2.  test h->k->l
        tst = (ushort)2*_gsmi(gcut2ext, b1l2, kpts[i].crds[0]) + (ushort)1;
        hklsv[0] = MAX(tst, hklsv[0]);
        tst = (ushort)2*_gsmi(gcut2ext, b2l2, kpts[i].crds[1]) + (ushort)1;
        hklsv[1] = MAX(tst, hklsv[1]);
        tst = (ushort)2*_gsmi(gcut2ext, b3l2, kpts[i].crds[2]) + (ushort)1;
        hklsv[2] = MAX(tst, hklsv[2]);
    }
    ///bring to the nearest power of 2 so that fft algos work.  set log_2(...)
    for(i = (ushort)0; i < (ushort)3; ++i){
        hklsv[i] = _npo2(hklsv[i]);
        hklsv[i] = MAX((ushort)2, hklsv[i]);
        hkll2sv[i] = _lo2(hklsv[i]);
    }

    //figure out the size of the max # of pws that will ever be needed
    *sizet = ((ushort)2*hklmt[0] + (ushort)1) *
             ((ushort)2*hklmt[1] + (ushort)1) *
             ((ushort)2*hklmt[2] + (ushort)1);
#if(!SAFE_SIZE_ALLOC)
    ///we can save a bit of space by taking into account that the parallelpiped 
    ///volume > the energy cutoff sphere's volume
    ///magic number is 4/3*pi / 8 ... the 8 comes from doubling the B vec lens
    *sizet = (uint)((fp)(*sizet) * ((fp)0.523598776*gcut2*gcut / 
                                    (hklmt[0]*SQRT(b1l2)*hklmt[1]*SQRT(b2l2)*
                                     hklmt[2]*SQRT(b3l2))
                                   )
                   );
#endif


}

/*
*  Sets the k-dependent kinetic energy array of the hamiltonian.  Also sets the
*  miller indices for the basis, i.e. the set (h, k, l) that obey |k+G| < Gcut.
*/
void sethamkin(const fp gcut2, const ushort absmillmax[restrict 3],
               const fp B[restrict 3][3], const fp meff,
               hamil*restrict ham){
    short h, k, l, hm, km, lm;
    uint npw, ind;
    fp tmul;
    fp hpk, kpk, lpk;
    fp kpc[3], kpg[3], kpg2;

    //Modifies ke operator with the effective mass, expected in units of m0
    tmul = HBARSQD_OVER_TWOM / meff;

    //Check each possible h, k, l and add to basis if it has a low enough ke
    //note that ~half of the G vecs are included, so doing this some fancy way
    //e.g. a KD tree would actually be worse
    hm = (short)absmillmax[0]; 
    km = (short)absmillmax[1]; 
    lm = (short)absmillmax[2];
    npw = 0u;
    ind = 0u;
    kpc[0] = ham->kp->crds[0];
    kpc[1] = ham->kp->crds[1];
    kpc[2] = ham->kp->crds[2];
    for(h = -hm; h <= hm; ++h){
    hpk = (fp)h + kpc[0];
    for(k = -km; k <= km; ++k){
    kpk = (fp)k + kpc[1];
    for(l = -lm; l <= lm; ++l){
    lpk = (fp)l + kpc[2];

        ///create the vector k + G
        ///k+G = kpt + (h*b1 + k*b2 + l*b3) ... if kpt is in direct coords,
        ///then its 3 coords in fractional units of b1, b2, b3.  Then k+G 
        ///becomes k+G = (h+kpt[0])*b1 + (k+kpt[1])*b2 + (l+kpt[2])*b3
        kpg[0] = hpk*B[0][0] + kpk*B[1][0] + lpk*B[2][0];
        kpg[1] = hpk*B[0][1] + kpk*B[1][1] + lpk*B[2][1];
        kpg[2] = hpk*B[0][2] + kpk*B[1][2] + lpk*B[2][2];
        
        ///test for |k+G| < Gcut.  If it is true, add to basis
        kpg2 = kpg[0]*kpg[0] + kpg[1]*kpg[1] + kpg[2]*kpg[2];
        if(kpg2 < gcut2){
            //ind = ACC2(3u, size, 0u);
            ham->mills[ind]      = h;
            ham->mills[ind + 1u] = k;
            ham->mills[ind + 2u] = l;
            ham->ke[npw] = tmul*kpg2;
            
            npw++;
            ind += 3u;
        }

    }
    }
    }

    ham->npw = npw;
}


/*
*  Sets the local potential V(G) on the reciprocal space grid
*  MIND: replace lattice->A with B before using, set atom cartesian crds,
*        and make sure vloc is memset to zero
*/
void sethamvloc(const lat lattice, const pot*restrict pots,
                hamil*restrict ham){
#define B lattice.A
    ushort i, j;
    short h, k, l;
    ushort hh, hk, hl, hi, ki, li;
    fp* invcounts = malloc((lattice.nspecs)*sizeof(fp)); 
    fp* invdqs = malloc((lattice.nspecs)*sizeof(fp));    
    fp g[3], g2;
    ushort acc; uint ind;
    cfp s; fp va, gdt;
    fp sam; ushort loi, hii;

    //some useful quantities
    for(i = (ushort)0; i < (ushort)lattice.nspecs; ++i){
        invcounts[i] = ONE / (fp)lattice.speccounts[i];
        invdqs[i] = (fp)(pots[i].nsamples - (ushort)1) / 
                    (pots[i].q2hi - pots[i].q2lo); 
    }
    hh = ham->dims[0] >> (ushort)1; ///you'd better hope that dims are still
    hk = ham->dims[1] >> (ushort)1; ///powers of two or this will break 
    hl = ham->dims[2] >> (ushort)1; ///everything

    //now actually set V(G)
    for(hi = (ushort)0; hi < ham->dims[0]; ++hi){
    h = (short)(hi - hh);
    for(ki = (ushort)0; ki < ham->dims[1]; ++ki){
    k = (short)(ki - hk);
    for(li = (ushort)0; li < ham->dims[2]; ++li){
    l = (short)(li - hl);
        
        ///g = h*b1 + k*b2 + l*b3
        g[0] = (fp)h*B[0][0] + (fp)k*B[1][0] + (fp)l*B[2][0];
        g[1] = (fp)h*B[0][1] + (fp)k*B[1][1] + (fp)l*B[2][1];
        g[2] = (fp)h*B[0][2] + (fp)k*B[1][2] + (fp)l*B[2][2];
        g2 = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];

        ///V(G) = S(G) * v_atom(|G|)
        ind = ACC3((uint)ham->dims[1], (uint)ham->dims[2], 
                   (uint)hi, (uint)ki, (uint)li);
        ///both S and v_atom need to be calculated species-by-species
        for(i = (ushort)0, acc=(ushort)0; i < (ushort)(lattice.nspecs); ++i){
            ///compute v_atom(G) w/ linear interpolation of pot samples
            sam = (g2 - pots[i].q2lo)*invdqs[i];
            loi = (ushort)sam;
            hii = loi + (ushort)1;
            if(loi >= hii || hii > pots[i].nsamples){ ///tst 0 catches underflow
                acc += lattice.speccounts[i];
                continue; ////no need to compute S if v is zero ... move on   
            }
            sam -= (fp)loi;
            va = pots[i].samples[loi]*(ONE-sam) + pots[i].samples[hii]*(sam);

            ///compute S(G) = 1/natoms * sum_natoms e^(-i G dot t), t=atom crds
            s.re = ZERO; s.im = ZERO;
            for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc){
                gdt = g[0]*lattice.sites[acc].crdsc[0] +
                      g[1]*lattice.sites[acc].crdsc[1] +
                      g[2]*lattice.sites[acc].crdsc[2];
                s.re += COS(gdt); s.im -= SIN(gdt);
            }
            s.re *= invcounts[i]; s.im *= invcounts[i];

            ///finally, we can add to the element-wise sum V(G)
            ham->vloc[ind].re += s.re*va; ham->vloc[ind].im += s.im*va;
        }
    }
    }   
    }

    //set V(G=0) = 0 for a consistent rigid band shift ... may want to change
    //this later to match exp. work functions, but 0 is fine for now
    //whatever you do, DONT leave G = 0 untreated unless you want numeric issues
    ind = ACC3((uint)ham->dims[1], (uint)ham->dims[2], 
               (uint)hh, (uint)hk, (uint)hl); 
    ham->vloc[ind].re = ZERO; ham->vloc[ind].im = ZERO;


    free(invcounts);
    free(invdqs);
#undef B
}


/*
*  Builds an explicit hamiltonian (i.e. the full basis! careful!) and fills the 
*  pre-allocated array H with its values
*/
void buildexphamil(const hamil haminfo, cfp*restrict H){
    uint i, j, i3, j3, indh, indv;
    short ht, kt, lt, shh, shk, shl;

    shh = (short)(haminfo.dims[0] >> (ushort)1);
    shk = (short)(haminfo.dims[1] >> (ushort)1);
    shl = (short)(haminfo.dims[2] >> (ushort)1);

    //H is hermitian, and it might be tempting to write the loop as:
    //0 < i < npw -> i+1 < j < npw, then set H[i][j] and H[j][i] at once ...
    // ... but this would ruin cache performance.  The 'dumb' way is actually
    //better for most sizes of H that pop up in physics
    indh = 0u;
    for(i = 0u, i3 = 0u; i < haminfo.npw; ++i, i3 += 3u){

        ///lower portion
        for(j = 0u, j3 = 0u; j < i; ++j, ++indh, j3 += 3u){
#if(SAFE_EXPL_HAM)
            ///make sure we're not about to oob V.  we need both checks since
            ///this is hermitian
            ht = haminfo.mills[i3]    - haminfo.mills[j3];
            if((ht <= -shh) || (shh <= ht)) goto panic0;
            kt = haminfo.mills[i3+1u] - haminfo.mills[j3+1u];
            if((kt <= -shk) || (shk <= kt)) goto panic0;
            lt = haminfo.mills[i3+2u] - haminfo.mills[j3+2u];
            if((lt <= -shl) || (shl <= lt)) goto panic0;
            ///grab the index of mills_i - mills_j in the vloc array
            ht += shh; kt += shk; lt += shl;
#else
            ///grab the index of mills_i - mills_j in the vloc array
            ht = haminfo.mills[i3]    - haminfo.mills[j3]    + shh;
            kt = haminfo.mills[i3+1u] - haminfo.mills[j3+1u] + shk;
            lt = haminfo.mills[i3+2u] - haminfo.mills[j3+2u] + shl;
#endif

            ///potential energy term: sample vloc at gi - gj
            indv = ACC3((uint)haminfo.dims[1], (uint)haminfo.dims[2], 
                        (uint)ht, (uint)kt, (uint)lt);
            H[indh].re = haminfo.vloc[indv].re;
            H[indh].im = haminfo.vloc[indv].im;
#if(SAFE_EXPL_HAM)
            continue;

            panic0:
            H[indh].re = ZERO;
            H[indh].im = ZERO;
#endif          
        }

        ///kinetic energy term: this is slightly easier than V
        H[indh].re = haminfo.ke[i]; 
        H[indh].im = ZERO; ///dont fucking delete this
        indh++;

        ///upper portion
        for(j = i+1u, j3 = 3u*j; j < haminfo.npw; ++j, ++indh, j3 += 3u){
#if(SAFE_EXPL_HAM)
            ///make sure we're not about to oob V.  we need both checks since
            ///this is hermitian
            ht = haminfo.mills[i3]    - haminfo.mills[j3];
            if((ht <= -shh) || (shh <= ht)) goto panic1;
            kt = haminfo.mills[i3+1u] - haminfo.mills[j3+1u];
            if((kt <= -shk) || (shk <= kt)) goto panic1;
            lt = haminfo.mills[i3+2u] - haminfo.mills[j3+2u];
            if((lt <= -shl) || (shl <= lt)) goto panic1;
            ///grab the index of mills_i - mills_j in the vloc array
            ht += shh; kt += shk; lt += shl;
#else
            ///grab the index of mills_i - mills_j in the vloc array
            ht = haminfo.mills[i3]    - haminfo.mills[j3]    + shh;
            kt = haminfo.mills[i3+1u] - haminfo.mills[j3+1u] + shk;
            lt = haminfo.mills[i3+2u] - haminfo.mills[j3+2u] + shl;
#endif

            ///potential energy term: sample vloc at gi - gj
            indv = ACC3((uint)haminfo.dims[1], (uint)haminfo.dims[2], 
                        (uint)ht, (uint)kt, (uint)lt);
            H[indh].re = haminfo.vloc[indv].re;
            H[indh].im = haminfo.vloc[indv].im;
#if(SAFE_EXPL_HAM)
            continue;

            panic1:
            H[indh].re = ZERO;
            H[indh].im = ZERO;
#endif             
        }
    }
/*    int co = 0;
    for(uint v = 0; v < haminfo.npw; ++v){
        for(uint b = 0; b < haminfo.npw; ++b){
            fp q = (*H)[ACC2(haminfo.npw, v, b)].re;
            fp w = (*H)[ACC2(haminfo.npw, v, b)].im;
            fp e = (*H)[ACC2(haminfo.npw, b, v)].re;
            fp r = (*H)[ACC2(haminfo.npw, b, v)].im;
            if(ABSF(q - e) > 1.0e-5 || ABSF(w + r) > 1.0e-5){
                //printf("ruh roh raggy ... %i %i\n", v, b);
                //printf(FPFORMAT" "FPFORMAT" vs ... "FPFORMAT" "FPFORMAT"\n",
                //        q, w, e, r);
                co++;
            }
        }
    }
    printf("------->>>%i\n", co);*/
}

