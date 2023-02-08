//Implementation of convolutions based on fft
#include "def.h"

//sin, cos tables for trig rec. initialization.  Covers 1d ffts up to size 2^32  
//!NOTE! check if flushing denormals -> zero ruins this (it shouldn't for 
//double prec, but maybe does for single prec)
//!NOTE2! we may want to set this as a lfp instead of fp ... the fft is v. 
//sensitive to numerical error
#if(PREC == 1)
static const fp CT[32] = {
0.0000000e+00F, 0.0000000e+00F, 6.1232340e-17F, 7.0710678e-01F, 9.2387953e-01F,
9.8078528e-01F, 9.9518473e-01F, 9.9879546e-01F, 9.9969882e-01F, 9.9992470e-01F,
9.9998118e-01F, 9.9999529e-01F, 9.9999882e-01F, 9.9999971e-01F, 9.9999993e-01F,
9.9999998e-01F, 1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F,
1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F,
1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F, 1.0000000e+00F,
1.0000000e+00F, 1.0000000e+00F
};
static const fp ST[32] = {
0.0000000e+00F, 0.0000000e+00F, 1.0000000e+00F, 7.0710678e-01F, 3.8268343e-01F,
1.9509032e-01F, 9.8017140e-02F, 4.9067674e-02F, 2.4541229e-02F, 1.2271538e-02F,
6.1358846e-03F, 3.0679568e-03F, 1.5339802e-03F, 7.6699032e-04F, 3.8349519e-04F,
1.9174760e-04F, 9.5873799e-05F, 4.7936900e-05F, 2.3968450e-05F, 1.1984225e-05F,
5.9921125e-06F, 2.9960562e-06F, 1.4980281e-06F, 7.4901406e-07F, 3.7450703e-07F,
1.8725351e-07F, 9.3626757e-08F, 4.6813379e-08F, 2.3406689e-08F, 1.1703345e-08F,
5.8516723e-09F, 2.9258362e-09F
};
#elif(PREC == 2)
static const fp CT[32] = {
0.0000000000000000e+00, 0.0000000000000000e+00, 6.1232339957367660e-17,
7.0710678118654757e-01, 9.2387953251128674e-01, 9.8078528040323043e-01,
9.9518472667219693e-01, 9.9879545620517241e-01, 9.9969881869620425e-01,
9.9992470183914450e-01, 9.9998117528260111e-01, 9.9999529380957619e-01,
9.9999882345170188e-01, 9.9999970586288223e-01, 9.9999992646571789e-01,
9.9999998161642933e-01, 9.9999999540410733e-01, 9.9999999885102686e-01,
9.9999999971275666e-01, 9.9999999992818922e-01, 9.9999999998204725e-01,
9.9999999999551181e-01, 9.9999999999887801e-01, 9.9999999999971945e-01,
9.9999999999992983e-01, 9.9999999999998246e-01, 9.9999999999999567e-01,
9.9999999999999889e-01, 9.9999999999999978e-01, 9.9999999999999989e-01,
1.0000000000000000e+00, 1.0000000000000000e+00
};
static const fp ST[32] = {
0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
7.0710678118654746e-01, 3.8268343236508978e-01, 1.9509032201612825e-01,
9.8017140329560604e-02, 4.9067674327418015e-02, 2.4541228522912288e-02,
1.2271538285719925e-02, 6.1358846491544753e-03, 3.0679567629659761e-03,
1.5339801862847655e-03, 7.6699031874270449e-04, 3.8349518757139556e-04,
1.9174759731070329e-04, 9.5873799095977345e-05, 4.7936899603066881e-05,
2.3968449808418219e-05, 1.1984224905069705e-05, 5.9921124526424275e-06,
2.9960562263346608e-06, 1.4980281131690111e-06, 7.4901405658471574e-07,
3.7450702829238413e-07, 1.8725351414619535e-07, 9.3626757073098084e-08,
4.6813378536549088e-08, 2.3406689268274551e-08, 1.1703344634137277e-08,
5.8516723170686385e-09, 2.9258361585343192e-09
};
#endif


//Applies a 1D fft to arr, whose (pow of 2) size is n and ln2 is log_n(n),
//Forwards or backwards fft is controlled by setting sgn to +/- 1
//This is the 'decimation in frequency' version with no bit reversal section
//Based on 'Matters Computational' - Arndt  (main routine similar to 21.2.1.2)
//         'Numerical Recipes in C' - Press (trig rec. portion)
static forceinline void _fft1difbase(const ushort n, const uchar l2n,
                                     cfp*restrict*restrict arr, const char sgn){
    ushort i, j, ij, ld, m, mh, ijmh;
    fp wr, wi, cn, sn, cntmp;
    cfp u, v, umv;

#define arr (*arr)
    //most of the transform
    for(ld = (ushort)l2n; ld >= (ushort)2; --ld){
        m = (ushort)1 << ld;
        mh = m >> (ushort)1;

        for(i = (ushort)0; i < n; i += m){
            ///some explicit 1+0i multiplications
            ijmh = i + mh;
            u.re = arr[i].re + arr[ijmh].re; u.im = arr[i].im + arr[ijmh].im;
            v.re = arr[i].re - arr[ijmh].re; v.im = arr[i].im - arr[ijmh].im;
            arr[i].re    = u.re; arr[i].im    = u.im;
            arr[ijmh].re = v.re; arr[ijmh].im = v.im;

            ///init trig rec
            wr = CT[ld]; 
            wi = (fp)sgn*ST[ld];
            cn = wr;    
            sn = wi;

            ///apply butterfly
            for(j = (ushort)1; j < mh; ++j){
                ij = i + j;
                ijmh = ij + mh;
                ////this part is u := arr[ij], v := arr[ijmh]
                ////followed by  arr[ij] := u+v, arr[ijmh] := (u-v)*(cn + isn)
                u.re = arr[ij].re;   u.im = arr[ij].im;
                v.re = arr[ijmh].re; v.im = arr[ijmh].im;
                umv.re = u.re - v.re; umv.im = u.im - v.im;
                arr[ij].re   = u.re + v.re; 
                arr[ij].im   = u.im + v.im;
                arr[ijmh].re = umv.re*cn - umv.im*sn;
                arr[ijmh].im = umv.re*sn + umv.im*cn;

                ////update trig rec
                cntmp = cn;
                cn = wr*cn    - wi*sn;
                sn = wi*cntmp + wr*sn;
            }
        }
    }

    //the rest of the explicit 1+0i multiplications
    for(i = (ushort)0; i < n; i += (ushort)2){
        ij = i + (ushort)1;
        u.re = arr[i].re + arr[ij].re; u.im = arr[i].im + arr[ij].im;
        v.re = arr[i].re - arr[ij].re; v.im = arr[i].im - arr[ij].im;
        arr[i].re  = u.re; arr[i].im  = u.im;
        arr[ij].re = v.re; arr[ij].im = v.im;
    }
#undef arr
}


//Applies a 1D fft to arr, whose (pow of 2) size is n and ln2 is log_n(n),
//Forwards or backwards fft is controlled by setting sgn to +/- 1
//This is the 'decimation in time' version with no bit reversal section
//Based on 'Matters Computational' - Arndt  (main routine similar to 21.2.1.2)
//         'Numerical Recipes in C' - Press (trig rec. portion)
static forceinline void _fft1ditbase(const ushort n, const uchar l2n,
                                     cfp*restrict*restrict arr, const char sgn){
    ushort i, j, ij, ld, m, mh, ijmh;
    fp wr, wi, cn, sn, cntmp;
    cfp u, v;

#define arr (*arr)
    //explicit 1+0i multiplications
    for(i = (ushort)0; i < n; i += (ushort)2){
        ij = i + (ushort)1;
        u.re = arr[i].re + arr[ij].re; u.im = arr[i].im + arr[ij].im;
        v.re = arr[i].re - arr[ij].re; v.im = arr[i].im - arr[ij].im;
        arr[i].re  = u.re; arr[i].im  = u.im;
        arr[ij].re = v.re; arr[ij].im = v.im;
    }

    //remaining portions of the transform
    for(ld = (ushort)2; ld < (ushort)(l2n + (uchar)1); ++ld){
        m = (ushort)1 << ld;
        mh = m >> (ushort)1;

        for(i = (ushort)0; i < n; i += m){
            ///again, some explicit 1+0i multiplications
            ijmh = i + mh;
            u.re = arr[i].re + arr[ijmh].re; u.im = arr[i].im + arr[ijmh].im;
            v.re = arr[i].re - arr[ijmh].re; v.im = arr[i].im - arr[ijmh].im;
            arr[i].re    = u.re; arr[i].im    = u.im;
            arr[ijmh].re = v.re; arr[ijmh].im = v.im;

            ///init trig rec
            wr = CT[ld]; 
            wi = (fp)sgn*ST[ld];
            cn = wr;    
            sn = wi;

            ///apply butterfly
            for(j = (ushort)1; j < mh; ++j){
                ij = i + j;
                ijmh = ij + mh;
                ////this part is u := arr[ij], v := arr[ijmh]*(cn+isn)
                ////followed by  arr[ij] := u+v, arr[ijmh] := u-v
                u.re = arr[ij].re; 
                u.im = arr[ij].im;
                v.re = arr[ijmh].re*cn - arr[ijmh].im*sn;
                v.im = arr[ijmh].re*sn + arr[ijmh].im*cn;
                arr[ij].re   = u.re + v.re; arr[ij].im   = u.im + v.im;
                arr[ijmh].re = u.re - v.re; arr[ijmh].im = u.im - v.im;

                ////update trig rec
                cntmp = cn;
                cn = wr*cn    - wi*sn;
                sn = wi*cntmp + wr*sn;
            }
        }
    }
#undef arr
}


//The first part of a 3d fft convolution.  Acts on 1D input, transforming it 
//w/o bit reversal.  arr and buff are both pow-of-2 dim lens[0]xlens[1]xlens[2],
//and l2lens are log_2(lens)
//This portion is the forward dif fft and permutes arr from ijk -> kij order
void fftconv3dif(const ushort lens[restrict 3], const uchar l2lens[restrict 3],
                 cfp*restrict*restrict arr, cfp*restrict*restrict buf){
    uint i, j, k;
    uint lo, la, lb;
    const uint llens[3] = {(uint)lens[0], (uint)lens[1], (uint)lens[2]};

    //1d ffts on rows of contigious memory (index k)
    for(i = 0u; i < llens[0]; ++i){
        for(j = 0u; j < llens[1]; ++j){
            lo = ACC3(llens[1], llens[2], i, j, 0u);
            _fft1difbase(lens[2], l2lens[2], arr + lo, (char)+1);
            ///permute ijk -> jki
            for(k = 0u; k < llens[2]; ++k){
                lb = ACC3(llens[2], llens[0], j, k, i);
                la = lo + k;   
                (*buf)[lb].re = (*arr)[la].re; (*buf)[lb].im = (*arr)[la].im;
            }
        }
    }
    //1d ffts on rows of contigious memory (index i)
    for(j = 0u; j < llens[1]; ++j){
        for(k = 0u; k < llens[2]; ++k){
            lo = ACC3(llens[2], llens[0], j, k, 0u);
            _fft1difbase(lens[0], l2lens[0], buf + lo, (char)+1);
            ///permute jki -> kij
            for(i = 0u; i < llens[0]; ++i){
                la = ACC3(llens[0], llens[1], k, i, j);
                lb = lo + i;   
                (*arr)[la].re = (*buf)[lb].re; (*arr)[la].im = (*buf)[lb].im;
            }
        }
    }
    //1d ffts on rows of contigious memory (index j)
    for(k = 0u; j < llens[2]; ++k){
        for(i = 0u; i < llens[0]; ++i){
            lo = ACC3(llens[0], llens[1], k, i, 0u);
            _fft1difbase(lens[1], l2lens[1], arr + lo, (char)+1);
            ///no need to permute here - just do the next transform backwards
        }
    }
}   



//The second part of a 3d fft convolution.  Acts on 1D input, transforming it 
//w/o bit reversal.  arr and buff are both pow-of-2 dim lens[0]xlens[1]xlens[2],
//and l2lens are log_2(lens)
//This portion is the backwards dit fft and permutes arr from kij -> ijk order
void fftconv3dit(const ushort lens[restrict 3], const uchar l2lens[restrict 3],
                 cfp*restrict*restrict arr, cfp*restrict*restrict buf){
    uint i, j, k;
    uint lo, la, lb;
    const uint llens[3] = {(uint)lens[0], (uint)lens[1], (uint)lens[2]};

    //1d ffts on rows of contigious memory (index j)
    for(k = 0u; k < llens[2]; ++k){
        for(i = 0u; i < llens[0]; ++i){
            lo = ACC3(llens[0], llens[1], k, i, 0u);
            _fft1ditbase(lens[1], l2lens[1], arr + lo, (char)-1);
            ///permute kij -> jki
            for(j = 0u; j < llens[1]; ++j){
                lb = ACC3(llens[2], llens[0], j, k, i);
                la = lo + j;   
                (*buf)[lb].re = (*arr)[la].re; (*buf)[lb].im = (*arr)[la].im;
            }
        }
    }
    //1d ffts on rows of contigious memory (index i)
    for(j = 0u; j < llens[1]; ++j){
        for(k = 0u; k < llens[2]; ++k){
            lo = ACC3(llens[2], llens[0], j, k, 0u);
            _fft1ditbase(lens[0], l2lens[0], buf + lo, (char)-1);
            ///permute jki -> ijk
            for(i = 0u; i < llens[0]; ++i){
                la = ACC3(llens[1], llens[2], i, j, k);
                lb = lo + i;   
                (*arr)[la].re = (*buf)[lb].re; (*arr)[la].im = (*buf)[lb].im;
            }
        }
    }
    //1d ffts on rows of contigious memory (index k)
    for(i = 0u; i < llens[0]; ++i){
        for(j = 0u; j < llens[1]; ++j){
            lo = ACC3(llens[1], llens[2], i, j, 0u);
            _fft1ditbase(lens[2], l2lens[2], arr + lo, (char)-1);
            ///no need to permute here - we're back in the original (ijk) order
        }
    }
}






