//Includes structure / function definitions, typedefs, and defines
#ifndef DEF_H
#define DEF_H

#include <math.h>


#define COMPILE_TEST_ROUTINES 0
#define PARALLEL 0 ///TODO: figure out pthreads so that this is useful

//---Program control------------------------------------------------------------
#define PREC 1 ///1 = single point (~7 digits) 2 = double point (~16 digits)

#define INFILE_LATTICE   "lattice"
#define INFILE_POTENTIAL "potential"
#define INFILE_RUNPARAMS "runparams"
#define INFILE_KPOINTS   "kpoints"

#define OUTFILE_EIGENVALS "bands"


//---Function optimization / control--------------------------------------------
#define SAFE_SIZE_ALLOC 0 ///overestimate memory needed for pws
#define SAFE_QR 0 ///stop overflow in qr eigensolver
#define SAFE_EXPL_HAM 1 ///stop oob in V(G) for explicit H if fftgmul < 2.0

//eigensolver defaults
#define QR_ITR_LIM 25u
#if(PREC == 1)
#define DEF_QR_EPS 1.0000000e-6F
#define DEF_JD_EPS 1.0000000e-4F
#define DEF_DA_EPS 1.0000000e-3F
#elif(PREC == 2)
#define DEF_QR_EPS 1.0000000000000000e-12
#define DEF_JD_EPS 1.0000000000000000e-08
#define DEF_DA_EPS 1.0000000000000000e-06
#endif

//---Physics, math constants----------------------------------------------------
///units are eV, angstroms, etc. unless otherwise noted
#if(PREC == 1)
#define ZERO              0.0000000F
#define HALF              0.5000000F
#define ONE               1.0000000F
#define PI                3.1415927F
#define TWOPI             6.2831853F
#define HBARSQD_OVER_TWOM 3.8099821F
#elif(PREC == 2)
#define ZERO              0.0000000000000000
#define HALF              0.5000000000000000
#define ONE               1.0000000000000000
#define PI                3.1415926535897932
#define TWOPI             6.2831853071795865
#define HBARSQD_OVER_TWOM 3.8099821161431488
#endif

//---Typedefs-------------------------------------------------------------------
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
#if(PREC == 1)
typedef float fp;
typedef double lfp;
#elif(PREC == 2)
typedef double fp;
typedef long double lfp;
#endif

typedef struct cfp cfp;
typedef struct kpt kpt;
typedef struct site site;
typedef struct lat lat;
typedef struct pot pot;
typedef struct job job;
typedef struct hamil hamil;

//---Structures-----------------------------------------------------------------
///for complex numbers without the insane complexity of the stdlib
struct cfp{
    fp re;
    fp im;
};

///k-points
struct kpt{
    fp wgt; ///k-point weight, = 1/(num kpts) unless symmetry is applied
    fp crds[3]; ///coords in recip space, usually direct but may be 2pi/angstrm

    ///info about this k-point's eigenvalues / eigenvectors
    fp*restrict evals;
    cfp*restrict*restrict evecs; ///dim neigs x npw
};

///data for sites / atoms.  Most stuff is for force field calcs
struct site{
    ///stuff for this specific site
    uchar* spec; ///id # -> element type.  pointer b/c of the s.a. atom hopping
    fp crdsf[3]; ///fractional (crystal) coords
    fp crdsc[3]; ///cartesian coords

    ///stuff for this site's neighbors
    uchar npairs;
    site**restrict pairs; ///points to nbrs 
    uchar*restrict ntrips; ///ntrips[i] = length of ith trip
    site**restrict*restrict trips; ///trips[i][j] = i'th pair's j'th nbr

    ///stuff for this site's periodic images
    ushort nimgs;
    site**restrict imgs;

    ///may be useful later?  points to parent atom
    //site* self;
};

///data describing real-space the lattice
struct lat{
    fp A[3][3]; ///{{a1}, {a2}, {a3}] (i.e. a* is the *th row of A)

    uchar nspecs;
    ushort*restrict speccounts;

    site*restrict sites; ///nsites = sum(speccounts)
};

///stuff for psueudopotentials.  a struct is a bit overkill here, but may need
///to be extended later (esp. for nonlocal contributions, etc)
///for now, has elemental local potentials as a function of q^2 = |G_i - G_j|^2
struct pot{
    uchar nvel; ///number of explicit electrons for this potential
    //fp encut; ///not total cutoff, just cutoff that this potential was fit to

    ///stuff needed for linear interpolation between the known data points 
    ushort nsamples;
    fp*restrict samples; ///equispaced points { v_elem(q^2) }

    fp q2lo;
    fp q2hi;
};


struct job{
    uchar runmode; ///elec struct, relaxation, both? ... etc
    uchar diagmode; ///how to diagonalize H ... direct, dav, folded spectrum?
    ushort nthreads; ///number of threads that work on a job at once
    ushort nbands; ///how many eigpairs to keep.  lowest or closest to fstarget?
    fp mbsmul; ///multiplies nbands for max basis size before reset in dav
    
    fp meff; ///the effective mass m* mentioned in the 'hamil' struct 
    fp fstarget; ///reference energy for folded spectrum method

    fp entol; ///diff between eigvals of two itrns before converg. achieved
    fp encut; ///pw energy cutoff (eV) -> to gcut^2 (2pi/angstrm)^2
    fp fftgmul; ///multiplies cutoff for width of fft mesh (between 1.0 and 2.0)

    //uchar readpsi;
    //uchar writepsi;
};

///data needed for the hamiltonian
struct hamil{
    kpt*restrict kp;

    ///basis
    uint npw; ///number of plane waves
    short*restrict mills; ///dim npw x 3, holds all (h, k, l) w/ |k+G| < Gcut

    ///kinetic energy
    fp*restrict ke; ///inline w/ mills' dim 0: hbar^2/2m* * |k+G|^2 in basis

    ///potential energy (local)
    ushort dims[3]; ///length in dirs h, k, l
    uchar l2dims[3]; ///log_2(dims)
    cfp*restrict vloc; ///the local portion of V(G) or V(r) on a grid
};


//---Macros (oOoOoO ... scary!)-------------------------------------------------
//---force gcc to do things ... I shouldn't have to even do this
#define forceinline __attribute__((always_inline)) inline

//---access matrices stored as arrays (row-major order, as god intended)
#define ACC2(d1,i,j) ((i)*(d1)+(j)) ///i,j of d0xd1 matrix
#define ACC3(d1,d2,i,j,k) (((i)*(d1)+(j))*(d2)+(k)) ///i,j,k of d0xd1xd2 tensor

#define MIN(a,b) (((a)<(b)) ? (a):(b))
#define MAX(a,b) (((a)>(b)) ? (a):(b))


//---Program-wide functions-----------------------------------------------------
//---eigen.c: implementation of eigensolvers (qr, jacobi, dav, ...)
uint qrh(const uint n, cfp*restrict*restrict A, 
		 fp*restrict*restrict e, cfp*restrict*restrict*restrict V,
		 const uint itrLim, const fp eps, const uchar vecs);

//---fft.c: implementation of fft-like functions (a "real" fft isn't necessary)
void fftconv3dif(const ushort lens[restrict 3], const uchar l2lens[restrict 3],
                 cfp*restrict*restrict arr, cfp*restrict*restrict buf);
void fftconv3dit(const ushort lens[restrict 3], const uchar l2lens[restrict 3],
                 cfp*restrict*restrict arr, cfp*restrict*restrict buf);

//---hamil.c: dealing with init / manip of the hamiltonian op
void setmaxdims(const fp gcut, const fp gcutvmul, const fp B[restrict 3][3],
                const ushort nkpts, const kpt*restrict kpts,                  
                ushort hklmt[restrict 3], ushort hklsv[restrict 3],
                uchar hkll2sv[restrict 3], uint*restrict sizet);
void sethamkin(const fp gcut2, const ushort absmillmax[restrict 3],
               const fp B[restrict 3][3], const fp meff,
               hamil*restrict ham); 
void sethamvloc(const lat lattice, const pot*restrict pots,
                hamil*restrict ham);
void buildexphamil(const hamil haminfo, cfp*restrict*restrict H);

//---io.c: input / output
uchar readjob(job*restrict runparams);
uchar readlattice(lat*restrict lattice);
uchar readpotential(const uchar nspecs, pot*restrict*restrict potentials);
uchar readkpoints(ushort*restrict nkpoints, kpt*restrict*restrict kpoints);
void reportjob(const job* runparams);
void reportlattice(const lat* lattice);
void reportpotential(const uchar nspecs, const pot* potentials);
void reportkpoints(const ushort nkpoints, const kpt* kpoints);
void writeheader();
void writebands(const fp B[restrict 3][3], const ushort nbands, 
                const ushort nkpts, const kpt*restrict kpts);

//---lattice.c: working on real/recip lattice
void tocrdsc(lat* lattice);
void torecipbasis(fp A[restrict 3][3]);

//---routines.c: different runmodes called from main
void runstaticbs(job* runparams);
void testio(job runparams);
void testqrh();
void testexphamilgen(job runparams);

//---"Overloaded" common functions that need to be called differently depending
// on the precision level-------------------------------------------------------
#if(PREC == 1)
//---io
#define FPFORMAT "%f"
//---math
#define ABSF(x) fabsf(x)
#define SQRT(x) sqrtf(x)
#define COS(x) cosf(x)
#define ACOS(x) acosf(x)
#define SIN(x) sinf(x)
#define ASIN(x) asinf(x)
#elif(PREC == 2)
//---io
#define FPFORMAT "%lf" 
//---math
#define ABSF(x) fabs(x) 
#define SQRT(x) sqrt(x)
#define COS(x) cos(x)
#define ACOS(x) acos(x)
#define SIN(x) sin(x)
#define ASIN(x) asin(x)
#endif

#endif
