//Includes structure / function definitions, typedefs, and defines
#ifndef DEF_H
#define DEF_H

#include <math.h>


#define COMPILE_TEST_ROUTINES 0
#define COMPILE_PSO_FILE 1
#define PARALLEL 0 ///TODO: figure out pthreads so that this is useful

//---Program control------------------------------------------------------------
#define PREC 1 ///1 = single point (~7 digits) 2 = double point (~16 digits)

#define INFILE_LATTICE   "lattice"
#define INFILE_POTENTIAL "potential"
#define INFILE_RUNPARAMS "runparams"
#define INFILE_KPOINTS   "kpoints"
#define INFILE_KEATING   "keating"

#define OUTFILE_EIGENVAL_BASE "bands-\0\0" ///DO NOT REMOVE null terminator!
#define OUTFILE_LATTICE  "flattice"


//---Function optimization / control--------------------------------------------
#define USE_BAD_PRECONDITIONER 0 ///use Zunger's fsm precndtnr for testing

#define SAFE_SIZE_ALLOC 1 ///1 = overestimate memory needed for pws
#define SAFE_NBR_ALLOC 1 ///1 = don't attempt to shrink extra periodic image mem
#define SAFE_QR 0 ///1 = stop overflow in qr eigensolver
#define SAFE_EXPL_HAM 1 ///1 = stop oob in V(G) for explicit H if fftgmul < 2.0
#define SAFE_LCG 1 ///normalize curr (1) / curr+prev (2) |psi> in lopcg's ogs

//eigensolver defaults
#define DEF_QR_ITRLIM 25u
#define DEF_JD_ITRLIM 40u
#define DEF_CG_ITRLIM 300u
#if(PREC == 1)
#define DEF_QR_EPS 1.0000000e-6F
#define DEF_JD_EPS 1.0000000e-4F
#define DEF_DA_EPS 1.0000000e-3F
#define DEF_CG_EPS 1.0000000e-3F
#elif(PREC == 2)
#define DEF_QR_EPS 1.0000000000000000e-12
#define DEF_JD_EPS 1.0000000000000000e-08
#define DEF_CG_EPS 1.0000000000000000e-06
#endif
///conditional normalization of residual vectors in OGS ... lo, hi bounds
#if(PREC == 1)
#define DEF_CG_CNLIMLO 1.0000000e-24F
#define DEF_CG_CNLIMHI 1.0000000e+00F ///smaller = more stability
#elif(PREC == 2)
#define DEF_CG_CNLIMLO 1.0000000000000000e-36
#define DEF_CG_CNLIMHI 1.0000000000000000e+00 ///smaller = more stability
#endif

//md defaults
#define DEF_MD_ITRLIM 100u
#define MD_QUENCH_STYLE 1 ///0 = no quench, 1 = vasp style, 2 = Mattoni style
#if(PREC == 1)
#define DEF_MD_EPS              1.0000000e-3F
#define DEF_MD_FTOL             1.0000000e-2F
#define DEF_MD_VASP_QUENCH_EPS  1.0000000e-2F
#elif(PREC == 2)
#define DEF_MD_EPS              1.0000000000000000e-6
#define DEF_MD_FTOL             1.0000000000000000e-2
#define DEF_MD_VASP_QUENCH_EPS  1.0000000000000000e-2
#endif

//---Physics, math constants----------------------------------------------------
///units are eV, angstroms, amu, fs, etc. unless otherwise noted
#if(PREC == 1)
#define SMALL             1.0000000e-6F
#define ZERO              0.0000000F
#define HALF              0.5000000F
#define ONE               1.0000000F
#define PI                3.1415927F
#define TWOPI             6.2831853F
#define HBARSQD_OVER_TWOM 3.8099821F
#define VERLET_ACC_CONV   4.8242666e-3F ///1/2 (eV/A/dal*fs^2)^-1 -> A 
#elif(PREC == 2)
#define SMALL             1.0000000000000000e-15
#define ZERO              0.0000000000000000
#define HALF              0.5000000000000000
#define ONE               1.0000000000000000
#define PI                3.1415926535897932
#define TWOPI             6.2831853071795865
#define HBARSQD_OVER_TWOM 3.8099821161431488
#define VERLET_ACC_CONV   4.8242666078326639e-3 ///1/2 (eV/A/dal*fs^2)^-1 -> A
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
typedef struct keat keat;
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

    ///points to parent atom
    site* self;
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
    fp mass; ///in daltons
    //fp encut; ///not total cutoff, just cutoff that this potential was fit to

    ///stuff needed for linear interpolation between the known data points 
    ushort nsamples;
    fp*restrict samples; ///equispaced points { v_elem(q^2) }

    fp q2lo;
    fp q2hi;
};


struct job{
    uchar runmode; ///elec struct, relaxation, both? ... etc
    uchar diagmode; ///how to diagonalize H ... direct, lobcg?
    ushort nthreads; ///number of threads that work on a job at once
    ushort nbands; ///how many eigpairs to keep.  lowest or closest to fstarget?
    
    ushort itrlim; ///number of unconverged eigensolver steps before aborting
    uchar stab; ///0, 1, 2 - stability/speed tradeoff for lopcg (higher=slower)    

    fp meff; ///the effective mass m* mentioned in the 'hamil' struct 
    uchar nfstargets; ///number of reference energies for folded spectrum
    fp*restrict fstargets; ///reference energy for folded spectrum method

    fp entol; ///diff between eigvals of two itrns before converg. achieved
    fp encut; ///pw energy cutoff (eV) -> to gcut^2 (2pi/angstrm)^2
    fp fftgmul; ///multiplies cutoff for width of fft mesh (between 1.0 and 2.0)

    fp timestep; ///vel. verlet time step in fs
    ushort mditrlim; ///number of unconverged md steps before aborting
    fp rcut; ///cutoff radius for pairs, triplets (angstrm)
    fp acut; ///cutoff angle for triplets (degrees)
    fp ftol; ///magnitude of maximum force on an atom before converg. achieved

    //uchar readpsi;
    //uchar writepsi;
};

//params needed for keating potential dynamics
struct keat{
    //two-body terms
    fp*restrict r0inv; ///inverse of ideal bond distd ij, 1/agstrms.  len nspecs
    fp*restrict al; ///alpha param, eV/angstrm^2.  len nspecs

    ///three-body terms
    fp*restrict c0; ///cos(ideal angle jik).  len nspecs^2
    fp*restrict sqrtbe; ///sqrt(beta param), sqrt(eV/angstrm^2).  len nspecs
};

///data needed for the hamiltonian
struct hamil{
    kpt*restrict kp;

    ///basis
    uint npw; ///number of plane waves
    short*restrict mills; ///dim npw x 3, holds all (h, k, l) w/ |k+G| < Gcut

    ///kinetic energy, g vec magnitudes
#if USE_BAD_PRECONDITIONER
    fp*restrict la; ///inline w/ mills' dim 0: hbar^2/2m* * |G|^2 in basis
#endif
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
#define SWAP(T, a, b) do{T SWAP = a; a = b; b = SWAP;} while(0)

//---Program-wide functions-----------------------------------------------------
//---eigen.c: implementation of eigensolvers (qr, jacobi, dav, ...)
uchar qrh(const uint n, cfp*restrict A, fp*restrict e, cfp*restrict*restrict V,
		  const uint itrlim, const fp eps, const uchar vecs);
uchar lcg(const uint neig, const uint n, const hamil ham,
		  //workspace params:
		  uchar*restrict IWS,
		  fp*restrict RWS, cfp*restrict CWS, cfp*restrict*restrict CCWS,
		  //initial guess params: V0 = initial guess for eigvec, dim v0rs x v0cs
		  const uint v0rs, const uint v0cs, const cfp*restrict*restrict V0,
		  //return params: e = neig x 1 = vals, X = neig x npw = vecs                                    
		  fp*restrict e, cfp*restrict*restrict X,	  
		  //control params
		  const uchar fsm, const fp eref, const uint itrlim, const fp eps, 
		  const uchar stab);

//---fft.c: implementation of fft-like functions (a "real" fft isn't necessary)
void fftconv3dif(const ushort lens[restrict 3], const uchar l2lens[restrict 3],
                 cfp*restrict arr, cfp*restrict buf);
void fftconv3dit(const ushort lens[restrict 3], const uchar l2lens[restrict 3],
                 cfp*restrict arr, cfp*restrict buf);

//---hamil.c: dealing with init / manip of the hamiltonian op
void setmaxdims(const fp gcut, const fp gcutvmul, const fp B[restrict 3][3],
                const ushort nkpts, const kpt*restrict kpts,                  
                ushort hklmt[restrict 3], ushort hklsv[restrict 3],
                uchar hkll2sv[restrict 3], uint*restrict sizet);
void sethamkin(const fp gcut2, const ushort absmillmax[restrict 3],
               const fp B[restrict 3][3], const fp meff,
               hamil*restrict ham);
void sethamlap(const fp B[restrict 3][3], const fp meff, hamil*restrict ham); 
void sethamvloc(const lat lattice, const fp B[restrict 3][3],
                const pot*restrict pots, hamil*restrict ham);
void buildexphamil(const hamil haminfo, cfp*restrict H);

//---io.c: input / output
uchar readjob(job*restrict runparams);
uchar readlattice(lat*restrict lattice);
uchar readpotential(const uchar nspecs, pot*restrict*restrict potentials);
uchar readkpoints(const fp B[restrict 3][3],
                  ushort*restrict nkpoints, kpt*restrict*restrict kpoints);
uchar readkeating(const uchar nspecs, keat*restrict*restrict keatings);
void writeheader(const uint seed);
void writebands(const uchar id,
                const fp B[restrict 3][3], const ushort nbands, 
                const ushort nkpts, const kpt*restrict kpts);
void writelattice(const lat lattice, const char*restrict specmap);

//---lattice.c: working on real/recip lattice
void tocrdsc(const lat lattice);
void tocrdsf(const fp B[restrict 3][3], const lat lattice);
void backtocell(const lat lattice, const uchar imgs);
void torecipbasis(fp A[restrict 3][3]);
void toinverse(fp A[restrict 3][3]);

//---md.c: molecular dynamics
void setnbrs(const lat lattice, 
             ushort*restrict nxatoms, site*restrict*restrict xatoms,
             const fp rcut, const fp acut, const fp eps);
uchar relax(fp*restrict RWS, const lat lattice, const pot*restrict pots,
            const keat*restrict keatings, const fp ts,
            const fp maxforcetol, const uint itrlim);

//---monpac.c: k-point grids
void mpgridfi(const fp B[restrict 3][3], const ushort sds[restrict 3],
              kpt*restrict kpts);

//---pso.c: particle swarm optimization of keating params
void psomin(job* runparams);

//---(test)routines.c: different runmodes called from main
void runstaticbs(job* runparams);
void runmd(job* runparams);

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
#define DEG2RAD(x) (x*1.7453293e-2F)
#define RAD2DEG(x) (x*5.7295780e+1F)
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
#define DEG2RAD(x) (x*1.7453292519943296e-2)
#define RAD2DEG(x) (x*5.7295779513082321e+1)
#endif

#endif
//16, 7