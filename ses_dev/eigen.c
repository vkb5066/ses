#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "def.h"



//Aux function for calculating (a^2 + 1)^(1/2) w/o under/overflow
//If you aren't a coward, you'll turn off SAFE_QR and calculate it w/o
//any fancy manipulations
static forceinline fp _shypot(const fp a){
#if SAFE_QR
	fp absa;
	absa = ABSF(a);
	if(absa > ONE) return absa*SQRT(ONE + ONE/(absa*absa));
	return SQRT(ONE + absa*absa);
#else
	return SQRT(a*a + ONE);
#endif
}

//Sorts eigenvals / vecs low -> high according to eigval
//vecs are expected to be stored as rows
//Uses insertion sort
static forceinline void _epairsort(const uint n, fp*restrict*restrict va,
								   cfp*restrict*restrict*restrict ve,
								   const uchar vecs){
	uint i, j;
    fp* ea_; fp* eb_; fp et_;
    cfp*restrict* va_; cfp*restrict* vb_; cfp* vt_;
	
	if(vecs){ ///sort evals and evecs at the same time
		for(i = 0u; i < n; ++i){
			for(j = i; j > 0u;){
					eb_ = *va + j;
					vb_ = *ve + j;
					ea_ = *va + --j;
					va_ = *ve + j;
				if(*eb_ < *ea_){
					et_ = *ea_;
					*ea_ = *eb_;
					*eb_ = et_;
					vt_ = *va_;
					*va_ = *vb_;
					*vb_ = vt_;
				}
			}
		}
	}
	else{ ///only sort evals
		for(i = 0u; i < n; ++i){
			for(j = i; j > 0u;){
					eb_ = *va + j;
					ea_ = *va + --j;
				if(*eb_ < *ea_){
					et_ = *ea_;
					*ea_ = *eb_;
					*eb_ = et_;
				}
			}
		}
	}
}

//QR Decomposition of a hermitian matrix into eigenvalues / eigenvectors
/* A:    input matrix in 1-D form: A = {a00, a01, a10, a11} for a 2x2 matrix
** n:    dimension of A
** e:    pre-allocated array of eigenvalues, size n (no need to initialize)
** V:    pre-allocated ragged matrix of eigenvectors (vectors held row-wise)
** vecs: compute eigvecs?  if not, V may be a dummy pointer - no malloc needed
** Returns: 0 on success, 1 if the itr limit is reached
**         the ith entry of (e, V) is the ith eigenval, eigenvec pair.
** Note that A will be destroyed in this function!
*/
uchar qrh(const uint n, cfp*restrict A, 
		  fp*restrict e, cfp*restrict*restrict V,
		  const uint itrlim, const fp eps, const uchar vecs){
//This is a re-implementation of an algorithm from 'mpmath':
//https://github.com/fredrik-johansson/mpmath
//which is itself based on a fortran implementation in 'EISPACK'.    
//This specific version has been edited to make the most accesses in row-order 
//instead of column order.
//#define A (*A)
//#define V (*V)
	uint i, j, k, l, m;
	uint i0, i1, i2;	

	fp fp1, fp2, fp3, fp4, fp5, fp6, fp7;
	fp dg; cfp cg, cfp1, f0;
	fp* d; cfp* t;

	memset(e, 0, n*sizeof(fp)); ///evals -> zero for iterative updates 
	d = calloc(n, sizeof(fp));
	t = calloc(n, sizeof(cfp));

	//--------------------------//
	//--- TRIDIAGONALIZATION ---//
	//--------------------------//
	t[n - 1u].re = ONE; t[n - 1u].im = ZERO;
	for(i = n - 1u; i > 0u; i--){
		l = i - 1u;

		///vector scaling
		fp1 = ZERO;
		for(i0 = ACC2(n, i, 0u); i0 < ACC2(n, i, i); ++i0){
			fp1 += ABSF(A[i0].re) + ABSF(A[i0].im);
		}

		///skipping, stopping criteria
		if(ABSF(fp1) < eps){
			d[i] = ZERO;
			e[i] = ZERO;
			t[l].re = ONE; t[l].im = ZERO;
			continue;
		}

		i0 = ACC2(n, i, l);
		if(i == 1u){
			f0.re = A[i0].re; f0.im = A[i0].im; 
			fp3 = SQRT(f0.re*f0.re + f0.im*f0.im);
			fp2 = ONE / fp3;
			d[i] = fp3;
			e[i] = ZERO;
			if(fp3 > eps){
				t[l].re = (t[i].re*f0.re - t[i].im*f0.im) * fp2;
				t[l].im = (t[i].re*f0.im + t[i].im*f0.re) * fp2;
			}
			else{
				t[l].re = t[i].re; t[l].im = t[i].im;
			}
			continue;
		}

		///setup for householder transformation
		fp2 = ONE / fp1; fp4= ZERO;
		for(i1 = ACC2(n, i, 0u); i1 < ACC2(n, i, i); ++i1){
			A[i1].re *= fp2; A[i1].im *= fp2;
			fp4 += A[i1].re*A[i1].re + A[i1].im*A[i1].im;
		}

		f0.re = A[i0].re; f0.im = A[i0].im;
		fp3 = SQRT(f0.re*f0.re + f0.im*f0.im);
		cg.re = SQRT(fp4); cg.im = ZERO;
		fp4+= cg.re*fp3;   ///at this point, g has no imaginary component
		d[i] = fp1*cg.re;  ///ditto ^
		if(fp3 > eps){
			f0.re /= fp3; f0.im /= fp3;
			cfp1.re =  f0.im*t[i].im - f0.re*t[i].re;
			cfp1.im = -f0.re*t[i].im - f0.im*t[i].re;
			fp5 = cg.re;
			cg.re = cg.re*f0.re; 
			cg.im = fp5*f0.im; 
		}
		else{
			cfp1.re = -t[i].re; cfp1.im = -t[i].im;
		}

		A[i0].re += cg.re; A[i0].im += cg.im;
		f0.re = ZERO; f0.im = ZERO;

		///apply householder transformation
		for(j = 0u, i0 = ACC2(n, i, 0u), k = ACC2(n, 0u, i); i0 < ACC2(n, i, i); /*ended here.  I've been getting rid of simple counters in place of i* counters.  i.e. if we only care about j < i, but we use i1 instead of j, we don't even need to increment or store j*/ 
			++j, ++i0, k += n){
			A[k].re = A[i0].re / fp4; A[k].im = A[i0].im / fp4;

			////A dot U
			cg.re = ZERO; cg.im = ZERO;
			for(m = 0u, i1 = ACC2(n, j, 0u), i2 = ACC2(n, i, 0u); m < j + 1u; 
				++m, ++i1, ++i2){
				/////g += A(j, k)*  dot  A(i, k)
				cg.re += A[i1].re*A[i2].re + A[i1].im*A[i2].im;
				cg.im += A[i1].re*A[i2].im - A[i1].im*A[i2].re;
			}
			for(m = j + 1u, i1 = ACC2(n, j+1u, j), i2 = ACC2(n, i, j+1u); m < i;
			 	++m, i1 += n, ++i2){
				/////g += A(k, j)  dot  A(i, k)
				cg.re += A[i1].re*A[i2].re - A[i1].im*A[i2].im;
				cg.im += A[i1].re*A[i2].im + A[i1].im*A[i2].re;
			}

			////P (f0 += t[j]*  dot  A(i, j) 
			t[j].re = cg.re / fp4; t[j].im = cg.im / fp4;
			f0.re += t[j].re*A[i0].re + t[j].im*A[i0].im;
			f0.im += t[j].re*A[i0].im - t[j].im*A[i0].re;
		}

		///reduce A, get Q
		fp5 = HALF*f0.re/fp4; ////by now, f0 is real (to numerical precision)
		for(j = 0u, i0 = ACC2(n, i, 0u); j < i; ++j, ++i0){
			f0.re = A[i0].re; f0.im = A[i0].im;
			cg.re = t[j].re - fp5*f0.re; cg.im = t[j].im - fp5*f0.im;
			t[j].re = cg.re; t[j].im = cg.im;

			////this part is expensive
			for(k = 0u, i1 = ACC2(n, j, 0u), i2 = ACC2(n, i, 0u); k < j + 1u; 
				++k, ++i1, ++i2){
				/////A(j, k) -= f0*  dot  t[k]    +    g*  dot  A(i, k)
				A[i1].re -= f0.re*t[k].re  + f0.im*t[k].im + 
							cg.re*A[i2].re + cg.im*A[i2].im;
				A[i1].im -= f0.re*t[k].im  - f0.im*t[k].re +
							cg.re*A[i2].im - cg.im*A[i2].re;
			}
		}

		t[l].re = cfp1.re; t[l].im = cfp1.im;
		e[i] = fp4;
	}

	//--------------------------//
	//------ INTERMISSION ------//
	//--------------------------//
	//Convient to shift off-diagonals by 1
	for(i = 1u; i < n; ++i) d[i - 1u] = d[i];
	d[n - 1u] = ZERO;
	//Also need to update current accum eigs and A
	e[0] = ZERO;
	for(i = 0u, i0 = 0u; i < n; ++i, i0 += (n + 1u)){
		fp1 = e[i];
		e[i] = A[i0].re;
		A[i0].re = fp1;
	}

	//Finally, set V to identity for eigenvectors
	if(vecs){
    	for(i = 0u; i < n; ++i){
			for(j = 0u; j < n; ++j){
				V[i][j].re = ZERO; V[i][j].im = ZERO;
			}
			V[i][i].re = ONE; ///shouldn't have to set the imag portion
		}
	}


	//--------------------------//
	//--- EIGENVALUES / VECS ---//
	//--------------------------//
	for(l = 0u; l < n; ++l){
		k = 0u; ///iteration counter: breaks if it reaches itrlim

		while(1){
			///Grab a small off-diag element
			m = l;
			while(1){
				if(m + 1u == n) break;
				if(ABSF(d[m]) < eps*(ABSF(e[m]) + ABSF(e[m+1u]))) break;
				m++;
			}
			if(m == l) break;
	
			if(k >= itrlim) return (uchar)1;  ////prevent hanging if QR fails		
			k++;

			///shift
			dg = HALF*(e[l+1u] - e[l])/d[l];
			fp2 = _shypot(dg);
			fp3 = (dg < ZERO)? dg - fp2 : dg + fp2;
			dg = e[m] - e[l] + d[l]/fp3;

			///plane->Givens rotations: get back to tridiagonal form
			fp3 = ONE; fp4 = ONE; fp5 = ZERO;
			for(i = m; i > l; --i){
				fp7 = fp3*d[i-1u]; 
                fp6 = fp4*d[i-1u];

				////improves convergence by choosing a large denom
				if(ABSF(fp7) > ABSF(dg)){
					fp4 = dg/fp7;
					fp2 = _shypot(fp4);
					d[i] = fp7*fp2;
					fp3 = ONE/fp2;
					fp4 *= fp3;
				}
				else{
					fp3 = fp7/dg;
					fp2 = _shypot(fp3);
					d[i] = dg*fp2;
					fp4 = ONE/fp2;
					fp3 *= fp4;
				}

				dg = e[i] - fp5;
				fp2 = (e[i-1u] - dg)*fp3 + ((fp)2.0)*fp4*fp6;
				fp5 = fp3*fp2;
				e[i] = dg + fp5;
				dg = fp4*fp2 - fp6;

				////accum eigenvectors
				///At this point, we can deal with all real arithmatic
				if(vecs){
                	for(j = 0u; j < n; ++j){
                    	fp7 = V[i][j].re;
                    	////V(i+1, j) = s*V(i, j) + c*f
                    	V[i][j].re = fp3*V[i-1u][j].re + fp4*fp7;
                    	////V(i, j) = c*V(i, j) - s*f
                    	V[i-1u][j].re = fp4*V[i-1u][j].re - fp3*fp7;
                	}
				}
			}

			///finish up updating
			e[l] -= fp5;
			d[l] = dg;
			d[m] = ZERO;
		}
	}
	free(d);


	//--------------------------//
	//---- FINISH EIGENVECS ----//
	//--------------------------//
	//Grab the imaginary components for the eigenvectors
	if(vecs){
		for(i = 0u; i < n; ++i){
			for(j = 0u; j < n; ++j){
				fp7 = V[i][j].re;
				V[i][j].re = V[i][j].re*t[j].re + V[i][j].im*t[j].im;
				V[i][j].im = V[i][j].im*t[j].re - fp7*t[j].im;
			}
		} 
		for(i = 0u, i0 = 0u; i < n; ++i, i0 += (n+1u)){
			if(ABSF(A[i0].re) < eps) continue;
			for(j = 0u; j < n; ++j){
				cfp1.re = ZERO; cfp1.im = ZERO;
				///Expensive part of the algo + complicated math = no functions
				///last time, I promise :)
				for(k = 0u, i1 = ACC2(n, 0, i); k < i; ++k, i1 += n){
					////ctmp += A(k, i)*  dot  V(j, k)
					cfp1.re += A[i1].re*V[j][k].re - A[i1].im*V[j][k].im;
					cfp1.im += A[i1].re*V[j][k].im + A[i1].im*V[j][k].re;
				}
				for(k = 0u, i1 = ACC2(n, i, 0u); k < i; ++k, ++i1){
					////V(j, k) -= ctmp  dot  A(i, k)
					V[j][k].re -= cfp1.re*A[i1].re + cfp1.im*A[i1].im;
					V[j][k].im += cfp1.re*A[i1].im - cfp1.im*A[i1].re;
				}
			}
		}
	}
	free(t); 

//#undef A
//#undef V		
	 
	_epairsort(n, &e, &V, vecs);
	return (uchar)0;
}


//places |psi> in vector form onto the grid
static forceinline void _psivtopsig(const uint npw, const short*restrict mills,
									const cfp*restrict psiv, 
									const ushort dims[restrict 3], 
									cfp*restrict psig){
	short hm, km, lm;
	uint i, i3, ind;
	uint hi, ki, li;

	hm = (short)dims[0] >> (short)1;
	km = (short)dims[1] >> (short)1;
	lm = (short)dims[2] >> (short)1;

	for(i = 0u, i3 = 0u; i < npw; ++i, i3 += 3u){
		hi = (uint)(mills[i3]      + hm);
		ki = (uint)(mills[i3 + 1u] + km);
		li = (uint)(mills[i3 + 2u] + lm);

		ind = ACC3((uint)dims[1], (uint)dims[2], hi, ki, li);
		psig[ind].re = psiv[i].re;
		psig[ind].im = psiv[i].im;
	}
}

//places |psi> on the grid into vector form
static forceinline void _psigtopsiv(const uint npw, const short*restrict mills,
									const cfp*restrict psig,
									const ushort dims[restrict 3],
									cfp*restrict psiv){
	short hm, km, lm;
	uint i, i3, ind;
	uint hi, ki, li;

	hm = dims[0] - (short)1;
	km = dims[1] - (short)1;
	lm = dims[2] - (short)1;

	for(i = 0u, i3 = 0u; i < npw; ++i, i3 += 3u){
		hi = (uint)((mills[i3]      + (short)dims[0]) & hm);
		ki = (uint)((mills[i3 + 1u] + (short)dims[1]) & km);
		li = (uint)((mills[i3 + 2u] + (short)dims[2]) & lm);

		ind = ACC3((uint)dims[1], (uint)dims[2], hi, ki, li);
		psiv[i].re = psig[ind].re;
		psiv[i].im = psig[ind].im;
	}
}


//Computes the action of a wavefunction on the hamiltonian and stores the result
//in res by ffting psiv into real space, multipling to vr, then ffting back into
//reciprocal space.  psiv, the original wavefunctions, are not updated (since we
//need them later in the dav code)
//dims are the dimensions of the real space potential vr and l2dims are their
//base-2 logs.  
static forceinline void _action(const ushort dims[restrict 3], /// = d 
								const uchar l2dims[restrict 3],
								const cfp*restrict vr, ///d[0] x d[1] x d[2]
								cfp*restrict buffft, ///ditto ^
								cfp*restrict psig, ///ditto ^ 
								const uint npw, 
								const short*restrict mills, ///npw x 3 
								const fp*restrict ke, ///npw x 1
								const cfp*restrict psiv, ///npw x 1
								cfp*restrict res){ ///npw x 1
	uint i, tlen;
	fp norm, vrnorm;

	tlen = (uint)dims[0] * (uint)dims[1] * (uint)dims[2];
	memset(psig, 0, tlen*sizeof(cfp)); ///needed since we reuse psig

	//|psi> coeffs -> real space grid
	_psivtopsig(npw, mills, psiv, dims, psig);
	fftconv3dif(dims, l2dims, psig, buffft);
	//action calculation in real space (note V(G) is hermitian so V(r) is real)
	norm = ONE / (fp)tlen; ///takes care of 1/(N1N2N3) for inv. fft
	for(i = 0u; i < tlen; ++i){
		vrnorm = vr[i].re * norm;
		psig[i].re *= vrnorm;
		psig[i].im *= vrnorm;
	}
	//|psi> on real space grid -> into result
	fftconv3dit(dims, l2dims, psig, buffft);
	_psigtopsiv(npw, mills, (const cfp*restrict)psig, dims, res);

	//finish up with the kinetic energy part
	for(i = 0u; i < npw; ++i){
		res[i].re += psiv[i].re * ke[i];
		res[i].im += psiv[i].im * ke[i];
	}
}


//modified grahm schmidt on V of size nvec x dimvec.  Q is V's buffer
//TODO: can we get rid of all of these damn calls to ACC2?? seems like we're only going over vectors, so this shouldn't be too hard
void _mgs(const uint nvec, const uint dimvec, cfp*restrict*restrict V, cfp*restrict*restrict Q){
#define V (*V)
#define Q (*Q)
	uint i, j, k, ij, ik, jk;
	fp norm, impt;
	cfp*restrict tmp;

	for(i = 0u; i < nvec; ++i){
		//normalize V into Q
		norm = ZERO;
		for(j = 0u; j < dimvec; ++j){
			ij = ACC2(dimvec, i, j);
			norm += V[ij].re*V[ij].re + V[ij].im*V[ij].im;  ///V(i,j)* @ V(i,j)
		}
		norm = ONE / SQRT(norm);
		for(j = 0u; j < dimvec; ++j){
			ij = ACC2(dimvec, i, j);
			Q[ij].re = V[ij].re * norm; Q[ij].im = V[ij].im * norm;
		}

		//orthogonalize all j > i to i in 0, ... , j
		for(j = i + 1u; j < nvec; ++j){
			norm = ZERO; impt = ZERO; ///treat as complex val = norm + i(impt)
			///val = Q(i,k)* @ V(j,k)
			for(k = 0u; k < dimvec; ++k){
				ik = ACC2(dimvec, i, k);
				jk = ACC2(dimvec, j, k);
				norm += (Q[ik].re*V[jk].re + Q[ik].im*V[jk].im);
				impt += (Q[ik].re*V[jk].im - Q[ik].im*V[jk].re);
			}
			///V(j, k) -= val @ Q(i,k)
			for(k = 0u; k < dimvec; ++k){
				ik = ACC2(dimvec, i, k);
				jk = ACC2(dimvec, j, k);
				V[jk].re -= (norm*Q[ik].re - impt*Q[ik].im);
				V[jk].im -= (norm*Q[ik].im + impt*Q[ik].re);
			}
		}
	}
#undef V
#undef Q
	//now we're done, but we've updated Q instead of V ... just flip some
	//pointers around and noone will know any better
	tmp = *V; *V = *Q; *Q = tmp;
}

//recovers the "correct" eigenpairs e,v from (lam-eref)^2
//fsm leavs v unaffected, so eigenvalues are the ones that best approximate
//Hv - ev = 0
//workspace WS must be at least neig*npw*sizeof(cfp) long, will be overwritten
//workspaces psig* must be fftdim[0]*[1]*[2]*sizeof(cfp) long
static void _fsmrecover(const hamil ham, const uint neig, const uint npw,
						const fp eref,
						cfp*restrict WS, cfp*restrict psig1, cfp*restrict psig2,
						fp*restrict e, const cfp*restrict*restrict v){
	uint i, j, ij;
	fp avg1, avg2, lam1, lam2, resre, resim;

	for(i = 0u; i < neig; ++i){
		ij = i*npw;
		//set product <psi|H* -> WS
		_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, npw, 
				ham.mills, ham.ke, v[i], WS + ij);		

		//determine which trial eigval gives the smallest error in H*v - lam*v
		avg1 = ZERO; avg2 = ZERO;
		resre = SQRT(ABSF(e[i])); lam1 = eref + resre; lam2 = eref - resre;
		for(j = 0u; j < npw; ++j){
			///test 1
			resre = WS[ij+j].re - lam1*v[i][j].re;
			resim = WS[ij+j].im - lam1*v[i][j].im;
			avg1 += resre*resre + resim*resim;
			///test 2
			resre = WS[ij+j].re - lam2*v[i][j].re;
			resim = WS[ij+j].im - lam2*v[i][j].im;
			avg2 += resre*resre + resim*resim;
		}
															//printf(FPFORMAT" vs "FPFORMAT" ("FPFORMAT" vs "FPFORMAT")\n", lam1, lam2, avg1, avg2);
		e[i] = (avg1 < avg2)? lam1 : lam2;
	}
}

//Davidson iterative eigenpair solver for diagonally dominant hermitian 
//matrices, specialized for electronic structure calculations (i.e. input is
//a hamiltonian constructor as opposed to an explicit matrix)
//too many parameters to write out here, read the comments in the source
//code
uchar dav(const uint n, const hamil ham,  	 
		  //buffer params: V, Q, W = mbs x npw, psig* = fftdim[0]*[1]*[2]
		  //               L = pointer-to-pointer mbs x mbs
		  const uint mbs, 
		  cfp*restrict V, cfp*restrict Q, cfp*restrict W, 
		  cfp*restrict*restrict L,
		  cfp*restrict psig1, cfp*restrict psig2,
		  //initial guess params: V0 = initial guess for eigvec, dim v0rs x v0cs
		  const uint v0rs, const uint v0cs, const cfp*restrict*restrict V0,
		  //return params: e = rbs x 1 = vals, X = rbs x npw = vecs
		  const uint rbs, ///(rbs is usually neigs)                                          
		  fp*restrict e, cfp*restrict*restrict X,	
		  //control params
		  const uint fsm, const fp eref, const uint itrlim, const fp eps){
//this is based on the version of the Davidson method given in section 3.1 of
//'The Davidson Method' - Crouzeix et. al. (doi 10.1137/0915004) 
//this version differs from the version in the text in the following ways:
//1. it works for complex problems
//2. step 1 of the algo changed from a matrix multiplication to a convolution
//3. explicit calc / storage of the residuals is ditched (conv. criteria only
//   based on eigenvalues)
//4. lots (an almost worrying amount) of re-use of pre-allocated buffer space
//   (i.e. the rayleigh matrix, the matrix of search directions, ... all use
//   the same memory - see comments within)
//5. it can handle the folded spectrum method ('Solving Schrodinger's Equation
//   Around a Desired Energy: Application to Silicon Quantum Dots' - Wang &
//   Zunger (doi 10.1063/1.466486))
//6. it has been re-arranged to make as many calculations row-order friendly
//   as I could manage.  the ritz vectors and search directions still have
//   bad access patterns, though
//that out of the way, the 'steps' that I refer to in the below comments will 
//correspond to the numbered steps in the aformentioned algo 3.1 of Crouzeix

	uint i, ii, j, ij, ji, ijp, k, ik, jk, kj, l; ///TODO: some of these arent needed ... actually a lot aren't needed
	uint cbs, cnt;
	fp norm, curre, laste, diffe;
	cfp keavg2; fp vavgsc, pre;

	//un-fixable errors (if any of these throw, work arrays are too small!)
	if(n < rbs) return (uchar)2;
	if(mbs < 2u*rbs) return(uchar)3;
	if(n <= mbs) return (uchar)4;


	//--------------------------//
	//----- INITIALIZATION -----//
	//--------------------------//
	//set V ...
	k = MIN(mbs, v0rs);
	l = MIN(n, v0cs);
	///... from V0 if possible ...
	for(i = 0u; i < k; ++i){
		for(j = 0u; j < l; ++j){
			ij = ACC2(n, i, j);
			V[ij].re = V0[i][j].re; V[ij].im = V0[i][j].im;
		}
	}
	///... and fill the rest randomly ...
	for(i = k; i < mbs; ++i){
		for(j = l; j < n; ++j){
			ij = ACC2(n, i, j);
			V[ij].re = (fp)rand() / (fp)RAND_MAX;
			V[ij].im = (fp)rand() / (fp)RAND_MAX;
		}
	}

	///... and make sure its normalized
	for(i = 0u; i < mbs; ++i){
		norm = ZERO;
		for(j = 0u; j < n; ++j){
			ij = ACC2(n, i, j);
			norm += V[ij].re*V[ij].re + V[ij].im*V[ij].im;
		}
		norm = ONE / SQRT(norm);
		for(j = 0u; j < n; ++j){
			ij = ACC2(n, i, j);
			V[ij].re *= norm; V[ij].im *= norm;
		}
	}
	/* technically, we should now guarentee that V is orthogonal. 
	 * realistically, any reasonable guess of V will satisfy this, and the 
	 * random init makes it _very_ likley that we'll satisfy the reqmnt.
	 * (if we do decide to call MGS(V), we can skip the above normalization)
	*/  

	//a bit more stuff is needed if we're doing the folded spectrum method:
	if(fsm){
		///operate I @ eref on hamiltonian
		for(i = 0u; i < n; ++i) ham.ke[i] -= eref;
		///and setup the average real-space potential for preconditioning
		///(scaled by eref)
		i = (uint)ham.dims[0] * (uint)ham.dims[1] * (uint)ham.dims[2];
		vavgsc = ZERO;
		for(j = 0u; j < i; ++j) vavgsc += ham.vloc[j].re;
		vavgsc = vavgsc/(fp)i - eref;
	}


	//--------------------------//
	//------- MAIN  LOOP -------//
	//--------------------------//
	laste = ZERO;
	for(cnt = 0u, cbs = rbs; cnt < itrlim; ++cnt, cbs += rbs){

		//step 1: calculate H|psi> -> W through fft-like convolutions
		if(fsm){
			///apply H operator twice to mimic squaring the hamiltonian
			///also use Q for the first calc since restricted pointers prohibit
			///double-using W for the second part
			for(i = 0u; i < cbs; ++i){
				ij = i*n;
				_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, n, 
						ham.mills, ham.ke, V + ij, Q + ij);
				_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, n, 
						ham.mills, ham.ke, Q + ij, W + ij);
			}
		}
		else{
			///this is the 'normal' way of applying <psi| to H, nothing fancy
			for(i = 0u; i < cbs; ++i){
				ij = i*n;
				_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, n, 
						ham.mills, ham.ke, V + ij, W + ij);
			}
		}

		//step 2: calculate the rayleigh matrix -> Q (cbs x cbs)
		//H is hermitian, so this is a bit more complicated than the normal
		//matrix product formula.  can this be done with FFTs?  maybe worth
		//checking out ...
		for(i = 0u; i < cbs; ++i){
			ii = ACC2(cbs, i, i);
			Q[ii].re = ZERO; Q[ii].im = ZERO;
			///init diags
			for(j = 0u; j < n; ++j){
				ij = ACC2(n, i, j);
				Q[ii].re += (V[ij].re*W[ij].re + V[ij].im*W[ij].im);
				Q[ii].im += (V[ij].re*W[ij].im - V[ij].im*W[ij].re);
			}
			///the rest of the matrix
			for(j = i+1u; j < cbs; ++j){
				ij = ACC2(cbs, i, j);
				Q[ij].re = ZERO; Q[ij].im = ZERO;
				for(k = 0u; k < n; ++k){
					ik = ACC2(n, i, k); jk = ACC2(n, j, k);
					Q[ij].re += (V[ik].re*W[jk].re + V[ik].im*W[jk].im);
					Q[ij].im += (V[ik].re*W[jk].im - V[ik].im*W[jk].re);
				}
				ji = ACC2(cbs, j, i);
				Q[ji].re = Q[ij].re; Q[ji].im = -Q[ij].im;
			}
		}

		//step 3: calculate the eigenpairs of the rayleigh matrix
		//qrh already sorts the vals and vecs
		qrh(cbs, Q, e, L, DEF_QR_ITRLIM, eps*(fp)0.001, (uchar)1);

		//step 4: get ritz vectors (approx. eigvecs of H, in-line w/ evalues)
		for(i = 0u; i < rbs; ++i){
			for(j = 0u; j < n; ++j){
				X[i][j].re = ZERO; X[i][j].im = ZERO;
				for(k = 0u; k < cbs; ++k){
					kj = ACC2(n, k, j);
					X[i][j].re += (L[i][k].re*V[kj].re - L[i][k].im*V[kj].im);
					X[i][j].im += (L[i][k].re*V[kj].im + L[i][k].im*V[kj].re);
				}
			}
		}

		//step 5: convergence check (w/o ritz vectors)
		curre = ZERO;
		for(i = 0u; i < rbs; ++i){
			curre += e[i];
		}
		curre /= (fp)rbs;
		if(cnt){
			diffe = ABSF(laste - curre);
			printf("(step = %04u  dE = %0.2e)\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", cnt, diffe);
			if(diffe < eps) break;
		}
		laste = curre;

		//step 5/6: calculate new search directions -> Q
		//in previous testing, i found that davidson's original preconditioner
		//was much better than, say, the diagonals + the row above/below them
		//(the inverse of which has a closed, relativly easy-to-compute form)
		//this may be different since previous testing was done in python which
		//is much slower than this and also only has double-point prec.  in case
		//of numerical problems, using the tridiagonal preconditioner may be
		//necessary
		///folded spectrum preconditioner (TODO: this is technically incorrect !!! see lobcg implementation)
		if(fsm){
			for(i = 0u; i < rbs; ++i){
				///calculate the this wavefunction's average kinetic energy
				keavg2.re = ZERO; keavg2.im = ZERO;
				for(j = 0u; j < n; ++j){
					keavg2.re += X[i][j].re * (ham.ke[j] + eref);
					keavg2.im += X[i][j].im * (ham.ke[j] + eref);
				}
				keavg2.re = (keavg2.re*keavg2.re + keavg2.im*keavg2.im) / (fp)n;
				keavg2.re = 3.5;
				for(j = 0u; j < n; ++j){
					ij = ACC2(n, i, j);
					Q[ij].re = e[i]*X[i][j].re; Q[ij].im = e[i]*X[i][j].im;
					for(k = 0u; k < cbs; ++k){
						kj = ACC2(n, k, j);
						Q[ij].re -= (L[i][k].re*W[kj].re - L[i][k].im*W[kj].im);
						Q[ij].im -= (L[i][k].re*W[kj].im + L[i][k].im*W[kj].re);
					}
					pre = keavg2.re / ((ham.la[j]+vavgsc)*(ham.la[j]+vavgsc) + 
									    keavg2.re); ///~ 1/(H-eref)^2
					Q[ij].re *= pre;
					Q[ij].im *= pre;
				}
			}
		}
		///normal preconditioner
		else{
			for(i = 0u; i < rbs; ++i){
				for(j = 0u; j < n; ++j){
					ij = ACC2(n, i, j);
					Q[ij].re = e[i]*X[i][j].re; Q[ij].im = e[i]*X[i][j].im;
					for(k = 0u; k < cbs; ++k){
						kj = ACC2(n, k, j);
						Q[ij].re -= (L[i][k].re*W[kj].re - L[i][k].im*W[kj].im);
						Q[ij].im -= (L[i][k].re*W[kj].im + L[i][k].im*W[kj].re);
					}
					Q[ij].re /= (e[i] - ham.ke[j]);
					Q[ij].im /= (e[i] - ham.ke[j]);
				}
			}
		}


		//step 7: increase basis or restart
		if(cbs + rbs < mbs){ 
			///increase basis size by rbs
			for(i = 0u; i < rbs; ++i){
				for(j = 0u; j < n; ++j){
					ijp = ACC2(n, cbs+i, j); ij = ACC2(n, i, j);
					V[ijp].re = Q[ij].re; V[ijp].im = Q[ij].im; ///TODO: might be able to memcpy Q -> V, but it'd be messy 
				}
			}
			_mgs(cbs+rbs, n, &V, &Q);	
		}
		else{ 
			///restart back to 2*rbs
			cbs = rbs; ///loop itr takes care of the rest of the cbs update
			for(i = 0u; i < rbs; ++i){
				for(j = 0u; j < n; ++j){
					ijp = ACC2(n, rbs+i, j); ij = ACC2(n, i, j);
					V[ij].re  = X[i][j].re; V[ij].im  = X[i][j].im;
					V[ijp].re = Q[ij].re;   V[ijp].im = Q[ij].im; ///TODO: might be able to memcpy Q -> V, but it'd be messy 
				}
			}
			_mgs(2u*rbs, n, &V, &Q);
		}
	}


	//recover actual eigenvalues from folded  ones if necessary
	if(fsm){
		for(i = 0u; i < n; ++i) ham.ke[i] += eref; ///undo the fsm		
		_fsmrecover(ham, rbs, n, eref, W, psig1, psig2, e, 
					(const cfp*restrict*restrict)X);
		_epairsort(rbs, &e, &X, (uchar)1);
	}

	if(cnt == itrlim) return (uchar)1;
	return (uchar)0; 
}






#define LOPCG_DEBUG 0  ////////////////////////////////////////////////////////////////////////////////////////

//normalizes imaginary vector v
static forceinline void _normalize(const uint n, cfp*restrict v){
	uint i; fp norm;

	norm = ZERO;
	for(i = 0u; i < n; ++i){
		norm += v[i].re*v[i].re + v[i].im*v[i].im;
	}

	norm = ONE / SQRT(norm);
	for(i = 0u; i < n; ++i){
		v[i].re *= norm; v[i].im *= norm;
	}
}

//conditionally normalizes imaginary vector v - normalization DOES NOT occur
//if limlo < norm(v) < limhi
static forceinline void _cnormalize(const uint n, cfp*restrict v,
									const fp limlo, const fp limhi){
	uint i; fp norm;

	norm = ZERO;
	for(i = 0u; i < n; ++i){
		norm += v[i].re*v[i].re + v[i].im*v[i].im;
	}

	norm = SQRT(norm);
	if(limlo < norm && norm < limhi) return;

	norm = ONE / norm;
	for(i = 0u; i < n; ++i){
		v[i].re *= norm; v[i].im *= norm;
	}
}

//deals with all necessary initialization for lopcg routine
static forceinline void _lcginit(const uchar stab, const uchar fsm,
								 const uint neig, const uint npw,
								 uint*restrict nlocks, uchar*restrict locked,
								 const uint v0r, const uint v0c, 
								 const cfp*restrict*restrict V0, 
								 cfp*restrict Xi, cfp*restrict Xm, 
								 cfp*restrict XmH, cfp*restrict BuH,
								 const hamil ham, 
								 cfp*restrict psig1, cfp*restrict psig2){
	uint i, j, ij;
	uint lim1, lim2;

	//unlock all eigenpairs
	*nlocks = 0u;
	memset(locked, 0, neig*sizeof(uchar));

	//set Xm randomly
	lim1 = neig*npw;
	for(i = 0u; i < lim1; ++i){
		Xm[i].re = (fp)rand() / (fp)RAND_MAX;
		Xm[i].re = (fp)rand() / (fp)RAND_MAX;
	}
	//and Xi from V0 if possible ...
	lim1 = MIN(neig, v0r);
	lim2 = MIN(npw,  v0c);
	for(i = 0u; i < lim1; ++i){
		for(j = 0u; j < lim2; ++j){
			ij = ACC2(npw, i, j);
			Xi[ij].re = V0[i][j].re;
			Xi[ij].im = V0[i][j].im;
		}
	}
	//... and fill the rest randomly
	for(i = lim1; i < neig; ++i){
		for(j = lim2; j < npw; ++j){
			ij = ACC2(npw, i, j);
			Xi[ij].re = (fp)rand() / (fp)RAND_MAX;
			Xi[ij].im = (fp)rand() / (fp)RAND_MAX;
		}
	}
	
#if LOPCG_DEBUG
	//consistent but 'random enough' values for testing
	for(i = 0u; i < neig; ++i){
		for(j = 0u; j < npw; ++j){
			fp fac = (fp)(npw*i) / (fp)(neig-1u);
			ij = ACC2(npw, i, j);
			Xi[ij].re = ONE / ABSF((fp)0.1 + (fp)j - fac);
			Xi[ij].im = fac;
			Xm[ij].re = ABSF((fp)(i+j) - fac);
			Xm[ij].im = ONE / (fac + (fp)1);
		}
	}
#endif

	lim1 = neig*npw;
	for(i = 0u; i < lim1; i += npw) _normalize(npw, Xm + i);
	for(i = 0u; i < lim1; i += npw) _normalize(npw, Xi + i);

	//now, if we're feeling dangerous, we can go with no orthogonalization
	//(stability = 0).  in this case, we have to initialize <psi|H*
	if(!stab){
		if(fsm) for(i = 0u; i < lim1; i += npw){
			_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, npw, 
			ham.mills, ham.ke, Xm + i, BuH + i);
			_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, npw, 
			ham.mills, ham.ke, BuH + i, XmH + i);
		}
		else for(i = 0u; i < lim1; i += npw){
			_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, npw, 
			ham.mills, ham.ke, Xm + i, XmH + i);
		} 
	}	
}

//set the 'normal' preconditioner for LCG.  this specific form is from
//Teter, Payne, Allan - "Solution of Schrodinger's Equation for large systems"
// - PRB 1989
static forceinline void _lcgpdia(const uint npw, 
								 const fp*restrict ke, const cfp*restrict psi, 
								 fp*restrict T){
	uint i;
	fp ttot, x, x3, fac;

	//get total wavefunction kinetic energy
	ttot = ZERO;
	for(i = 0u; i < npw; ++i){
		ttot += ke[i] * (psi[i].re*psi[i].re + psi[i].im*psi[i].im);
	}
	ttot = ONE / ttot;

	//T[i] = 27+18x+12x^2+8x^3 / 27+18x+12x^2+8x^3+16x^4 ; x = ke_i / ke_tot
	for(i = 0u; i < npw; ++i){
		x = ke[i]*ttot;
		x3 = x*x*x;
		fac = (fp)27 + ((fp)18 + (fp)12*x)*x + (fp)8*x3;
		T[i] = fac / (fac + (fp)16*x*x3);
	}
}

//set the 'fsm' preconditioner for LCG
//the kinetic energy input should already be scaled by eref 
//(kesc = ke_orig - eref), as should vavgsc.  la are hbar^2/2m* * |G|^2
//in addition to the preconditioner of TPA (see above), the other option (that 
//i'll probably get rid of b/c it seems to suck) is from 
//Wang, Zunger = "Solving Schrodinger's Equation around a desired energy: 
//application to Si QDs" - Jrn. of Chem Phys, 1994 
//with alpha = <T>
#define PRECTEST 1
#if PRECTEST
static forceinline void _lcgpfsm(const uint npw, const fp*restrict kesc, 
								 const cfp*restrict psi, fp*restrict T){	
#else
static forceinline void _lcgpfsm(const uint npw, const fp eref, 
								 const fp*restrict kesc, const fp*restrict la,
								 const fp vavgsc, const cfp*restrict psi,
								 fp*restrict T){
#endif
#if PRECTEST ///why does this work so well vs. the original??
	uint i;
	fp ttot, x, x3, fac;

	//get total wavefunction kinetic energy
	ttot = ZERO;
	for(i = 0u; i < npw; ++i){
		ttot += kesc[i]*kesc[i] * (psi[i].re*psi[i].re + psi[i].im*psi[i].im);
	}

	ttot = ONE / ttot;
	//T[i] = 27+18x+12x^2+8x^3 / 27+18x+12x^2+8x^3+16x^4 ; x = ke_i / ke_tot
	for(i = 0u; i < npw; ++i){
		x = kesc[i]*kesc[i]*ttot;
		x3 = x*x*x;
		fac = (fp)27 + ((fp)18 + (fp)12*x)*x + (fp)8*x3;
		T[i] = fac / (fac + (fp)16*x*x3);
	}
#else
	uint i;
	fp tavg2;

	//get average wavefunction kinetic energy (squared)
	tavg2 = ZERO;
	for(i = 0u; i < npw; ++i){
		tavg2 += (kesc[i] + eref) * (psi[i].re*psi[i].re + psi[i].im*psi[i].im);
	}
	tavg2 /= (fp)npw; 
	tavg2 *= tavg2;

	//T[i] ~= K^2 / ((G^2 + v0 - eref)^2 + K^2)
	for(i = 0u; i < npw; ++i){
		T[i] = tavg2 / ((la[i]+vavgsc)*(la[i]+vavgsc) + tavg2);
	}
#endif
}

//grahm schmidt orthogonalization for the triplet of basis vectors in LOPCG
//the new search directions (Wi) are NOT normalized, the current
//eigenvector approximations (Xi) are OPTIONALLY normalized (SAFE_LCG = 1 to do
//so), and the previous approximations are OPTIONALLY normalized as well 
//(SAFE_LCG = 2 for both normalizations)
static forceinline void _lcgogsh(const uint n, const cfp*restrict orig, 
								 const cfp*restrict agst, cfp*restrict res){
	uint i;
	fp facre, facim;

	facre = ZERO; facim = ZERO;
	for(i = 0u; i < n; ++i){ //fac = (agst)* @ orig
		facre += agst[i].re*orig[i].re + agst[i].im*orig[i].im; 
		facim += agst[i].re*orig[i].im - agst[i].im*orig[i].re;
	}
	for(i = 0u; i < n; ++i){ //res -= fac*agst
		res[i].re -= (facre*agst[i].re - facim*agst[i].im);
		res[i].im -= (facre*agst[i].im + facim*agst[i].re);
	}
}
static void _lcgogs(const uint neig, const uint npw, 
					cfp*restrict wi, cfp*restrict xi, cfp*restrict xm,
					cfp*restrict wib, cfp*restrict xib, cfp*restrict xmb,
					const uchar*restrict locked){
	uint i, j, ii, jj;

	memcpy(wib, wi, neig*npw*sizeof(cfp));
	memcpy(xib, xi, neig*npw*sizeof(cfp));
	memcpy(xmb, xm, neig*npw*sizeof(cfp));

	//ortho Wi
	for(i = 0u, ii = 0u; i < neig; ++i, ii += npw){ 
		if(locked[i]) continue;
		///against locked vectors first
		for(j = 0u, jj = 0u; j < neig; ++j, jj += npw){
			if(!locked[j]) continue;
			_lcgogsh(npw, wib+ii, wib+jj, wi+ii);
 			_lcgogsh(npw, wib+ii, xib+jj, wi+ii);
			_lcgogsh(npw, wib+ii, xmb+jj, wi+ii);
		}
		///now against previous wi
		for(j = 0u, jj = 0u; j < i; ++j, jj += npw){
			if(locked[j]) continue;
			_lcgogsh(npw, wib+ii, wi+jj, wi+ii);
		}
#if (SAFE_LCG == 1)
		_cnormalize(npw, wi+ii, DEF_CG_CNLIMLO, DEF_CG_CNLIMHI);
#elif (SAFE_LCG > 1)
		_normalize(npw, wi+ii);
#endif
	}

	//ortho Xi
	for(i = 0u, ii = 0u; i < neig; ++i, ii += npw){ 
		if(locked[i]) continue;
		///against locked vectors first
		for(j = 0u, jj = 0u; j < neig; ++j, jj += npw){
			if(!locked[j]) continue;
			_lcgogsh(npw, xib+ii, wib+jj, xi+ii);
 			_lcgogsh(npw, xib+ii, xib+jj, xi+ii);
			_lcgogsh(npw, xib+ii, xmb+jj, xi+ii);
		}
		///now against previous wi, xi
		for(j = 0u, jj = 0u; j < neig; ++j, jj += npw){
			if(locked[j]) continue;
			_lcgogsh(npw, xib+ii, wi+jj, xi+ii);
		}
		for(j = 0u, jj = 0u; j < i; ++j, jj += npw){
			if(locked[j]) continue;
			_lcgogsh(npw, xib+ii, xi+jj, xi+ii);
		}
#if SAFE_LCG
		_normalize(npw, xi+ii); 
#endif		
	}

	//ortho Xm
	for(i = 0u, ii = 0u; i < neig; ++i, ii += npw){ 
		if(locked[i]) continue;
		///against locked vectors first
		for(j = 0u, jj = 0u; j < neig; ++j, jj += npw){
			if(!locked[j]) continue;
			_lcgogsh(npw, xmb+ii, wib+jj, xm+ii);
 			_lcgogsh(npw, xmb+ii, xib+jj, xm+ii);
			_lcgogsh(npw, xmb+ii, xmb+jj, xm+ii);
		}
		///now against previous wi, xi, xm
		for(j = 0u, jj = 0u; j < neig; ++j, jj += npw){
			if(locked[j]) continue;
			_lcgogsh(npw, xmb+ii, wi+jj, xm+ii);
			_lcgogsh(npw, xmb+ii, xi+jj, xm+ii);
		}
		for(j = 0u, jj = 0u; j < i; ++j, jj += npw){
			if(locked[j]) continue;
			_lcgogsh(npw, xmb+ii, xm+jj, xm+ii);
		}
#if (SAFE_LCG > 1)
		_normalize(npw, xm+ii); 
#endif		
	}
}

//rayleigh-ritz procedure on {Wi, Xi, Xm}
static forceinline void _lcgritzh(const uint n, 
								  const cfp*restrict a, const cfp*restrict b,
								  cfp*restrict res){
	uint i;
	fp rre, rim; ///gcc 4.x WON'T VECTORIZE THIS LOOP if 'res' is used or 
				 ///these local variables are replaced with a single 'cfp' var
	rre = ZERO; rim = ZERO;
	for(i = 0u; i < n; ++i){ ///a* * b
		rre += (a[i].re*b[i].re + a[i].im*b[i].im);
		rim += (a[i].re*b[i].im - a[i].im*b[i].re);
	}

	res->re = rre;
	res->im = rim;
}
static void _lcgritz(const uint neig, const uint npw, 
					 const cfp*restrict Wi, const cfp*restrict WiH, 
					 const cfp*restrict Xi, const cfp*restrict XiH,
					 const cfp*restrict Xm, const cfp*restrict XmH,
					 const uint nlocked, uchar*restrict lockedu, 
					 fp*restrict lockedf,
					 cfp*restrict rH, cfp*restrict rS, fp*restrict re,
					 cfp*restrict*restrict rX, 
					 fp*restrict e, cfp*restrict*restrict X, 
					 const fp qreps){
	uint i, j, k, ii, jj, ij, ik, jk;
	uint c0, c1;
	uint n1, n2, n3;
	cfp su;	

	//create H*, S.  this has a fairly awful access pattern thanks to the 
	//seperation of the basis into Wi, Xi, and Xm.  this is also fairly 
	//complex thanks to the locks and my refusal to not take advantage of the
	//fact that H* and S are hermitian
	n1 = neig - nlocked; n2 = 2u*n1; n3 = 3u*n1;
	c0 = 0u;
	for(i = 0u, ii = 0u; i < neig; ++i, ii += npw){ 
		if(lockedu[i]) continue;
		///WiH* @ (Wi, Xi, Xm) -> H* & Wi* @ (Wi, Xi, Xm) -> S
		for(j = 0u, jj = 0u, c1 = 0u; j < neig; ++j, jj += npw){ 
			if(lockedu[j]) continue;
			if(c1    >= c0   ){
				ij = ACC2(n3, c0, c1   );
				_lcgritzh(npw, WiH+ii, Wi+jj, rH+ij);
				_lcgritzh(npw, Wi+ii , Wi+jj, rS+ij);
			}
			if(c1+n1 >= c0   ){
				ij = ACC2(n3, c0, c1+n1);
				_lcgritzh(npw, WiH+ii, Xi+jj, rH+ij);
				_lcgritzh(npw, Wi+ii , Xi+jj, rS+ij);
			}
			if(c1+n2 >= c0   ){
				ij = ACC2(n3, c0, c1+n2);
				_lcgritzh(npw, WiH+ii, Xm+jj, rH+ij);
				_lcgritzh(npw, Wi+ii , Xm+jj, rS+ij);
			}
			c1++;
		}
		///XiH* @ (Wi, Xi, Xm) -> H* & Xi* @ (Wi, Xi, Xm) -> S
		for(j = 0u, jj = 0u, c1 = 0u; j < neig; ++j, jj += npw){ 
			if(lockedu[j]) continue;
			if(c1    >= c0+n1){
				ij = ACC2(n3, c0+n1, c1   );
				_lcgritzh(npw, XiH+ii, Wi+jj, rH+ij);
				_lcgritzh(npw, Xi+ii , Wi+jj, rS+ij);
			}
			if(c1+n1 >= c0+n1){
				ij = ACC2(n3, c0+n1, c1+n1);
				_lcgritzh(npw, XiH+ii, Xi+jj, rH+ij);
				_lcgritzh(npw, Xi+ii , Xi+jj, rS+ij);
			}
			if(c1+n2 >= c0+n1){
				ij = ACC2(n3, c0+n1, c1+n2);
				_lcgritzh(npw, XiH+ii, Xm+jj, rH+ij);
				_lcgritzh(npw, Xi+ii , Xm+jj, rS+ij);
			}
			c1++;
		}
		///XmH* @ (Wi, Xi, Xm) -> H* & Xm* @ (Wi, Xi, Xm) -> S
		for(j = 0u, jj = 0u, c1 = 0u; j < neig; ++j, jj += npw){ 
			if(lockedu[j]) continue;
			if(c1    >= c0+n2){
				ij = ACC2(n3, c0+n2, c1   );
				_lcgritzh(npw, XmH+ii, Wi+jj, rH+ij);
				_lcgritzh(npw, Xm+ii , Wi+jj, rS+ij);
			}
			if(c1+n1 >= c0+n2){
				ij = ACC2(n3, c0+n2, c1+n1);
				_lcgritzh(npw, XmH+ii, Xi+jj, rH+ij);
				_lcgritzh(npw, Xm+ii , Xi+jj, rS+ij);
			}
			if(c1+n2 >= c0+n2){
				ij = ACC2(n3, c0+n2, c1+n2);
				_lcgritzh(npw, XmH+ii, Xm+jj, rH+ij);
				_lcgritzh(npw, Xm+ii , Xm+jj, rS+ij);
			}
			c1++;
		}
		c0++;
	}
	///finish the unfilled halfs
	for(i = 0u, ii = 0u; i < n3; ++i, ii += n3){
		for(j = 0u; j < i; ++j){
			ij = ii + j;
			jj = ACC2(n3, j, i);
			rH[ij].re = rH[jj].re; rH[ij].im = -rH[jj].im;
			rS[ij].re = rS[jj].re; rS[ij].im = -rS[jj].im;
		}
	}

//    for(int h = 0; h < neig*npw; ++h) printf("%i "FPFORMAT" "FPFORMAT"\n",h, Wi[h].re, Wi[h].im);
//    for(int h = 0; h < n3*n3; ++h) printf("%i "FPFORMAT" "FPFORMAT"\n",h, rH[h].re, rH[h].im);
	//transform the generalized eigenproblem Rh*rX = re*rS*rX into a simplified
	//eigenproblem.  See netlib.org/lapack/lug/node54.html for details
	///rS is herm. pos definite: rS - >LL+ (cholesky)
	for(i = 0u, ii = 0u; i < n3; ++i, ii += n3+1u){
		for(j = 0u, jj = 0u; j < i; ++j, jj += n3+1u){
			su.re = ZERO; su.im = ZERO;
			for(k = 0u; k < j; ++k){
				ik = ACC2(n3, i, k); jk = ACC2(n3, j, k);
				su.re += (rS[ik].re*rS[jk].re + rS[ik].im*rS[jk].im);
				su.im += (rS[ik].im*rS[jk].re - rS[ik].re*rS[jk].im); 
			}
			ij = ACC2(n3, i, j);
			rS[ij].re = (rS[ij].re - su.re) / rS[jj].re;
			rS[ij].im = (rS[ij].im - su.im) / rS[jj].re;
		}
		///this next section is pure real
		su.re = ZERO;
		for(j = 0u; j < i; ++j){
			ij = ACC2(n3, i, j);
			su.re += rS[ij].re*rS[ij].re + rS[ij].im*rS[ij].im;
		}
if(rS[ii].re - su.re < ZERO) printf("god damn it "FPFORMAT" "FPFORMAT"\n", rS[ii].re , su.re);
		rS[ii].re = SQRT(rS[ii].re - su.re);
		rS[ii].im = ZERO; ///may be unnecessary
	}
	///rH -> X, but only 1/2 of X is needed (solve L@X = rH for X)
	for(i = 0u, ii = 0u; i < n3; ++i, ii += n3+1u){
		for(j = i; j < n3; ++j){
			ij = ACC2(n3, i, j);
			for(k = 0u; k < i; ++k){ ///S(ik) * H(kj)
				ik = ACC2(n3, i, k); jk = ACC2(n3, k, j); 
				rH[ij].re -= (rS[ik].re*rH[jk].re - rS[ik].im*rH[jk].im);
				rH[ij].im -= (rS[ik].im*rH[jk].re + rS[ik].re*rH[jk].im);
			}
			rH[ij].re /= rS[ii].re;
			rH[ij].im /= rS[ii].re;
		}
	}
	///C is hermitian: X -> C (solve C@L+ = X for C)
	for(i = 0u, ii = 0u; i < n3; ++i, ii += n3+1u){
		///this first section is pure real
		for(j = 0u; j < i; ++j){
			ij = ACC2(n3, i, j);
			rH[ii].re -= rH[ij].re*rS[ij].re + rH[ij].im*rS[ij].im;
		}
		rH[ii].re /= rS[ii].re;
		rH[ii].im = ZERO; ///may be unnecessary
		///and now the general case
		for(j = i+1u; j < n3; ++j){
			ij = ACC2(n3, i, j);
			for(k = 0u; k < j; ++k){ ///H(ij) -= H(ik) * S*(jk)
				ik = ACC2(n3, i, k); jk = ACC2(n3, j, k);
				rH[ij].re -= (rH[ik].re*rS[jk].re + rH[ik].im*rS[jk].im);
				rH[ij].im -= (rH[ik].im*rS[jk].re - rH[ik].re*rS[jk].im);
			}
			jj = ACC2(n3, j, j);
			rH[ij].re /= rS[jj].re; rH[ij].im /= rS[jj].re;
			ik = ACC2(n3, j, i);
			rH[ik].re = rH[ij].re; rH[ik].im = -rH[ij].im; ///hermitian property
		}
	}
	///now, the eigenvalue problem is transformed.  solve it ... 
//for(int h = 0; h < n3*n3; ++h) printf("%i "FPFORMAT" "FPFORMAT"\n",h, rH[h].re, rH[h].im);
	uchar ghahah = qrh(n3, rH, re, rX, DEF_QR_ITRLIM, qreps, (uchar)1);
	///... but the vectors are still transformed.  need to get them back:
	for(i = n3; i --> 0u ;){ 	  ///unsigned ints and counting backwards,
		for(j = n3; j --> 0u ;){  ///what could go wrong?
			for(k = j+1u; k < n3; ++k){ ///X(ij) -= S*(kj) * X(ik)
				ik = ACC2(n3, i, k); jk = ACC2(n3, k, j);
				rX[i][j].re -= (rS[jk].re*rX[i][k].re + rS[jk].im*rX[i][k].im);
				rX[i][j].im -= (rS[jk].re*rX[i][k].im - rS[jk].im*rX[i][k].re);
			}
			jj = ACC2(n3, j, j);
			rX[i][j].re /= rS[jj].re; rX[i][j].im /= rS[jj].re;
		}							  
	}

	//update the eigenvalues that we're interested in
	for(i = 0u, c0 = 0u; i < neig; ++i){
		if(lockedu[i]) continue;
		e[i] = re[c0];
		c0++;
	}

	//compute X = rX @ (Wi, Xi, Xm) to recover the eigenvector approximations 
	//from ritz vectors ... but theres a problem:  we have to move down the COLS
	//of (.), but the ROWS are locked i.e. we'd normally need an if statement in 
	//the innermost loop.  Instead, doing the dot product as if there were no 
	//locks, but multiplying zero to the locked vectors is more efficient
	for(i = 0u; i < neig; ++i) lockedu[i] = (uchar)(!lockedu[i]); ///lck->unlck
	for(i = 0u; i < neig; ++i) lockedf[i] = (fp)(lockedu[i]); ///1 or 0
	for(i = 0u, c0 = 0u; i < neig; ++i){
		if(!lockedu[i]) continue;
		for(j = 0u; j < npw; ++j){
			c1 = 0u;
			X[i][j].re = ZERO; X[i][j].im = ZERO;
			for(k = 0u; c1 < n1; ++k){
				ii = ACC2(npw, k, j);
				//X[i][j] += rX[c0][c1] * lockedf[k]*Wi[ACC2(npw, k, j)]
				X[i][j].re += lockedf[k]*(rX[c0][c1].re*Wi[ii].re - 
										  rX[c0][c1].im*Wi[ii].im);
				X[i][j].im += lockedf[k]*(rX[c0][c1].re*Wi[ii].im + 
										  rX[c0][c1].im*Wi[ii].re);
				c1 += (uint)lockedu[k];
			}
			for(k = 0u; c1 < n2; ++k){
				ii = ACC2(npw, k, j);
				//X[i][j] += rX[c0][c1] * lockedf[k]*Xi[ACC2(npw, k, j)]
				X[i][j].re += lockedf[k]*(rX[c0][c1].re*Xi[ii].re - 
										  rX[c0][c1].im*Xi[ii].im);
				X[i][j].im += lockedf[k]*(rX[c0][c1].re*Xi[ii].im + 
										  rX[c0][c1].im*Xi[ii].re);
				c1 += (uint)lockedu[k];
			}
			for(k = 0u; c1 < n3; ++k){
				ii = ACC2(npw, k, j);
				//X[i][j] += rX[c0][c1] * lockedf[k]*Xm[ACC2(npw, k, j)]
				X[i][j].re += lockedf[k]*(rX[c0][c1].re*Xm[ii].re - 
										  rX[c0][c1].im*Xm[ii].im);
				X[i][j].im += lockedf[k]*(rX[c0][c1].re*Xm[ii].im + 
										  rX[c0][c1].im*Xm[ii].re);
				c1 += (uint)lockedu[k];
			}
		}
		c0++;
	}
	for(i = 0u; i < neig; ++i) lockedu[i] = (uchar)(!lockedu[i]); ///unlck->lck
}

//update convergence info for lopcg.  returns 1 if all eigenvalues are converged
//to tol
static int _lcgupdateconv(const uint neigs, const fp*restrict ecurr,
						   fp*restrict eprev, fp*restrict ediff, 
						   uint*restrict nlocked, uchar*restrict locked,
						   const uchar stab, const fp tol){
	uint i;
	int ret;
	fp avgdiffs;

	//update lockes
	avgdiffs = ZERO;
	for(i = 0u; i < neigs; ++i){
		if(locked[i]) continue;
		ediff[i] = ABSF(ecurr[i] - eprev[i]);
		if(ediff[i] < tol){
			(*nlocked)++;
			locked[i] = (uchar)1;
		}
		else{
			avgdiffs += ediff[i];
		}
	}

	//update previous eigenvalue approximations
	for(i = 0u; i < neigs; ++i) eprev[i] = ecurr[i];

	printf(" ncon: %u / %u ... dE(uncon) = "FPFORMAT" eV\n", *nlocked, neigs, avgdiffs/(fp)(neigs - *nlocked));
	ret = (neigs == *nlocked);

	//if we're trying to be very careful, we leave all vectors unlocked.  here,
	//just reset them all to zero
	if(stab == (uchar)2){
		*nlocked = 0u;
		memset(locked, 0, neigs*sizeof(uchar));
	}

	return ret;
}

//locally optimal PCG eigenpair solver for general hermitian matrices, 
//specialized for electronic structure calculations (i.e. input is a hamiltonian 
//constructor as opposed to an explicit matrix)
//too many parameters to write out here, read the comments in the source
//code
//the following must be allocated prior to use:
// ---workspace---
//  (neigs) * sizeof(uchar)                                           b for IWS
//  (6*neig + npw) * sizeof(fp)                                       b for RWS
//  (18*neig*neig + 6*neig*npw + npw + 2*fftgridsize) * sizeof(cfp)   b for CWS
//  (3*neigs) * sizeof(cfp) x (3*neigs) * sizeof(cfp)                 b for CCWS
// ---return values---
//  (neigs) * sizeof(fp)                                              b for e
//  (neigs) * sizeof(cfp) x (npw) * sizeof(cfp)                       b for X
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
		  const uchar stab){

	uint i, j, ii, ij, cnt;
	fp vavgsc, lam;
	cfp*restrict Wi,  *restrict Xi,  *restrict Xm;
	cfp*restrict WiH, *restrict XiH, *restrict XmH, *restrict BuH;
	cfp*restrict psigrid1, *restrict psigrid2;
	fp*restrict T;
	cfp*restrict HR, *restrict SR; fp*restrict er; cfp*restrict*restrict vr;
	uint nlocks; uchar*restrict lockedu; fp*restrict lockedf;
	fp*restrict eprev, *restrict ediff; 

#define ACTIONNRM(PSI, PSIH) do{_action(ham.dims, ham.l2dims, ham.vloc,\
										psigrid1, psigrid2, n, ham.mills,\
										ham.ke, PSI, PSIH);} while(0)
#define ACTIONFSM(PSI, PSIH) do{_action(ham.dims, ham.l2dims, ham.vloc,\
										psigrid1, psigrid2, n, ham.mills,\
										ham.ke, PSI, BuH);\
								_action(ham.dims, ham.l2dims, ham.vloc,\
										psigrid1, psigrid2, n, ham.mills,\
										ham.ke, BuH, PSIH);} while(0)

	//workspace setup
	ij = neig*n;
	///subspace basis vectors: cplx wkspce 0 -> 3*neig*npw
	Wi = CWS + 0; Xi = Wi + ij; Xm = Xi + ij;
	///action of hamiltonian on wfnctns: += 3*neig*npw + npw
	WiH = Xm + ij; XiH = WiH + ij; XmH = XiH + ij; BuH = XmH + ij;
	///fft grids: cplx wkspce += 2*fftgridsize[0]*[1]*[2]
	ij = (uint)(ham.dims[0]) * (uint)(ham.dims[1]) * (uint)(ham.dims[2]);
	psigrid1 = BuH + n; psigrid2 = psigrid1 + ij;
	///ritz eigenpair stuff: cplx wkspce += 2*9*neig*neig
	HR = psigrid2 + ij; SR = HR + 9u*neig*neig;
	///ritz eigenpair stuff: real wkspce 0 -> 3*neig
	er = RWS + 0; vr = CCWS;
	///preconditioner: real workspace += npw
	T = er + 3u*neig;
	///convergence checking: real wkspce += 2*neig
	eprev = T + n; ediff = eprev + neig;	
	///locking: real wkspce += neig, uchar wkspce 0 -> neig
	lockedf = ediff + neig; lockedu = IWS + 0;


	//init (set locks, basis vectors, apply hamiltonian when needed)
	//plus a bit of exra setup is required if we're doing fsm
	if(fsm){
		for(i = 0u; i < n; ++i) ham.ke[i] -= eref;
#if (!PRECTEST)
		vavgsc = ZERO;
		ij = (uint)(ham.dims[0]) * (uint)(ham.dims[1]) * (uint)(ham.dims[2]);
		for(i = 0u; i < ij; ++i) vavgsc += ham.vloc[i].re;
		vavgsc = vavgsc/(fp)ij - eref;
#endif
	}
	_lcginit(stab, fsm, neig, n, &nlocks, lockedu, v0rs, v0cs, V0, 
			 Xi, Xm, XmH, BuH, ham, psigrid1, psigrid2);


	//main loop
	for(cnt = 0u; cnt < itrlim; ++cnt){
		
		//new search directions (aka residuals) Wi
		for(i = 0u, ii = 0u; i < neig; ++i, ii += n){
			if(lockedu[i]) continue;
			if(fsm){
				ACTIONFSM(Xi+ii, XiH+ii);
#if PRECTEST
				_lcgpfsm(n, ham.ke, Xi+ii, T);
#else
				_lcgpfsm(n, eref, ham.ke, ham.la, vavgsc, Xi+ii, T);
#endif
			} 
			else{
				ACTIONNRM(Xi+ii, XiH+ii);
				_lcgpdia(n, ham.ke, Xi+ii, T);
			}

			lam = ZERO; /// = <psi|H|psi>
			for(ij = ii; ij < ii+n; ++ij){
				lam += XiH[ij].re*Xi[ij].re + XiH[ij].im*Xi[ij].im;
			}
			
			lam = ONE / lam;
			for(j = 0u, ij = ii; j < n; ++j, ++ij){
				Wi[ij].re = T[j] * (Xi[ij].re - XiH[ij].re*lam);
				Wi[ij].im = T[j] * (Xi[ij].im - XiH[ij].im*lam);
			}
		}

		//init action on trial psi after orthogonalization (if necessary)
		if(stab){
			_lcgogs(neig, n, Wi, Xi, Xm, WiH, XiH, XmH, lockedu);
			for(i = 0u, ii = 0u; i < neig; ++i, ii += n){
				if(lockedu[i]) continue;
				if(fsm){
					ACTIONFSM(Wi+ii, WiH+ii); 
					ACTIONFSM(Xi+ii, XiH+ii);
					ACTIONFSM(Xm+ii, XmH+ii);	
				}
				else{
					ACTIONNRM(Wi+ii, WiH+ii); 
					ACTIONNRM(Xi+ii, XiH+ii);
					ACTIONNRM(Xm+ii, XmH+ii);
				}
			}
		}
		else{ ///unstable, only need action on W
			for(i = 0u, ii = 0u; i < neig; ++i, ii += n){
				if(lockedu[i]) continue;
				if(fsm) ACTIONFSM(Wi+ii, WiH+ii); 
				else    ACTIONNRM(Wi+ii, WiH+ii); 
			}
		}

		//ray-ritz on the trial subspaces
		_lcgritz(neig, n, Wi, WiH, Xi, XiH, Xm, XmH, nlocks, lockedu, lockedf,
				 HR, SR, er, vr, e, X, /*DEF_QR_EPS*/(fp)0.001*eps);

		//with no orthog, we can skip an application of <psi|H*
		if(!stab){
			for(i = 0u, ii = 0u; i < neig; ++i, ii += n){
				if(lockedu[i]) continue;
				memcpy(XmH+ii, XiH+ii, n*sizeof(cfp));
			}
		}

		printf("step %u / %u :", cnt, itrlim);
		//termination conditions, update locks
		if(_lcgupdateconv(neig, e, eprev, ediff, &nlocks, lockedu, stab, eps))
			break;

		//prep for next iteration
		for(i = 0u, ii = 0u; i < neig; ++i, ii += n){
			if(lockedu[i]) continue;
			memcpy(Xm+ii, Xi+ii, n*sizeof(cfp));
			for(j = 0u, ij = ii; j < n; ++j, ++ij){
				Xi[ij].re = X[i][j].re; Xi[ij].im = X[i][j].im;
			}
		}
	}

	//last bit: if we did fsm, we need to recover the original eigenpairs
	if(fsm){
		for(i = 0u; i < n; ++i) ham.ke[i] += eref; ///undo the fsm		
		_fsmrecover(ham, neig, n, eref, Wi, psigrid1, psigrid2, e, 
					(const cfp*restrict*restrict)X);
		_epairsort(neig, &e, &X, (uchar)1);
	}	

#undef ACTIONNRM
#undef ACTIONFSM
	return (uchar)(cnt == itrlim);
}







