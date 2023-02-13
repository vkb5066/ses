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
		  const uint itrlim, const fp eps, uchar vecs){
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
	//memset((cfp**)psig, 0, tlen*sizeof(cfp*)); ///needed since we reuse psig figure out how to use this
	for(i = 0u; i < tlen; ++i){
		psig[i].re = ZERO; psig[i].im = ZERO;
buffft[i].re = ZERO; buffft[i].im = ZERO;
	}

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

//llots of restricts should be able to be taken out of here ... what the fuck is even going on
uchar dav(const uint n, const hamil ham,  	 
		  //buffer params: V, Q, W = mbs x npw, psig* = fftdim[0]*[1]*[2]
		  //               L = pointer-to-pointer mbs x npw
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


	//un-fixable errors (if any of these throw, work arrays are too small!)
	if(n < rbs) return (uchar)1;
	if(mbs < 2u*rbs) return(uchar)2;
	if(n <= mbs) return (uchar)3;


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
			//V[ij].re = (fp)0;  //                          make sure you delete these two lines lololololololololololololloololololollololo
			//V[ij].im = (fp)0; //
		}
		//V[ACC2(n, i, i)].re = ONE; V[ACC2(n, i, i)].im = ZERO; ////and meeeeeeeeeeeeeeeeeeeeeeeee
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
	//operate I @ eref on hamiltonian if we're doing folded spectrum
	if(fsm){
		for(i = 0u; i < n; ++i) ham.ke[i] -= eref;
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
						ham.mills, ham.ke, 
						V + ij, Q + ij);
				_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, n, 
						ham.mills, ham.ke, 
						Q + ij, W + ij);
			}
		}
		else{
			///this is the 'normal' way of applying <psi| to H, nothing fancy
			for(i = 0u; i < cbs; ++i){
				ij = i*n;
				_action(ham.dims, ham.l2dims, ham.vloc, psig1, psig2, n, 
						ham.mills, ham.ke, 
						V + ij, W + ij);
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
		qrh(cbs, Q, e, L, QR_ITR_LIM, 1.0e-12, (uchar)1); ///!!!swap def_qr_eps to use defined eps / 100 or so after testing finished 

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
			diffe = laste - curre;
//			printf("%u  "FPFORMAT"\n", cnt, diffe);
			if(ABSF(diffe) < eps) break;
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

	if(cnt == itrlim) return (uchar)1;

	//--------------------------//
	//----- FSM  TREATMENT -----//
	//--------------------------//	
	//if we're doing the FSM, the eigenvalues we have now obey e = (lam - ref)^2
	//and we want to find lam (note the eigenvectors are unaffected).  
	//i have yet to come up with a more efficient solutions, so we just check
	//H|psi> - lam_trial|psi> for both solutions to the above equation.
	//the lam_trial that gives the closest approximation to the zero vector is
	//what we're after.  
	if(fsm){
		for(i = 0u; i < n; ++i) ham.ke[i] += eref; ///undo the fsm
		
	}

	return (uchar)0;
}





