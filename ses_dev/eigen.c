#include <stdlib.h>
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

//TODO: QRHERM CAN BE CLEANED UP A BIT:
//GETTING RID OF A FEW INITIAL VARIABLES

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
uint qrh(const uint n, cfp*restrict*restrict A, 
		 fp*restrict*restrict e, cfp*restrict*restrict*restrict V,
		 const uint itrLim, const fp eps, uchar vecs){
//This is a re-implementation of an algorithm from 'mpmath':
//https://github.com/fredrik-johansson/mpmath
//which is itself based on a fortran implementation in 'EISPACK'.    
//This specific version has been edited to make the most accesses in row-order 
//instead of column order.
#define A (*A)
#define V (*V)
	uint i, j, k, l, m;
	uint i0, i1, i2;	

	fp fp1, fp2, fp3, fp4, fp5, fp6, fp7;
	fp dg; cfp cg, cfp1, f0;
	fp* d; cfp* t;

	memset(*e, 0, n*sizeof(fp)); ///evals -> zero for iterative updates 
	d = calloc(n, sizeof(fp));
	t = calloc(n, sizeof(cfp));

	//--------------------------//
	//--- TRIDIAGONALIZATION ---//
	//--------------------------//
	t[n - 1u].re = ONE; t[n - 1u].im = ZERO;
	for(i = n - 1; i > 0; i--){
		l = i - 1;

		///vector scaling
		fp1 = ZERO;
		for(j = 0; j < i; ++j){
			i0 = ACC2(n, i, j);
			fp1 += ABSF(A[i0].re) + ABSF(A[i0].im);
		}

		///skipping, stopping criteria
		if(ABSF(fp1) < eps){
			d[i] = ZERO;
			(*e)[i] = ZERO;
			t[l].re = ONE; t[l].im = ZERO;
			continue;
		}

		i0 = ACC2(n, i, l);
		if(i == 1){
			f0.re = A[i0].re; f0.im = A[i0].im; 
			fp3 = SQRT(f0.re*f0.re + f0.im*f0.im);
			fp2 = ONE / fp3;
			d[i] = fp3;
			(*e)[i] = ZERO;
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
		fp2 = ONE/fp1; fp4= ZERO;
		for(j = 0; j < i; ++j){
			i1 = ACC2(n, i, j);
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
		for(j = 0; j < i; ++j){
			i0 = ACC2(n, i, j);
			k = ACC2(n, j, i);
			A[k].re = A[i0].re / fp4; A[k].im = A[i0].im / fp4;

			////A dot U
			cg.re = ZERO; cg.im = ZERO;
			for(k = 0; k < j + 1; ++k){
				/////g += A(j, k)*  dot  A(i, k)
				i1 = ACC2(n, j, k); i2 = ACC2(n, i, k);
				cg.re += A[i1].re*A[i2].re + A[i1].im*A[i2].im;
				cg.im += A[i1].re*A[i2].im - A[i1].im*A[i2].re;
			}
			for(k = j + 1; k < i; ++k){
				/////g += A(k, j)  dot  A(i, k)
				i1 = ACC2(n, k, j); i2 = ACC2(n, i, k);
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
		for(j = 0; j < i; ++j){
			i0 = ACC2(n, i, j);
			f0.re = A[i0].re; f0.im = A[i0].im;
			cg.re = t[j].re - fp5*f0.re; cg.im = t[j].im - fp5*f0.im;
			t[j].re = cg.re; t[j].im = cg.im;

			////this part is expensive
			for(k = 0; k < j + 1; k++){
				/////A(j, k) -= f0*  dot  t[k]    +    g*  dot  A(i, k)
				i0 = ACC2(n, j, k); i1 = ACC2(n, i, k);
				A[i0].re -= f0.re*t[k].re  + f0.im*t[k].im + 
							cg.re*A[i1].re + cg.im*A[i1].im;
				A[i0].im -= f0.re*t[k].im  - f0.im*t[k].re +
							cg.re*A[i1].im - cg.im*A[i1].re;
			}
		}

		t[l].re = cfp1.re; t[l].im = cfp1.im;
		(*e)[i] = fp4;
	}

	//--------------------------//
	//------ INTERMISSION ------//
	//--------------------------//
	//Convient to shift off-diagonals by 1
	for(i = 1; i < n; ++i) d[i - 1] = d[i];
	d[n - 1] = ZERO;
	//Also need to update current accum eigs and A
	(*e)[0] = ZERO;
	for(i = 0; i < n; ++i){
		i0 = ACC2(n, i, i);
		fp1 = (*e)[i];
		(*e)[i] = A[i0].re;
		A[i0].re = fp1;
	}

	//Finally, set V to identity for eigenvectors
	if(vecs){
    	for(i = 0; i < n; ++i){
			for(j = 0; j < n; ++j){
				V[i][j].re = ZERO; V[i][j].im = ZERO;
			}
			V[i][i].re = ONE; ///shouldn't have to set the imag portion
		}
	}


	//--------------------------//
	//--- EIGENVALUES / VECS ---//
	//--------------------------//
	for(l = 0; l < n; ++l){
		k = 0; ///iteration counter: breaks if it reaches itrlim

		while(1){
			///Grab a small off-diag element
			m = l;
			while(1){
				if(m + 1 == n) break;
				if(ABSF(d[m]) < eps*(ABSF((*e)[m]) + ABSF((*e)[m+1]))) break;
				m += 1;
			}
			if(m == l) break;
	
			if(k >= itrLim) return 1u;  ////prevent hanging if QR fails		
			k++;

			///shift
			dg = HALF*((*e)[l+1] - (*e)[l])/d[l];
			fp2 = _shypot(dg);
			fp3 = (dg < ZERO)? dg - fp2 : dg + fp2;
			dg = (*e)[m] - (*e)[l] + d[l]/fp3;

			///plane->Givens rotations: get back to tridiagonal form
			fp3 = ONE; fp4 = ONE; fp5 = ZERO;
			for(i = m; i > l; --i){
				fp7 = fp3*d[i-1]; 
                fp6 = fp4*d[i-1];

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

				dg = (*e)[i] - fp5;
				fp2 = ((*e)[i-1] - dg)*fp3 + ((fp)2.0)*fp4*fp6;
				fp5 = fp3*fp2;
				(*e)[i] = dg + fp5;
				dg = fp4*fp2 - fp6;

				////accum eigenvectors
				///At this point, we can deal with all real arithmatic
				if(vecs){
                	for(j = 0; j < n; ++j){
                    	fp7 = V[i][j].re;
                    	////V(i+1, j) = s*V(i, j) + c*f
                    	V[i][j].re = fp3*V[i-1][j].re + fp4*fp7;
                    	////V(i, j) = c*V(i, j) - s*f
                    	V[i-1][j].re = fp4*V[i-1][j].re - fp3*fp7;
                	}
				}
			}

			///finish up updating
			(*e)[l] -= fp5;
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
		for(i = 0; i < n; ++i){
			for(j = 0; j < n; ++j){
				fp7 = V[i][j].re;
				V[i][j].re = V[i][j].re*t[j].re + V[i][j].im*t[j].im;
				V[i][j].im = V[i][j].im*t[j].re - fp7*t[j].im;
			}
		} 
		for(i = 0; i < n; ++i){
			if(ABSF(A[ACC2(n, i, i)].re) < eps) continue;
			for(j = 0; j < n; ++j){
				cfp1.re = ZERO; cfp1.im = ZERO;
				///Expensive part of the algo + complicated math = no functions
				///last time, I promise :)
				for(k = 0; k < i; ++k){
					////ctmp += A(k, i)*  dot  V(j, k)
					i0 = ACC2(n, k, i);
					cfp1.re += A[i0].re*V[j][k].re - A[i0].im*V[j][k].im;
					cfp1.im += A[i0].re*V[j][k].im + A[i0].im*V[j][k].re;
				}
				for(k = 0; k < i; ++k){
					////V(j, k) -= ctmp  dot  A(i, k)
					i0 = ACC2(n, i, k);
					V[j][k].re -= cfp1.re*A[i0].re + cfp1.im*A[i0].im;
					V[j][k].im += cfp1.re*A[i0].im - cfp1.im*A[i0].re;
				}
			}
		}
	}
	free(t); 

#undef A
#undef V		
	 
	_epairsort(n, e, V, vecs);
	return 0u;
}









