// This file is part of the imaging2 class library. 
// 
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.

//! Header file for extern LAPACK routines. 


#ifndef IMAGING_LAPACK_H
#define IMAGING_LAPACK_H


#ifdef LAPACK_ACML

#include <acml.h>

#define IMAGING_DGETRI(n, a, lda, ipiv, work, lwork, info) dgetri(n, a, lda, ipiv, info)
#define IMAGING_SGETRI(n, a, lda, ipiv, work, lwork, info) sgetri(n, a, lda, ipiv, info)

#define IMAGING_DGETRF(m, n, a, lda, ipiv, info) dgetrf(m, n, a, lda, ipiv, info)
#define IMAGING_SGETRF(m, n, a, lda, ipiv, info) sgetrf(m, n, a, lda, ipiv, info)

#define IMAGING_DGESV(n, nrhs, a, lda, ipiv, b, ldb, info) dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
#define IMAGING_SGESV(n, nrhs, a, lda, ipiv, b, ldb, info) sgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

#define IMAGING_DORGTR(uplo, n, a, lda, tau, work, lwork, info) dorgtr(uplo, n, a, lda, tau, info)
#define IMAGING_SORGTR(uplo, n, a, lda, tau, work, lwork, info) sorgtr(uplo, n, a, lda, tau, info)

#define IMAGING_DSYTRD(uplo, n, a, lda, d, e, tau, work, lwork, info) dsytrd(uplo, n, a, lda, d, e, tau, info)
#define IMAGING_SSYTRD(uplo, n, a, lda, d, e, tau, work, lwork, info) ssytrd(uplo, n, a, lda, d, e, tau, info)

#define IMAGING_DSTEQR(compz, n, d, e, z, ldz, work, info) dsteqr(compz, n, d, e, z, ldz, info)
#define IMAGING_SSTEQR(compz, n, d, e, z, ldz, work, info) ssteqr(compz, n, d, e, z, ldz, info)

#else

#define IMAGING_DGETRI(n, a, lda, ipiv, work, lwork, info) dgetri_(&n, a, &lda, ipiv, work, &lwork, info)
#define IMAGING_SGETRI(n, a, lda, ipiv, work, lwork, info) sgetri_(&n, a, &lda, ipiv, work, &lwork, info)

#define IMAGING_DGETRF(m, n, a, lda, ipiv, info) dgetrf_(&m, &n, a, &lda, ipiv, info)
#define IMAGING_SGETRF(m, n, a, lda, ipiv, info) sgetrf_(&m, &n, a, &lda, ipiv, info)

#define IMAGING_DGESV(n, nrhs, a, lda, ipiv, b, ldb, info) dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info)
#define IMAGING_SGESV(n, nrhs, a, lda, ipiv, b, ldb, info) sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info)

#define IMAGING_DORGTR(uplo, n, a, lda, tau, work, lwork, info) dorgtr_(&uplo, &n, a, &lda, tau, work, &lwork, info)
#define IMAGING_SORGTR(uplo, n, a, lda, tau, work, lwork, info) sorgtr_(&uplo, &n, a, &lda, tau, work, &lwork, info)

#define IMAGING_DSYTRD(uplo, n, a, lda, d, e, tau, work, lwork, info) dsytrd_(&uplo, &n, a, &lda, d, e, tau, work, &lwork, info)
#define IMAGING_SSYTRD(uplo, n, a, lda, d, e, tau, work, lwork, info) ssytrd_(&uplo, &n, a, &lda, d, e, tau, work, &lwork, info)

#define IMAGING_DSTEQR(compz, n, d, e, z, ldz, work, info) dsteqr_(&compz, &n, d, e, z, &ldz, work, info)
#define IMAGING_SSTEQR(compz, n, d, e, z, ldz, work, info) ssteqr_(&compz, &n, d, e, z, &ldz, work, info)

#endif /* LAPACK_ACML */



extern "C" {
void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
void sgetri_(int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);

void dgetrf_(int *m,int *n, double *a, int *lda, int *ipiv, int *info);
void sgetrf_(int *m,int *n, float *a, int *lda, int *ipiv, int *info);

void dgesv_(int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);
void sgesv_(int *n,int *nrhs,float *a,int *lda,int *ipiv,float *b,int *ldb,int *info);

void dorgtr_(char *uplo, int *n, double *a, int *lda, double *tau,
	    double *work, int *lwork, int *info);
void sorgtr_(char *uplo, int *n, float *a, int *lda, float *tau,
		float *work, int *lwork, int *info);

void dsytrd_(char *uplo, int *n, double *a, int *lda, double *d, double *e, 
	       double *tau, double *work, int *lwork, int *info);
void ssytrd_(char *uplo, int *n, float *a, int *lda, float *d, float *e, 
	       float *tau, float *work, int *lwork, int *info);

void dsteqr_(char *compz, int *n,double *d, double *e, 
	       double *z, int *ldz, double *work, int *info); 
void ssteqr_(char *compz,int *n, float *d, float *e, 
	       float *z, int *ldz, float *work, int *info); 
}

#endif
