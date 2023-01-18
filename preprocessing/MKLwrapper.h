// WRAPPER.h
//
// C++ wrapper for BLAS and LAPACK fortran routines
//

#ifndef __WRAPPER
#define __WRAPPER

#define ADD_UNDERSCORE 

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define THIS_SOL2    
#define THIS_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define THIS_SGI
#define THIS_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define THIS_LINUX
#define THIS_ARCHITECTURE "Linux"

#elif defined (__APPLE__)
#define THIS_MAC
#define THIS_ARCHITECTURE "Mac"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define THIS_AIX
#define THIS_ARCHITECTURE "IBM AIX"

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define THIS_ALPHA
#define THIS_ARCHITECTURE "Compaq Alpha"

#elif defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#if defined (__MINGW32__) || defined (__MINGW32__)
#define THIS_MINGW
#elif defined (__CYGWIN32__) || defined (__CYGWIN32__)
#define THIS_CYGWIN
#else
#define THIS_WINDOWS
#undefine ADD_UNDERSCORE 
#endif
#define THIS_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define THIS_HP
#define THIS_ARCHITECTURE "HP Unix"
#undefine ADD_UNDERSCORE 

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define THIS_HP
#define THIS_ARCHITECTURE "HP 700 Unix"
#undefine ADD_UNDERSCORE 

#else
#define THIS_UNKNOWN
#define THIS_ARCHITECTURE "unknown"
#endif


#if defined(THIS_SOL2) && !defined(NSUNPERF) && defined(BLAS64)
#define SUN64
#endif


#if defined SUN64

//#define XERBLA xerbla_64_
#define DGESV dgesv_64_
#define DGETRF dgetrf_64_
#define DGETRS dgetrs_64_
#define DGETRI dgetri_64_
#define DPOTRF dpotrf_64_
#define DPOTRS dpotrs_64_
#define DPOTRI dpotri_64_
#define DGELS dgels_64_
#define DGEEV dgeev_64_
#define DGEMM dgemm_64_
#define SGEMM sgemm_64_
#define DGEMV dgemv_64_
#define SGEMV sgemv_64_
#define DGER dger_64_
#define DAXPY daxpy_64_
#define SAXPY saxpy_64_
#define DSCAL dscal_64_
#define SSCAL sscal_64_
#define DCOPY dcopy_64_
#define SCOPY scopy_64_
#define DNRM2 dnrm2_64_
#define SNRM2 snrm2_64_
#define DDOT ddot_64_
#define SDOT sdot_64_
#define DSYEV dsyev_64_

#define ZGESV zgesv_64_
#define ZGETRF zgetrf_64_
#define ZGETRS zgetrs_64_
#define ZGETRI zgetri_64_
#define ZPOTRF zpotrf_64_
#define ZPOTRS zpotrs_64_
#define ZPOTRI zpotri_64_
#define ZGELS zgels_64_
#define ZGEEV zgeev_64_
#define ZGEMM zgemm_64_
#define ZGEMV zgemv_64_
#define ZGERU zgeru_64_
#define ZAXPY zaxpy_64_
#define ZSCAL zscal_64_
#define ZCOPY zcopy_64_
#define ZNRM2 znrm2_64_
#define ZDOT zdot_64_

#elif defined ADD_UNDERSCORE

//#define XERBLA xerbla_
#define DGESV dgesv_
#define DGETRF dgetrf_
#define DGETRS dgetrs_
#define DGETRI dgetri_
#define DPOTRF dpotrf_
#define DPOTRS dpotrs_
#define DPOTRI dpotri_
#define DGELS dgels_
#define DGEEV dgeev_
#define DGEMM dgemm_
#define SGEMM sgemm_
#define DGEMV dgemv_
#define SGEMV sgemv_
#define DGER dger_
#define DAXPY daxpy_
#define SAXPY saxpy_
#define DSCAL dscal_
#define SSCAL sscal_
#define DCOPY dcopy_
#define SCOPY scopy_
#define DNRM2 dnrm2_
#define SNRM2 snrm2_
#define DDOT ddot_
#define SDOT sdot_
#define DSYEV dsyev_

#define ZGESV zgesv_
#define ZGETRF zgetrf_
#define ZGETRS zgetrs_
#define ZGETRI zgetri_
#define ZPOTRF zpotrf_
#define ZPOTRS zpotrs_
#define ZPOTRI zpotri_
#define ZGELS zgels_
#define ZGEEV zgeev_
#define ZGEMM zgemm_
#define ZGEMV zgemv_
#define ZGERU zgeru_
#define ZAXPY zaxpy_
#define ZSCAL zscal_
#define ZCOPY zcopy_
#define ZNRM2 znrm2_
#define ZDOT zdot_

#else

//#define XERBLA xerbla
#define DGESV dgesv
#define DGETRF dgetrf
#define DGETRS dgetrs
#define DGETRI dgetri
#define DPOTRF dpotrf
#define DPOTRS dpotrs
#define DPOTRI dpotri
#define DGELS dgels
#define DGEEV dgeev
#define DGEMM dgemm
#define SGEMM sgemm
#define DGEMV dgemv
#define SGEMV sgemv
#define DGER dger
#define DAXPY daxpy
#define SAXPY saxpy
#define DSCAL dscal
#define SSCAL sscal
#define DCOPY dcopy
#define SCOPY scopy
#define DNRM2 dnrm2
#define SNRM2 snrm2
#define DDOT ddot
#define SDOT sdot
#define DSYEV dsyev

#define ZGESV zgesv
#define ZGETRF zgetrf
#define ZGETRS zgetrs
#define ZGETRI zgetri
#define ZPOTRF zpotrf
#define ZPOTRS zpotrs
#define ZPOTRI zpotri
#define ZGELS zgels
#define ZGEEV zgeev
#define ZGEMM zgemm
#define ZGEMV zgemv
#define ZGERU zgeru
#define ZAXPY zaxpy
#define ZSCAL zscal
#define ZCOPY zcopy
#define ZNRM2 znrm2
#define ZDOT zdot

#endif

typedef long long int Int;

/*typedef struct{ double r, i; } Complex;*/
#include <complex>
typedef std::complex<double> Complex ;

#endif
