#ifndef _COMPLEX_H
#define _COMPLEX_H

/* Shut up MSVC " possible loss of data" warnings */ 
#		pragma warning( disable : 4244 )
#		pragma warning( disable : 4273 )	/* the bloody "inconsistent dll linkage" -- DRASTIC */
#		pragma warning( disable : 4305 )

#include "../../../gmtsar/declspec.h"
EXTERN_MSC fcomplex Cmul(fcomplex x, fcomplex y);
EXTERN_MSC fcomplex Cexp(float theta);
EXTERN_MSC fcomplex Conjg(fcomplex z);
EXTERN_MSC fcomplex RCmul(float a, fcomplex z);
EXTERN_MSC fcomplex Cadd(fcomplex x, fcomplex y);
EXTERN_MSC float Cabs(fcomplex z);

#endif /* _COMPLEX_H */
