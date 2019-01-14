#ifndef __EIGENINCLUDES_H
#define __EIGENINCLUDES_H

#include <platformDependent.h>

#ifdef __GNUC__
WARNING_OFF(old-style-cast)
#pragma GCC system_header
#endif

// Eigen has an enum that clashes with the X11 preprocessor define 'Success'
#ifdef X_H
#ifdef Success
#undef Success
#endif
#endif

#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#ifdef USE_EXTERNAL_SUITESPARSE
#include <Eigen/SPQRSupport>
#endif

#ifdef __GNUC__
/* WARNING_ON(old-style-cast) */
#endif

#endif // __EIGENINCLUDES_H
