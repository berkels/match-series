#ifndef __CIMGINCLUDES_H
#define __CIMGINCLUDES_H

#include <platformDependent.h>

// QCustomPlot must be included before CImg, since CImg defines things like 'Bool' and 'Status'
#ifdef USE_MODULES_QT
#include <qcustomplotIncludes.h>
#include <qtIncludes.h>
#endif

#ifdef __GNUC__
WARNING_OFF(old-style-cast)
#pragma GCC system_header
#endif

#if defined (__clang__)
#pragma clang system_header
#endif

#ifdef USE_LIB_PNG
// Enable LibPNG support
#define cimg_use_png 1
#endif
#include <CImg.h>

#ifdef PI
#undef PI
#endif
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#ifdef __GNUC__
/* WARNING_ON(old-style-cast) */
#endif

#endif // __CIMGINCLUDES_H
