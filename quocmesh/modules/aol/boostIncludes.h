#ifndef __BOOSTINCLUDES_H
#define __BOOSTINCLUDES_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

#ifdef USE_BOOST
#ifndef Q_MOC_RUN // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/gamma.hpp>
#endif
#endif

#endif
