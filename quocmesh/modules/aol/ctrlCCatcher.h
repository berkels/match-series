#ifndef __CTRLCCATCHER_H
#define __CTRLCCATCHER_H

#include <platformDependent.h>

typedef void (*sigfunc)(int);

// Since old GCC versions don't support WARNING_OFF, turn off all warnings in here
// to get rid of the old style cast warning caused by the definition of SIG_ERR and SIG_DFL
// in the system header files.
#if ( defined ( __GNUC__ ) ) && ( GCC_VERSION < 40200 )
#pragma GCC system_header
#endif

WARNING_OFF(old-style-cast)

#include <csignal>

// These signal value SIG_ERR causes some trouble when you try to avoid
// old-style casts: it is always implemented as simple integers --
// sometimes c-casted to function pointer type "sigfunc", sometimes
// not, depending on the compiler's C library. So we have to cast
// it for ourselves. And because we cannot be sure that it
// is not c-casted yet, we have to enclose its usage
// (every usage, as it is a #define, what is why we use
// our own value "ErrorSignal") with WARNING_ON.
// Thus, it would be useless to use static_cast<sigfunc> here.
const sigfunc ErrorSignal = reinterpret_cast<sigfunc> (SIG_ERR);
const int InterruptSignal = SIGINT;

const sigfunc DefaultHandler = reinterpret_cast<sigfunc> (SIG_DFL);

WARNING_ON(old-style-cast)

using namespace std;

namespace aol {

void ctrlCHandler ( int sig );
bool getCtrlCState();
void resetCtrlCState();

//! Prints "[action (y/n, Ctrl-C to terminate)? " in the given color
bool askForAction ( string action, string color = aol::color::error );

} // end of namespace aol.

#endif
