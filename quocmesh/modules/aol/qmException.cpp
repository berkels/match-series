#include <qmException.h>
#include <aol.h>

namespace aol {
  bool Exception::consumed = false;
  bool Exception::wasCopied = false;
  int Exception::copyNum = 0;

//! When throwing exceptions in openmp threads and not catching them in the same thread,
//! the terminate handler is called. When using openmp, this is redirected to catchGlobal,
//! which can on most platforms re-throw and then catch the exception, or at least print a general error message.
//! std::set_terminate ( abort ) restores the standard behavior.
void catchGlobal() {
#ifdef _OPENMP
#pragma omp critical (catchGlobal)
#endif
  {
    // The try-throw-catch trick below doesn't seem to work under VC++, so
    // print something to make sure that the user at least sees that
    // aol::catchGlobal is called..
#ifdef  _MSC_VER
    cerr << "Terminating with aol::catchGlobal()\n";
#endif
    cerr << aol::color::error << endl;
    try {
      throw;
    }
    catch ( aol::Exception &el ) {
      cerr << "\nUnhandled aol::Exception caught!\n";
      el.dump();
    }
    catch ( std::bad_alloc& ba ) {
      cerr << "\nstd::bad_alloc caught: " << ba.what() << endl;
      cerr << "Current memusage according to aol::memusage: " << aol::memusage() / ( 1024 * 1024 ) << " MB\n";
    }
    catch ( std::exception &ex ) {
      cerr << "\nUnhandled std::exception caught!\n";
      cerr << ex.what () << endl;
    }
    catch (...) {
      cerr << "\nUnhandled unknown exception caught!\n";
    }
    cerr << aol::color::reset;
    aol::callSystemPauseIfNecessaryOnPlatform();
    abort();
  }
}

// Since exceptions thrown inside OpenMP blocks can't be caught outside the block,
// we use our custom terminate handler to have a chance to dump the exception message.
#ifdef _OPENMP
std::terminate_handler oldTerminateHandler = std::set_terminate(catchGlobal);
#endif

} // namespace aol
