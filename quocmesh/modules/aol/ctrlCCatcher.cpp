#include <aol.h>
#include <ctrlCCatcher.h>

namespace aol {

static bool ctrlC_pressed = false;

void ctrlCHandler ( int sig ) {
  signal(sig, ctrlCHandler);
  ctrlC_pressed = true;
}

bool getCtrlCState() {
  return ctrlC_pressed;
}

void resetCtrlCState() {
  ctrlC_pressed = false;
}

bool askForAction ( string action, string color ) {

  cerr << color << action << " (y/n, Ctrl-C to terminate)? " << aol::color::reset;
  string choice;
  sigfunc currentHandler = signal ( InterruptSignal, DefaultHandler );
  cin >> choice;
  signal ( InterruptSignal, currentHandler );

  return ( choice == "y" || choice == "yes" || choice == "j" || choice == "ja" );
}

} // end of namespace aol.
