#include <qmException.h>

void foo () {
  throw aol::OutOfBoundsException ( "foo", __FILE__, __LINE__ );
}

void bar () {
  throw aol::UnimplementedCodeException ( "bar", __FILE__, __LINE__ );
}

int main ()
{
  // If the exception has been handled, a call to consume
  // avoids printing an error message
  try { foo (); }
  catch (aol::OutOfBoundsException& e ) { e.consume (); };

  // Catch by value creates an additional copy of the exception
  // but still consuming it once is sufficient
  try { foo (); }
  catch (aol::OutOfBoundsException e ) { e.consume (); };

  // It the exception is caught but not consumed,
  // the error message will be printed automatically
  try { bar (); }
  catch (aol::UnimplementedCodeException& /*e*/ ) { };

  // this is a somewhat special selfTest that should not print a success message
}
