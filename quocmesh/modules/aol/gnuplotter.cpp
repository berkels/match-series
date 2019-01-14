#include <gnuplotter.h>

namespace aol {

bool runGnuplot ( const char *GnuplotCommandFileName ) {
  string systemCommand;
  // Since the Windows version of gnuplot finally also has an exe called "gnuplot" the ifdef
  // shouldn't be necessary anymore, it won't hurt to keep it here for now though.
#ifdef WIN32
  systemCommand = "gnuplot.exe ";
#else
  systemCommand = "gnuplot ";
#endif
  systemCommand += GnuplotCommandFileName;
  const bool failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
  if ( failed )
    cerr << "aol::runGnuplot: Calling gnuplot returned an error.\n";
  return !failed;
}
  
void plotPrecisionMatches ( const string GnuplotdatOutputFile,
                            const string &BaseOutName,
                            const string &BackgroundImageFile,
                            const string &PrecXMatchFile,
                            const string &PrecYMatchFile ) {
  std::ofstream gnuplotdat ( GnuplotdatOutputFile.c_str() );
  gnuplotdat << "set terminal postscript eps color\n";
  gnuplotdat << "set output \"" << BaseOutName << "\"\n";
  gnuplotdat << "unset xtics\n";
  gnuplotdat << "unset ytics\n";
  gnuplotdat << "set size square 1, 1.4285715\n";
  gnuplotdat << "plot \"" << BackgroundImageFile << "\" binary filetype=png flipy w rgbimage notitle, \""
             << PrecXMatchFile << "\" with vectors title \"X\", \""
             << PrecYMatchFile << "\" with vectors title \"Y\"";
  gnuplotdat.close();
  aol::runGnuplot ( GnuplotdatOutputFile.c_str() );
}
  
}
