#include "multiStreambuf.h"

namespace aol {

//---------------------------------------------------------------------------

MultiStreambuf::MultiStreambuf() {
  setp(0, 0);
}
//---------------------------------------------------------------------------

size_t MultiStreambuf::addStreambuf(streambuf * buffer) {
  _streambufPtrList.push_back(buffer);
  return _streambufPtrList.size();
}
//---------------------------------------------------------------------------

size_t MultiStreambuf::removeStreambuf(streambuf * buffer) {
  _streambufPtrList.remove(buffer);
  return _streambufPtrList.size();
}
//---------------------------------------------------------------------------

MultiStreambuf::int_type MultiStreambuf::overflow(MultiStreambuf::int_type c) {
  if (!traits_type::eq_int_type(c, traits_type::eof())) {
    StreambufList::iterator iter = _streambufPtrList.begin();
    for (; iter != _streambufPtrList.end(); ++iter)
      (*iter)->sputc(c);
  }
  return traits_type::not_eof(c);
}
//---------------------------------------------------------------------------

//! synchronizes all buffers. Returns for failure if at least
//! one participation buffer has failed and returned -1.
MultiStreambuf::int_type MultiStreambuf::sync() {
  int ret = 0;
  StreambufList::iterator iter = _streambufPtrList.begin();
  for (; iter != _streambufPtrList.end(); ++iter)
    // if ret has already set to value "-1" (failed),
    // do not overwrite this, but keep.
    if ( ret == -1 )
      (*iter)->pubsync();
    // otherwise, give *iter a chance to indicate
    // failure.
    else
      ret = (*iter)->pubsync();
  return ret;
}


//===========================================================================


AdditionalOutputToFile::AdditionalOutputToFile ( string Filename )
  : _filename ( Filename )
  , _filestream ( Filename.c_str() )
  , _previousCoutStreambuf ( cout.rdbuf() )
  , _previousCerrStreambuf ( cerr.rdbuf() )
  , _previousClogStreambuf ( clog.rdbuf() ) {

  _multiCoutStreambuf.addStreambuf ( _previousCoutStreambuf );
  _multiCerrStreambuf.addStreambuf ( _previousCerrStreambuf );
  _multiClogStreambuf.addStreambuf ( _previousClogStreambuf );

  _multiCoutStreambuf.addStreambuf ( _filestream.rdbuf() );
  _multiCerrStreambuf.addStreambuf ( _filestream.rdbuf() );
  _multiClogStreambuf.addStreambuf ( _filestream.rdbuf() );

  cout.rdbuf ( &_multiCoutStreambuf );
  cerr.rdbuf ( &_multiCerrStreambuf );
  clog.rdbuf ( &_multiClogStreambuf );
}
//---------------------------------------------------------------------------

AdditionalOutputToFile::~AdditionalOutputToFile () {
  cout.flush();
  cerr.flush();
  clog.flush();

  cout.rdbuf ( _previousCoutStreambuf );
  cerr.rdbuf ( _previousCerrStreambuf );
  clog.rdbuf ( _previousCerrStreambuf );

  clog << "All output has been written to file " << _filename << endl;
}

//---------------------------------------------------------------------------

} // end of namespace aol.
