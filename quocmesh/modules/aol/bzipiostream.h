#ifndef __BZIPIOSTREAM_H
#define __BZIPIOSTREAM_H

#include "aol.h"

#ifdef USE_LIB_BZ2
#include <bzlib.h>
// Unfortunately bzlib.h defines a lot of crap under some platforms.
// We need to kill the defines here.
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef PACKED
#undef PACKED
#endif
#ifdef rad1
#undef rad1
#endif
#ifdef rad2
#undef rad2
#endif
#ifdef rad3
#undef rad3
#endif
#endif

#ifdef USE_ZLIB
#include <zlib.h>
#endif

namespace aol {

// VC++ 2012 always produces warning C4250 when deriving from the std stream classes, but this warning should be ignored.
#if defined ( _MSC_VER ) && ( _MSC_VER == 1700 )
#pragma warning ( push )
#pragma warning ( disable : 4250 )
#endif

/**
 * \brief Can be used like an ifstream, but supports on the fly bzip2 decompression using libbz2.
 *
 * Automatically decides whether to decompress the input file or not based on the suffix of
 * the constructor argument FileName.
 *
 * \author Berkels
 */
class Bzipifstream : public stringstream{
public:
  Bzipifstream ( const char *FileName ){
    FILE *tempf;
    // copied from http://www.zlib.net/zlib_how.html
    // chunk is simply the buffer size for feeding data to and pulling data from the zlib routines.
    // Larger buffer sizes would be more efficient, especially for inflate(). If the memory is available,
    // buffers sizes on the order of 128K or 256K bytes should be used.
    const int chunk = 262144;
    char buf[chunk];
    tempf = fopen ( FileName, "rb" );
    if ( !tempf ){
      string errorMessage = "Could not open \"";
      errorMessage += FileName;
      errorMessage += "\" for reading.\n";
      throw Exception ( errorMessage.c_str(), __FILE__, __LINE__ );
    }

    if ( hasBzipSuffix ( FileName ) ) {
#ifdef USE_LIB_BZ2
      BZFILE *tempb;
      int bzerror;

      tempb = BZ2_bzReadOpen ( &bzerror, tempf, 0, 0, NULL, 0 );
      while ( bzerror == BZ_OK ) {
        const int bytesRead = BZ2_bzRead ( &bzerror, tempb, buf, chunk );
        (*this).write ( buf, bytesRead );
      }
      BZ2_bzReadClose ( &bzerror, tempb );
#else
      throw Exception ( "Reading bz2 compressed files with Bzipifstream without using bzlib is impossible.\nDefine USE_LIB_BZ2 and link bzlib, e.g by using CFLAG += -DUSE_LIB_BZ2 and LFLAGS += -lbz2, to remedy this.\n", __FILE__, __LINE__ );
#endif // USE_LIB_BZ2
    } else {
      if ( hasZlibSuffix ( FileName ) || hasGzipSuffix ( FileName ) )  {
#ifdef USE_ZLIB
        // ret will be used for zlib return codes.
        // have is the amount of data returned from inflate().
        // The strm structure is used to pass information to and from the zlib routines, and to maintain the inflate() state.
        // buf and out are the input and output buffers for inflate().
        int ret;
        unsigned have;
        z_stream strm;
        unsigned char out[chunk];

        // allocate inflate state
        strm.zalloc = Z_NULL;
        strm.zfree = Z_NULL;
        strm.opaque = Z_NULL;
        strm.avail_in = 0;
        strm.next_in = Z_NULL;

        // window bits for inflate: logBase2 of window size + automatic header detection
        const int windowBits = 15 | 32;

        // intialize
// Only recent GCC versions allow to the diagnostic pragma inside functions.
#if ( defined ( __GNUC__ ) ) && ( ( GCC_VERSION >= 40600 ) || ( defined (__clang__) ) )
        WARNING_OFF ( old-style-cast )
#endif
        ret = inflateInit2 (&strm, windowBits);
#if ( defined ( __GNUC__ ) ) && ( ( GCC_VERSION >= 40600 ) || ( defined (__clang__) ) )
        WARNING_ON ( old-style-cast )
#endif
        if (ret != Z_OK)  {
          cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
          return;
        }

        // decompress until deflate stream ends or end of file
        // We read input data and set the strm structure accordingly. If we've reached the end of the input file,
        // then we leave the outer loop and report an error, since the compressed data is incomplete.
        do {
          strm.avail_in = fread(buf, 1, chunk, tempf);
          if (ferror(tempf)) {
            (void)inflateEnd(&strm);
            cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
            return;
          }
          if (strm.avail_in == 0)
            break;
          strm.next_in = reinterpret_cast<unsigned char *> ( buf );

          // run inflate() on input until output buffer not full
          do {
            strm.avail_out = chunk;
            strm.next_out = out;

            ret = inflate(&strm, Z_NO_FLUSH);
            if ( ( ret == Z_STREAM_ERROR ) || ( ret == Z_NEED_DICT ) || ( ret == Z_DATA_ERROR ) || ( ret == Z_MEM_ERROR ) )  {
                (void)inflateEnd(&strm);
                cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
                return;
            }

            have = chunk - strm.avail_out;
            (*this).write ( reinterpret_cast<char *> ( out ), have );
            if ( this->bad() ) {
              (void)inflateEnd(&strm);
              cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
              return;
            }
          } while (strm.avail_out == 0);
          // done when inflate() says it's done
        } while (ret != Z_STREAM_END);

        // clean up
        (void)inflateEnd(&strm);
        if ( ret != Z_STREAM_END ) cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
#else
        throw Exception ( "Reading zlib compressed files with Bzipifstream without using zlib is impossible.\n", __FILE__, __LINE__ );
#endif // USE_ZLIB
      }
      else {
        while ( !feof(tempf) ) {
          const int bytesRead = static_cast<int>(fread ( buf, 1, chunk, tempf ));
          (*this).write ( buf, bytesRead );
        }
      }
    }
    fclose ( tempf );
  }

private:
  template< typename AnyThing >
  Bzipifstream& operator<< ( const AnyThing& ); // do not implement
};

/**
 * \brief Can be used like an ofstream, but supports on the fly bzip2 compression using libbz2.
 *
 * \author Berkels
 */
class Bzipofstream : public stringstream{
  FILE *outFile;
  const bool compressOutputBzip;
  const bool compressOutputZlib;
  const bool compressOutputGzip;
public:
  Bzipofstream ( const char *FileName )
    : outFile ( NULL ),
      compressOutputBzip ( hasBzipSuffix ( FileName ) ),
      compressOutputZlib ( hasZlibSuffix ( FileName ) ),
      compressOutputGzip ( hasGzipSuffix ( FileName ) )
  {
    outFile = fopen ( FileName, "wb" );
    if ( !outFile ){
      string errorMessage = "Could not open \"";
      errorMessage += FileName;
      errorMessage += "\" for output.\n";
      throw Exception ( errorMessage.c_str(), __FILE__, __LINE__ );
    }
  }
  ~Bzipofstream(){
    close();
  }
  void close() {
    // copied from http://www.zlib.net/zlib_how.html
    // chunk is simply the buffer size for feeding data to and pulling data from the zlib routines.
    // Larger buffer sizes would be more efficient, especially for inflate(). If the memory is available,
    // buffers sizes on the order of 128K or 256K bytes should be used.
    const int chunk = 262144;
    char buf[chunk];
    // We can't do anything when there is no outfile (probably close was already called).
    if ( outFile == NULL )
      return;

    if ( compressOutputBzip ){
#ifdef USE_LIB_BZ2
      BZFILE *b;
      int bzerror;
      b = BZ2_bzWriteOpen ( &bzerror, outFile, 9, 0, 0 );
      while ( !this->eof() ) {
        (*this).read ( buf, chunk );
        const int bytesRead = static_cast<int>(this->gcount());
        BZ2_bzWrite ( &bzerror, b, buf, bytesRead );
      }
      unsigned int temp;
      BZ2_bzWriteClose ( &bzerror, b, 0, &temp, &temp );
#else
      throw Exception ( "Writing bz2 compressed files with Bzipofstream without using bzlib is impossible.\nDefine USE_LIB_BZ2 and link bzlib, e.g by using CFLAG += -DUSE_LIB_BZ2 and LFLAGS += -lbz2, to remedy this.\n", __FILE__, __LINE__ );
#endif
    }
    else{
      if ( compressOutputZlib || compressOutputGzip )  {
#ifdef USE_ZLIB
        // ret will be used for zlib return codes. flush will keep track of the current flushing state for deflate(),
        // which is either no flushing, or flush to completion after the end of the input file is reached.
        // have is the amount of data returned from deflate().
        // The strm structure is used to pass information to and from the zlib routines, and to maintain the deflate() state.
        // buf and out are the input and output buffers for deflate().
        int ret, flush;
        unsigned have;
        z_stream strm;
        unsigned char out[chunk];

        // allocate deflate state
        strm.zalloc = Z_NULL;
        strm.zfree = Z_NULL;
        strm.opaque = Z_NULL;

        // window bits for deflate: logBase2 of window size
        int windowBits = 15; // maximum
        if ( compressOutputGzip ) windowBits |= 16; // Gzip encoding

        // intialize
#if ( defined ( __GNUC__ ) ) && ( ( GCC_VERSION >= 40600 ) || ( defined (__clang__) ) )
        WARNING_OFF ( old-style-cast )
#endif
        ret = deflateInit2 (&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits, 9 /* maximum memory level*/, Z_DEFAULT_STRATEGY);
#if ( defined ( __GNUC__ ) ) && ( ( GCC_VERSION >= 40600 ) || ( defined (__clang__) ) )
        WARNING_ON ( old-style-cast )
#endif
        if (ret != Z_OK)  {
          cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
          return;
        }

        // compress until end of file
        // We start off by reading data from the input file. The number of bytes read is put directly into avail_in,
        // and a pointer to those bytes is put into next_in. We also check to see if end-of-file on the input has been
        // reached using feof(). If we are at the end of file, then flush is set to the zlib constant Z_FINISH, which is
        // later passed to deflate() to indicate that this is the last chunk of input data to compress. If we are not yet
        // at the end of the input, then the zlib constant Z_NO_FLUSH will be passed to deflate to indicate that we are still in the middle of the uncompressed data.
        do {
          // fread(in, 1, chunk, source);
          (*this).read ( buf, chunk );
          strm.avail_in = static_cast<int>(this->gcount());
          if ( this->bad() ) {
            (void)deflateEnd(&strm);
            cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
            return;
          }
          flush = this->eof() ? Z_FINISH : Z_NO_FLUSH;
          strm.next_in = reinterpret_cast<unsigned char*>( buf );

          // The inner do-loop passes our chunk of input data to deflate(), and then keeps calling deflate() until it is done producing output.
          // Once there is no more new output, deflate() is guaranteed to have consumed all of the input, i.e., avail_in will be zero.
          // run deflate() on input until output buffer not full, finish compression if all of source has been read in */
          do {
            strm.avail_out = chunk;
            strm.next_out = out;
            ret = deflate(&strm, flush);
            QUOC_ASSERT(ret != Z_STREAM_ERROR);
            have = chunk - strm.avail_out;
            if (fwrite(out, 1, have, outFile) != have || ferror(outFile)) {
              (void)deflateEnd(&strm);
              cerr << aol::color::error << "zlib returned an error! Aborting..." << aol::color::reset << endl;
              return;
            }
          } while (strm.avail_out == 0);
          QUOC_ASSERT(strm.avail_in == 0);

          // done when last data in file processed
        } while (flush != Z_FINISH);
        QUOC_ASSERT(ret == Z_STREAM_END);

        // clean up
        (void)deflateEnd(&strm);
#else
        throw Exception ( "Writing bz2 compressed files with Bzipofstream without using bzlib is impossible.\nDefine USE_LIB_BZ2 and link bzlib, e.g by using CFLAG += -DUSE_LIB_BZ2 and LFLAGS += -lbz2, to remedy this.\n", __FILE__, __LINE__ );
#endif

      } else {
        while ( !this->eof() ) {
          (*this).read ( buf, chunk );
          const int bytesRead = static_cast<int>(this->gcount());
          fwrite( buf, 1, bytesRead, outFile );
        }
      }
    }
    fclose( outFile );
    outFile = NULL;
  }

private:
  template< typename AnyThing >
  Bzipofstream& operator>> ( const AnyThing& ); // do not implement
};

// Turn on the warning C4250 again.
#if defined ( _MSC_VER ) && ( _MSC_VER == 1700 )
#pragma warning ( pop )
#endif

}// end namespace aol

#endif //__BZIPIOSTREAM_H

