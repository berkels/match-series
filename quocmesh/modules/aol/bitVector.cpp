#include <bitVector.h>
#include <bzipiostream.h>

namespace aol {

ostream& BitVector::dump ( ostream & out ) const {
  for ( int i = 0 ; i < this->_size; i++ ) {
    if ( i % 64 == 0 ) {
      out << endl;
    } else if ( i % 8 == 0 ) {
      out << " ";
    }
    out << static_cast<int> ( ( *this ) [i] );
  }
  if ( this->_size % 64 != 0 ) {
    out << endl;
  }
  return ( out );
}

void BitVector::saveToFile ( const char *filename ) const {
  aol::Bzipofstream out ( filename );
  out << aol::VectorFileMagicChar::BitVector << FileFormatMagicNumber<bool>::FFType << endl;
  out << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::BitVector << FileFormatMagicNumber<bool>::FFType << " storing an aol::BitVector" << endl;
  out << this->size() << endl;
  const char* buffer = reinterpret_cast<char*> ( _pFieldData );
  out.write ( buffer, convertBitVectorLengthToCharLength ( this->size() ) * sizeof ( bool ) );
}

void BitVector::loadFromFile ( const char *filename ) {
  aol::Bzipifstream reader ( filename );
  char L = 0;
  int ident;
  reader >> L;
  reader >> ident;
  if ( ( L != aol::VectorFileMagicChar::BitVector ) || ( ident != FileFormatMagicNumber<bool>::FFType ) ) {
    cerr << L << ident << ", should be " << aol::VectorFileMagicChar::BitVector << FileFormatMagicNumber<bool>::FFType << " (" << FileFormatMagicNumber<bool>::FFContainedName << ")" << endl;
    throw aol::Exception ( "Illegal magic number for aol::BitVector", __FILE__, __LINE__ );
  }
  aol::READ_COMMENTS ( reader );
  int size;
  reader >> size;
  char buffer[1024];
  reader.getline ( buffer, 1 );

  this->reallocate ( size );
  reader.read ( reinterpret_cast<char*>( _pFieldData ), convertBitVectorLengthToCharLength ( size ) * sizeof( bool ) );
}

ostream &operator<< ( ostream &os, const BitVector &BVec ) {
  return ( BVec.dump ( os ) );
}

}
