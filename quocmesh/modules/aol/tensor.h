#ifndef __TENSOR_H
#define __TENSOR_H

#include <smallMat.h>

namespace aol {

//! Tensor of order n+1 if VecType is of order n
template <class VecType, int dim>
class Tensor {
public:
  Tensor () {}

  Tensor ( const Tensor& tensor ) {
    for ( int i = 0; i < dim; ++i )
      this->_vec [i] = tensor._vec [i];
  }

  Tensor& operator = ( const Tensor& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;

    for ( int i = 0; i < dim; ++i )
      this->_vec [i] = tensor._vec [i];

    return *this;
  }

~Tensor () {}

  const VecType& operator [] ( const int i ) const {
    return this->_vec [i];
  }

  VecType& operator [] ( const int i ) {
    return this->_vec [i];
  }

  ostream& print ( ostream& os ) const {
    for ( int i = 0; i < dim; ++i )
      os << this->_vec [i] << endl;
    return os;
  }

  void setZero () {
    for ( int i = 0; i < dim; ++i )
      this->_vec [i].setZero ();
  }
  
  template < typename RealType >
  Tensor &operator*= ( const RealType Alpha ) {
    for ( int i = 0; i < dim; ++i )
      this->_vec[i] *= Alpha;
    return *this; 
  }
  

protected:
  // Allow dummy tensors with dim == 0, but avoid arrays of size 0.
  VecType _vec [dim > 0 ? dim : 1];
};


template <class DataType>
class Tensor222 : public Tensor<Matrix22<DataType>, 2 > {
public:
  Tensor222 ()
      : Tensor<Matrix22<DataType>, 2 > () {}

  Tensor222 ( const Tensor222<DataType>& tensor )
      : Tensor<Matrix22<DataType>, 2 > ( tensor ) {}

  Tensor222<DataType>& operator = ( const Tensor222<DataType>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix22<DataType>, 2 >::operator = ( tensor );
    return *this;
  }

~Tensor222 () {}

  const DataType& operator () ( const int i, const int j, const int k ) const {
    return this->_vec [i][j][k];
  }

  DataType& operator () ( const int i, const int j, const int k ) {
    return this->_vec [i][j][k];
  }

  DataType get ( const int i, const int j, const int k ) const {
      return this->_vec [i][j][k];
    }

  void set ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] = value;
  }

  void add ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] += value;
  }
};

template <class DataType>
class Tensor333 : public Tensor<Matrix33<DataType>, 3 > {
public:
  Tensor333 ()
      : Tensor<Matrix33<DataType>, 3 > () {}

  Tensor333 ( const Tensor333<DataType>& tensor )
      : Tensor<Matrix33<DataType>, 3 > ( tensor ) {}

  Tensor333<DataType>& operator = ( const Tensor333<DataType>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix33<DataType>, 3 >::operator = ( tensor );
    return *this;
  }

~Tensor333 () {}

  const DataType& operator () ( const int i, const int j, const int k ) const {
    return this->_vec [i][j][k];
  }

  DataType& operator () ( const int i, const int j, const int k ) {
    return this->_vec [i][j][k];
  }

  DataType get ( const int i, const int j, const int k ) const {
      return this->_vec [i][j][k];
    }

  void set ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] = value;
  }

  void add ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] += value;
  }
};

template <class DataType>
ostream& operator << ( ostream& os, const Tensor222<DataType>& t ) {
  return t.print ( os );
}

template <class DataType>
class Tensor322 : public Tensor<Matrix22<DataType>, 3 > {
public:
  Tensor322 ()
      : Tensor<Matrix22<DataType>, 3 > () {}

  Tensor322 ( const Tensor322<DataType>& tensor )
      : Tensor<Matrix22<DataType>, 3 > ( tensor ) {}

  Tensor322<DataType>& operator = ( const Tensor322<DataType>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Matrix22<DataType>, 3 >::operator = ( tensor );
    return *this;
  }

~Tensor322 () {}

  const DataType& operator () ( const int i, const int j, const int k ) const {
    return this->_vec [i][j][k];
  }

  DataType& operator () ( const int i, const int j, const int k ) {
    return this->_vec [i][j][k];
  }

  DataType get ( const int i, const int j, const int k ) const {
      return this->_vec [i][j][k];
    }

  void set ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] = value;
  }

  void add ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] += value;
  }
  
  aol::Vec3<DataType> getVector(const int j, const int k){
    aol::Vec3<DataType> vector;
    vector.set(this->_vec[0][j][k],this->_vec[1][j][k],this->_vec[2][j][k]);
    return vector;
  }
  
  void getVector(aol::Vec3<DataType>& vector , const int j, const int k){
	vector.set(this->_vec[0][j][k],this->_vec[1][j][k],this->_vec[2][j][k]);
  }
  
  void setVector(const aol::Vec3<DataType>& vector , const int j, const int k){
	 for (int i =0;i<3;i++){
	  this->_vec[i][j][k] = vector.get(i);
    }
  }
  
  DataType normSqr ( ) const{
    DataType result = 0.0;
    for( int i=0; i<3; ++i )
      for( int j=0; j<2; ++j )
        for( int k = 0; k <2; ++k )
          result += aol::Sqr( this->_vec[i][j][k] );
        
    return result;
  }
};

template <typename DataType>
class Tensor2222 : public Tensor<Tensor222<DataType>, 2 > {
public:
  Tensor2222 ()
      : Tensor<Tensor222<DataType>, 2 > () {}

  Tensor2222 ( const Tensor2222<DataType>& tensor )
      : Tensor<Tensor222<DataType>, 2 > ( tensor ) {}

  Tensor2222<DataType>& operator = ( const Tensor2222<DataType>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Tensor222<DataType>, 2 >::operator = ( tensor );
    return *this;
  }

~Tensor2222 () {}

  void apply ( const Mat<2, 2, DataType>& arg, Mat<2, 2, DataType>& dest ) const {
    dest.setZero ();
    applyAdd( arg, dest );
  }

  void applyAdd ( const Mat<2, 2, DataType>& arg, Mat<2, 2, DataType>& dest ) const {
    for ( int i = 0; i < 2; i++ )
      for ( int j = 0; j < 2; j++ )
        for ( int k = 0; k < 2; k++ )
          for ( int l = 0; l < 2; l++ )
            dest [i][j] += this->_vec [i][j][k][l] * arg [k][l];
  }

  const DataType& operator () ( const int i, const int j, const int k, const int l ) const {
    return this->_vec [i][j][k][l];
  }

  DataType& operator () ( const int i, const int j, const int k, const int l ) {
    return this->_vec [i][j][k][l];
  }

  DataType get ( const int i, const int j, const int k, const int l ) const {
      return this->_vec [i][j][k][l];
    }

  void set ( const int i, const int j, const int k, const int l, const DataType value ) {
    this->_vec [i][j][k][l] = value;
  }

  void add ( const int i, const int j, const int k, const int l, const DataType value ) {
    this->_vec [i][j][k][l] += value;
  }
};

template <class DataType>
ostream& operator << ( ostream& os, const Tensor2222<DataType>& t ) {
  return t.print ( os );
}

template <typename DataType>
class Tensor3322 : public Tensor<Tensor322<DataType>, 3 > {
public:
  Tensor3322 ()
      : Tensor<Tensor322<DataType>, 3 > () {}

  Tensor3322 ( const Tensor3322<DataType>& tensor )
      : Tensor<Tensor322<DataType>, 3 > ( tensor ) {}

  Tensor3322<DataType>& operator = ( const Tensor3322<DataType>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Tensor322<DataType>, 3 >::operator = ( tensor );
    return *this;
  }

~Tensor3322 () {}

  void apply ( const Mat<2, 2, DataType>& arg, Mat<2, 2, DataType>& dest ) const {
    dest.setZero ();
    applyAdd( arg, dest );
  }

  void applyAdd ( const Mat<2, 2, DataType>& arg, Mat<2, 2, DataType>& dest ) const {
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        for ( int k = 0; k < 2; k++ )
          for ( int l = 0; l < 2; l++ )
            dest [i][j] += this->_vec [i][j][k][l] * arg [k][l];
  }

  const DataType& operator () ( const int i, const int j, const int k, const int l ) const {
    return this->_vec [i][j][k][l];
  }

  DataType& operator () ( const int i, const int j, const int k, const int l ) {
    return this->_vec [i][j][k][l];
  }

  DataType get ( const int i, const int j, const int k, const int l ) const {
      return this->_vec [i][j][k][l];
    }

  void set ( const int i, const int j, const int k, const int l, const DataType value ) {
    this->_vec [i][j][k][l] = value;
  }

  void add ( const int i, const int j, const int k, const int l, const DataType value ) {
    this->_vec [i][j][k][l] += value;
  }
};


template <class DataType>
ostream& operator << ( ostream& os, const Tensor3322<DataType>& t ) {
  return t.print ( os );
}

template <int K, int L, int M, class DataType>
class Tensor3rdOrder : public Tensor<Mat<L, M, DataType>, K > {
public:
  Tensor3rdOrder () : Tensor<Mat<L, M, DataType>, K > () {}

  Tensor3rdOrder ( const Tensor3rdOrder<K, L, M, DataType>& tensor ) : Tensor<Mat<L, M, DataType>, K > ( tensor ) {}

  Tensor3rdOrder<K, L, M, DataType>& operator = ( const Tensor3rdOrder<K, L, M, DataType>& tensor ) {
    // Beware of self-assignment
    if ( this == &tensor ) return *this;
    Tensor<Mat<L, M, DataType>, K >::operator = ( tensor );
    return *this;
  }

  ~Tensor3rdOrder () {}

  const DataType& operator () ( const int i, const int j, const int k ) const {
    return this->_vec [i][j][k];
  }

  DataType& operator () ( const int i, const int j, const int k ) {
    return this->_vec [i][j][k];
  }

  DataType get ( const int i, const int j, const int k ) const {
    return this->_vec [i][j][k];
  }

  void set ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] = value;
  }

  void add ( const int i, const int j, const int k, const DataType value ) {
    this->_vec [i][j][k] += value;
  }
};

//! Elastic tensor
template <class DataType>
struct ElasticTensor : aol::Tensor2222<DataType> {
  enum Type {LAMENAVIER, YOUNGPOISSON};
  ElasticTensor ( DataType c11, DataType c22, DataType c12, DataType c44 );
  ElasticTensor ( DataType c11, DataType c12, DataType c44 );
  ElasticTensor ( DataType c1 = 1, DataType c2 = 1, Type t = LAMENAVIER );
  istream& read ( istream& is );
  DataType C11, C22, C12, C44;
  aol::Matrix22<DataType> C;
  void isoOffset ( DataType off, DataType eps = 0 );
  void rotate ( DataType angle );
  void symmetrize ();
};

// Elastic Tensor
template <class DataType>
void ElasticTensor<DataType>::isoOffset ( DataType off, DataType eps ) {

  // off = abs ( off ); // should be positive for stability

  do {
    C44 *= 1 + off;
    this->set ( 0, 1, 0, 1, C44 );
    this->set ( 1, 0, 1, 0, C44 );
    this->set ( 0, 1, 1, 0, C44 );
    this->set ( 1, 0, 0, 1, C44 );
  } while ( abs (C11 - C12 - 2 * C44) < eps * C11 );
}

template <class DataType>
ElasticTensor<DataType>::ElasticTensor ( DataType c1, DataType c2, Type t ) {
  DataType lambda, mu;
  if ( t == LAMENAVIER ) {
    lambda = c1;
    mu = c2;
  } else {
    DataType E = c1, nu = c2;
    lambda = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
    mu = E / ( 2 * ( 1 + nu ) );
  }
  ElasticTensor temp ( lambda + 2*mu, lambda, mu );
  this->operator= ( temp );
}

template <class DataType>
ElasticTensor<DataType>::ElasticTensor ( DataType c11, DataType c12, DataType c44 )  {
  ElasticTensor temp ( c11, c11, c12, c44 );
  this->operator= ( temp );
}

template <class DataType>
ElasticTensor<DataType>::ElasticTensor ( DataType c11, DataType c22, DataType c12, DataType c44 )
: C11 ( c11 ), C22 ( c22 ), C12 ( c12 ), C44 ( c44 ), C ( C11, C12, C12, C22 ) {
  int i, j, k, l;

  for ( i = 0; i < 2; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < 2; l++ ) {
          this->set ( i, j, k, l, 0 );
          if ( ( i == j ) && ( k == l ) )
            this->set ( i, j, k, l, C [i][k] );
          if ( ( ( ( i == k ) && ( j == l ) ) || ( ( i == l ) && ( j == k ) ) ) && ( i != j ) )
            this->set ( i, j, k, l, C44 );
        }
}

template <class DataType>
istream& ElasticTensor<DataType>::read ( istream& is ) {
  aol::Vec3<DataType> data;
  is >> data;
  ElasticTensor<DataType> temp ( data [0], data [1], data [2] );
  ( *this ) = temp;
  return is;
}

//! change to coordinate system rotated by angle ccw
template <class DataType>
void ElasticTensor<DataType>::rotate ( DataType angle ) {

  // setup rotation matrix
  aol::Matrix22<DataType> rotation;
  rotation.makeRotation( angle );

  // save current entries
  aol::ElasticTensor<DataType> tmp ( *this );
  this->setZero();

  // do multiplication
  for ( int i = 0; i < 2; ++i )
    for ( int j = 0; j < 2; ++j )
      for ( int k = 0; k < 2; ++k )
        for ( int l = 0; l < 2; ++l )
          for ( int m = 0; m < 2; ++m )
            for ( int n = 0; n < 2; ++n )
              for ( int o = 0; o < 2; ++o )
                for ( int p = 0; p < 2; ++p )
                  (*this)[i][j][k][l] += rotation[i][m] * rotation[j][n] * rotation[k][o] * rotation[l][p] * tmp[m][n][o][p];
}

//! symmetrize ( \f$ C_{ijkl}=C_{jikl}=C_{ijlk}=C_{klij} \f$ ) by averaging
template <class DataType>
void ElasticTensor<DataType>::symmetrize() {

  (*this)[0][0][1][1] = (*this)[1][1][0][0] =                                              .5 * ( (*this)[0][0][1][1] + (*this)[1][1][0][0] );
  (*this)[0][1][0][1] = (*this)[1][0][0][1] = (*this)[0][1][1][0] = (*this)[1][0][1][0] = .25 * ( (*this)[0][1][0][1] + (*this)[1][0][0][1] + (*this)[0][1][1][0] + (*this)[1][0][1][0] );
  (*this)[0][0][0][1] = (*this)[0][0][1][0] = (*this)[0][1][0][0] = (*this)[1][0][0][0] = .25 * ( (*this)[0][0][0][1] + (*this)[0][0][1][0] + (*this)[0][1][0][0] + (*this)[1][0][0][0] );
  (*this)[1][1][1][0] = (*this)[1][1][0][1] = (*this)[1][0][1][1] = (*this)[0][1][1][1] = .25 * ( (*this)[1][1][1][0] + (*this)[1][1][0][1] + (*this)[1][0][1][1] + (*this)[0][1][1][1] );
}

template <class DataType>
inline istream& operator >> ( istream& is, ElasticTensor<DataType>& green ) {
  return green.read ( is );
}

} // namespace aol

#endif // __TENSOR_H
