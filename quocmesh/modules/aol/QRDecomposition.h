#ifndef __QRDECOMPOSITION_H
#define __QRDECOMPOSITION_H

#include <aol.h>
#include <matrix.h>

namespace aol {

/*!
 * Class for managing the concatenation of elementary rotation matrices.
 * \author Droske
 */
template <typename RealType>
class ConcatJacobiMatrices {
public:
  ConcatJacobiMatrices() {}

  /*! Append an elementary Jacobi-Rotation matrix to the whole operator.
   * \f[ Q = U_{1}\cdot\ldots\cdot U_{n-1} \mapsto Q U(p,q,\phi) \f]
   */
  void appendJacobiMatrix ( int P, int Q, RealType SinPhi, RealType CosPhi ) {
    matrices.push_back ( jacobi_rot_mat ( P, Q, SinPhi, CosPhi ) );
  }

  //! clears all matrices
  void clear() {
    matrices.erase ( matrices.begin(), matrices.end() );
  }

  /*! Apply the matrices iteratively to the given matrix from the left, i.e.,
   * \f[ M \mapsto U_{n}^T\ldots U_{1}^T M \f]
   */
  void apply ( FullMatrix<RealType> &Mat ) const {
    for ( typename vector<jacobi_rot_mat>::const_iterator it = matrices.begin(); it != matrices.end(); ++it ) {
      applyJacobi ( it->p, it->q, it->s, it->c, Mat );
    }
  }

  void applyTranspose ( FullMatrix<RealType> &Mat ) const {
    //for ( vector<jacobi_rot_mat>::const_iterator it = matrices.end(); it != matrices.begin(); it-- ) {
    for ( int i = static_cast<int> ( matrices.size() - 1 ); i >= 0; i-- ) {
      const jacobi_rot_mat *it = &matrices[i];
      applyTransposeJacobi ( it->p, it->q, it->s, it->c, Mat );
    }
  }

  void applyJacobi ( int P, int Q, RealType Phi, FullMatrix<RealType> &Mat ) const {
    applyJacobi ( P, Q, sin ( Phi ), cos ( Phi ), Mat );
  }

  /*!
   * \f[ a'_{pj} = a_{pj} \cos \phi - a_{qj} \sin \phi \f]
   * \f[ a'_{qj} = a_{pj} \sin \phi + a_{qj} \cos \phi \f]
   * \f[ a'_{ij} = a_{ij} \mbox{ for } i\neq p, q \f]
   */
  void applyJacobi ( int P, int Q, RealType SinPhi, RealType CosPhi, FullMatrix<RealType> &Mat ) const {
    const RealType s = SinPhi;
    const RealType c = CosPhi;

    if ( P >= static_cast<int> ( Mat.getNumRows() ) || Q >= static_cast<int> ( Mat.getNumRows() ) ) {
      throw Exception ( "ConcatJacobiMatrices::applyJacobi: indices P or Q out of bounds for given matrix Mat.", __FILE__, __LINE__ );
    }

    for ( int j = 0; j < Mat.getNumCols(); j++ ) {
      RealType ap = Mat.get ( P, j ) * c - Mat.get ( Q, j ) * s;
      RealType aq = Mat.get ( P, j ) * s + Mat.get ( Q, j ) * c;
      Mat.set ( P, j, ap );
      Mat.set ( Q, j, aq );
    }
  }

  /*!
   * \f[ a'_{pj} = a_{pj} \cos \phi + a_{qj} \sin \phi \f]
   * \f[ a'_{qj} = a_{pj} -\sin \phi + a_{qj} \cos \phi \f]
   * \f[ a'_{ij} = a_{ij} \mbox{ for } i\neq p, q \f]
   */
  void applyTransposeJacobi ( int P, int Q, RealType SinPhi, RealType CosPhi, FullMatrix<RealType> &Mat ) const {
    const RealType s = SinPhi;
    const RealType c = CosPhi;

    if ( P >= static_cast<int> ( Mat.getNumRows() ) || Q >= static_cast<int> ( Mat.getNumRows() ) ) {
      throw Exception ( "ConcatJacobiMatrices::applyJacobi: indices P or Q out of bounds for given matrix Mat.", __FILE__, __LINE__ );
    }

    for ( int j = 0; j < Mat.getNumCols(); j++ ) {
      RealType ap = Mat.get ( P, j ) * c + Mat.get ( Q, j ) * s;
      RealType aq = -Mat.get ( P, j ) * s + Mat.get ( Q, j ) * c;
      Mat.set ( P, j, ap );
      Mat.set ( Q, j, aq );
    }
  }

protected:
  struct jacobi_rot_mat {
    jacobi_rot_mat ( int P, int Q, RealType Phi ) :
        p ( P ), q ( Q ) {
      s = sin ( Phi );
      c = cos ( Phi );
    }

    jacobi_rot_mat ( int P, int Q, RealType SinPhi, RealType CosPhi )
        : p ( P ), q ( Q ),
        s ( SinPhi ), c ( CosPhi ) {}

    int p, q;
    RealType s, c;
  };

  vector<jacobi_rot_mat> matrices;
};

/*!
 * \brief Base class for QR-decomposition transforms based on Given-Rotations.
 *        Provides member function for the elimination of a matrix element.
 * \author Droske
 * \ingroup directSolver
 */
template <typename RealType>
class QRDecomposeGivensBase {
public:
  virtual void transform ( const FullMatrix<RealType> &H, FullMatrix<RealType> &R, FullMatrix<RealType> &Q ) = 0;

  virtual ~QRDecomposeGivensBase() {}

protected:
  void eliminate ( FullMatrix<RealType> &R, int Column, int ToEliminate, int EliminateFrom ) {
    double v1 = R.get ( EliminateFrom, Column );
    double v2 = R.get ( ToEliminate, Column );

    double c = fabs ( v1 ) / sqrt ( Sqr ( v2 ) + Sqr ( v1 ) );
    double s = v2 / sqrt ( Sqr ( v2 ) + Sqr ( v1 ) );

    if ( v1 > 0 ) {
      s = -s;
    }
    jacobi.appendJacobiMatrix ( EliminateFrom, ToEliminate, s, c );
    jacobi.applyJacobi ( EliminateFrom, ToEliminate, s, c, R );
  }


  ConcatJacobiMatrices<RealType> jacobi;
};

//! \brief QR decomposition for a upper triangular matrix
//!        with an additional non-zero secondary diagonal,
//!        NOT a full Householder decomposition!
//! \ingroup directSolver
template <typename RealType>
class QRDecomposeHouseholderMatrix : public QRDecomposeGivensBase<RealType> {
public:
  virtual void transform ( const FullMatrix<RealType> &H, FullMatrix<RealType> &R, FullMatrix<RealType> &Q ) {
    if ( Q.getNumRows() != H.getNumRows() || Q.getNumCols() != H.getNumRows() ) {
      throw Exception ( "Size of Q has to be number of rows of H.", __FILE__, __LINE__ );
    }
    R = H;

    this->jacobi.clear();
    for ( int i = 0; i < H.getNumCols(); i++ ) {
      this->eliminate ( R, i, i + 1, i );
    }

    Q.setIdentity();
    this->jacobi.applyTranspose ( Q );
  }
};

/*!
 * \brief QR-decomposition based on the modified Gram-Schmidt (MGS) process described in
 *        "Efficient Implementation of Minimal Polynomial and Reduced Rank Extrapolation Methods"
 *        by Avram Sidi, August 1990
 *
 * There is not point to derive this from QRDecomposeGivensBase.
 *
 * \author Berkels
 * \ingroup directSolver
 */
template <typename RealType>
class QRDecomposeModifiedGramSchmidt {
public:
  //! Performs a QR decomposition of the Matrix H, the columns of Q are stored in the components of QVectors.
  void transform ( const FullMatrix<RealType> &H, FullMatrix<RealType> &R, aol::MultiVector<RealType> &QVectors ) {
    // For easier handling, we copy the columns of H into the components of the MultiVector QVectors.
    for ( int i = 0; i < H.getNumCols(); i++ )
      H.getColumn ( i, QVectors[i] );
    R.setZero();
    // Step one in Si90
    R.set ( 0, 0, QVectors[0].norm() );
    QVectors[0] /= R.get ( 0, 0 );
    // Step two in Si90
    for ( int k = 1; k < H.getNumCols(); k++ ) {
      for ( int j = 0; j < k; j++ ) {
        R.set ( j, k, QVectors[j] * QVectors[k] );
        QVectors[k].addMultiple ( QVectors[j], ( -1. ) *R.get ( j, k ) );
      }
      R.set ( k, k, QVectors[k].norm() );
      QVectors[k] /= R.get ( k, k );
    }
  }
  //! Performs a QR decomposition of the Matrix H
  void transform ( const FullMatrix<RealType> &H, FullMatrix<RealType> &R, FullMatrix<RealType> &Q ) {
    aol::MultiVector<RealType> qVectors ( H.getNumCols(), H.getNumRows() );
    transform ( H, R, qVectors );
    // Copy the calculated orthonormal vectors into the Matrix Q.
    for ( int i = 0; i < H.getNumCols(); i++ )
      Q.setColumn ( i, qVectors[i] );
  }
  //! Performs a QR decomposition of the Matrix H, but does not return Q.
  //! If you don't need Q this is a little more efficient than using transform( H, R, Q ).
  void transform ( const FullMatrix<RealType> &H, FullMatrix<RealType> &R ) {
    aol::MultiVector<RealType> qVectors ( H.getNumCols(), H.getNumRows() );
    transform ( H, R, qVectors );
  }
};
  
/*!
 * Performs ortho-normalization of a given basis (FullMatrix: each column is a basis vector)
 *
 * \author Mevenkamp
 */
template <typename RealType>
void orthoNormalizeBasis ( aol::FullMatrix<RealType> &Basis ) {
  const int m = Basis.getNumRows ( ), n = Basis.getNumCols ( );
  aol::FullMatrix<RealType> Q ( m, n ), R ( n, n );
  aol::QRDecomposeModifiedGramSchmidt<RealType> qrGramSchmidt;
  qrGramSchmidt.transform ( Basis, R, Q );
  Basis = Q;
}

/*!
 * Performs ortho-normalization of a given basis (MultiVector: each component is a basis vector)
 *
 * \author Mevenkamp
 */
template <typename RealType>
void orthoNormalizeBasis ( aol::MultiVector<RealType> &Basis ) {
  const int n = Basis.numComponents ( );
  if ( n == 0 ) throw aol::Exception ( "Empty basis!", __FILE__, __LINE__ );
  const int m = Basis[0].size ( );
  aol::FullMatrix<RealType> A ( m, n );
  for ( int i = 0; i < m ; ++i )
    for ( int j = 0; j < n ; ++j )
      A.set ( i, j, Basis[j][i] );
  orthoNormalizeBasis ( A );
  for ( int i = 0; i < m ; ++i )
    for ( int j = 0; j < n ; ++j )
      Basis[j][i] = A.get ( i, j );
}

} // namespace aol

#endif
