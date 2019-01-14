#ifndef __RESTRICTION_H
#define __RESTRICTION_H

#include <gridBase.h>
#include <sparseMatrices.h>
#include <scalarArray.h>

namespace qc {

enum RestrictType {
  STD_MG_RESTRICT,         // standard multigrid restriction: unscaled transpose of prolongation
  STD_QUOC_RESTRICT,       // former "standard" restriction: row-wise scaled transpose of prolongation, maps all-1 vector to all-1 vector, can be used for downsampling images
  THROW_AWAY_RESTRICT,     // restriction by throwing away values at the intermediate grid points
  PERIODIC_RESTRICT        // restricted image is periodic
};

template <typename RealType, RestrictType RestrType >
class RestrictOp : public aol::BiOp< aol::Vector<RealType> > {

protected:
  const GridDefinition &_coarseGrid;  //!< Coarse grid
  const GridDefinition &_fineGrid;    //!< Fine grid
  aol::OperatorType _opType;
  mutable aol::SparseMatrix<RealType>* mat; // may want to use specially designed matrix for this purpose to reduce computational and memory overhead

public:
  /** Construct the operator to restrict between the given grids.
   *  There are various approaches one can use for restriction, this behavior can be controlled.
   *  As the standard behavior for quoc grids implemented earlier is not the standard behavior for multigrid methods, there is no default here.
   *  \author Schwen (based on code by others)
   */
  RestrictOp ( const GridDefinition &Coarse,
               const GridDefinition &Fine,
               aol::OperatorType OpType = aol::ONTHEFLY ) :
      _coarseGrid ( Coarse ), _fineGrid ( Fine ), _opType ( OpType ), mat ( NULL ) {

    if ( _fineGrid.getGridDepth() != ( _coarseGrid.getGridDepth() + 1 ) ) {
      cerr << "Incompatible grid depths! Fine level = " << _fineGrid.getGridDepth() << ", coarse level = " << _coarseGrid.getGridDepth() << endl;
      throw aol::Exception ( "qc::RestrictOp: Incompatible grid levels", __FILE__, __LINE__ );
    }

    if ( _fineGrid.getDimOfWorld() != _coarseGrid.getDimOfWorld() ) {
      throw aol::Exception ( "qc::RestrictOp: Incompatible dimensions of coarse grid", __FILE__, __LINE__ );
    }

  }

  ~RestrictOp() {
    if ( mat )
      delete ( mat );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    Dest.setZero(); // might be unnecessary
    if ( _opType == aol::ONTHEFLY ) {
      switch ( RestrType ) {
      case STD_MG_RESTRICT :
        std_mg_restrict ( Arg, Dest );
        break;
      case STD_QUOC_RESTRICT:
        std_quoc_restrict ( Arg, Dest );
        break;
      case THROW_AWAY_RESTRICT:
        throw_away_restrict ( Arg, Dest );
        break;
      case PERIODIC_RESTRICT:
        periodic_restrict ( Arg, Dest );
        break;
      default:
        throw aol::Exception ( "qc::RestrictOp::apply: Invalid RestrictType", __FILE__, __LINE__ );
      }
    }

    if ( _opType == aol::ASSEMBLED ) {
      if ( !mat ) {
        assembleMatrix();
      }
      mat->apply ( Arg, Dest );
    }

  }

  using aol::BiOp<aol::Vector<RealType> >::apply;

  void assembleAddMatrix ( aol::SparseMatrix<RealType> &matrix ) const {

    aol::SparseMatrix<RealType> matToAdd ( matrix.getNumRows(), matrix.getNumCols() );

    if ( matToAdd.getNumRows() != _coarseGrid.getNumberOfNodes() || matToAdd.getNumCols() != _fineGrid.getNumberOfNodes() ) {
      throw aol::Exception ( "qc::RestrictOp::assembleMatrix: Incompatible matrix", __FILE__, __LINE__ );
    }

    switch ( RestrType ) {
    case STD_MG_RESTRICT :
      assemble_std_mg_restrict_matrix ( matToAdd );
      break;
    case STD_QUOC_RESTRICT:
      assemble_std_quoc_restrict_matrix ( matToAdd );
      break;
    case THROW_AWAY_RESTRICT:
      assemble_throw_away_restrict_matrix ( matToAdd );
      break;
    default:
      throw aol::Exception ( "qc::RestrictOp::assembleAddMatrix: Invalid RestrictType", __FILE__, __LINE__ );
    }

    matrix += matToAdd;

  }


protected:
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp<aol::Vector<RealType> >::applyAdd;

  void assembleMatrix() const {
    // const is complete nonsense here. unfortunately, apply and apply add are const, thus require constness here.
    if ( mat != NULL )
      delete ( mat );
    mat = new aol::SparseMatrix<RealType> ( _coarseGrid.getNumberOfNodes(), _fineGrid.getNumberOfNodes() );
    assembleAddMatrix ( *mat );
  }

  void assemble_std_mg_restrict_matrix ( aol::SparseMatrix<RealType> &mat ) const ;

  void assemble_std_quoc_restrict_matrix ( aol::SparseMatrix<RealType> &mat ) const;

  void assemble_throw_away_restrict_matrix ( aol::SparseMatrix<RealType> &mat ) const ;


  void std_mg_restrict ( const aol::Vector<RealType> &fineVector, aol::Vector<RealType> &coarseVector ) const ;

  void std_quoc_restrict ( const aol::Vector<RealType> &fineVector, aol::Vector<RealType> &coarseVector ) const ;

  void throw_away_restrict ( const aol::Vector<RealType> &fineVector, aol::Vector<RealType> &coarseVector ) const ;
  
  void periodic_restrict ( const aol::Vector<RealType> &fineVector, aol::Vector<RealType> &coarseVector ) const ;

  RestrictOp ( const RestrictOp< RealType, RestrType > & ); // do not implement
  RestrictOp< RealType, RestrType >& operator= ( const RestrictOp< RealType, RestrType > & );

 // end of class qc::RestrictOp
};

}

#endif
