#ifndef __PROLONGATION_H
#define __PROLONGATION_H

#include <gridBase.h>
#include <sparseMatrices.h>
#include <scalarArray.h>

namespace qc {

template <typename RealType>
class ProlongOp : public aol::BiOp< aol::Vector<RealType> > {
protected:
  const GridDefinition &_coarseGrid;  //!< Coarse grid
  const GridDefinition &_fineGrid;    //!< Fine grid
  aol::OperatorType _opType;
  mutable aol::SparseMatrix< RealType > *mat; // may want to use specially designed matrix for this purpose to reduce computational and memory overhead

public:
  /** Construct the operator to prolongate between the given grids
   */
  ProlongOp ( const GridDefinition &Coarse,
              const GridDefinition &Fine,
              aol::OperatorType OpType = aol::ONTHEFLY ) :
      _coarseGrid ( Coarse ), _fineGrid ( Fine ), _opType ( OpType ), mat ( NULL ) {}

  ~ProlongOp() {
    if ( mat )
      delete ( mat );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( _opType == aol::ONTHEFLY ) {
      std_mg_prolongate ( _fineGrid, Dest, _coarseGrid, Arg );
    }

    if ( _opType == aol::ASSEMBLED ) {
      if ( !mat ) {
        assembleMatrix();
      }
      mat->apply ( Arg, Dest );
    }
  }

  using aol::BiOp<aol::Vector< RealType> >::apply;

  virtual void assembleAddMatrix ( aol::SparseMatrix<RealType> &Mat ) const {
    // const is complete nonsense here. unfortunately, apply and apply add are const, thus require constness here.
    // could allow more general matrix, but want consistency with restriction operator

    assemble_std_prolong_matrix ( Mat );

  }

protected:
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp<aol::Vector< RealType> >::applyAdd;

  void assembleMatrix() const {
    // const is complete nonsense here. unfortunately, apply and apply add are const, thus require constness here.
    if ( mat != NULL )
      delete ( mat );
    mat = new aol::SparseMatrix<RealType> ( _fineGrid.getNumberOfNodes(), _coarseGrid.getNumberOfNodes() );
    assembleAddMatrix ( *mat );
  }

  void assemble_std_prolong_matrix ( aol::SparseMatrix<RealType> &mat ) const ;

  static void std_mg_prolongate ( const GridDefinition &FineGrid,
                                  aol::Vector<RealType> &FineVector,
                                  const GridDefinition &CoarseGrid,
                                  const aol::Vector<RealType> &CoarseVector );

  ProlongOp ( const ProlongOp<RealType>& ); // do not implement
  ProlongOp<RealType>& operator= ( const ProlongOp<RealType>& );
 // end class ProlongOp
};


namespace simplex {

//! Class for prolongation of data on a simplicial grid consistent with qc::simplex::TopologyLookup.
template <typename RealType, qc::Dimension Dim>
class ProlongOp :
  public aol::Op< aol::Vector<RealType> > {

protected:
  const GridDefinition &_coarseGrid;
  const GridDefinition &_fineGrid;

public:
  ProlongOp ( const GridDefinition &CoarseGrid,
              const GridDefinition &FineGrid ) :
    _coarseGrid ( CoarseGrid ),
    _fineGrid ( FineGrid ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    qc::ScalarArray<RealType,Dim> arg( Arg, _coarseGrid, aol::FLAT_COPY ), dest( Dest, _fineGrid, aol::FLAT_COPY );

    for ( int z1 = 0, z2 = 0; z1 < _coarseGrid.getNumZ(); z1++, z2 += 2 )
      for ( int y1 = 0, y2 = 0; y1 < _coarseGrid.getNumY(); y1++, y2 += 2 )
        for ( int x1 = 0, x2 = 0; x1 < _coarseGrid.getNumX(); x1++, x2 += 2 ) {
          qc::CoordType pos1( x1, y1, z1 ), auxPos1, pos2( x2, y2, z2 ), auxPos2, aux;
          // copy values at old nodes
          dest.add( pos2, arg.get( pos1 ) );
          // interpolate values at new nodes on edges
          for ( int i = 0; i < Dim; ++i ) {
            auxPos1 = pos1;
            auxPos1[i] += 1;
            auxPos2 = pos2;
            auxPos2[i] += 1;
            if ( auxPos2[i] < _fineGrid.getSize()[i] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( pos1 ) ) );
          }
          // interpolate values at new nodes on xy-faces
          auxPos1 = pos1;
          auxPos1[0] += 1;
          aux = pos1;
          aux[1] += 1;
          auxPos2 = pos2;
          auxPos2[0] += 1;
          auxPos2[1] += 1;
          if ( auxPos2[0] < _fineGrid.getSize()[0] && auxPos2[1] < _fineGrid.getSize()[1] )
            dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( aux ) ) );
          // in 3d, interpolate values at new nodes on xz- and yz-faces and in the cube-middle
          if ( Dim == qc::QC_3D ) {
            auxPos1 = pos1;
            auxPos1[0] += 1;
            auxPos1[2] += 1;
            auxPos2 = pos2;
            auxPos2[0] += 1;
            auxPos2[2] += 1;
            if ( auxPos2[0] < _fineGrid.getSize()[0] && auxPos2[2] < _fineGrid.getSize()[2] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( pos1 ) ) );
            auxPos1 = pos1;
            auxPos1[1] += 1;
            auxPos1[2] += 1;
            auxPos2 = pos2;
            auxPos2[1] += 1;
            auxPos2[2] += 1;
            if ( auxPos2[1] < _fineGrid.getSize()[1] && auxPos2[2] < _fineGrid.getSize()[2] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( pos1 ) ) );
            auxPos1 = pos1;
            auxPos1[0] += 1;
            aux = pos1;
            aux[1] += 1;
            aux[2] += 1;
            auxPos2 = pos2;
            for ( int i = 0; i < Dim; ++i )
              auxPos2[i] += 1;
            if ( auxPos2[0] < _fineGrid.getSize()[0] && auxPos2[1] < _fineGrid.getSize()[1] && auxPos2[2] < _fineGrid.getSize()[2] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( aux ) ) );
          }
        }
  }
}; // end class ProlongOpSimplex

} // end of namespace simplex

namespace splines {

template < typename RealType >
class ProlongOp : public aol::Op< aol::Vector<RealType> > {
protected:
  const GridDefinition &_coarseGrid;  //!< Coarse grid
  const GridDefinition &_fineGrid;    //!< Fine grid
  
public:
  ProlongOp ( const GridDefinition &Coarse, const GridDefinition &Fine ) : _coarseGrid ( Coarse ), _fineGrid ( Fine ) {
    if ( _fineGrid.getDimOfWorld() != qc::QC_2D )
      throw aol::Exception ( "Unsupported dimension!", __FILE__, __LINE__ );
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const int coarseWidth = _coarseGrid.getWidth() - 1;
    const int fineWidth = _fineGrid.getWidth() - 1;

    int x = 0;
    int y = 0;
    
    qc::ScalarArray < RealType, qc::QC_2D > ArgRef ( Arg, _coarseGrid, aol::FLAT_COPY );
    qc::ScalarArray < RealType, qc::QC_2D > DestRef ( Dest, _fineGrid, aol::FLAT_COPY );
       
    // Interpolate corners
    DestRef.add ( 0, 0, ArgRef.get ( 0, 0 ) );
    DestRef.add ( fineWidth, 0, ArgRef.get ( coarseWidth, 0 ) ); 
    DestRef.add ( 0, fineWidth, ArgRef.get ( 0, coarseWidth ) );
    DestRef.add ( fineWidth, fineWidth, ArgRef.get ( coarseWidth, coarseWidth ) );
    
    
    // Interpolate left boundary
    for ( y = 1; y < fineWidth; ++y ) {
      if( (y % 2) == 0 ){
        DestRef.add ( 0, y, 0.75 * ArgRef.get ( 0, y / 2 ) );
        DestRef.add ( 0, y, 0.125 * ArgRef.get ( 0, y / 2 + 1 ) );
        DestRef.add ( 0, y, 0.125 * ArgRef.get ( 0, y / 2 - 1 ) );
      }else{
        DestRef.add ( 0, y, 0.5 * ArgRef.get ( 0, y / 2 ) );
        DestRef.add ( 0, y, 0.5 * ArgRef.get ( 0, y / 2 + 1 ) );
      }
    }
    // Interpolate right boundary
    for ( y = 1; y < fineWidth; ++y ) {
      if( (y % 2) == 0 ){
        DestRef.add ( fineWidth, y, 0.75 * ArgRef.get ( coarseWidth, y / 2 ) );
        DestRef.add ( fineWidth, y, 0.125 * ArgRef.get ( coarseWidth, y / 2 + 1 ) );
        DestRef.add ( fineWidth, y, 0.125 * ArgRef.get ( coarseWidth, y / 2 - 1 ) );
      }else{
        DestRef.add ( fineWidth, y, 0.5 * ArgRef.get ( coarseWidth, y / 2 ) );
        DestRef.add ( fineWidth, y, 0.5 * ArgRef.get ( coarseWidth, y / 2 + 1 ) );
      }
    }    
    // Interpolate bottom boundary
    for ( x = 1; x < fineWidth; ++x ) {
      if( (x % 2) == 0 ){
        DestRef.add ( x, 0, 0.75 * ArgRef.get ( x/2, 0 ) );
        DestRef.add ( x, 0, 0.125 * ArgRef.get ( x / 2 + 1, 0 ) );
        DestRef.add ( x, 0, 0.125 * ArgRef.get ( x / 2 - 1, 0 ) );
      }else{
        DestRef.add ( x, 0, 0.5 * ArgRef.get ( x / 2, 0 ) );
        DestRef.add ( x, 0, 0.5 * ArgRef.get ( x / 2 + 1, 0 ) );
      }
    }    
    // Interpolate top boundary
    for ( x = 1; x < fineWidth; ++x ) {
      if( (x % 2) == 0 ){
        DestRef.add ( x, fineWidth, 0.75 * ArgRef.get ( x/2, coarseWidth ) );
        DestRef.add ( x, fineWidth, 0.125 * ArgRef.get ( x / 2 + 1, coarseWidth ) );
        DestRef.add ( x, fineWidth, 0.125 * ArgRef.get ( x / 2 - 1, coarseWidth ) );
      }else{
        DestRef.add ( x, fineWidth, 0.5 * ArgRef.get ( x / 2, coarseWidth ) );
        DestRef.add ( x, fineWidth, 0.5 * ArgRef.get ( x / 2 + 1, coarseWidth ) );
      }
    }
    
    const RealType F1O4 = aol::ZOTrait<RealType>::one / static_cast<RealType> ( 4 );
    const RealType F1O16 = aol::ZOTrait<RealType>::one / static_cast<RealType> ( 16 );
    const RealType F3O8 = static_cast<RealType> ( 3 ) / static_cast<RealType> ( 8 );
    const RealType F1O64 = aol::ZOTrait<RealType>::one / static_cast<RealType> ( 64 );
    const RealType F3O32 = static_cast<RealType> ( 3 ) / static_cast<RealType> ( 32 );
    const RealType F9O16 = static_cast<RealType> ( 9 ) / static_cast<RealType> ( 16 );
    
    // Interpolate interior points
    for ( x = 1; x < fineWidth; ++x ) {
      for ( y = 1; y < fineWidth; ++y ) {
        const bool xx = ( x % 2 == 0 );
        const bool yy = ( y % 2 == 0 );
        if ( xx ) {
          if ( yy ) {
            DestRef.add ( x, y, F1O64 * ( ArgRef.get ( x/2 - 1, y/2 - 1 ) + ArgRef.get ( x/2 - 1, y/2 + 1 ) + ArgRef.get ( x/2 + 1, y/2 - 1 ) + ArgRef.get ( x/2 + 1, y/2 + 1 ) ) 
              + F3O32 * ( ArgRef.get ( x/2 , y/2 - 1 ) + ArgRef.get ( x/2, y/2 + 1 ) + ArgRef.get ( x/2 - 1, y/2 ) + ArgRef.get ( x/2 + 1, y/2 ) ) 
              + F9O16 * ArgRef.get ( x/2, y/2 ) );
          }
          else {
            DestRef.add ( x, y, F1O16 * ( ArgRef.get ( x/2 - 1, y/2 ) + ArgRef.get ( x/2 + 1, y/2 ) + ArgRef.get ( x/2 - 1, y/2 + 1 ) + ArgRef.get ( x/2 + 1, y/2 + 1 ) ) 
              + F3O8 * ( ArgRef.get ( x/2 , y/2 ) + ArgRef.get ( x/2, y/2 + 1 ) ) ); 
          }
        }
        else {
          if ( yy ) {
            DestRef.add ( x, y, F1O16 * ( ArgRef.get ( x/2, y/2 - 1 ) + ArgRef.get ( x/2, y/2 + 1 ) + ArgRef.get ( x/2 + 1, y/2 - 1 ) + ArgRef.get ( x/2 + 1, y/2 + 1 ) ) 
              + F3O8 * ( ArgRef.get ( x/2 , y/2 ) + ArgRef.get ( x/2 + 1, y/2 ) ) );            
          }
          else {
            DestRef.add ( x, y, F1O4 * ( ArgRef.get ( x/2, y/2 ) + ArgRef.get ( x/2, y/2 + 1 ) + ArgRef.get ( x/2 + 1, y/2 ) + ArgRef.get ( x/2 + 1, y/2 + 1 ) ) );
          }
        }
      }
    }
  }
};

// in 2D the same as above, in 3D approximate prolongation op
template < typename CoarseConfiguratorType, typename FineConfiguratorType, qc::Dimension Dim >
class GeneralProlongOp : public aol::Op< aol::Vector<typename FineConfiguratorType::RealType> > {
};

template < typename CoarseConfiguratorType, typename FineConfiguratorType  >
class GeneralProlongOp < CoarseConfiguratorType, FineConfiguratorType, qc::QC_2D > : public aol::Op< aol::Vector<typename FineConfiguratorType::RealType> > {
  typedef typename FineConfiguratorType::RealType RealType;
  ProlongOp < typename FineConfiguratorType::RealType > _prolongOp;

public:
  GeneralProlongOp ( const typename CoarseConfiguratorType::InitType& CoarseGrid, const typename FineConfiguratorType::InitType& FineGrid )
: _prolongOp ( CoarseGrid, FineGrid ) { }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _prolongOp.applyAdd ( Arg, Dest );
  }
};

template < typename CoarseConfiguratorType, typename FineConfiguratorType >
class GeneralProlongOp < CoarseConfiguratorType, FineConfiguratorType, qc::QC_3D > : public aol::Op< aol::Vector<typename FineConfiguratorType::RealType> > {
  typedef typename FineConfiguratorType::RealType RealType;
  const CoarseConfiguratorType _coarseConfig;
  const FineConfiguratorType _fineConfig;

public:
  GeneralProlongOp ( const typename CoarseConfiguratorType::InitType& CoarseGrid, const typename FineConfiguratorType::InitType& FineGrid )
: _coarseConfig ( CoarseGrid ), _fineConfig ( FineGrid )  { }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector < RealType > argNodal ( Arg, aol::STRUCT_COPY );
    _coarseConfig.convertSplineCoefficientsToNodalValues ( Arg, argNodal );
    aol::Vector < RealType > destNodal ( Dest, aol::STRUCT_COPY );
    qc::ProlongOp < RealType > ( _coarseConfig.getInitializer (), _fineConfig.getInitializer () ).apply ( argNodal, destNodal );
    _fineConfig.convertSplineCoefficientsToNodalValuesAdd ( destNodal, Dest );
  }

};
  
} // end of namespace splines
  
} // end of namespace qc

#endif
