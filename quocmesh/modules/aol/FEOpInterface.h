#ifndef __FEOPINTERFACE_H
#define __FEOPINTERFACE_H

#include <aol.h>
#include <op.h>
#include <baseFunctionSet.h>
#include <sparseMatrices.h>
#include <fastUniformGridMatrix.h>
#include <discreteFunction.h>
#include <maskedVector.h>
#include <array.h>
#include <smallVec.h>
#include <pointerClasses.h>
#include <elementMask.h>

namespace aol {

template <typename ConfiguratorType, typename _DomainType, typename _RangeType=_DomainType>
class FEOpInterface : public Op<_DomainType, _RangeType> {

public:
  FEOpInterface ( const typename ConfiguratorType::InitType & Grid )
    : _config ( new ConfiguratorType ( Grid ), true )
  {}

  FEOpInterface ( const ConfiguratorType & Configurator )
    : _config ( &Configurator, false )
  {}

  const ConfiguratorType & getConfigurator() const {
    return *_config;
  }

  const typename ConfiguratorType::BaseFuncSetType &getBaseFunctionSet ( const typename ConfiguratorType::ElementType &El ) const {
    return getConfigurator().getBaseFunctionSet ( El );
  }

  int getNumGlobalDofs( ) const {
    return getConfigurator().getNumGlobalDofs( );
  }

protected:
  DeleteFlagPointer<const ConfiguratorType> _config;
};


//! Helper class for local assembly making partial specialization possible
template <typename FEOpType, typename MatrixType, GridGlobalIndexMode IndexMode>
struct LocalAssemblyHelper
{};

//! Local assembly helper class: partial specialization for QUOC_GRID_INDEX_MODE
template <typename FEOpType, typename MatrixType>
struct LocalAssemblyHelper<FEOpType, MatrixType, QUOC_GRID_INDEX_MODE>
{
public:
  typedef typename FEOpType::RealType                                                                    RealType;
  typedef typename FEOpType::ConfiguratorType                                                            ConfiguratorType;
  typedef typename FEOpType::VectorType                                                                  VectorType;
  typedef typename ConfiguratorType::ElementIteratorType                                                 IteratorType;
  typedef typename IteratorType::EndType                                                                 IteratorEndType;
  typedef aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType>         MatType;
  typedef qc::ElementMask<typename ConfiguratorType::InitType, ConfiguratorType, ConfiguratorType::Dim>  ElementMaskType;

  static void doLocalAssembly ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                MatrixType &Mat, const RealType Factor )  {

    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix for the current element
      feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, i );
      }

      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[ j ];

          Mat.add ( glob_i, glob_j, Factor * localMatrix [ i ][ j ] );
        }
      }
    }
  }

  static void doLocalAssemblyDirichlet ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                         MatrixType &Mat, const aol::BitVector * DirichletMaskRow, bool setDirichletNodes,
                                         const ElementMaskType * elMask, const aol::BitVector * DirichletMaskCol )  {

    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      //skip elements that are not in the desired region determined by elementMapper, if elementMapper != NULL
      if( elMask != NULL )
        if( !(elMask->isInDesiredRegion(*it)) )
          continue;

      // assemble the local matrix for the current element
        feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );

        const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, i );
      }

      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        // write the ith row only if the ith node is not a Dirichlet node
        if ( ! (*DirichletMaskRow)[glob_i] ) {
          for ( int j = 0; j < numLocalDofs; ++j ) {
            int glob_j = globalDofs[ j ];
            // write the jth column only if the jth node is not a Dirichlet node
            if ( ! (*DirichletMaskCol)[glob_j] )
              Mat.add ( glob_i, glob_j, localMatrix [ i ][ j ] );
          }
        }
      }
    }

    // set ones on the diagonal for Dirichlet nodes
    if ( setDirichletNodes )
      for ( int i = 0; i < DirichletMaskRow->size(); ++i )
        if ( ( (*DirichletMaskRow)[i] ) && ( (*DirichletMaskCol)[i] ) )
          Mat.add( i, i, aol::NumberTrait<RealType>::one );
  }

  // valid for QUOC_GRID_INDEX_MODE and DT_GRID_INDEX_MODE, specialization for QUOC_ADAPTIVEGRID_INDEX_MODE below
  static void doMultiplyOnTheFly ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                   const VectorType &Arg, VectorType &Dest )  {

    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix belonging to the current element
      feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices to the current Dofs
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, i );
      }
      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[ j ];

          Dest[ glob_i ] += localMatrix[i][j] * Arg[ glob_j ] ;
        }
      }
    }
  }
};

//! Local assembly helper class: partial specialization for DT_GRID_INDEX_MODE
template <typename FEOpType, typename MatrixType>
struct LocalAssemblyHelper<FEOpType, MatrixType, DT_GRID_INDEX_MODE>
{
public:
  typedef typename FEOpType::RealType                                                             RealType;
  typedef typename FEOpType::ConfiguratorType                                                     ConfiguratorType;
  typedef typename ConfiguratorType::ElementIteratorType                                          IteratorType;
  typedef typename IteratorType::EndType                                                          IteratorEndType;
  typedef aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType>  MatType;
  typedef qc::ElementMask<typename ConfiguratorType::InitType, ConfiguratorType, ConfiguratorType::Dim>  ElementMaskType;

  static void doLocalAssembly ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                MatrixType &Mat, const RealType Factor )  {

    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix for the current element
      feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );

      Vec2<int> glob;
      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i ) {
        for ( int j = 0; j < numLocalDofs; ++j ) {
          feopPtr->getConfigurator().localToGlobal ( *it, i, j, glob );

          Mat.add ( glob[0], glob[1], Factor * localMatrix [ i ][ j ] );
        }
      }
    }
  }

  static void doLocalAssemblyDirichlet ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                         MatrixType &Mat, const aol::BitVector * DirichletMaskRow, bool setDirichletNodes,
                                         const ElementMaskType * /*elMask*/, const aol::BitVector * DirichletMaskCol )  {

    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix for the current element
      feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices of the local Dofs of the current element
      // Here we can use the usual localToGlobal, because i describes the index
      // of a row.
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, i );
      }

#if 0
      // finally add the locally computed values to the matrix
      Vec2<int> glob;
      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        // write the ith row only if the ith node is not a Dirichlet node
        if ( !(*DirichletMask)[glob_i] ) {
          for ( int j = 0; j < numLocalDofs; ++j ) {
            feopPtr->getConfigurator().localToGlobal ( *it, i, j, glob );

            // write the jth column only if the jth node is not a Dirichlet node
            //               if( !DirichletMask[ globalDofs[j] ] )
            if ( aol::isNaN ( localMatrix[i][j] ) )
              throw aol::Exception ( "FELinOpInterface: localMatrix-entry is nan!", __FILE__, __LINE__ );
            Mat.add ( glob[0], glob[1], localMatrix[i][j] );
          }
        }
      }
      // set 1.0 on the diagonal if ith node is a Dirichlet node
      if ( setDirichletNodes )
        for ( int i = 0; i < numLocalDofs; ++i ) {
          int glob_i = globalDofs[ i ];
          if ( (*DirichletMask)[glob_i] )
            Mat.add ( glob_i, glob_i, 1.0 );
        }
#else
      // this block doesn't handle Dirichlet values at all, because the maskedDTMatrix
      // already handles this.

      // finally add the locally computed values to the matrix
      Vec2<int> glob;
      for ( int i = 0; i < numLocalDofs; ++i ) {
        // write the ith row only if the ith node is not a Dirichlet node
        for ( int j = 0; j < numLocalDofs; ++j ) {
          feopPtr->getConfigurator().localToGlobal ( *it, i, j, glob );
          Mat.add ( glob[0], glob[1], localMatrix[i][j] );
        }
      }
#endif
    }
  }

  // valid for QUOC_GRID_INDEX_MODE and DT_GRID_INDEX_MODE, specialization for QUOC_ADAPTIVEGRID_INDEX_MODE below
  static void doMultiplyOnTheFly ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                   const Vector<RealType> &Arg, Vector<RealType> &Dest )  {
    LocalAssemblyHelper<FEOpType, MatrixType, QUOC_GRID_INDEX_MODE>::doMultiplyOnTheFly( feopPtr, localMatrix, end_it, Arg, Dest );
  }
};

//! Local assembly helper class: partial specialization for QUOC_ADAPTIVEGRID_INDEX_MODE
template <typename FEOpType, typename MatrixType>
struct LocalAssemblyHelper<FEOpType, MatrixType, QUOC_ADAPTIVEGRID_INDEX_MODE>
{
public:
  typedef typename FEOpType::RealType                                                             RealType;
  typedef typename FEOpType::ConfiguratorType                                                     ConfiguratorType;
  typedef typename ConfiguratorType::ElementIteratorType                                          IteratorType;
  typedef typename IteratorType::EndType                                                          IteratorEndType;
  typedef aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType>  MatType;
  typedef qc::ElementMask<typename ConfiguratorType::InitType, ConfiguratorType, ConfiguratorType::Dim>  ElementMaskType;

  static void doLocalAssembly ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                MatrixType &Mat, const RealType Factor )  {

    if( ConfiguratorType::DimOfWorld != qc::QC_2D){
      throw Exception ( "FEOpInterface:assembleAddMatrix adaptive grids are implemented for 2d only.\n", __FILE__, __LINE__ );
    }
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    int globalDofsConstraint[ ConfiguratorType::maxNumLocalDofs ];
    int globalDofsConstraintOtherElem[ ConfiguratorType::maxNumLocalDofs ];

    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix for the current element
      feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );
      const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );
      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, i );
        globalDofsConstraintOtherElem[i] = globalDofs[i];
        globalDofsConstraint[i] = globalDofs[i];
        if(feopPtr->getConfigurator().getInitializer().isAdaptive()){
          const int elSize = 1 << ( feopPtr->getConfigurator().getInitializer().getGridDepth() - (*it).level() );
          int hanging = feopPtr->getConfigurator().getInitializer().checkForHangingNode(*it,i);
          if(hanging == i){
            globalDofsConstraintOtherElem[i] = globalDofs[i];
          }
          else{
            // constraining nodes in x direction
            if( ( ( (*it).type() == 0) && (i == 1) ) || (((*it).type() == 1) && (i == 0)) || (((*it).type() == 2) && (i == 3)) || (((*it).type() == 3) && (i == 2)) ){
              globalDofsConstraint[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, hanging );
              if( ((*it).type() == 0) || ((*it).type() == 2)){
                globalDofsConstraintOtherElem[i] = globalDofs[i] + 1 * elSize;
              }
              else{
                globalDofsConstraintOtherElem[i] = globalDofs[i] - 1 * elSize;
              }
            }
            else if((((*it).type() == 0) && (i == 2)) || (((*it).type() == 1) && (i == 3)) || (((*it).type() == 2) && (i == 0)) || (((*it).type() == 3) && (i == 1))){
              globalDofsConstraint[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, hanging );
              if(  (*it).type() == 0 || (*it).type() == 1){
                globalDofsConstraintOtherElem[i] = globalDofs[i] + 1 * feopPtr->getConfigurator().getInitializer().getNumX() * elSize;
              }
              else{
                globalDofsConstraintOtherElem[i] = globalDofs[i] - 1 * feopPtr->getConfigurator().getInitializer().getNumX() * elSize;
              }
            }
            else{
              throw Exception ( "FEOpInterface::assembleAddMatrix unknown hanging node configuration.\n", __FILE__, __LINE__ );
            }
          }
        }
      }

      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i ) {
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofsConstraint[ j ];
          int glob_i = globalDofsConstraint[ i ];

          Mat.add ( glob_i, glob_j, 0.25 * Factor * localMatrix [ i ][ j ] );
          Mat.add ( globalDofsConstraintOtherElem[i], glob_j, 0.25 * Factor * localMatrix [ i ][ j ] );
          Mat.add ( glob_i, globalDofsConstraintOtherElem[j], 0.25 * Factor * localMatrix [ i ][ j ] );
          Mat.add ( globalDofsConstraintOtherElem[i], globalDofsConstraintOtherElem[j], 0.25 * Factor * localMatrix [ i ][ j ] );
        }
      }
    }
  }

  static void doLocalAssemblyDirichlet ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                         MatrixType &Mat, const aol::BitVector * DirichletMaskRow, bool setDirichletNodes,
                                         const ElementMaskType * /*elMask*/, const aol::BitVector * DirichletMaskCol )  {
    throw aol::UnimplementedCodeException( "FELinOpInterface::assembleAddMAtrix: specialization for QUOC_ADAPTIVEGRID_INDEX_MODE missing.\n", __FILE__, __LINE__ );
  }

  static void doMultiplyOnTheFly ( const FEOpType * feopPtr, MatType & localMatrix, const IteratorEndType & end_it,
                                   const Vector<RealType> &Arg, Vector<RealType> &Dest )  {

    if( ConfiguratorType::DimOfWorld != qc::QC_2D){
      throw Exception ( "FEOpInterface:multiplyOnTheFly adaptive grids are implemented for 2d only.\n", __FILE__, __LINE__ );
    }
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    int globalDofsConstraint[ ConfiguratorType::maxNumLocalDofs ];
    int globalDofsConstraintOtherElem[ ConfiguratorType::maxNumLocalDofs ];

    for ( IteratorType it = feopPtr->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix belonging to the current element
      feopPtr->asImp().prepareLocalMatrix ( *it, localMatrix );
      const int numLocalDofs = feopPtr->getConfigurator().getNumLocalDofs ( *it );
      // get the global indices to the current Dofs
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, i );
        globalDofsConstraintOtherElem[i] = globalDofs[ i ];
        globalDofsConstraint[i] = globalDofs[ i ];
        if(feopPtr->getConfigurator().getInitializer().isAdaptive()){
          const int elSize = 1 << ( feopPtr->getConfigurator().getInitializer().getGridDepth() - (*it).level() );
          int hanging = feopPtr->getConfigurator().getInitializer().checkForHangingNode(*it,i);
          if(hanging == i){
            globalDofsConstraintOtherElem[i] = globalDofs[i];
          }
          else{
            // contraining nodes in x direction
            if(((*it).type() == 0 && i == 1) || ((*it).type() == 1 && i == 0) || ((*it).type() == 2 && i == 3) || ((*it).type() == 3 && i == 2) ){
              globalDofsConstraint[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, hanging );
              if( (*it).type() == 0 || (*it).type() == 2){
                globalDofsConstraintOtherElem[i] = globalDofs[i] + 1 * elSize;
              }
              else{
                globalDofsConstraintOtherElem[i] = globalDofs[i] - 1 * elSize;
              }
            }
            else if(((*it).type() == 0 && i == 2) || ((*it).type() == 1 && i == 3) || ((*it).type() == 2 && i == 0) || ((*it).type() == 3 && i == 1)){
              globalDofsConstraint[ i ] = feopPtr->getConfigurator().localToGlobal ( *it, hanging );
              if(  (*it).type() == 0 || (*it).type() == 1){
                globalDofsConstraintOtherElem[i] = globalDofs[i] + 1 * feopPtr->getConfigurator().getInitializer().getNumX() * elSize;
              }
              else{
                globalDofsConstraintOtherElem[i] = globalDofs[i] - 1 * feopPtr->getConfigurator().getInitializer().getNumX() * elSize;
              }
            }
            else{
              throw Exception ( "FEOpInterface:multiplyOnTheFly unknown hanging node configuration.\n", __FILE__, __LINE__ );
            }
          }
        }
      }
      for ( int i = 0; i < numLocalDofs; ++i ) {
        for ( int j = 0; j < numLocalDofs; ++j ) {
          Dest[ globalDofsConstraintOtherElem[i] ] += localMatrix[i][j] * (Arg[ globalDofsConstraint[j] ] + Arg[globalDofsConstraintOtherElem[j]]) / 2.0 / 2.0 ;
          Dest[ globalDofsConstraint[i] ] += localMatrix[i][j] * (Arg[ globalDofsConstraint[j] ] + Arg[globalDofsConstraintOtherElem[j]]) / 2.0 / 2.0 ;
        }
      }
    }
  }
};


//! General Interface for efficient Finite Element operators for all types of meshes and basis functions.
/*!
 * \author Droske
 */
template <typename _RealType, typename _ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinOpInterface : public FEOpInterface<_ConfiguratorType, typename _ConfiguratorType::VectorType> {
public:
  typedef _RealType         RealType;
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::VectorType VectorType;

  explicit FELinOpInterface ( const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FEOpInterface<ConfiguratorType, VectorType> ( Grid ), _mat ( NULL ), _opType ( OpType ) {}

  explicit FELinOpInterface ( const ConfiguratorType &Config, OperatorType OpType = ONTHEFLY )
      : FEOpInterface<ConfiguratorType, VectorType> ( Config ), _mat ( NULL ), _opType ( OpType ) {}

  virtual ~FELinOpInterface( ) {
    delete _mat;
  }

  //! clears the assembled matrix
  void reset( ) {
    if ( _mat ) {
      delete _mat;
    }
    _mat = NULL;
  }

  OperatorType getOpType() const {
    return _opType;
  }

  template<typename BitMaskFunctorType>
  void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    switch ( _opType ) {
    case ONTHEFLY:
      multiplyOnTheFly<BitMaskFunctorType> ( Arg, Dest );
      break;
    case ASSEMBLED:
#ifdef _OPENMP
#pragma omp critical (aol_FEOpInterface_callAssembleMatrix)
#endif
      if ( !_mat ) {
        assembleMatrix( );
      }
#if defined (__GNUC__)
      _mat->template applyAdd<BitMaskFunctorType> ( Arg, Dest );
#else
      _mat->applyAdd<BitMaskFunctorType> ( Arg, Dest );
#endif
      break;
    default:
      throw aol::UnimplementedCodeException ( "FELinOpInterface::applyAdd: unsupported opType", __FILE__, __LINE__ );
    }
  }


  void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    switch ( _opType ) {
    case ONTHEFLY:
      multiplyOnTheFly ( Arg, Dest );
      break;
    case ASSEMBLED:
#ifdef _OPENMP
#pragma omp critical (aol_FEOpInterface_callAssembleMatrix)
#endif
      if ( !_mat ) {
        assembleMatrix( );
      }
      _mat->applyAdd ( Arg, Dest );
      break;
    default:
      throw aol::UnimplementedCodeException ( "FELinOpInterface::applyAdd: unsupported opType", __FILE__, __LINE__ );
    }
  }

  typename ConfiguratorType::MatrixType& getMatrix( ) const {
    if ( !_mat ) {
      assembleMatrix( );
    }
    //return dynamic_cast<typename ConfiguratorType::MatrixType&>(*_mat);
    return *_mat;
  }

  void makeDiagonal ( DiagonalMatrix<RealType> &Mat ) {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    Mat.setZero();

    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    const typename IteratorType::EndType end_it = asImp().end();

    for ( IteratorType it = asImp().begin(); it != end_it; ++it ) {
      this->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = this->asImp().getNumLocalDofs ( *it );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int global_i = this->asImp().localToGlobal ( *it, i );
        Mat.add ( global_i, global_i, localMatrix [ i ][ i ] );
      }
    }
  }

protected:
  void multiplyOnTheFly ( const VectorType &Arg, VectorType &Dest ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    // traverse the elements of the grid
    // Note that the cases QUOC_GRID_INDEX_MODE and DT_GRID_INDEX_MODE should be treated identically
    // in the multiplyOnTheFly method
    LocalAssemblyHelper<FELinOpInterface<RealType,ConfiguratorType,Imp,IndexMode>,typename ConfiguratorType::MatrixType,IndexMode>::doMultiplyOnTheFly( this, localMatrix, end_it, Arg, Dest );
  }

  template<typename BitMaskFunctorType>
  void multiplyOnTheFly ( const VectorType &Arg, VectorType &Dest ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    typedef aol::DofMask::iterator *DofMaskIteratorPtr;

    BitMaskFunctorType bitMaskFunctor;
    const aol::DofMask& dofMask = this->getConfigurator().getNodeBoundaryMask();
    DofMaskIteratorPtr *dofMaskStencil;

    dofMaskStencil = new DofMaskIteratorPtr[ ConfiguratorType::maxNumLocalDofs ];
    for ( int j = 0; j < ConfiguratorType::maxNumLocalDofs; ++j ) {
      dofMaskStencil[j] = new aol::DofMask::iterator ( dofMask.begin() );
    }

    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    // traverse the elements of the grid
    // Note that the cases QUOC_GRID_INDEX_MODE and DT_GRID_INDEX_MODE should be treated identically
    // in the multiplyOnTheFly method
    switch ( IndexMode ) {
      case QUOC_GRID_INDEX_MODE:
      case DT_GRID_INDEX_MODE:
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
          // assemble the local matrix belonging to the current element
          this->asImp().prepareLocalMatrix ( *it, localMatrix );

          const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

          // get the global indices to the current Dofs
          for ( int i = 0; i < numLocalDofs; ++i ) {
            globalDofs[ i ] = this->getConfigurator().localToGlobal ( *it, i );
            dofMaskStencil[i]->incrementUntilIndex ( globalDofs[i] );
          }
          for ( int i = 0; i < numLocalDofs; ++i ) {
            int glob_i = globalDofs[ i ];
            // The bitMaskFunctor is evaluated on the expression:
            // **(dofMaskStencil[i]) == globalDof )
            // This expression is true in the case of an internal node
            //OLD (ERROR): if ( bitMaskFunctor ( **(dofMaskStencil[i]) == glob_i ) )
            // A node is set only if it is an internal node:
            if ( ** ( dofMaskStencil[i] ) == glob_i ) {
              for ( int j = 0; j < numLocalDofs; ++j ) {
                int glob_j = globalDofs[ j ];

                // A node is included in the matrix product only if the bitMaskFunctor applied to the bitMask is true
                // i.e. we apply the matrix only to either boundary or inner nodes.
                if ( bitMaskFunctor ( ** ( dofMaskStencil[j] ) == glob_j ) ) {
                  Dest[ glob_i ] += localMatrix[i][j] * Arg[ glob_j ] ;
                }


              }
            }
          }
        }

        for ( int j = 0; j < ConfiguratorType::maxNumLocalDofs; ++j ) {
          delete dofMaskStencil[j];
        }
        delete[] dofMaskStencil;
        break;

      case QUOC_ADAPTIVEGRID_INDEX_MODE:
        throw aol::Exception ( "FELinOpInterface::multiplyOnTheFly<BitMaskFunctorType>: Adaptive grid mode not implemented!", __FILE__, __LINE__ );
        break;

      default:
        throw aol::Exception ( "FELinOpInterface: Unknown global index mode!", __FILE__, __LINE__ );
    }


  }

  void assembleMatrix( ) const {
    if ( _mat ) delete _mat;
    _mat = this->getConfigurator().createNewMatrix( );
    assembleAddMatrix ( *_mat );
  }

public:
  /** (this assembled matrix * Factor) is added to Mat  */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    // this iterator traverses the elements of the grid (defined in the configurator)
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

    LocalAssemblyHelper<FELinOpInterface<RealType,ConfiguratorType,Imp,IndexMode>,MatrixType,IndexMode>::doLocalAssembly( this, localMatrix, end_it, Mat, Factor );
  }


  /** assemble Matrix with respect to homogeneous Dirichlet boundary condition */
  /* WARNING: element mask is only implemented for QUOC_GRID_INDEX_MODE */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const aol::BitVector * DirichletMask, bool setDirichletNodes = false,
                           qc::ElementMask<typename ConfiguratorType::InitType, ConfiguratorType, ConfiguratorType::Dim>* elMask = NULL,
                           const aol::BitVector * DirichletMaskCol = NULL ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    QUOC_ASSERT ( DirichletMask != NULL );
    if ( DirichletMaskCol == NULL ) DirichletMaskCol = DirichletMask;
    QUOC_ASSERT ( DirichletMask->size() == DirichletMaskCol->size() );

    // this iterator traverses the elements of the grid (defined in the configurator)
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

    LocalAssemblyHelper<FELinOpInterface<RealType,ConfiguratorType,Imp,IndexMode>,MatrixType,IndexMode>::doLocalAssemblyDirichlet
      ( this, localMatrix, end_it, Mat, DirichletMask, setDirichletNodes, elMask, DirichletMaskCol );
  }

    /*! Assembling routine with support for homogeneous Dirichlet boundary condition
    * and additional support nodes that should not belong to the simulation domain.
    * This assembling routine is only implemented for QUOC_GRID_INDEX_MODE and will raise an exception otherwise.
    * \author meier
    */
    template <typename MatrixType>
    void assembleAddMatrix ( MatrixType &Mat, const aol::BitVector * DirichletMask, bool setDirichletNodes, const aol::BitVector* nonDomainMask) const {
      typedef typename ConfiguratorType::ElementIteratorType IteratorType;

      QUOC_ASSERT ( DirichletMask != NULL );
      QUOC_ASSERT ( nonDomainMask != NULL );

      // this iterator traverses the elements of the grid (defined in the configurator)
      const typename IteratorType::EndType end_it = this->getConfigurator().end();
      int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
      aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

      switch ( IndexMode ) {

      case QUOC_GRID_INDEX_MODE:
        for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {

            const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

            // skip elements that have at least one local DOF that does not belong to the domain
            bool skipElement = false;
            for (int i=0; i<numLocalDofs; i++) {
                skipElement |= (*nonDomainMask)[this->getConfigurator().localToGlobal(*it, i)];
            }
            if (skipElement) {
                continue;
            }

          // assemble the local matrix for the current element
          this->asImp().prepareLocalMatrix ( *it, localMatrix );

          // get the global indices of the local Dofs of the current element
          for ( int i = 0; i < numLocalDofs; ++i ) {
            globalDofs[ i ] = this->getConfigurator().localToGlobal ( *it, i );
          }

          // finally add the locally computed values to the matrix
          for ( int i = 0; i < numLocalDofs; ++i ) {
            int glob_i = globalDofs[ i ];
            // write the ith row only if the ith node is not a Dirichlet node
            if ( ! (*DirichletMask)[glob_i] ) {
              for ( int j = 0; j < numLocalDofs; ++j ) {
                int glob_j = globalDofs[ j ];
                // write the jth column only if the jth node is not a Dirichlet node
                if ( ! (*DirichletMask)[glob_j] )
                  Mat.add ( glob_i, glob_j, localMatrix [ i ][ j ] );
              }
            }
          }
        }

        // set ones on the diagonal for Dirichlet nodes
        if ( setDirichletNodes )
          for ( int i = 0; i < DirichletMask->size(); ++i )
            if ( (*DirichletMask)[i] )
              Mat.add( i, i, aol::NumberTrait<RealType>::one );

        break;

      default:
          throw aol::Exception ( "FELinOpInterface::assembleAddMAtrix(): not implemented for this grid index mode!", __FILE__, __LINE__ );
      }
  }


protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

  mutable typename ConfiguratorType::MatrixType *_mat;
  OperatorType _opType;

  template <typename FEOpType, typename MatrixType, GridGlobalIndexMode indexMode> friend struct LocalAssemblyHelper;

private:
  FELinOpInterface ( const FELinOpInterface < RealType, ConfiguratorType, Imp, IndexMode >& ); // do not implement
  FELinOpInterface < RealType, ConfiguratorType, Imp, IndexMode >& operator= ( const FELinOpInterface < RealType, ConfiguratorType, Imp, IndexMode >& );
};


//! General Interface for efficient Finite Element operators for all types of meshes and vector-valued basis functions.
/*!
 * \author Benedikt Wirth
 */
template <typename ConfiguratorType, typename Imp, int NumVecCompsArg = ConfiguratorType::Dim, int NumVecCompsDest = NumVecCompsArg>
class FELinVectorOpInterface :
  public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

  mutable aol::BlockMatrix<typename ConfiguratorType::MatrixType> *_mat;
  const OperatorType _opType;

public:
  explicit FELinVectorOpInterface ( const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY ) :
    FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Grid ), _mat ( NULL ), _opType ( OpType ) {}

  explicit FELinVectorOpInterface ( const ConfiguratorType &Config, OperatorType OpType = ONTHEFLY ) :
    FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Config ), _mat ( NULL ), _opType ( OpType ) {}

  virtual ~FELinVectorOpInterface( ) {
    delete _mat;
  }

  //! clears the assembled matrix
  void reset( ) {
    if ( _mat )
      delete _mat;
    _mat = NULL;
  }

  void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    switch ( _opType ) {
    case ONTHEFLY:
      multiplyOnTheFly ( Arg, Dest );
      break;
    case ASSEMBLED:
#ifdef _OPENMP
#pragma omp critical (aol_FEOpInterface_callAssembleMatrix)
#endif
      if ( !_mat )
        assembleMatrix( );
      _mat->applyAdd ( Arg, Dest );
      break;
    default:
      throw aol::UnimplementedCodeException ( "FELinVectorOpInterface::applyAdd: unsupported opType", __FILE__, __LINE__ );
    }
  }

  /** (this assembled matrix * Factor) is added to Mat  */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    // this iterator traverses the elements of the grid (defined in the configurator)
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix[NumVecCompsArg][NumVecCompsDest];

    // traverse the elements of the grid
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = this->getConfigurator().localToGlobal ( *it, i );

      // finally add the locally computed values to the matrix
      for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
        for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
          for ( int i = 0; i < numLocalDofs; ++i ) {
            int glob_i = globalDofs[ i ];
            for ( int j = 0; j < numLocalDofs; ++j ) {
              int glob_j = globalDofs[ j ];
              Mat.getReference( destComp, argComp ).add( glob_i, glob_j, Factor * localMatrix[argComp][destComp][i][j] );
            }
          }
    }
  }

  /** assemble Matrix with respect to homogeneous Dirichlet boundary condition */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const aol::BitVector * DirichletMask, bool setDirichletNodes = false ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    QUOC_ASSERT ( DirichletMask != NULL );

    // this iterator traverses the elements of the grid (defined in the configurator)
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix[NumVecCompsArg][NumVecCompsDest];

    // traverse the elements of the grid
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = this->getConfigurator().localToGlobal ( *it, i );

      // finally add the locally computed values to the matrix
      for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
        for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
          for ( int i = 0; i < numLocalDofs; ++i ) {
            int glob_i = globalDofs[ i ];
            // write the ith row only if the ith node is not a Dirichlet node
            if ( ! (*DirichletMask)[glob_i] )
              for ( int j = 0; j < numLocalDofs; ++j ) {
                int glob_j = globalDofs[ j ];
                // write the jth column only if the jth node is not a Dirichlet node
                if ( ! (*DirichletMask)[glob_j] )
                  Mat.getReference( destComp, argComp ).add( glob_i, glob_j, localMatrix[argComp][destComp][i][j] );
              }
          }
    }

    // set ones on the diagonal for Dirichlet nodes
    if ( setDirichletNodes )
      for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
	for ( int i = 0; i < DirichletMask->size(); ++i )
	  if ( (*DirichletMask)[i] )
	    Mat.getReference( argComp, argComp ).add( i, i, aol::NumberTrait<RealType>::one );

  }

protected:
  void multiplyOnTheFly ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    // this iterator traverses the elements of the grid (defined in the configurator)
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix[NumVecCompsArg][NumVecCompsDest];

    // traverse the elements of the grid
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( *it, localMatrix );

      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = this->getConfigurator().localToGlobal ( *it, i );

      // add quadrature summands to correct positions in result vector
      for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
        for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
          for ( int i = 0; i < numLocalDofs; ++i ) {
            int glob_i = globalDofs[ i ];
            for ( int j = 0; j < numLocalDofs; ++j ) {
              int glob_j = globalDofs[ j ];
              Dest[destComp][ glob_i ] += localMatrix[argComp][destComp][i][j] * Arg[argComp][ glob_j ] ;
            }
          }
    }
  }

  void assembleMatrix( ) const {
    if ( _mat ) delete _mat;
    _mat = new aol::BlockMatrix<typename ConfiguratorType::MatrixType>( NumVecCompsDest, NumVecCompsArg, this->getNumGlobalDofs(), this->getNumGlobalDofs() );
    assembleAddMatrix ( *_mat );
  }

private:
  FELinVectorOpInterface ( const FELinVectorOpInterface < ConfiguratorType, Imp, NumVecCompsArg, NumVecCompsDest >& ); // do not implement
  FELinVectorOpInterface < ConfiguratorType, Imp, NumVecCompsArg, NumVecCompsDest >& operator= ( const FELinVectorOpInterface < ConfiguratorType, Imp, NumVecCompsArg, NumVecCompsDest >& );
};



enum LUMPED_MASS_OP_MODE {
  DO_NOT_INVERT,
  INVERT
};

/** General Interface for efficient Finite-Element operators for all types of meshes and basis functions.
 *  \author Droske
 */
template <typename RealType, typename ConfiguratorType, typename Imp>
class LumpedFELinOpInterface : public FEOpInterface<ConfiguratorType, aol::Vector<RealType> > {
protected:
  mutable aol::DiagonalMatrix<typename ConfiguratorType::RealType> *_mat;
  bool _invert;
public:

  LumpedFELinOpInterface ( const typename ConfiguratorType::InitType &Grid, LUMPED_MASS_OP_MODE Invert )
  : FEOpInterface<ConfiguratorType, aol::Vector<RealType> > ( Grid ), _mat ( NULL ) {
    if ( Invert == INVERT ) _invert = true;
    else _invert = false;
  }

  LumpedFELinOpInterface ( const ConfiguratorType &Config, LUMPED_MASS_OP_MODE Invert )
  : FEOpInterface<ConfiguratorType, aol::Vector<RealType> > ( Config ), _mat ( NULL ) {
    if ( Invert == INVERT ) _invert = true;
    else _invert = false;
  }

  virtual ~LumpedFELinOpInterface( ) {
    delete _mat;
  }

  //! clears the assembled matrix
  void reset( ) {
    if ( _mat ) {
      delete _mat;
    }
    _mat = NULL;
  }

  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    if ( !_mat ) {
      assembleMatrix( );
    }
    _mat->applyAdd ( Arg, Dest );
  }

  void applySingle ( Vector<RealType> &Arg ) const {
    if ( !_mat )
      assembleMatrix( );
    _mat->applySingle ( Arg );
  }

  aol::DiagonalMatrix<typename ConfiguratorType::RealType>& getMatrix( ) {
    if ( !_mat ) {
      assembleMatrix( );
    }
    return *_mat;
  }

  const aol::DiagonalMatrix<typename ConfiguratorType::RealType>& getMatrix( ) const {
    if ( !_mat ) {
      assembleMatrix( );
    }
    return *_mat;
  }

protected:
  void assembleMatrix( ) const {
    if ( !_mat ||  _mat->getNumRows() != this->getNumGlobalDofs()  ) {
      delete _mat;
      _mat = new aol::DiagonalMatrix<typename ConfiguratorType::RealType> ( this->getNumGlobalDofs() );
    }
    _mat->setZero();
    assembleAddMatrix ( *_mat );
  }

public:
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    aol::Vec<ConfiguratorType::maxNumLocalDofs,RealType> localDiagMatrix;

    if ( !_invert ) {
      for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
        this->asImp().prepareLocalMatrix ( *it, localDiagMatrix );

        const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

        for ( int i = 0; i < numLocalDofs; ++i ) {
          const int glob_i = this->getConfigurator().localToGlobal ( *it, i );
          Mat.add ( glob_i, glob_i, localDiagMatrix [ i ] );
        }
      }
    } else {
      // save the computed values in a vector and then invert the entries
      // that are not zero. It's NOT possible just to invert all matrix entries,
      // because some rows of the matrix might not exist (e.g. in a subGridSparseMatrix).
      aol::Vector< typename ConfiguratorType::RealType > tmp ( Mat.getNumRows () );
      tmp.setZero();

      for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
        this->asImp().prepareLocalMatrix ( *it, localDiagMatrix );

        const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

        for ( int i = 0; i < numLocalDofs; ++i ) {
          const int glob_i = this->getConfigurator().localToGlobal ( *it, i );
          tmp.add ( glob_i, localDiagMatrix [ i ] );
        }
      }

      // now invert the entries
      for ( int i = 0; i < tmp.size(); ++i )
        if ( tmp[i] != 0 )  Mat.add ( i, i, 1. / tmp[i] );
    }
  }


protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

private:
  LumpedFELinOpInterface ( const LumpedFELinOpInterface < RealType, ConfiguratorType, Imp >& ); // do not implement
  LumpedFELinOpInterface < RealType, ConfiguratorType, Imp >& operator= ( const LumpedFELinOpInterface < RealType, ConfiguratorType, Imp >& );
};


//! class to provide a virtual function v_getCoeff. This class has no template parameter for the derived class,
//! hence it can be used as class type in a list of operators which all provide the getCoeff-function.
//! (As for example used in the simultaneouslyAssembledLinCombOp).
//!  \author nemitz
template <typename ConfiguratorType>
class provideVirtualGetCoeffClass {
  public:
  virtual typename ConfiguratorType::RealType v_getCoeff ( const typename ConfiguratorType::ElementType &El,
                                                  int QuadPoint,
                                                  const typename ConfiguratorType::DomVecType& RefCoord ) const = 0;
  virtual ~provideVirtualGetCoeffClass() {}
};

//! class to provide a virtual function v_getCoeffMatrix. This class has no template parameter for the derived class,
//! hence it can be used as class type in a list of operators which all provide the getCoeffMatrix-function.
//! (As for example used in the simultaneouslyAssembledLinCombOp).
//!  \author nemitz
template <typename ConfiguratorType>
class provideVirtualGetCoeffMatrixClass {
  public:
    virtual void getCoeffMatrix (  const typename ConfiguratorType::ElementType &El, int QuadPoint,
                              const typename ConfiguratorType::DomVecType& DomRefCoord,
                              typename ConfiguratorType::MatType &Matrix ) const = 0;

    virtual ~provideVirtualGetCoeffMatrixClass() {}
};


/**
 * \brief Provides a Finite-Element operators of the form \f$ \int_\Omega w(x)\partial_s \phi_i \partial_t \phi_j dx \f$.
 * This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
 * \author Droske
 * \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinScalarWeightedMixedDiffInterface :
    public aol::FELinOpInterface<typename ConfiguratorType::RealType, ConfiguratorType, Imp, IndexMode> {
protected:
  const typename ConfiguratorType::InitType &_grid;
  const int _s, _t;
public:
  FELinScalarWeightedMixedDiffInterface ( const typename ConfiguratorType::InitType &Grid,
                                int s, int t, aol::OperatorType OpType = aol::ONTHEFLY ) :
      aol::FELinOpInterface<RealType, ConfiguratorType, Imp, IndexMode > ( Grid, OpType ),
      _grid ( Grid ), _s ( s ), _t ( t ) {}

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType  VecType;

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    LocalMatrix.setZero();

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      RealType coeff = this->asImp().getCoeff ( El, q, bfs.getRefCoord ( q ) );
      for ( int i = 0; i < numDofs; ++i ) {
        const VecType &basisveci  = bfs.evaluateGradient ( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          const VecType &basisvecj  = bfs.evaluateGradient ( j, q );
          LocalMatrix[i][j] += coeff * basisveci[_s] * basisvecj[_t] * bfs.getWeight ( q );
        }
      }
    }
    LocalMatrix *= this->getConfigurator().vol ( El );
  }

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }
};

/**
 * \brief \f$ \left(\int_\Omega \partial_s \phi_i \partial_t \phi_j dx\right)_{ij} \f$.
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class FEOpMixedDerivative : public FELinScalarWeightedMixedDiffInterface<ConfiguratorType, FEOpMixedDerivative<ConfiguratorType> > {
public:
  FEOpMixedDerivative ( const typename ConfiguratorType::InitType &Grid,
                        int s, int t, aol::OperatorType OpType = aol::ONTHEFLY )
      : FELinScalarWeightedMixedDiffInterface<ConfiguratorType, FEOpMixedDerivative<ConfiguratorType> > ( Grid, s, t, OpType ) {}

  typedef typename ConfiguratorType::RealType RealType;

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return 1.;
  }
};


/** \brief Operator for mixed derivatives weighted by a scalar.
 *  \author Schwen (MEVIS)
 *  \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class WeightedFEOpMixedDerivative : public FELinScalarWeightedMixedDiffInterface<ConfiguratorType, WeightedFEOpMixedDerivative<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrWeight;

  WeightedFEOpMixedDerivative ( const typename ConfiguratorType::InitType &Grid,
                                const aol::Vector<RealType> &Weight,
                                int s, int t, aol::OperatorType OpType = aol::ONTHEFLY )
    : FELinScalarWeightedMixedDiffInterface<ConfiguratorType, WeightedFEOpMixedDerivative<ConfiguratorType> > ( Grid, s, t, OpType ),
      _discrWeight ( Grid, Weight ) {
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El,
                             int QuadPoint,
                             const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return _discrWeight.evaluateAtQuadPoint(El, QuadPoint);
  }
};


//! \brief provides an easy interface to Finite Element operators of the form \f$ \mbox{div}(A(x)\nabla u)\f$, where \f$A\f$
//!        is a SYMMETRIC coefficient matrix.
/*!
 * This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
 * \author Droske
 * \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinMatrixWeightedStiffInterface :
    public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                          FELinMatrixWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode >,
    public provideVirtualGetCoeffMatrixClass<ConfiguratorType> {
public:
  explicit FELinMatrixWeightedStiffInterface ( const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FELinOpInterface<RealType, ConfiguratorType, FELinMatrixWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Grid, OpType ),
      _grid ( Grid ) {
  }

  FELinMatrixWeightedStiffInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FELinOpInterface<RealType, ConfiguratorType, FELinMatrixWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Config, OpType ),
      _grid ( Grid ) {
  }

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;


  //! Implementation of the virtual function v_getCoeff, which just calls the getCoeff-method
  void v_getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                     const DomVecType& DomRefCoord, MatType &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, DomRefCoord, Matrix );
  }

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const DomVecType& DomRefCoord, typename ConfiguratorType::MatType &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, DomRefCoord, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    MatType mat;
    VecType matgrad1;

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix ( El, q, bfs.getRefCoord ( q ), mat );
      for ( int i = 0; i < numDofs; ++i ) {
        mat.mult ( bfs.evaluateGradient ( i, q ), matgrad1 );
        for ( int j = i; j < numDofs; ++j ) {
          LocalMatrix[i][j] += ( matgrad1 * bfs.evaluateGradient ( j, q ) ) * bfs.getWeight ( q );
        }
      }
    }
    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = i; j < numDofs; ++j ) {
        LocalMatrix[i][j] *= this->getConfigurator().vol ( El );
      }
    }

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = i + 1; j < numDofs; ++j ) {
        LocalMatrix[j][i] = LocalMatrix[i][j];
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  const typename ConfiguratorType::InitType &_grid;
};


//! \brief provides an easy interface to Finite Element operators of the form \f$ \mbox{div}(A(x)\nabla u)\f$, where \f$A\f$
//!        is a ASYMMETRIC coefficient matrix.
/*!
 * This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
 * The corresponding matrix assembly yields \f$ \left(\int_\Omega \nabla\phi_i\cdot A(x)\nabla\phi_j dx\right)_{ij} \f$
 * for FE basis functions \f$ \phi_i,\phi_j \f$.
 *
 * \author Droske
 * \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinAsymMatrixWeightedStiffInterface :
      public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                            FELinAsymMatrixWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > {
protected:
public:
  explicit FELinAsymMatrixWeightedStiffInterface ( const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FELinOpInterface<RealType, ConfiguratorType, FELinAsymMatrixWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Grid, OpType ),
      _grid ( Grid ) {
  }

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType VecType;

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType& DomRefCoord, typename ConfiguratorType::MatType &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, DomRefCoord, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    MatType mat;
    VecType matgrad1;

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      getCoeffMatrix ( El, q, bfs.getRefCoord ( q ), mat );
      for ( int i = 0; i < numDofs; ++i ) {
        mat.mult ( bfs.evaluateGradient ( i, q ), matgrad1 );
        for ( int j = 0; j < numDofs; ++j ) {
          LocalMatrix[j][i] += ( matgrad1 * bfs.evaluateGradient ( j, q ) ) * bfs.getWeight ( q );
        }
      }
    }
    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] *= this->getConfigurator().vol ( El );
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  const typename ConfiguratorType::InitType &_grid;
};


//! \brief provides an easy interface to Finite Element operators of the form \f$ \mbox{div}(a(x)\nabla u) \f$, where \f$a\f$
//!        is a positive scalar function.
/*!
 * This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
 * \author Droske
 * \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinScalarWeightedStiffInterface :
    public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                          FELinScalarWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode >,
    public provideVirtualGetCoeffClass<ConfiguratorType> {
public:
  explicit FELinScalarWeightedStiffInterface ( const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY ) :
      FELinOpInterface<RealType, ConfiguratorType, FELinScalarWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Grid, OpType ),
      _grid ( Grid ) {
  }

  FELinScalarWeightedStiffInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY ) :
      FELinOpInterface<RealType, ConfiguratorType, FELinScalarWeightedStiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Config, OpType ),
      _grid ( Grid ) {
  }

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType  VecType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;


  //! Implementation of the virtual function v_getCoeff, which just calls the getCoeff-method
  RealType v_getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType& RefCoord ) const {
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType& RefCoord ) const {
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    // loop over the quadrature points
    for ( int q = 0; q < numQuadPoints; ++q ) {
      // get the coefficient from the derived class
      RealType coeff = getCoeff ( El, q, bfs.getRefCoord ( q ) );
      // now add the contribution of each pair of basisfunctions to the local matrix
      for ( int i = 0; i < numDofs; ++i ) {
        const VecType &basisveci  = bfs.evaluateGradient ( i, q );
        for ( int j = i; j < numDofs; ++j ) {
          LocalMatrix[i][j] += ( basisveci * bfs.evaluateGradient ( j, q ) ) * coeff * bfs.getWeight ( q );
        }
      }
    }

    // finally multiply the local matrix with the volume of the element (H^2 in most cases)
    // and make use of symmetry to compute the lower triangular entries.
    const RealType vol = this->getConfigurator().vol ( El );
    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = i; j < numDofs; ++j ) {
        LocalMatrix[j][i] = LocalMatrix[i][j] = LocalMatrix[i][j] * vol;
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  const typename ConfiguratorType::InitType &_grid;
};


//! \brief provides an easy interface to scaled mass matrices of the form \f$ \int_\Omega \gamma(x)\phi_i \phi_j\f$, where \f$\gamma\f$
//!        is a positive scalar function.
/*!
 * This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
 * \author Droske
 * \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinScalarWeightedMassInterface :
      public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                            FELinScalarWeightedMassInterface<ConfiguratorType, Imp, IndexMode>, IndexMode >,
      public provideVirtualGetCoeffClass<ConfiguratorType> {
public:
  explicit FELinScalarWeightedMassInterface (  const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
      FELinScalarWeightedMassInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Grid, OpType ),
      _grid ( Grid ) {
  }
  FELinScalarWeightedMassInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
      FELinScalarWeightedMassInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Config, OpType ),
      _grid ( Grid ) {
  }
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::QuadType QuadType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  //! Implementation of the virtual function v_getCoeff, which just calls the getCoeff-method
  RealType v_getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType& RefCoord ) const {
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType& RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      RealType coeff = this->asImp().getCoeff ( El, q, bfs.getRefCoord ( q ) );
      for ( int i = 0; i < numDofs; ++i ) {
        RealType basisi = bfs.evaluate ( i, q );
        for ( int j = i; j < numDofs; ++j ) {
          LocalMatrix[i][j] +=  basisi * bfs.evaluate ( j, q ) * coeff * bfs.getWeight ( q ) ;
        }
      }
    }

    RealType vol = this->getConfigurator().vol ( El );
    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = i; j < numDofs; ++j ) {
        LocalMatrix[j][i] = LocalMatrix[i][j] = LocalMatrix[i][j] * vol;
      }
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  const typename ConfiguratorType::InitType &_grid;
};

//! \brief Class to assemble matrices of the form \f$ \left(\int_\Omega A(x)\phi_i\cdot\phi_j\mathrm{d}x\right)_{ij} \f$,
//!        where \f$ A \f$ is a matrix-valued weight function and the \f$ \phi_k \f$ are the finite element basis functions.
/*!
 * \author Benedikt Wirth
 * \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, int NumVecCompsArg = ConfiguratorType::Dim, int NumVecCompsDest = NumVecCompsArg>
class FELinMatrixWeightedVectorMassInterface :
  public FELinVectorOpInterface<ConfiguratorType,FELinMatrixWeightedVectorMassInterface<ConfiguratorType, Imp, NumVecCompsArg, NumVecCompsDest> > {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::QuadType   QuadType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

public:
  explicit FELinMatrixWeightedVectorMassInterface (  const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY ) :
    FELinVectorOpInterface<ConfiguratorType,FELinMatrixWeightedVectorMassInterface<ConfiguratorType, Imp, NumVecCompsArg, NumVecCompsDest> > ( Grid, OpType ) {}

  explicit FELinMatrixWeightedVectorMassInterface ( const ConfiguratorType & Config, OperatorType OpType = ONTHEFLY ) :
    FELinVectorOpInterface<ConfiguratorType,FELinMatrixWeightedVectorMassInterface<ConfiguratorType, Imp, NumVecCompsArg, NumVecCompsDest> > ( Config, OpType ) {}

  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const DomVecType &/*RefCoord*/, MatType &/*Matrix*/ ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El,
                            Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> LocalMatrix[NumVecCompsArg][NumVecCompsDest] ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
      for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
        for ( int i = 0; i < numDofs; ++i )
          for ( int j = 0; j < numDofs; ++j )
            LocalMatrix[argComp][destComp][i][j] = 0.;

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    MatType coeffMatrix;
    for ( int q = 0; q < numQuadPoints; ++q ) {
      this->asImp().getCoeffMatrix( El, q, bfs.getRefCoord ( q ), coeffMatrix );
      for ( int i = 0; i < numDofs; ++i ) {
        RealType basisi = bfs.evaluate ( i, q );
        for ( int j = i; j < numDofs; ++j )
          for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
            for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
              LocalMatrix[argComp][destComp][i][j] += basisi * bfs.evaluate ( j, q ) * coeffMatrix[destComp][argComp] * bfs.getWeight ( q );
      }
    }

    RealType vol = this->getConfigurator().vol ( El );
    for ( int argComp = 0; argComp < NumVecCompsArg; ++argComp )
      for ( int destComp = 0; destComp < NumVecCompsDest; ++destComp )
        for ( int i = 0; i < numDofs; ++i )
          for ( int j = i; j < numDofs; ++j )
            LocalMatrix[argComp][destComp][j][i] = LocalMatrix[argComp][destComp][i][j] = LocalMatrix[argComp][destComp][i][j] * vol;
  }
};


//! provides an easy interface to lumped mass matrices. Attention: works only for nodal basis Finite Elements!
/*!
 * This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
 * \author Droske
 */
template <typename ConfiguratorType, typename Imp>
class LumpedMassOpInterface :
      public LumpedFELinOpInterface<typename ConfiguratorType::RealType, ConfiguratorType, LumpedMassOpInterface<ConfiguratorType, Imp> > {
public:

  LumpedMassOpInterface (  const typename ConfiguratorType::InitType &Initializer, LUMPED_MASS_OP_MODE Invert )
    : LumpedFELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
                            LumpedMassOpInterface<ConfiguratorType, Imp> > ( Initializer, Invert ) {}

  LumpedMassOpInterface (  const ConfiguratorType &Config, LUMPED_MASS_OP_MODE Invert )
    : LumpedFELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
                            LumpedMassOpInterface<ConfiguratorType, Imp> > ( Config, Invert ) {}

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::QuadType QuadType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType& RefCoord ) const {
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Vec<ConfiguratorType::maxNumLocalDofs, RealType> &LocalDiagMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      LocalDiagMatrix[i] = 0.;
    }

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      RealType coeff = getCoeff ( El, q, bfs.getRefCoord ( q ) );
      for ( int i = 0; i < numDofs; ++i ) {
        LocalDiagMatrix[i] +=  coeff * bfs.evaluate ( i, q ) * bfs.getWeight ( q ) ;
      }
    }

    for ( int i = 0; i < numDofs; ++i ) {
      LocalDiagMatrix[i] *= this->getConfigurator().vol ( El );
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};


//! \brief General Interface for nonlinear FE-operators depending on x, which are locally assembled.
//!        only getNonlinearity has to be provided.
//!
//! not for operators depending on derivatives of
//! the function, use FENonlinDiffOpInterface instead.
/*!
 * \author Droske
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename Imp>
class FENonlinOpInterface : public FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinOpInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  FENonlinOpInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~FENonlinOpInterface( ) {}

  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc ( this->getConfigurator(), Arg );

    typedef RealType NLTYPE;

    NLTYPE *nl_cache = new NLTYPE[ this->getConfigurator().maxNumQuadPoints() ];

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( discrFunc, *it, q, bfs.getRefCoord ( q ), nl_cache[q] );
      }

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        RealType a = 0.;

        for ( int q = 0; q < numQuadPoints; ++q ) {
          a += ( nl_cache[q] * bfs.evaluate ( dof, q ) ) * bfs.getWeight ( q );
        }

        a *= this->getConfigurator().vol ( *it );

        Dest[ this->getConfigurator().localToGlobal ( *it, dof ) ] += a;
      }
    }
    delete[] nl_cache;
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         typename ConfiguratorType::RealType &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFunc, El, QuadPoint, RefCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//!Operator for evaluating a function on each element.
/**Computes \f$ \left(\frac{\int_\Omega \vec{f}(\vec\phi(x),\nabla\vec\phi(x),x)\varphi_i(x)dx}{\int_\Omega\varphi_i(x)dx}\right)_i \f$,
 * where \f$\varphi_i\f$ is the \f$i\f$th finite element basis function, \f$\vec f\f$ is specified by "getNonlinearity",
 * and \f$\vec\phi\f$ is passed to "apply" as argument. The size of \f$\vec f\f$ and \f$\vec\phi\f$ may be chosen freely.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int DimArg, int DimDest, typename Imp>
class FENonlinEvaluationOpInterface : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinEvaluationOpInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  FENonlinEvaluationOpInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~FENonlinEvaluationOpInterface( ) {}

  void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    MultiVector<RealType> dest( Dest, aol::STRUCT_COPY );
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,DimArg> discrFunc ( _initializer, Arg );

    aol::Vector<RealType> locWeight( Arg[0], aol::STRUCT_COPY );

    typedef aol::Vec<DimDest,RealType> NLTYPE;

    NLTYPE *nl_cache = new NLTYPE[ this->getConfigurator().maxNumQuadPoints() ];

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( discrFunc, *it, q, bfs.getRefCoord ( q ), nl_cache[q] );
      }

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        NLTYPE a;
        RealType b = 0.;

        for ( int q = 0; q < numQuadPoints; ++q ) {
          RealType fac = bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
          a += nl_cache[q] * fac;
          b += fac;
        }

        a *= this->getConfigurator().vol ( *it );
        b *= this->getConfigurator().vol ( *it );

        for ( int j = 0; j < DimDest; ++j )
          dest[j][ this->getConfigurator().localToGlobal ( *it, dof ) ] += a[j];
        locWeight[ this->getConfigurator().localToGlobal ( *it, dof ) ] += b;
      }
    }
    delete[] nl_cache;

    for ( int i = 0; i < locWeight.size(); ++i )
      for ( int j = 0; j < DimDest; ++j )
        dest[j][i] /= locWeight[i];
    Dest += dest;
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,DimArg> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         typename aol::Vec<DimDest,RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFunc, El, QuadPoint, RefCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! \brief General Interface for nonlinear FE-operators depending on x, which are locally assembled.
//!        Computes \f$ \left(\int_\Omega \vec{f}\left(\vec\phi(x),\nabla\vec\phi(x),x\right)\cdot\vec\varphi_i(x) dx\right)_i \f$,
//!        only getNonlinearity has to be provided.
//!
//! not for operators depending on derivatives of
//! the function, use FENonlinDiffOpInterface instead.
/*!
 * \author Droske
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinVectorOpInterface : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinVectorOpInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  FENonlinVectorOpInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~FENonlinVectorOpInterface( ) {}

  void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {

    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " "
      << Dest.numComponents() << " "
      << NumCompArg << " "
      << NumCompDest << " " << endl;

      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumCompArg; c++ ) {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );
    }

    typedef typename aol::Vec<NumCompDest, RealType> NLTYPE;

    NLTYPE *nl_cache = new NLTYPE[ this->getConfigurator().maxNumQuadPoints() ];

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const RealType vol = this->getConfigurator().vol ( *it );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( discrFuncsArg, *it, q, bfs.getRefCoord ( q ), nl_cache[q] );
      }

      aol::Vec<NumCompDest, RealType> a;
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        a.setZero();

        for ( int q = 0; q < numQuadPoints; ++q ) {
          a.addMultiple ( nl_cache[q], bfs.evaluate ( dof, q ) * bfs.getWeight ( q ) );
        }

        for ( int d = 0; d < NumCompDest; ++d ) {
          Dest[d][ this->getConfigurator().localToGlobal ( *it, dof ) ] += vol * a[d];
        }
      }
    }
    delete[] nl_cache;
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         aol::Vec<NumCompDest, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, QuadPoint, RefCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! General Interface
/*!
 * \author Droske
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename Imp>
class FENonlinDiffOpInterface : public FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinDiffOpInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  FENonlinDiffOpInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~FENonlinDiffOpInterface( ) {}

  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc ( this->getConfigurator(), Arg );

    typedef typename ConfiguratorType::VecType NLTYPE;

    NLTYPE *nl_cache = new NLTYPE[ this->getConfigurator().maxNumQuadPoints() ];

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( discrFunc, *it, q, bfs.getRefCoord ( q ), nl_cache[q] );
      }

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        RealType a = 0.;
        for ( int q = 0; q < numQuadPoints; ++q ) {
          a += ( nl_cache[q] * bfs.evaluateGradient ( dof, q ) ) * bfs.getWeight ( q );
        }
        a *= this->getConfigurator().vol ( *it ) ;
        Dest[ this->getConfigurator().localToGlobal ( *it, dof ) ] += a;
      }
    }
    delete[] nl_cache;
  }


  template<typename BitMaskFunctorType>
  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {

    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    BitMaskFunctorType bitMaskFunctor;
    const aol::DofMask& dofMask = this->getConfigurator().getNodeBoundaryMask();
    typedef aol::DofMask::iterator *DofMaskIteratorPtr;
    DofMaskIteratorPtr *dofMaskStencil;
    dofMaskStencil = new DofMaskIteratorPtr[ ConfiguratorType::maxNumLocalDofs ];
    for (int j=0; j<ConfiguratorType::maxNumLocalDofs; ++j) {
      dofMaskStencil[j] = new aol::DofMask::iterator(dofMask.begin());
    }

    // to implement the extension of the inner nodes or the masking of the inner nodes
    // respectively, we have to copy the argument and mask either the inner or the
    // boundary nodes with 0. depending on the bitMaskFunctor.
    aol::DofMask::iterator maskingIterator( dofMask.begin() );
    Vector<RealType> maskedArg( Arg.size() );
    for ( int i=0; i<Arg.size(); ++i ) {
      if ( bitMaskFunctor ( *(maskingIterator) == i ) )
      {
        maskedArg[i] = Arg[i];
      }

      if (*(maskingIterator) == i)
        ++maskingIterator;
    }

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc ( this->getConfigurator(), Arg );


    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {

        // Get the global indices of all nodes of the current element
        globalDofs[ dof ] = this->getConfigurator().localToGlobal ( *it, dof );
        // increment will jump over boundary indices
        dofMaskStencil[dof]->incrementUntilIndex( globalDofs[dof] );

        // Include contribution from node if the bitMaskFunctor evaluates to true at that node
        if ( **(dofMaskStencil[dof]) == globalDofs[dof] )
        {

          RealType a = 0.;
          for ( int q = 0; q < numQuadPoints; ++q ) {
            typename ConfiguratorType::VecType nl;

            this->asImp().getNonlinearity ( discrFunc, *it, q, bfs.getRefCoord ( q ), nl );

            a += ( nl * bfs.evaluateGradient ( dof, q ) ) * bfs.getWeight ( q );
          }
          a *= this->getConfigurator().vol ( *it ) ;
          Dest[ globalDofs[dof] ] += a;

        }

      }
    }


    for (int j=0; j<ConfiguratorType::maxNumLocalDofs; ++j) {
      delete dofMaskStencil[j];
    }
    delete[] dofMaskStencil;

  }



  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         typename ConfiguratorType::VecType &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFunc, El, QuadPoint, RefCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

/*!
 * \brief General Interface for nonlinear first order FE-operators depending on x, which are locally assembled.
 * only getNonlinearity has to be provided. not for operators depending on derivatives of
 * the function, use FENonlinDiffOpInterface instead.
 * The discrete operator is of the form \f[ [F(\Phi)]_i = \int_\Omega tr\left( f( x, \nabla \phi )^T \cdot \nabla \psi_i \right) dx \f]
 * where \f$ \psi_i \f$ corresponds to the i-th basis function, \f$f:R^d\times R^{d,d} \to R^{d,d} \f$, which has to be implemented in getNonlinearity.
 * \author Droske
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp, typename CompType = aol::auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > >
class FENonlinVectorDiffOpInterface : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinVectorDiffOpInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  FENonlinVectorDiffOpInterface ( const ConfiguratorType &Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~FENonlinVectorDiffOpInterface( ) {}

  void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {

    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " "
      << Dest.numComponents() << " "
      << NumCompArg << " "
      << NumCompDest << " " << endl;

      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumCompArg; c++ ) {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );
    }

    applyAdd( discrFuncsArg, Dest);
  }

  void applyAdd ( const CompType &Arg, MultiVector<RealType> &Dest ) const {

    typedef aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> NLTYPE;
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    NLTYPE *nl_cache = new NLTYPE[ this->getConfigurator().maxNumQuadPoints() ];

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( Arg, *it, q, bfs.getRefCoord ( q ), nl_cache[q] );
      }

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        aol::Vec<NumCompDest, typename ConfiguratorType::RealType> a;
        a.setZero();

        for ( int q = 0; q < numQuadPoints; ++q ) {
          typename aol::Mat<NumCompDest, ConfiguratorType::Dim, typename ConfiguratorType::RealType> nl;
          typename ConfiguratorType::VecType grad;
          typename aol::Vec<NumCompDest, typename ConfiguratorType::RealType> tmp;

          grad = bfs.evaluateGradient ( dof, q );
          nl = nl_cache[q];

          nl *= bfs.getWeight ( q );

          nl.mult ( grad, tmp );
          a += tmp;
        }

        a *= this->getConfigurator().vol ( *it );

        for ( int d = 0; d < NumCompDest; d++ ) {
          Dest[d][ this->getConfigurator().localToGlobal ( *it, dof ) ] += a[d];
        }
      }
    }
    delete[] nl_cache;
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const CompType &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         aol::Mat<NumCompDest, ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, QuadPoint, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


/*!
 * \brief General Interface for nonlinear first order FE-operators depending on x, which are locally assembled.
 * only getNonlinearity has to be provided. not for operators depending on derivatives of
 * the function, use FENonlinDiffOpInterface instead.
 * The discrete operator is of the form \f[ [F(\Phi)]_i = \int_\Omega tr\left( f( x, \nabla \phi )^T \cdot \nabla \psi_i \right) dx \f]
 * where \f$ \psi_i \f$ corresponds to the i-th basis function, \f$ :R^d\times R^{d,d} \to R^{d,d} \f$, which has to be implemented in getNonlinearity.
 * \author Droske
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp>
class NonlinearCombinedDiffVectorOpInterface : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit NonlinearCombinedDiffVectorOpInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  NonlinearCombinedDiffVectorOpInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~NonlinearCombinedDiffVectorOpInterface( ) {}

  void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {

    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " "
      << Dest.numComponents() << " "
      << NumCompArg << " "
      << NumCompDest << " " << endl;

      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumCompArg; c++ ) {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );
    }


    typedef typename ConfiguratorType::MatType NL_MAT_TYPE;
    typedef typename ConfiguratorType::VecType NL_VEC_TYPE;

    NL_MAT_TYPE nl_mat_cache = new NL_MAT_TYPE[ this->getConfigurator().maxNumQuadPoints() ];
    NL_VEC_TYPE nl_vec_cache = new NL_VEC_TYPE[ this->getConfigurator().maxNumQuadPoints() ];


    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( discrFuncsArg, *it, q, bfs.getRefCoord ( q ), nl_mat_cache[q], nl_vec_cache[q] );
      }

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        aol::Vec<NumCompDest, typename ConfiguratorType::RealType> a;
        a.setZero();

        for ( int q = 0; q < numQuadPoints; ++q ) {
          NL_MAT_TYPE nl_mat;
          NL_VEC_TYPE grad, nl_vec, tmp;

          grad = bfs.evaluateGradient ( dof, q );

          nl_mat = nl_mat_cache[q];
          nl_vec = nl_vec_cache[q];

          nl_mat *= bfs.getWeight ( q );

          nl_mat.mult ( grad, tmp );
          a += tmp;

          nl_vec *= bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
          a += nl_vec;
        }
        a *= this->getConfigurator().vol ( *it );

        for ( int d = 0; d < NumCompDest; d++ ) {
          Dest[d][ this->getConfigurator().localToGlobal ( *it, dof ) ] += a[d];
        }
      }
    }
    delete [] nl_mat_cache;
    delete [] nl_vec_cache;
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         aol::Mat<NumCompDest, NumCompDest, typename ConfiguratorType::RealType> &NL_mat,
                         aol::Vec<NumCompDest, typename ConfiguratorType::RealType> &NL_vec ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, QuadPoint, NL_mat, NL_vec );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


/**
 * \brief Provides an easy interface to operators of the form \f$ \int_\Omega f(x,p(x))[\phi]\varphi_i(p(x))\f$,
 *        where \f$\phi=\f$Arg, \f$f\f$ supplied by getNonlinearity and \f$p\f$ by getTransformedPosition.
 *
 * \author Berkels
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename Imp>
class FENonlinOpInterfaceWithMovedTestFunctions : public FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
public:
  explicit FENonlinOpInterfaceWithMovedTestFunctions ( const typename ConfiguratorType::InitType &Grid )
    : FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > ( Grid ),
      _grid ( Grid ) {
  }

  FENonlinOpInterfaceWithMovedTestFunctions ( const ConfiguratorType &Config, const typename ConfiguratorType::InitType &Grid )
    : FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > ( Config ),
      _grid ( Grid ) {
  }

  virtual ~FENonlinOpInterfaceWithMovedTestFunctions( ) {}

  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    const aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc ( this->getConfigurator(), Arg );

    typename ConfiguratorType::RealType valueAtQuadPoint;
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType transformedCoord;

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord") and f(\phi,\Phi(x),x) ("valueAtQuadPoint")
        if ( this->asImp().getTransformedPosition ( discrFunc, *it, q, bfs.getRefCoord ( q ), transformedEl, transformedCoord ) ){
          this->asImp().getNonlinearity ( discrFunc, *it, q, bfs.getRefCoord ( q ), transformedEl, transformedCoord, valueAtQuadPoint );
          // scale f(\phi,\Phi(x),x) with the volume of the current element and the correct quadrature weight
          valueAtQuadPoint *= this->getConfigurator().vol( *it ) * bfs.getWeight ( q );

          // for each \varphi_i with support of \varphi_i\circ\Phi in the current element add its contribution
          const typename ConfiguratorType::BaseFuncSetType &tbfs = this->getConfigurator().getBaseFunctionSet ( transformedEl );
          const int numTransformedLocalDofs = this->getConfigurator().getNumLocalDofs ( transformedEl );
          for ( int dof = 0; dof < numTransformedLocalDofs; dof++ )
            Dest[ this->getConfigurator().localToGlobal ( transformedEl, dof ) ] += valueAtQuadPoint * tbfs.evaluate ( dof, transformedCoord );
        }
      }
    }
  }

  //! Represents \f$f(\phi,\Phi(x),x)\f$ and has to be implemented in the derived class.
  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         const typename ConfiguratorType::ElementType &TransformedEl,
                         const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         typename ConfiguratorType::RealType &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFunc, El, QuadPoint, RefCoord, TransformedEl, TransformedLocalCoord, NL );
  }

  //! interface function, has to be provided in derived classes.
  bool getTransformedPosition ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                                const typename ConfiguratorType::ElementType &El,
                                int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                                typename ConfiguratorType::ElementType &TransformedEl,
                                typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getTransformedPosition ( DiscFunc, El, QuadPoint, RefCoord, TransformedEl, TransformedLocalCoord );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};


/**
 * \brief Standard stiffness matrix \f$ \left(\int_\Omega\nabla\varphi_i(x)\cdot\nabla\varphi_j(x)dx\right)_{ij}\f$.
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class StiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, StiffOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  explicit StiffOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = ONTHEFLY )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, StiffOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType )  {}

  StiffOp ( const ConfiguratorType & Configurator, const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = ONTHEFLY )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, StiffOp<ConfiguratorType, IndexMode>, IndexMode > ( Configurator, Initializer, OpType )  {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return 1.;
  }
};


/**
 * \brief provides an easy interface to Finite Element operators of the form \f$ \mbox{div}(a(x)\nabla u) \f$, where \f$a\f$
 *        is a Finite Element function given by the aol::Vector Weight.
 *
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class WeightedStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, WeightedStiffOp<ConfiguratorType, IndexMode>, IndexMode > {
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrWeight;
public:
  WeightedStiffOp ( const typename ConfiguratorType::InitType &Initializer,
                    const aol::Vector<RealType> &Weight,
                    aol::OperatorType OpType = ONTHEFLY )
    : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, WeightedStiffOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ),
      _discrWeight( Initializer, Weight )
  {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El,
                             int QuadPoint,
                             const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return _discrWeight.evaluateAtQuadPoint(El, QuadPoint);
  }
};


/**
 * \brief Standard mass matrix \f$ \left(\int_\Omega\varphi_i(x)\varphi_j(x)dx\right)_{ij} \f$.
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class MassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, MassOp<ConfiguratorType, IndexMode>, IndexMode > {
protected:
public:
  typedef typename ConfiguratorType::RealType RealType;

  explicit MassOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = ONTHEFLY )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, MassOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ) {}


  MassOp ( const ConfiguratorType &Config, const typename ConfiguratorType::InitType &Initializer, OperatorType OpType = ONTHEFLY )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, MassOp<ConfiguratorType, IndexMode>, IndexMode > ( Config, Initializer, OpType ) {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType& , int, const typename ConfiguratorType::DomVecType& ) const {
    return 1.;
  }
};


/**
 * \brief provides an easy interface to scaled mass matrices of the form \f$ \int_\Omega \gamma(x)\phi_i \phi_j\f$, where \f$\gamma\f$
 *        is a Finite Element function given by the aol::Vector Weight.
 *
 * \attention FELinScalarWeightedMassInterface says the weight has to be positive. I'm not sure that this is really needed.
 *
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class WeightedMassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, WeightedMassOp<ConfiguratorType, IndexMode>, IndexMode > {
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrWeight;
public:
  WeightedMassOp ( const typename ConfiguratorType::InitType &Initializer,
                   const aol::Vector<RealType> &Weight,
                   aol::OperatorType OpType = ONTHEFLY )
    : aol::FELinScalarWeightedMassInterface<ConfiguratorType, WeightedMassOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ),
      _discrWeight( Initializer, Weight )
  {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El,
                             int QuadPoint,
                             const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return _discrWeight.evaluateAtQuadPoint(El, QuadPoint);
  }
};

/**
 * \brief This class represents the weighted mass matrix \f$ \left(\int_\Omega w(x)^2\varphi_i(x)\varphi_j(x)dx\right)_{ij} \f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ w \f$ is passed to the constructor as aol::Vector.
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class SquaredWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredWeightMassOp<ConfiguratorType, IndexMode>, IndexMode > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the function $w$, which is to be squared
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;

public:
  SquaredWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                         const aol::Vector<RealType> &W,
                         aol::OperatorType OpType = aol::ONTHEFLY ) :
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredWeightMassOp<ConfiguratorType, IndexMode>, IndexMode >( Grid, OpType ),
    _w( Grid, W ) {}

  /**
   * Returns \f$ w(x)^2 \f$ at the point $x$ specified by element, quadrature point, and local coordinates.
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return aol::Sqr( _w.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

/**
 * \brief This class represents a scaled mass operator, assembling the matrix \f$ \left(\int_\Omega \|\nabla w\|^2\psi_i\psi_j\right)_{ij} \f$,
 * where the real-valued function \f$ w \f$ is passed to the constructor. \f$\psi_{i,j}\f$ represent the \f$i\f$-th and \f$j\f$-th FE basis function.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class SquaredDiffWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredDiffWeightMassOp<ConfiguratorType, IndexMode>, IndexMode> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the weight to be differentiated and then squared
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;

public:
  SquaredDiffWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                           const aol::Vector<RealType> &W,
                           aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredDiffWeightMassOp<ConfiguratorType, IndexMode>, IndexMode>( Grid, OpType ),
    // load the weight
    _w( Grid, W ) {}

  /**
   * \brief Returns \f$\|\nabla w\|^2\f$ evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    _w.evaluateGradientAtQuadPoint( El, QuadPoint, grad );
    return grad.normSqr();
  }
};

/**
 * \brief This class represents the weighted stiffness matrix \f$ \left(\int_\Omega w(x)^2\nabla\varphi_i(x)\cdot\nabla\varphi_j(x)dx\right)_{ij}\f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ w \f$ is passed to the constructor as aol::Vector.
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class SquaredWeightStiffOp :
  public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, SquaredWeightStiffOp<ConfiguratorType, IndexMode>, IndexMode > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the function $w$, which is to be squared
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;

public:
  SquaredWeightStiffOp( const typename ConfiguratorType::InitType &Grid,
                        const aol::Vector<RealType> &W,
                        aol::OperatorType OpType = aol::ONTHEFLY ) :
    aol::FELinScalarWeightedStiffInterface<ConfiguratorType, SquaredWeightStiffOp<ConfiguratorType, IndexMode>, IndexMode >( Grid, OpType ),
    _w( Grid, W ) {}

  /**
   * Returns \f$ w(x)^2 \f$ at the point $x$ specified by element, quadrature point, and local coordinates.
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return aol::Sqr( _w.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

/**
 * \brief Standard lumped mass matrix \f$\left(\int_\Omega \mathcal{I}_h(\varphi_i \varphi_j)) dx\right)_{ij}\f$
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class LumpedMassOp : public aol::LumpedMassOpInterface<ConfiguratorType, LumpedMassOp<ConfiguratorType> > {
protected:
public:
  typedef typename ConfiguratorType::RealType RealType;

  LumpedMassOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType ); // prevent accidental call, do not implement

  LumpedMassOp ( const typename ConfiguratorType::InitType &Initializer, LUMPED_MASS_OP_MODE Invert )
    : aol::LumpedMassOpInterface<ConfiguratorType, LumpedMassOp<ConfiguratorType> > ( Initializer, Invert ) {}

  LumpedMassOp ( const ConfiguratorType &Config, LUMPED_MASS_OP_MODE Invert )
    : aol::LumpedMassOpInterface<ConfiguratorType, LumpedMassOp<ConfiguratorType> > ( Config, Invert ) {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType& /*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return 1.;
  }
};



//! \brief General Interface for nonlinear FE-operators depending on x, which are locally assembled.
//!        only getNonlinearity has to be provided.
//!
//! not for operators depending on derivatives of
//! the function, use FENonlinDiffOpInterface instead.
/*!
 * \author Droske
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename Imp, int NumComponents = ConfiguratorType::Dim >
class FENonlinIntegrationVectorInterface
      : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  static const int NumOfComponents = NumComponents;
  typedef typename ConfiguratorType::RealType RealType;

protected:
  const typename ConfiguratorType::InitType &_initializer;

public:
  explicit FENonlinIntegrationVectorInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {}

  FENonlinIntegrationVectorInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {}

  virtual ~FENonlinIntegrationVectorInterface() {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    // initialize discrete functions from argument MultiVector
    auto_container<NumComponents, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumComponents; c++ ) {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );
    }

    RealType res = 0.;

    for ( IteratorType it = this->getConfigurator().begin(); it != this->getConfigurator().end(); ++it ) {
      // Necessary because the number of quadpoints may vary
      // getBaseFunctionSet() implicitly calls initializeFromElement
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      RealType a = 0.;
      for ( int q = 0; q < numQuadPoints; ++q )
        a += this->asImp().evaluateIntegrand ( discrFuncsArg, *it, q, bfs.getRefCoord ( q ) )
             * bfs.getWeight ( q );

      a *= this->getConfigurator().vol ( *it );
      res += a;
    }
    Dest += res;
  }

  void applyAddIntegrand ( const aol::MultiVector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> dest( Dest, aol::STRUCT_COPY );
    auto_container<NumComponents, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumComponents; c++ )
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );

    aol::Vector<RealType> locWeight( Dest, aol::STRUCT_COPY );

    RealType *nl_cache = new RealType[ this->getConfigurator().maxNumQuadPoints() ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      for ( int q = 0; q < numQuadPoints; ++q )
        nl_cache[q] = this->asImp().evaluateIntegrand ( discrFuncsArg, *it, q, bfs.getRefCoord ( q ) );

      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        RealType a = 0., b = 0.;

        for ( int q = 0; q < numQuadPoints; ++q ) {
          RealType fac = bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
          a += nl_cache[q] * fac;
          b += fac;
        }

        a *= this->getConfigurator().vol ( *it );
        b *= this->getConfigurator().vol ( *it );

        dest[ this->getConfigurator().localToGlobal ( *it, dof ) ] += a;
        locWeight[ this->getConfigurator().localToGlobal ( *it, dof ) ] += b;
      }
    }
    delete[] nl_cache;

    for ( int i = 0; i < locWeight.size(); ++i )
      dest[i] /= locWeight[i];
    Dest += dest;
  }

  void getIntegrand ( const aol::MultiVector < RealType >& Arg, aol::Vector < RealType >& Dest ) const {
    Dest.setZero ();
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    // initialize discrete functions from argument MultiVector
    aol::auto_container< NumComponents, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumComponents; ++c ) {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );
    }
    for ( IteratorType it = this->getConfigurator().begin(); it != this->getConfigurator().end(); ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType integrand = aol::ZTrait<RealType>::zero;
      for ( int q = 0; q < numQuadPoints; ++q )
        integrand += this->asImp().evaluateIntegrand ( discrFuncsArg, *it, q, bfs.getRefCoord ( q ) ) * bfs.getWeight ( q );
      integrand *= this->getConfigurator().vol ( *it );
      Dest[this->getConfigurator().getConsecutiveElementNumber(*it)] += integrand;
    }
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const auto_container<NumComponents, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};


//! \brief Generalization of the class FENonlinIntegrationVectorInterface: the number of components and the size of components in
//!        MultiVector Arg in applyAdd may differ
//! \author on(?)
//! \ingroup FENonlinOpInt
template <typename ConfiguratorType, typename Imp>
class FENonlinIntegrationVectorGeneralInterface
      : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinIntegrationVectorGeneralInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {}

  FENonlinIntegrationVectorGeneralInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {}

  virtual ~FENonlinIntegrationVectorGeneralInterface() {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {


    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;

    //! for every Arg[c] a descrete function is created
    for ( int c = 0; c < Arg.numComponents(); c++ )   //this line has changed from FENonlinIntegrationVectorInterface
    {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Arg[c] ) );
    }

    RealType res = 0.;

    for ( IteratorType it = this->getConfigurator().begin(); it != this->getConfigurator().end(); ++it ) {
      typedef typename ConfiguratorType::QuadType QType;
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );

      RealType a = 0.;
      for ( int q = 0; q < QType::numQuadPoints; ++q ) {
        a += this->asImp().evaluateIntegrand ( discrFuncsArg, *it, q, bfs.getRefCoord ( q ) ) * bfs.getWeight ( q );
      }
      a *= this->getConfigurator().vol( *it );
      res += a;
    }
    Dest += res;
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};


/**
 * \brief Like FENonlinIntegrationVectorInterface but for vector-valued integrands.
 *
 * \todo Discuss whether the class name is an acceptable extension of the FEOp naming scheme.
 *
 * \author Berkels, Linkmann
 * \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename Imp, int NumComponents = ConfiguratorType::Dim, int NumberOfIntegrandComponents = ConfiguratorType::Dim >
class VectorFENonlinIntegrationVectorInterface
      : public Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Vector<typename ConfiguratorType::RealType> > {
public:
  static const int NumOfComponents = NumComponents;
  typedef typename ConfiguratorType::RealType RealType;
protected:
  mutable ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_initializer;
public:
  VectorFENonlinIntegrationVectorInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      _config ( Initializer ), _initializer ( Initializer ) {}

  virtual ~VectorFENonlinIntegrationVectorInterface() {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {


    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    auto_container<NumComponents, aol::DiscreteFunctionDefault<ConfiguratorType> > discrFuncsArg;
    for ( int c = 0; c < NumComponents; c++ ) {
      discrFuncsArg.set_copy ( c, aol::DiscreteFunctionDefault<ConfiguratorType> ( _initializer, Arg[c] ) );
    }

    aol::Vec<NumberOfIntegrandComponents, RealType> res;

    for ( IteratorType it = _config.begin(); it != _config.end(); ++it ) {
      typedef typename ConfiguratorType::QuadType QType;
      aol::Vec<NumberOfIntegrandComponents, RealType> a;
      for ( int q = 0; q < QType::numQuadPoints; q++ ) {
        aol::Vec<NumberOfIntegrandComponents, RealType> temp;
        this->asImp().evaluateIntegrand ( discrFuncsArg, *it, q, _config.getBaseFunctionSet ( *it ).getRefCoord ( q ), temp );
        a += temp * _config.getBaseFunctionSet ( *it ).getWeight ( q );
      }

      a *= _config.vol ( *it );

      res += a;
    }
    for ( int i = 0; i < NumberOfIntegrandComponents; i++ ) {
      Dest[i] += res[i];
    }
    // Dest += res;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::MultiVector<RealType> argMVec ( 0, 0 );
    argMVec.appendReference ( Arg );
    applyAdd ( argMVec, Dest );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::MultiVector<RealType> argMVec ( 0, 0 );
    argMVec.appendReference ( Arg );
    Op<aol::MultiVector<RealType>, aol::Vector<RealType> >::apply ( argMVec, Dest );
  }

  using Op<aol::MultiVector<RealType>, aol::Vector<RealType> >::apply;

  //! interface function, has to be provided in derived classes.
  void evaluateIntegrand ( const auto_container<NumComponents, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                           const typename ConfiguratorType::ElementType &El,
                           int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,  aol::Vec<NumberOfIntegrandComponents, RealType> &Integrand ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().evaluateIntegrand ( DiscFuncs, El, QuadPoint, RefCoord, Integrand );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};



//! Interface to compute \f$\int_\Omega f(\phi,x) dx\f$,
/** where \f$\phi\f$ is the argument of the operator.
 *  The integrand can be computed using "applyAddIntegrand"
 *
 *  \author Droske, Wirth
 *  \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename Imp>
class FENonlinIntegrationScalarInterface
      : public FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FENonlinIntegrationScalarInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  FENonlinIntegrationScalarInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer ) :
      FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~FENonlinIntegrationScalarInterface( ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    const aol::DiscreteFunctionDefault<ConfiguratorType> discFunc ( this->getConfigurator(), Arg );

    RealType res = 0.;

    const typename IteratorType::EndType end = this->getConfigurator().end();
    for ( IteratorType it = this->getConfigurator().begin(); it != end; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      RealType a = 0.;
      for ( int q = 0; q < numQuadPoints; ++q ) {
        a += this->asImp().evaluateIntegrand ( discFunc, *it, q, bfs.getRefCoord ( q ) ) * bfs.getWeight ( q );
      }

      a *= this->getConfigurator().vol ( *it );
      res += a;
    }
    Dest += res;
  }

  void applyAddIntegrand ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> dest( Dest, aol::STRUCT_COPY );
    const aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc ( _initializer, Arg );

    aol::Vector<RealType> locWeight( Arg, aol::STRUCT_COPY );

    RealType *nl_cache = new RealType[ this->getConfigurator().maxNumQuadPoints() ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      for ( int q = 0; q < numQuadPoints; ++q )
        nl_cache[q] = this->asImp().evaluateIntegrand ( discrFunc, *it, q, bfs.getRefCoord ( q ) );

      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        RealType a = 0., b = 0.;

        for ( int q = 0; q < numQuadPoints; ++q ) {
          RealType fac = bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
          a += nl_cache[q] * fac;
          b += fac;
        }

        a *= this->getConfigurator().vol ( *it );
        b *= this->getConfigurator().vol ( *it );

        dest[ this->getConfigurator().localToGlobal ( *it, dof ) ] += a;
        locWeight[ this->getConfigurator().localToGlobal ( *it, dof ) ] += b;
      }
    }
    delete[] nl_cache;

    for ( int i = 0; i < locWeight.size(); ++i )
      dest[i] /= locWeight[i];
    Dest += dest;
  }

  void getIntegrand ( const aol::Vector< RealType > &Arg, aol::Vector < RealType > &Dest ) const {
    if ( Dest.size() != this->getConfigurator().getInitializer().getNumberOfElements () )
      throw aol::Exception ( "Size mismatch!", __FILE__, __LINE__ );
    Dest.setZero ();
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const aol::DiscreteFunctionDefault<ConfiguratorType> discFunc ( this->getConfigurator(), Arg );
    const typename IteratorType::EndType end = this->getConfigurator().end();
    for ( IteratorType it = this->getConfigurator().begin(); it != end; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType integrand = aol::ZTrait<RealType>::zero;
      for ( int q = 0; q < numQuadPoints; ++q ) {
        integrand += this->asImp().evaluateIntegrand ( discFunc, *it, q, bfs.getRefCoord ( q ) ) * bfs.getWeight ( q );
      }
      integrand *= this->getConfigurator().vol ( *it );
      Dest[this->getConfigurator().getConsecutiveElementNumber(*it)] += integrand;
    }
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



template <typename ConfiguratorType>
class IntegrateFEFunction
      : public FENonlinIntegrationScalarInterface<ConfiguratorType, IntegrateFEFunction<ConfiguratorType> > {
public:

  typedef typename ConfiguratorType::RealType RealType;

  explicit IntegrateFEFunction ( const typename ConfiguratorType::InitType &Initializer )
      :  FENonlinIntegrationScalarInterface<ConfiguratorType, IntegrateFEFunction<ConfiguratorType> > ( Initializer ) { }

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return DiscFunc.evaluateAtQuadPoint ( El, QuadPoint );
  }
};



//! \brief provides an easy interface to matrices of the form \f$ \int_\Omega \phi_i (\nabla \phi_j \cdot v(x)) \f$, where \f$v\f$
//!        is a vector valued function, which has to be provided by the member function getCoefficientVector.
/*!
* This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
* \author Berkels
* \ingroup FELinOpInt
*/
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinVectorWeightedSemiDiffInterface :
      public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                            FELinVectorWeightedSemiDiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > {
protected:
  const typename ConfiguratorType::InitType &_grid;
  const bool _transpose;
public:
  explicit FELinVectorWeightedSemiDiffInterface ( const typename ConfiguratorType::InitType &Grid,
                           const bool Transpose = false,
                           OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
      FELinVectorWeightedSemiDiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Grid, OpType ),
      _grid ( Grid ),
      _transpose ( Transpose ) {  }

  FELinVectorWeightedSemiDiffInterface ( const ConfiguratorType & Config,
                           const typename ConfiguratorType::InitType &Grid,
                           const bool Transpose = false,
                           OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
      FELinVectorWeightedSemiDiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Config, OpType ),
      _grid ( Grid ),
      _transpose ( Transpose ) {  }

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::QuadType QuadType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType  VecType;

  //! this function has to be provided in the implementation (derived class) of the interface
  void getCoefficientVector ( const typename ConfiguratorType::ElementType &El,
                              int QuadPoint,
                              const typename ConfiguratorType::VecType &/*RefCoord*/,
                              typename ConfiguratorType::VecType &Vector ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getCoefficientVector ( El, QuadPoint, Vector );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    for ( int q = 0; q < QuadType::numQuadPoints; ++q ) {
      typename ConfiguratorType::VecType vector;

      this->asImp().getCoefficientVector ( El, q, bfs.getRefCoord ( q ), vector );
      for ( int i = 0; i < numDofs; ++i ) {
        RealType basisi = bfs.evaluate ( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          if ( !_transpose )
            LocalMatrix[i][j] +=  basisi * ( vector * bfs.evaluateGradient ( j, q ) ) * bfs.getWeight ( q ) ;
          else
            LocalMatrix[j][i] +=  basisi * ( vector * bfs.evaluateGradient ( j, q ) ) * bfs.getWeight ( q ) ;
        }
      }
    }

    const RealType vol = this->getConfigurator().vol ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] *= vol;
      }
    }
  }

protected:
  inline Imp &asImp() {
    return static_cast<Imp&> ( *this );
  }
  const Imp &asImp() const {
    return static_cast<const Imp&> ( *this );
  }
};

/**
 * \brief calculates \f$ \int_\Omega\left(\frac{\nabla d}{|\nabla d|_\epsilon}\cdot\nabla \phi_j\right)\phi_i dx \f$.
 *
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class ImageNormalSemiDiffOp : public aol::FELinVectorWeightedSemiDiffInterface<ConfiguratorType, ImageNormalSemiDiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrD;
  const RealType _epsilonSqr;
public:
  ImageNormalSemiDiffOp( const typename ConfiguratorType::InitType &Initializer,
                         const aol::Vector<RealType> &DDofs,
                         const RealType Epsilon,
                         const bool Transponse = false,
                         aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinVectorWeightedSemiDiffInterface<ConfiguratorType,ImageNormalSemiDiffOp<ConfiguratorType> >( Initializer, Transponse, OpType ),
      _discrD( Initializer, DDofs ),
      _epsilonSqr(aol::Sqr(Epsilon)){
  }

  inline void getCoefficientVector( const typename ConfiguratorType::ElementType &El,
                                    int QuadPoint,
                                    const typename ConfiguratorType::VecType &/*RefCoord*/,
                                    typename ConfiguratorType::VecType &Vector ) const {
    _discrD.evaluateGradientAtQuadPoint( El, QuadPoint, Vector );
    Vector *= 1./sqrt ( Vector.normSqr() + _epsilonSqr );
  }
};

/**
 * \brief calculates \f$ \int_\Omega\left(\nabla d\cdot\nabla \phi_j\right)\phi_i dx \f$.
 *
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class ImageGradientSemiDiffOp : public aol::FELinVectorWeightedSemiDiffInterface<ConfiguratorType, ImageGradientSemiDiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrD;
public:
  ImageGradientSemiDiffOp( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Vector<RealType> &DDofs,
                           const bool Transponse = false,
                           aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinVectorWeightedSemiDiffInterface<ConfiguratorType,ImageGradientSemiDiffOp<ConfiguratorType> >( Initializer, Transponse, OpType ),
      _discrD( Initializer, DDofs ){
  }

  inline void getCoefficientVector( const typename ConfiguratorType::ElementType &El,
                                    int QuadPoint,
                                    const typename ConfiguratorType::VecType &/*RefCoord*/,
                                    typename ConfiguratorType::VecType &Vector ) const {
    _discrD.evaluateGradientAtQuadPoint( El, QuadPoint, Vector );
  }
};

//! \brief provides an easy interface to matrices of the form \f$ \int_\Omega \phi_i \partial_k\phi_j v(x) \f$, where \f$v\f$
//!        is a scalar valued function, which has to be provided by the member function getCoeff.
/*!
* This class works in 2D and 3D and can be endowed with arbitrary quadrature rules.
* \author Berkels
* \ingroup FELinOpInt
*/
template <typename ConfiguratorType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class FELinScalarWeightedSemiDiffInterface :
      public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                            FELinScalarWeightedSemiDiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > {
protected:
  const typename ConfiguratorType::InitType &_grid;
  const int _gradientIndex;
  const bool _transpose;
public:
  FELinScalarWeightedSemiDiffInterface ( const typename ConfiguratorType::InitType &Grid,
                                       const int IndexComponent,
                                       const bool Transpose = false,
                                       OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
      FELinScalarWeightedSemiDiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Grid, OpType ),
      _grid ( Grid ),
      _gradientIndex(IndexComponent),
      _transpose ( Transpose ) {  }

  FELinScalarWeightedSemiDiffInterface ( const ConfiguratorType & Config,
                                       const typename ConfiguratorType::InitType &Grid,
                                       const int IndexComponent,
                                       const bool Transpose = false,
                                       OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
      FELinScalarWeightedSemiDiffInterface<ConfiguratorType, Imp, IndexMode>, IndexMode > ( Config, OpType ),
      _grid ( Grid ),
      _gradientIndex(IndexComponent),
      _transpose ( Transpose ) {  }

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::QuadType QuadType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType  VecType;

  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El,
                             int QuadPoint,
                             const typename ConfiguratorType::DomVecType& RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    for ( int q = 0; q < QuadType::numQuadPoints; ++q ) {
      typename ConfiguratorType::RealType coeff;
      coeff = this->asImp().getCoeff( El, q, bfs.getRefCoord( q ) );
      for ( int i = 0; i < numDofs; ++i ) {
        RealType basisi = bfs.evaluate ( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          if ( !_transpose )
            LocalMatrix[i][j] +=  basisi * ( coeff * bfs.evaluateGradient ( j, q )[_gradientIndex] ) * bfs.getWeight ( q );
          else
            LocalMatrix[j][i] +=  basisi * ( coeff * bfs.evaluateGradient ( j, q )[_gradientIndex] ) * bfs.getWeight ( q );
        }
      }
    }

    const RealType vol = this->getConfigurator().vol ( El );

    for ( int i = 0; i < numDofs; ++i ) {
      for ( int j = 0; j < numDofs; ++j ) {
        LocalMatrix[i][j] *= vol;
      }
    }
  }

protected:
  inline Imp &asImp() {
    return static_cast<Imp&> ( *this );
  }
  const Imp &asImp() const {
    return static_cast<const Imp&> ( *this );
  }
};

/**
 * Provides an easy interface to mixed mass-stiffness matrices of the form
 * \f$ \left(\int_\Omega \gamma(x)\phi_i\nabla\phi_j dx\right)_{ij} \f$, where \f$ \gamma \f$
 * is a Finite Element function given by the aol::Vector Weight.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class WeightedSemiDiffOp :
  public aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, WeightedSemiDiffOp< ConfiguratorType > > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrWeight;

public:
  WeightedSemiDiffOp( const typename ConfiguratorType::InitType &Initializer,
                      const aol::Vector<RealType> &Weight,
                      const int IndexComponent,
                      const bool Transpose = false,
                      aol::OperatorType OpType = ONTHEFLY ) :
    aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, WeightedSemiDiffOp< ConfiguratorType > > ( Initializer, IndexComponent, Transpose, OpType ),
    _discrWeight( Initializer, Weight ) {
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El,
                             int QuadPoint,
                             const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return _discrWeight.evaluateAtQuadPoint(El, QuadPoint);
  }
};

/**
 * As aol::WeightedSemiDiffOp<> with constant weight function, i.e. \f$ \gamma \equiv 1 \f$
 *
 * \author Heeren
 */
template <typename ConfiguratorType>
class SemiDiffOp :
  public aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, SemiDiffOp< ConfiguratorType > > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

public:
  SemiDiffOp( const typename ConfiguratorType::InitType &Initializer,
              const int IndexComponent,
              const bool Transpose = false,
              aol::OperatorType OpType = ONTHEFLY ) :
    aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, SemiDiffOp< ConfiguratorType > > ( Initializer, IndexComponent, Transpose, OpType ){ }

  SemiDiffOp( const ConfiguratorType& Config,
              const typename ConfiguratorType::InitType &Initializer,
              const int IndexComponent,
              const bool Transpose = false,
              aol::OperatorType OpType = ONTHEFLY ) :
      aol::FELinScalarWeightedSemiDiffInterface< ConfiguratorType, SemiDiffOp< ConfiguratorType > > ( Config, Initializer, IndexComponent, Transpose, OpType ){ }


 inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
   return 1.;
 }
};

/*!
 * \brief Operator for stiffness matrices of the form \f[ \int_\Omega w_1 \times w_2 \nabla u \cdot \nabla w dx \f]
 * \author Droske
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class TensorStiffOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, TensorStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> _w1;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> _w2;

public:

  TensorStiffOp ( const typename ConfiguratorType::InitType &Initializer,
                  const aol::MultiVector<RealType> &W1Dofs,
                  const aol::MultiVector<RealType> &W2Dofs,
                  OperatorType OpType = ONTHEFLY )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, TensorStiffOp<ConfiguratorType> > ( Initializer,
                                                                          OpType ),
      _w1 ( Initializer, W1Dofs ), _w2 ( Initializer, W2Dofs ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::VecType&,
                               typename ConfiguratorType::MatType &Matrix ) const {
    typename ConfiguratorType::DomVecType w1, w2;

    _w1.evaluateAtQuadPoint ( El, QuadPoint, w1 );
    _w2.evaluateAtQuadPoint ( El, QuadPoint, w2 );

    Matrix.makeTensorProduct ( w1, w2 );
  }
};

/*!
 * \brief Interface that returns the summands inside the integrand implemented by overloading evaluateIntegrand as a vector for usage in least squares regression
 * \author Tatano
 * \ingroup FENonlinOpInt
 * \todo Make the class name more consistent with the existing FE interface classes.
 * \todo This class didn't initially take into account the element volumes. Since existing code depends
 *       on this old behavior, it's still the default behavior (see template scaleWithElVol).
 *       The existing code should be changed so that the old behavior can be dropped.
 */
template <typename ConfiguratorType, typename Imp, typename MatrixType, bool scaleWithElVol = false>
class FELeastSquaresFunctionalInterface : public FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> >{
public:
  typedef typename ConfiguratorType::RealType RealType;

  explicit FELeastSquaresFunctionalInterface ( const typename ConfiguratorType::InitType &Grid )
    : FEOpInterface<ConfiguratorType, aol::Vector<RealType> > ( Grid ) { }

  explicit FELeastSquaresFunctionalInterface ( const ConfiguratorType & Config )
    : FEOpInterface<ConfiguratorType, aol::Vector<RealType> > ( Config ) { }

  virtual ~FELeastSquaresFunctionalInterface( ) {}

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    const aol::DiscreteFunctionDefault<ConfiguratorType> discFunc ( this->getConfigurator(), Arg );

    const typename IteratorType::EndType end = this->getConfigurator().end();
    int i = 0;
    for ( IteratorType it = this->getConfigurator().begin(); it != end; ++it ) {
      const RealType scale = scaleWithElVol ? this->getConfigurator().vol ( *it ) : aol::ZOTrait<RealType>::one;
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      for ( int q = 0; q < numQuadPoints; ++q ) {
        //i+q = numQuadPoints * this->getConfigurator().getConsecutive
        Dest[i+q] = this->asImp().evaluateIntegrand( discFunc, *it, q, bfs.getRefCoord ( q ) ) * sqrt( bfs.getWeight ( q ) * scale );
      }
      i+=numQuadPoints;
    }
  }

  void applyAdd ( const aol::Vector<RealType> &, aol::Vector<RealType> & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }

  void applyMultBasisFunction ( const aol::Vector<RealType> &Arg, MatrixType &Dest ) const {
    aol::Vector<RealType> integrandWithoutBasisFunction ( Dest.getNumRows ( ) );
    apply ( Arg, integrandWithoutBasisFunction );

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    int i=0;
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );
      for ( int dof = 0; dof < numLocalDofs; dof++ )
        for ( int q = 0; q < numQuadPoints; ++q )
          Dest.set(i+q, this->getConfigurator().localToGlobal ( *it, dof ), integrandWithoutBasisFunction[i+q] * bfs.evaluate ( dof, q ) );
      i+=numQuadPoints;
    }
  }

  void applyMultDerBasisFunction ( const aol::Vector<RealType> &Arg, int gradComp, MatrixType &Dest ) const {
    aol::Vector<RealType> integrandWithoutBasisFunction ( Dest.getNumRows ( ) );
    apply ( Arg, integrandWithoutBasisFunction );

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = this->getConfigurator().end();
    int i=0;
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );
      for ( int dof = 0; dof < numLocalDofs; dof++ )
        for ( int q = 0; q < numQuadPoints; ++q )
          Dest.set(i+q, this->getConfigurator().localToGlobal ( *it, dof ),
                   integrandWithoutBasisFunction[i+q] * (bfs.evaluateGradient (dof, q))[gradComp] );
      i+=numQuadPoints;
    }
  }


  //! interface function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

/**
 * \author Berkels
 */
template < typename ConfiguratorType >
int totalNumberOfQuadPoints ( const ConfiguratorType &Config ) {
  typedef typename ConfiguratorType::ElementIteratorType IteratorType;
  const typename IteratorType::EndType end = Config.end();
  int i = 0;
  for ( IteratorType it = Config.begin(); it != end; ++it ) {
    i += Config.getBaseFunctionSet ( *it ).numQuadPoints( );
  }
  return i;
}

/**
 * \brief Interface for derivatives of aol::FELeastSquaresFunctionalInterface
 * \author Berkels
 * \ingroup FENonlinOpInt
 * \todo Make the class name more consistent with the existing FE interface classes.
 */
template <typename ConfiguratorType, typename Imp, typename MatrixType>
class FELeastSquaresFunctionalVectorDerivativeInterface : public FEOpInterface<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, MatrixType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
public:
  explicit FELeastSquaresFunctionalVectorDerivativeInterface ( const typename ConfiguratorType::InitType &Initializer )
    : FEOpInterface<ConfiguratorType, aol::MultiVector<RealType>, MatrixType> ( Initializer ), _initializer ( Initializer ) {}

  FELeastSquaresFunctionalVectorDerivativeInterface ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer )
    : FEOpInterface<ConfiguratorType, aol::MultiVector<RealType> > ( Config ), _initializer ( Initializer ) {}

  virtual ~FELeastSquaresFunctionalVectorDerivativeInterface( ) {}

  void apply ( const aol::MultiVector<RealType> &Arg, MatrixType &Dest ) const {

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> discrFunc ( _initializer, Arg );

    typename aol::Vec<ConfiguratorType::Dim, RealType> nl;

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    int i=0;
    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const RealType vol = this->getConfigurator().vol ( *it );

      for ( int q = 0; q < numQuadPoints; ++q ) {
        this->asImp().getNonlinearity ( discrFunc, *it, q, bfs.getRefCoord ( q ), nl );
        const RealType sqrtWeight = sqrt( bfs.getWeight ( q ) * vol );
        for ( int dof = 0; dof < numLocalDofs; ++dof ) {
          const RealType factor = bfs.evaluate ( dof, q ) * sqrtWeight;
          const int globalDofIdx = this->getConfigurator().localToGlobal ( *it, dof );
          for ( int d = 0; d < ConfiguratorType::Dim; ++d ) {
            Dest.getReference ( 0 , d ).set ( i+q, globalDofIdx, factor * nl[d] );
          }
        }
      }
      i+=numQuadPoints;
    }
  }

  void applyAdd ( const aol::MultiVector<RealType> &, MatrixType & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }

  int getDimRangeF ( ) const {
    return totalNumberOfQuadPoints ( this->getConfigurator() );
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         aol::Vec<ConfiguratorType::Dim, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, QuadPoint, RefCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

}  // end namespace





#endif
