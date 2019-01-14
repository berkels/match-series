#ifndef __SIMULTANEOUSLYASSEMBLED_H
#define __SIMULTANEOUSLYASSEMBLED_H

#include <aol.h>
#include <op.h>
#include <baseFunctionSet.h>
#include <sparseMatrices.h>
#include <fastUniformGridMatrix.h>
#include <discreteFunction.h>
#include <FEOpInterface.h>

namespace aol {

//! This operator stores a linear combination of different operators. The difference to
//! to the LinCombOp is that the mesh is traversed only once. For each element of the
//! mesh all operators are evaluated. Hence the SimultaneouslyAssembledLinCombOp shows
//! up with a lot of different possibilities for the evaluation of FE-basis-functions
//! (mass, stiffness, matrix-stiffness etc.).
//! This operator is useful if a traverse of the mesh is expensive (like in the RLE-
//! datastructure by Michael Nielsen).
/*!
 * \author Nemitz
 */


enum ASSEMBLE_GETCOEFF_OP_TYPE {
  SCALED_MASS_INT,
  SCALED_STIFF_INT
};

enum ASSEMBLE_GETCOEFFMATRIX_OP_TYPE {
  MATRIX_STIFF_INT
};


template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class SimultaneouslyAssembledLinCombOp :
  public FELinOpInterface< typename ConfiguratorType::RealType, ConfiguratorType,
                        SimultaneouslyAssembledLinCombOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::InitType                   InitType;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::QuadType                   QuadType;
  typedef typename ConfiguratorType::MatType                    MatType;
  typedef typename ConfiguratorType::VecType                    VecType;
  typedef Op<aol::Vector<RealType> >                            OpDR;
  typedef provideVirtualGetCoeffClass<ConfiguratorType>         getCoeffOpType;
  typedef provideVirtualGetCoeffMatrixClass<ConfiguratorType>   getCoeffMatrixOpType;

protected:
  ConfiguratorType              _config;
  const InitType                &_grid;

  // the single operator-lists
  list<const getCoeffOpType*>         _scaledMassIntOpList;
  list<const getCoeffOpType*>         _scaledStiffIntOpList;
  list<const getCoeffMatrixOpType*>   _matrixStiffIntOpList;

  // the list of the coefficients of the single operators
  list<RealType>                      _scaledMassIntCoeffList;
  list<RealType>                      _scaledStiffIntCoeffList;
  list<RealType>                      _matrixStiffIntCoeffList;

  // flags whether to process single types of operators or not
  bool                               _thereAreScaledMassIntOps;
  bool                               _thereAreScaledStiffIntOps;
  bool                               _thereAreMatrixStiffIntOps;


public:
  SimultaneouslyAssembledLinCombOp (  const typename ConfiguratorType::InitType &Grid, OperatorType OpType = ONTHEFLY )
      : FELinOpInterface < typename ConfiguratorType::RealType, ConfiguratorType,
        SimultaneouslyAssembledLinCombOp<ConfiguratorType, IndexMode>, IndexMode > ( _config, OpType ),
      _config ( Grid ),
      _grid ( Grid ) {
        _thereAreScaledMassIntOps = false;
        _thereAreScaledStiffIntOps = false;
        _thereAreMatrixStiffIntOps = false;
  }


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {
    return 1.;
    //     return this->asImp().getCoeff ( El, QuadPoint, RefCoord );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    const int numDofs = _config.getNumLocalDofs ( El );

    // clear the local matrix
    for ( int i = 0; i < numDofs; i++ ) {
      for ( int j = 0; j < numDofs; j++ ) {
        LocalMatrix[i][j] = 0.;
      }
    }

    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    // traverse the quadrature-points
    for ( int q = 0; q < numQuadPoints; q++ ) {

      // Here come the evaluations of the single interface-types (mass, stiffness, matrix-stiffness...)
      if ( _thereAreScaledMassIntOps ) {
        // traverse all the operators in the scaledMassIntList
        typename list<const getCoeffOpType* >::const_iterator    opIt;
        typename list<RealType>::const_iterator                  coeffIt;

        RealType coeff = 0;
        // call the getCoeff-method from each operator and multiply it with the belonging coefficient
        for ( opIt = _scaledMassIntOpList.begin(), coeffIt = _scaledMassIntCoeffList.begin();
                    opIt != _scaledMassIntOpList.end(); opIt++, coeffIt++ )
          coeff += ( *opIt )->v_getCoeff( El, q, bfs.getRefCoord ( q ) ) * ( *coeffIt );

        for ( int i = 0; i < numDofs; i++ ) {
          RealType basisi = bfs.evaluate ( i, q );
          for ( int j = i; j < numDofs; j++ ) {
            LocalMatrix[i][j] +=  basisi * bfs.evaluate ( j, q ) * coeff * bfs.getWeight ( q ) ;
          }
        }

      }   // end if there are scaledMassOps

      if ( _thereAreScaledStiffIntOps ) {
        // traverse all the operators in the scaledMassIntList
        typename list<const getCoeffOpType* >::const_iterator    opIt;
        typename list<RealType>::const_iterator                  coeffIt;

        RealType coeff = 0;
        // call the getCoeff-method from each operator and multiply it with the belonging coefficient
        for ( opIt = _scaledStiffIntOpList.begin(), coeffIt = _scaledStiffIntCoeffList.begin();
              opIt != _scaledStiffIntOpList.end(); opIt++, coeffIt++ )
          coeff += ( *opIt )->v_getCoeff( El, q, bfs.getRefCoord ( q ) ) * ( *coeffIt );

        // now add the contribution of each pair of basisfunctions to the local matrix
        for ( int i = 0; i < numDofs; i++ ) {
          const VecType &basisveci  = bfs.evaluateGradient ( i, q );
          for ( int j = i; j < numDofs; j++ ) {
            LocalMatrix[i][j] += ( basisveci * bfs.evaluateGradient ( j, q ) ) * coeff * bfs.getWeight ( q );
          }
        }

      }   // end if there are scaledStiffOps

      if ( _thereAreMatrixStiffIntOps ) {
        // traverse all the operators in the scaledMassIntList
        typename list<const getCoeffMatrixOpType* >::const_iterator    opIt;
        typename list<RealType>::const_iterator                        coeffIt;

        MatType mat;
        VecType matgrad1;

        for ( opIt = _matrixStiffIntOpList.begin(), coeffIt = _matrixStiffIntCoeffList.begin();
              opIt != _matrixStiffIntOpList.end(); opIt++, coeffIt++ )
        {
          // call the getCoeffMatrix-method from the current operator
          ( *opIt )->getCoeffMatrix ( El, q, bfs.getRefCoord ( q ), mat );

          for ( int i = 0; i < numDofs; i++ ) {
            mat.mult ( this->getBaseFunctionSet ( El ).evaluateGradient ( i, q ), matgrad1 );
            for ( int j = i; j < numDofs; j++ ) {
              LocalMatrix[i][j] += ( matgrad1 * this->getBaseFunctionSet ( El ).evaluateGradient ( j, q ) )
                                   * this->getBaseFunctionSet ( El ).getWeight ( q ) * ( *coeffIt );
            }
          }

        }   // end operator loop
      }   // end if there are scaledStiffOps

    }   // end quadrature-point-loop

    for ( int i = 0; i < numDofs; i++ ) {
      for ( int j = i; j < numDofs; j++ ) {
        LocalMatrix[j][i] = LocalMatrix[i][j] = LocalMatrix[i][j] * _config.vol ( El );
      }
    }
  }

  //! Append reference to a GetCoeff-operator
  void appendGetCoeffOpReference ( const getCoeffOpType &Op, ASSEMBLE_GETCOEFF_OP_TYPE assembleOpType,
                                   RealType coeff = static_cast<RealType> ( 1.0 ) ) {
    switch( assembleOpType )  {
      case SCALED_MASS_INT:
        _scaledMassIntOpList.push_back( &Op );
        _scaledMassIntCoeffList.push_back( coeff );
        _thereAreScaledMassIntOps = true;
        break;
      case SCALED_STIFF_INT:
        _scaledStiffIntOpList.push_back( &Op );
        _scaledStiffIntCoeffList.push_back( coeff );
        _thereAreScaledStiffIntOps = true;
        break;
      default:
        throw aol::UnimplementedCodeException ( "SimultaneouslyAssembledLinCombOp::appendGetCoeffOpReference: unknown ASSEMBLE_GETCOEFF_OP_TYPE", __FILE__, __LINE__ );
    }
  }

  //! Append reference to a GetCoeffMatrix-operator
  void appendGetCoeffMatrixOpReference ( const getCoeffMatrixOpType &Op, ASSEMBLE_GETCOEFFMATRIX_OP_TYPE assembleOpType,
                                         RealType coeff = static_cast<RealType> ( 1.0 ) ) {
    switch( assembleOpType )  {
      case MATRIX_STIFF_INT:
        _matrixStiffIntOpList.push_back( &Op );
        _matrixStiffIntCoeffList.push_back( coeff );
        _thereAreMatrixStiffIntOps = true;
        break;
      default:
        throw aol::UnimplementedCodeException ( "SimultaneouslyAssembledLinCombOp::appendGetCoeffOpReference: unknown ASSEMBLE_GETCOEFF_OP_TYPE", __FILE__, __LINE__ );
    }
  }

};



}

#endif
