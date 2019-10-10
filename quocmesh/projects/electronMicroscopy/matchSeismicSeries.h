#ifndef __MATCHSEISMICSERIES_H
#define __MATCHSEISMICSERIES_H

#include <matchSeries.h>

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ConfiguratorType1D, typename RegistrationConfiguratorType1D/*, typename GradientDescentType*/>
void interleaveImages ( const typename ConfiguratorType::ArrayType &Reference,
                        const aol::RandomAccessContainer<typename ConfiguratorType::ArrayType> &Templates,
                        const typename ConfiguratorType::RealType Lambda,
                        typename ConfiguratorType::ArrayType &InterleavedImage,
                        const bool Align = true,
                        const bool JustCopyLastRefLine = false ) {

  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::DirichletRegularizationConfigurator<ConfiguratorType1D> RegularizationConfiguratorType1D;
  typedef qc::StandardRegistration<ConfiguratorType1D, typename ConfiguratorType1D::ArrayType, RegistrationConfiguratorType1D, RegularizationConfiguratorType1D, aol::H1GradientDescent<ConfiguratorType1D, aol::MultiVector<RealType>, qc::CholeskyBasedInverseH1Metric<ConfiguratorType1D> > > STDRegisType;

  const int lineSize = Reference.getNumX();
  const int numLines = Reference.getNumY();

  InterleavedImage = Reference;

  const int refLineDistance = Templates.size() + 1;

  const typename ConfiguratorType1D::InitType grid1D ( aol::Vec3<int> ( lineSize, 1, 1 ) );

  int lastNonEmptyRefLine = -1;

  for ( int line = 0; line < numLines; ++line ) {
    qc::ScalarArray<RealType, qc::QC_1D> refLineA ( lineSize );
    Reference.getLine ( qc::QC_Y, line, refLineA );

    if ( refLineA.normSqr() == 0 ) {
      // If the reference is missing the first line, copy it from the template.
      // Could perhaps be improved by "half" aligning it to the next reference line.
      if ( line == 0 ) {
        qc::ScalarArray<RealType, qc::QC_1D> templLine ( lineSize );
        Templates[0].getLine ( qc::QC_Y, line, templLine );
        InterleavedImage.putLine ( qc::QC_Y, line, templLine );
      }
      else if ( JustCopyLastRefLine && ( lastNonEmptyRefLine >= 0 ) )
      {
        Reference.getLine ( qc::QC_Y, lastNonEmptyRefLine, refLineA );
        InterleavedImage.putLine ( qc::QC_Y, line, refLineA );
      }
      continue;
    }
    else
      lastNonEmptyRefLine = line;

    if ( JustCopyLastRefLine )
      continue;

    if ( line + refLineDistance >= numLines )
      continue;

    qc::ScalarArray<RealType, qc::QC_1D> refLineB ( lineSize );
    Reference.getLine ( qc::QC_Y, line + refLineDistance, refLineB );
    qc::ScalarArray<RealType, qc::QC_1D> templLine ( lineSize );

    if ( refLineB.normSqr() == 0 )
      throw aol::Exception ( aol::strprintf ( "Unexpected empty line in reference (current line index is %d)", line ).c_str(), __FILE__, __LINE__ );

    qc::MultiArray<RealType, qc::QC_1D> lineShiftA ( grid1D );
    qc::MultiArray<RealType, qc::QC_1D> lineShiftB ( grid1D );
    for ( int i = 0; i < refLineDistance - 1; ++i ){
      Templates[i].getLine ( qc::QC_Y, line + i + 1, templLine );

      if ( templLine.normSqr() == 0 )
        throw aol::Exception ( aol::strprintf ( "Unexpected empty line in template (current line index is %d)", line ).c_str(), __FILE__, __LINE__ );

      if ( Align ) {
        const RegistrationConfiguratorType1D regisConfig;
        const RegularizationConfiguratorType1D regulConfig ( grid1D );

        STDRegisType stdRegisA ( refLineA, templLine, regisConfig, regulConfig, Lambda );
        stdRegisA.findTransformation ( lineShiftA );

        STDRegisType stdRegisB ( refLineB, templLine, regisConfig, regulConfig, Lambda );
        stdRegisB.findTransformation ( lineShiftB );

        qc::MultiArray<RealType, qc::QC_1D> lineShift ( lineShiftA );
        lineShift.scaleAndAddMultiple ( ( refLineDistance - i - 1. ) / refLineDistance, lineShiftB, ( 1. + i ) / refLineDistance );

        qc::ScalarArray<RealType, qc::QC_1D> shiftedTemplLine ( lineSize );
        qc::DeformImage<ConfiguratorType1D> ( templLine, grid1D, shiftedTemplLine, lineShift );
        InterleavedImage.putLine ( qc::QC_Y, line + 1 + i, shiftedTemplLine );
      }
      else
        InterleavedImage.putLine ( qc::QC_Y, line + 1 + i, templLine );
    }
  }
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
void removeEmptyLinesFromImage ( const typename ConfiguratorType::ArrayType &Image, typename ConfiguratorType::ArrayType &CleanedImage, const bool KeepEmptyTopBottomBlocks = false ) {
  typedef typename ConfiguratorType::RealType RealType;
  const int lineSize = Image.getNumX();
  const int numLines = Image.getNumY();

  qc::ScalarArray<RealType, qc::QC_1D> imageLine ( lineSize );
  aol::Vector<int> numEmptyLineNums;
  aol::Vector<int> numNonEmptyLineNums;
  for ( int line = 0; line < numLines; ++line ) {
    Image.getLine ( qc::QC_Y, line, imageLine );

    // Line is empty or just contains NaNs.
    if ( ( imageLine.normSqr() == 0 ) || ( aol::VectorHelper<aol::Vector<RealType> >::numNonNaNs ( imageLine ) == 0 ) )
      numEmptyLineNums.pushBack ( line );
    else
      numNonEmptyLineNums.pushBack ( line );
  }

  if ( numNonEmptyLineNums.size() == 0 )
    throw aol::Exception ( "Input contains only empty lines", __FILE__, __LINE__ );

  const int firstNonEmptyLine = numNonEmptyLineNums[0];
  const int lastNonEmptyLine = numNonEmptyLineNums[numNonEmptyLineNums.size()-1];
  const int averageNonEmptyLineDist = aol::Rint ( ( lastNonEmptyLine - firstNonEmptyLine ) / static_cast<RealType> ( numNonEmptyLineNums.size() - 1 ) );

  const int numTopLinesToKeep = KeepEmptyTopBottomBlocks ? ( firstNonEmptyLine / averageNonEmptyLineDist ) : 0;
  const int numBottomLinesToKeep = KeepEmptyTopBottomBlocks ? ( ( numLines - 1 - lastNonEmptyLine ) / averageNonEmptyLineDist ) : 0;

  CleanedImage.reallocate ( lineSize, numNonEmptyLineNums.size() + numTopLinesToKeep + numBottomLinesToKeep );

  for ( int line = 0; line < numNonEmptyLineNums.size(); ++line ) {
    Image.getLine ( qc::QC_Y, numNonEmptyLineNums[line], imageLine );
    CleanedImage.putLine ( qc::QC_Y, line + numTopLinesToKeep, imageLine );
  }
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
void replaceEmptyLinesWithNaN ( typename ConfiguratorType::ArrayType &Image ) {
  const int lineSize = Image.getNumX();
  const int numLines = Image.getNumY();

  qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_1D> imageLine ( lineSize );
  qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_1D> nanLine ( lineSize );
  nanLine.setAll ( aol::NumberTrait<typename ConfiguratorType::RealType>::NaN );
  aol::Vector<int> numEmptyLineNums;
  for ( int line = 0; line < numLines; ++line ) {
    Image.getLine ( qc::QC_Y, line, imageLine );

    // Line is empty
    if ( imageLine.normSqr() == 0 )
      Image.putLine ( qc::QC_Y, line, nanLine );
  }
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, template<class> class RegistrationConfiguratorType,
          template<class> class RegularizationConfiguratorType = qc::DirichletRegularizationConfigurator,
          typename _DyadicConfiguratorType1D = qc::RectangularGridConfigurator<typename ConfiguratorType::RealType, qc::QC_1D, aol::GaussQuadrature<typename ConfiguratorType::RealType, qc::QC_1D, 3> >,
          typename DyadicConfiguratorType = qc::RectangularGridConfigurator<typename ConfiguratorType::RealType, qc::QC_2D, aol::GaussQuadrature<typename ConfiguratorType::RealType, qc::QC_2D, 3>, qc::CellCenteredCubicGrid<qc::QC_2D> >,
          typename _DyadicRegistrationMultilevelDescentType = XShiftRegistrationMultilevelDescent<DyadicConfiguratorType, _DyadicConfiguratorType1D, qc::NCCRegistrationConfigurator<_DyadicConfiguratorType1D> >,
          typename NonDyadicGradientDescentType = aol::H1GradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, qc::CholeskyBasedInverseH1Metric<ConfiguratorType> >
>
class XShiftNonDyadicRegistrationMultilevelDescent
  : public NonDyadicRegistrationMultilevelDescent
             <ConfiguratorType, RegistrationConfiguratorType, RegularizationConfiguratorType, DyadicConfiguratorType, _DyadicRegistrationMultilevelDescentType > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef _DyadicConfiguratorType1D DyadicConfiguratorType1D;
  typedef _DyadicRegistrationMultilevelDescentType DyadicRegistrationMultilevelDescentType;
  typedef typename qc::MultilevelArrayTrait<RealType, typename DyadicConfiguratorType::InitType>::MultilevelArrayType MultilevelArrayType;

  XShiftNonDyadicRegistrationMultilevelDescent ( const aol::ParameterParser &Parser )
    : NonDyadicRegistrationMultilevelDescent
        <ConfiguratorType, RegistrationConfiguratorType, RegularizationConfiguratorType, DyadicConfiguratorType, DyadicRegistrationMultilevelDescentType > ( Parser ) {
    initMLADroppingEmptyLines ( this->getTemplImageReference(), this->_pMld->getTemplImageMLAReference() );
    initMLADroppingEmptyLines ( this->getRefImageReference(), this->_pMld->getRefImageMLAReference()  );
    replaceEmptyLinesWithNaN<ConfiguratorType> ( this->getNonConstTemplImageReference() );
    replaceEmptyLinesWithNaN<ConfiguratorType> ( this->getNonConstRefImageReference() );
  }

  void initMLADroppingEmptyLines ( const ArrayType &Image, MultilevelArrayType &MLArray ) const {
    ArrayType cleanedImage;
    removeEmptyLinesFromImage<ConfiguratorType> ( Image, cleanedImage, this->_parser.getBool ( "keepEmptyTopBottomBlocks" ) );
    MLArray[this->_pMld->getMaxGridDepth()].resampleFrom ( cleanedImage );
    MLArray.levRestrict ( 0, this->_pMld->getMaxGridDepth() );
  }

  virtual RealType updateNonDyadicDeformation ( qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    typedef qc::RectangularGridConfigurator<RealType, qc::QC_1D, aol::GaussQuadrature<RealType, qc::QC_1D, 3> > DyConfType1D;
    if ( this->_parser.checkAndGetBool ( "forceNonDyadicXShiftYConst" ) == false )
      return findXShift<ConfiguratorType, DyConfType1D, qc::NCCRegistrationConfigurator<DyConfType1D>, NonDyadicGradientDescentType>
      ( this->_grid, this->getRefImageReference(), this->getTemplImageReference(), this->_lambda, this->_parser.getDouble ( "yRegFactor" ), this->_parser.getDouble ( "laplaceWeight" ), Phi[0], NULL, this->_parser.getIntOrDefault ( "maxGDIterations", 500 ) );
    else
      return findXShiftYConst<ConfiguratorType, DyConfType1D, qc::NCCRegistrationConfigurator<DyConfType1D> > ( this->_grid, this->getRefImageReference(), this->getTemplImageReference(), Phi[0], NULL, this->_parser.getIntOrDefault ( "maxGDIterations", 500 ) );
  }

  virtual void loadRefOrTemplate ( const char* Filename, const qc::REGISTRATION_INPUT_TYPE Input, const bool NoScaling ) {
    this->loadAndPrepareImage ( Filename, ( Input == qc::REFERENCE ) ? this->getNonConstRefImageReference() : this->getNonConstTemplImageReference(), NoScaling );
    initMLADroppingEmptyLines ( ( Input == qc::REFERENCE ) ? this->getRefImageReference() : this->getTemplImageReference(),
                               ( Input == qc::REFERENCE ) ? this->_pMld->getRefImageMLAReference() : this->_pMld->getTemplImageMLAReference ( ) );
    replaceEmptyLinesWithNaN<ConfiguratorType> ( ( Input == qc::REFERENCE ) ? this->getNonConstRefImageReference() : this->getNonConstTemplImageReference() );
  }

  virtual void setTemplate ( const ArrayType &Template ) {
    this->getNonConstTemplImageReference() = Template;
    initMLADroppingEmptyLines ( this->getTemplImageReference(), this->_pMld->getTemplImageMLAReference() );
    replaceEmptyLinesWithNaN<ConfiguratorType> ( this->getNonConstTemplImageReference() );
  }

  virtual void setReference ( const ArrayType &Reference ) {
    this->getNonConstRefImageReference() = Reference;
    initMLADroppingEmptyLines ( this->getRefImageReference(), this->_pMld->getRefImageMLAReference()  );
    replaceEmptyLinesWithNaN<ConfiguratorType> ( this->getNonConstRefImageReference() );
  }
};

#endif // MATCHSEISMICSERIES_H
