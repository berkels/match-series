#ifndef __REGISTRATION_H
#define __REGISTRATION_H

#include <parameterParser.h>
#include <AmbrosioTortorelli.h>
#include <hyperelastic.h>
#include <deformations.h>
#include <imageTools.h>
#include <gradientDescent.h>
#include <gradientflow.h>
#include <quocDescent.h>
#include <quocTimestepSaver.h>
#include <cellCenteredGrid.h>
#include <anisoStiffOps.h>
#include <quocOps.h>

namespace qc {

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ImageDOFType = typename ConfiguratorType::ArrayType>
class RegistrationInterface {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType InitType;

  static const int _dim = ConfiguratorType::Dim;
  RealType _tau;

  //! Parameter for the weight of the deformation regularization term
  const RealType _lambda;

private:
  //! Pointer to the reference image.
  aol::DeleteFlagPointer<ImageDOFType> _pReference0;

  //! Pointer to the template image.
  aol::DeleteFlagPointer<ImageDOFType> _pTemplate0;

public:

  const InitType _grid;

  RegistrationInterface ( const ImageDOFType &Reference0,
                          const ImageDOFType &Template0,
                          const RealType Lambda )
    : _tau ( aol::ZOTrait<RealType>::one ),
      _lambda ( Lambda ),
      _pReference0 ( ),
      _pTemplate0 ( ),
      _grid ( Reference0.getSize() ) {
    setRefImageReference ( Reference0 );
    setTemplImageReference ( Template0 );
  }

  RegistrationInterface ( const qc::GridSize<ConfiguratorType::Dim> &Size,
                          const RealType Lambda )
    : _tau ( aol::ZOTrait<RealType>::one ),
      _lambda ( Lambda ),
      _pReference0 ( ),
      _pTemplate0 ( ),
      _grid ( Size ) {}

  RegistrationInterface ( const char *ReferenceFileName,
                          const char *TemplateFileName,
                          const RealType Lambda,
                          const bool DontNormalizeInputImages = false )
    : _tau ( aol::ZOTrait<RealType>::one ),
      _lambda ( Lambda ),
      _pReference0 ( new ImageDOFType ( ReferenceFileName ), true ),
      _pTemplate0 ( new ImageDOFType ( TemplateFileName ), true ),
      _grid ( _pReference0->getSize() ) {
    if ( _pReference0->getSize() != _pTemplate0->getSize() )
      throw aol::Exception ( "qc::RegistrationInterface: Sizes of reference and template image don't match.\n", __FILE__, __LINE__ );

    // The const_casts aren't really nice, but since we just manually created these arrays, it's fine to rescale the values.
    if ( DontNormalizeInputImages ) {
      // Just scale PGM/PNG input data from [0,255] to [0,1] in this case.
      if ( aol::fileNameEndsWith ( ReferenceFileName, ".png" ) || aol::fileNameEndsWith ( ReferenceFileName, ".pgm" ) )
        const_cast<ImageDOFType&> ( getRefImageReference() ) /= 255;
      if ( aol::fileNameEndsWith ( TemplateFileName, ".png" ) || aol::fileNameEndsWith ( TemplateFileName, ".pgm" ) )
        const_cast<ImageDOFType&> ( getTemplImageReference() ) /= 255;
    }
    else {
      const_cast<ImageDOFType&> ( getRefImageReference() ).scaleValuesTo01();
      const_cast<ImageDOFType&> ( getTemplImageReference() ).scaleValuesTo01();
    }
  }

  const InitType& getGridReference ( ) const {
    return _grid;
  }

  const InitType& getInitializerRef( ) const {
    return _grid;
  }

  const ImageDOFType& getTemplImageReference ( ) const {
    return *_pTemplate0;
  }

  ImageDOFType& getNonConstTemplImageReference ( ) {
    if ( _pTemplate0->getDeleteFlag() == false )
      throw aol::Exception ( "qc::RegistrationInterface: Can't return non-const pointer to template, if we only store a reference to the data.\n", __FILE__, __LINE__ );
    return *_pTemplate0;
  }

  const ImageDOFType& getRefImageReference ( ) const {
    return *_pReference0;
  }

  ImageDOFType& getNonConstRefImageReference ( ) {
    if ( _pReference0->getDeleteFlag() == false )
      throw aol::Exception ( "qc::RegistrationInterface: Can't return non-const pointer to reference, if we only store a reference to the data.\n", __FILE__, __LINE__ );
    return *_pReference0;
  }

  void setTemplImageReference ( const ImageDOFType &Template ) {
    _pTemplate0.reset ( new ImageDOFType ( Template, aol::FLAT_COPY ), true );
  }

  void setTemplImageReference ( RealType* PTemplateData ) {
    _pTemplate0.reset ( new ImageDOFType ( PTemplateData, qc::GridSize<ConfiguratorType::Dim>::createFrom ( _grid ), aol::FLAT_COPY ), true );
  }

  void setRefImageReference ( const ImageDOFType &Reference ) {
    _pReference0.reset ( new ImageDOFType ( Reference, aol::FLAT_COPY ), true );
  }

  void setRefImageReference ( RealType* PReferenceData ) {
    _pReference0.reset ( new ImageDOFType ( PReferenceData, qc::GridSize<ConfiguratorType::Dim>::createFrom ( _grid ), aol::FLAT_COPY ), true );
  }

  void setTau ( const RealType Tau ) {
    _tau = Tau;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ImageDOFType = typename ConfiguratorType::ArrayType, bool OnlySaveCheckView = false, bool SaveDisplacement = true>
class RegistrationStepSaver : public aol::StepSaverBase<typename ConfiguratorType::RealType, aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_fineGrid;
  const typename ConfiguratorType::InitType &_coarseGrid;
  const ImageDOFType &_referenceFine;
  const ImageDOFType &_referenceCoarse;
  const ImageDOFType &_templateFine;
  const int _checkboxWidth;
  bool _onlySaveDisplacement;
  bool _saveImagesAsDouble;
public:
  RegistrationStepSaver ( const typename ConfiguratorType::InitType &FineGrid,
                          const typename ConfiguratorType::InitType &CoarseGrid,
                          const ImageDOFType &ReferenceFine,
                          const ImageDOFType &ReferenceCoarse,
                          const ImageDOFType &TemplateFine,
                          const int CheckboxWidth )
   : _fineGrid ( FineGrid ),
     _coarseGrid ( CoarseGrid ),
     _referenceFine ( ReferenceFine ),
     _referenceCoarse ( ReferenceCoarse ),
     _templateFine ( TemplateFine ),
     _checkboxWidth ( CheckboxWidth ),
     _onlySaveDisplacement ( false ),
     _saveImagesAsDouble ( false ) {}

  RegistrationStepSaver ( const typename ConfiguratorType::InitType &Grid,
                          const ImageDOFType &Reference,
                          const ImageDOFType &Template,
                          const int CheckboxWidth )
   : _fineGrid ( Grid ),
     _coarseGrid ( Grid ),
     _referenceFine ( Reference ),
     _referenceCoarse ( Reference ),
     _templateFine ( Template ),
     _checkboxWidth ( CheckboxWidth ),
     _onlySaveDisplacement ( false ),
     _saveImagesAsDouble ( false ) {}

  void setOnlySaveDisplacement ( const bool OnlySaveDisplacement ) {
    _onlySaveDisplacement = OnlySaveDisplacement;
  }

  void setSaveImagesAsDouble ( const bool SaveImagesAsDouble ) {
    _saveImagesAsDouble = SaveImagesAsDouble;
  }

private:
  void doSaveStep ( const aol::MultiVector<RealType> &SaveInput, const int Iteration, const char *OverrideBaseSaveName ) const {
    if ( _onlySaveDisplacement == false ) {
      ImageDOFType templateDefFine ( _fineGrid );
      qc::deformImageWithCoarseDeformation<ConfiguratorType> ( _templateFine, _fineGrid, _coarseGrid, templateDefFine, SaveInput );

      ImageDOFType tmpFine ( _fineGrid );
      const string checkViewBaseName = this->createSaveName ( "fine_checkview", "", Iteration, OverrideBaseSaveName );
      if ( ConfiguratorType::Dim != QC_1D ) {
        qc::DataGenerator<ConfiguratorType> generator ( _fineGrid );
        generator.generateCheckView ( tmpFine, _referenceFine, templateDefFine, _checkboxWidth, true );
        qc::writeImage<RealType> ( _fineGrid, tmpFine, checkViewBaseName.c_str(), false, qc::QC_X, _saveImagesAsDouble );
      }
      else {
        aol::PlotDataFileHandler<RealType> plotHandler;
        plotHandler.generateFunctionPlot ( _referenceFine );
        plotHandler.generateFunctionPlot ( templateDefFine );
        aol::Plotter<RealType> plotter;
        plotter.addPlotCommandsFromHandler ( plotHandler );
        plotter.set_outfile_base_name ( checkViewBaseName.c_str() );
        plotter.genPNG();
      }

      if ( OnlySaveCheckView == true )
        return;

      qc::writeImage<RealType> ( _fineGrid, templateDefFine, this->createSaveName ( "templ_def_fine", "", Iteration, OverrideBaseSaveName ).c_str(), false, qc::QC_X, _saveImagesAsDouble );

      //! \todo Implement qc::TransformFunction in 1D
      if ( ConfiguratorType::Dim != QC_1D ) {
        ImageDOFType tmpCoarse ( _coarseGrid );
        qc::TransformFunction<RealType, ConfiguratorType::Dim> transform ( _coarseGrid );
        transform.setDeformation ( SaveInput );
        transform.apply ( _referenceCoarse, tmpCoarse );
        qc::writeImage<RealType> ( _coarseGrid, tmpCoarse, this->createSaveName ( "ref_def", "", Iteration, OverrideBaseSaveName ).c_str(), false, qc::QC_X, _saveImagesAsDouble );
      }
    }

    if ( SaveDisplacement ) {
      qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( _coarseGrid, SaveInput, aol::FLAT_COPY );
      phi.save ( this->createSaveName ( "deformation", "_%d.dat.bz2", Iteration, OverrideBaseSaveName ).c_str(), qc::PGM_DOUBLE_BINARY );
    }
  };
};

enum REGISTRATION_INPUT_TYPE {
  REFERENCE,
  TEMPLATE
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ImageDOFType = typename ConfiguratorType::ArrayType>
class RegistrationMultilevelDescentInterfaceBase : public qc::MultilevelDescentInterface<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename qc::MultilevelArrayTrait<RealType, InitType>::ProlongOpType ProlongOpType;
  typedef typename qc::MultilevelArrayTrait<RealType, InitType>::RestrictOpType RestrictOpType;
  typedef qc::MultilevelArray<RealType, ImageDOFType, ProlongOpType, RestrictOpType, InitType> MultilevelArrayType;
protected:
  const aol::DeleteFlagPointer<const aol::ParameterParser> _pParser;
  // original image data
  MultilevelArrayType _org_template, _org_reference;
  string _saveDirectory;

public:
  const aol::ParameterParser& getParserReference ( ) const {
    return *_pParser;
  }

  RegistrationMultilevelDescentInterfaceBase ( const aol::ParameterParser &Parser )
    : qc::MultilevelDescentInterface<ConfiguratorType> ( Parser.getInt ( "precisionLevel" ) ),
      _pParser ( &Parser, false ),
      _org_template ( this->_grid ),
      _org_reference ( this->_grid ) {
    loadImageData( );

    // This creates the save directory and dumps the paramters to a file in that dir.
    qc::DefaultArraySaver<RealType, ConfiguratorType::Dim> saver;
    saver.initFromParser ( this->getParserReference(), true );
    _saveDirectory = saver.getSaveDirectory();
  }

  RegistrationMultilevelDescentInterfaceBase ( const int MaxDepth )
    : qc::MultilevelDescentInterface<ConfiguratorType> ( MaxDepth ),
      _pParser ( new aol::ParameterParser, true ),
      _org_template ( this->_grid ),
      _org_reference ( this->_grid ) {
  }

  void setLevel ( const int Level ) {
    qc::MultilevelDescentInterface<ConfiguratorType>::setLevel ( Level );
    _org_template.setCurLevel ( Level );
    _org_reference.setCurLevel ( Level );
  }

  const MultilevelArrayType& getTemplImageMLAReference ( ) const {
    return _org_template;
  }

  MultilevelArrayType& getTemplImageMLAReference ( ) {
    return _org_template;
  }

  const MultilevelArrayType& getRefImageMLAReference ( ) const {
    return _org_reference;
  }

  MultilevelArrayType& getRefImageMLAReference ( )  {
    return _org_reference;
  }

  const ImageDOFType& getTemplImageReference ( ) const {
    return _org_template[this->_maxDepth];
  }

  const ImageDOFType& getRefImageReference ( ) const {
    return _org_reference[this->_maxDepth];
  }

  const ImageDOFType& getCurrentTemplImageReference ( ) const {
    return _org_template[this->_curLevel];
  }

  const ImageDOFType& getCurrentRefImageReference ( ) const {
    return _org_reference[ this->_curLevel];
  }

  void setTemplate ( const ImageDOFType &Template ) {
    _org_template [this->_maxDepth] = Template;
    _org_template.levRestrict ( 0, this->_maxDepth );
  }

  void setReference ( const ImageDOFType &Reference ) {
    _org_reference [this->_maxDepth] = Reference;
    _org_reference.levRestrict ( 0, this->_maxDepth );
  }

  void prepareMultilevelArray ( MultilevelArrayType &Input, MultilevelArrayType &Dest, const bool NoScaling, const bool NoSmoothing ) const {
    if ( !NoScaling ) {
      if ( getParserReference().hasVariable ( "enhanceContrastSaturationPercentage" ) ) {
        const aol::Vec2<RealType> minMax = Input.current().getSaturatedMinMaxValue ( getParserReference().getDouble ( "enhanceContrastSaturationPercentage" ) );
        Input.current().clamp ( minMax[0], minMax[1] );
        Input.current().scaleValuesTo01();
      }
      else {
        RealType minValue = Input.current().getMinValue();
        if ( ( minValue < 0. ) || getParserReference().checkAndGetBool ( "normalizeMinToZero" ) )
          Input.current().addToAll ( -minValue );
        RealType maxValue = Input.current().getMaxValue();
        if ( maxValue != 0. )
          Input.current() /= maxValue;
      }
    }

    if ( this->_curLevel < Input.getCurLevel() ) {
      Input.levRestrict ( this->_curLevel, Input.getCurLevel() );
      Input.setCurLevel ( this->_curLevel );
    }

    if ( this->_curLevel > Input.getCurLevel() )
        throw aol::Exception ( "Input data size too small. Note: Unless \"resizeInput\" is used, \"precisionLevel\" "
                               "may not be bigger than the grid level of the input data.\n", __FILE__, __LINE__ );

    // presmooth the image if the users wants
    if ( !NoSmoothing && getParserReference().hasVariable ( "preSmoothSigma" ) ) {
      const InitType &inputGrid = *(this->_grids [Input.getCurLevel()]);
      typename qc::MultilevelArrayTrait<RealType, InitType>::LinSmoothType linSmooth ( inputGrid );
      linSmooth.setSigma( getParserReference().getDouble( "preSmoothSigma" ) * inputGrid.H() );
      linSmooth.applySingle( Input.current() );
    }

    Dest.current( ) = Input.current();
    Dest.levRestrict ( 0, this->_curLevel );
  }

  static void prepareImage ( const aol::ParameterParser &Parser, const bool NoResizeOrCrop, const ImageDOFType &InputImage, ImageDOFType &PreparedImage ) {
    if ( !NoResizeOrCrop && Parser.checkAndGetBool ( "resizeInput" ) ) {
      PreparedImage.resampleFrom ( InputImage );
    }
    else if ( !NoResizeOrCrop && Parser.checkAndGetBool ( "cropInput" ) ) {
      aol::Vec<ConfiguratorType::Dim, int> cropStart;
      cropStart[0] = Parser.getInt ( "cropStartX" );
      if ( ConfiguratorType::Dim >= 2 )
        cropStart[1] = Parser.getInt ( "cropStartY" );
      if ( ConfiguratorType::Dim >= 3 )
        cropStart[2] = Parser.getInt ( "cropStartZ" );
      InputImage.copyBlockTo ( cropStart, PreparedImage );
    }
    else {
      if ( PreparedImage.getSize() != InputImage.getSize() )
        throw aol::Exception ( "Input data size doesn't match the size of the computational grid. "
                               "Either use a different type of grid or change the size with the "
                               "input data with \"resizeInput\" or \"cropInput\".\n", __FILE__, __LINE__ );

      PreparedImage = InputImage;
    }
  }

  void initMultilevelArrayFromImage ( const ImageDOFType &InputData, MultilevelArrayType &Dest, const bool NoScaling = false, const bool NoResizeOrCrop = false, const bool NoSmoothing = false ) const {
    const int approximateInputLevel = qc::logBaseTwo ( InputData.getNumX() );
    const InitType fullGrid ( ( getParserReference().checkAndGetBool ( "cropInput" ) || ( getParserReference().checkAndGetBool ( "resizeInput" ) && this->_curLevel > approximateInputLevel ) )
                              ? this->_curLevel : approximateInputLevel, ConfiguratorType::Dim );
    MultilevelArrayType fullImage ( fullGrid );
    ImageDOFType fullImageArray ( fullImage.current( ), aol::FLAT_COPY );
    prepareImage ( getParserReference(), NoResizeOrCrop, InputData, fullImageArray );
    prepareMultilevelArray ( fullImage, Dest, NoScaling, NoSmoothing );
  }

  void loadAndPrepareImage ( const char* Filename, MultilevelArrayType &Dest, const bool NoScaling = false, const bool NoResizeOrCrop = false, const bool NoSmoothing = false ) const {
    ImageDOFType inputData ( Filename );
    initMultilevelArrayFromImage ( inputData, Dest, NoScaling, NoResizeOrCrop, NoSmoothing );
  }

  void loadAndPrepareImage ( const char* Filename, ImageDOFType &Dest, const bool NoScaling, const bool NoResizeOrCrop = false, const bool NoSmoothing = false ) const {
    // Using the qc::MultilevelArray version of loadAndPrepareImage requires one additional copy of the image, but allows to use the existing code.
    MultilevelArrayType tempMultilevelArray ( this->_grid );
    loadAndPrepareImage ( Filename, tempMultilevelArray, NoScaling, NoResizeOrCrop, NoSmoothing );
    Dest = tempMultilevelArray[this->_maxDepth];
  }

  void initRefOrTemplate ( const ImageDOFType &InputData, const REGISTRATION_INPUT_TYPE Input, const bool NoScaling ) {
    initMultilevelArrayFromImage ( InputData, ( Input == REFERENCE ) ? _org_reference : _org_template, NoScaling, ( Input == REFERENCE ) && getParserReference().checkAndGetBool ( "dontResizeOrCropReference" ) );
  }

  void loadRefOrTemplate ( const char* Filename, const REGISTRATION_INPUT_TYPE Input, const bool NoScaling ) {
    loadAndPrepareImage ( Filename, ( Input == REFERENCE ) ? _org_reference : _org_template, NoScaling, ( Input == REFERENCE ) && getParserReference().checkAndGetBool ( "dontResizeOrCropReference" ) );
  }

  void loadImageData( ) {
    if ( getParserReference().hasVariable ( "reference" ) ) {
      const string refFilename = getParserReference().getStringExpandTilde ( "reference" );
      cerr << "reading reference image from: " << refFilename << endl;
      loadRefOrTemplate ( refFilename.c_str(), REFERENCE, getParserReference().checkAndGetBool ( "dontNormalizeInputImages" ) );
    }
    else
      cerr << "No reference specified\n";

    if ( getParserReference().hasVariable ( "template" ) ) {
      const string templFilename = getParserReference().getStringExpandTilde ( "template" );
      cerr << "reading template image from: " << templFilename << endl;
      loadRefOrTemplate ( templFilename.c_str(), TEMPLATE, getParserReference().checkAndGetBool ( "dontNormalizeInputImages" ) );
    }
    else
      cerr << "No template specified\n";
  }

  const char* getSaveDirectory () const {
    return _saveDirectory.c_str();
  }

  void setSaveDirectory ( const char *SaveDirectory ) {
    _saveDirectory = SaveDirectory;
  }

  void makeAndSetSaveDirectory ( const char *SaveDirectory ) {
    aol::makeDirectory ( SaveDirectory );
    setSaveDirectory ( SaveDirectory );
  }

  void solve( const int StartLevel = -1, const int StopLevel = -1 ) {
    if ( getParserReference().checkAndGetBool ( "saveRefAndTempl" ) && ( getParserReference().checkAndGetBool ( "onlySaveDisplacement" ) == false ) ) {
      qc::writeImage<RealType> ( this->_grid, getRefImageReference(), ( string ( getSaveDirectory() ) + "reference" ).c_str(), false, qc::QC_X, getParserReference().checkAndGetBool ( "saveAsDouble" ) );
      qc::writeImage<RealType> ( this->_grid, getTemplImageReference(), ( string ( getSaveDirectory() ) + "template" ).c_str(), false, qc::QC_X, getParserReference().checkAndGetBool ( "saveAsDouble" ) );
    }
    qc::MultilevelDescentInterface<ConfiguratorType>::solve( ( StartLevel < 0 ) ? getParserReference().getInt ( "startLevel" ) : StartLevel,
                                                             ( StopLevel < 0 ) ? getParserReference().getInt ( "stopLevel" ) : StopLevel );
  }

  void solveAndProlongToMaxDepth( const int StartLevel = -1, const int StopLevel = -1 ) {
    solve( StartLevel, StopLevel );
    const int missingLevel = this->_maxDepth - this->_curLevel;
    for ( int j = 0; j < missingLevel; ++j )
      this->prolongate ( );
  }
};

/**
 * \author Berkels
 */
template <typename _ConfiguratorType, typename ImageDOFType = typename _ConfiguratorType::ArrayType>
class RegistrationMultilevelDescentInterface : public qc::RegistrationMultilevelDescentInterfaceBase<_ConfiguratorType, ImageDOFType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename qc::MultilevelArrayTrait<RealType, InitType>::MultilevelArrayType MultilevelArrayType;
  typedef typename qc::MultilevelArrayTrait<RealType, InitType>::MultiDimMultilevelArrayType MultiDimMultilevelArrayType;
  typedef typename MultilevelArrayType::ArrayType ArrayType;
  typedef qc::MultiArray<RealType, ConfiguratorType::Dim> TransformationDOFType;
public:
  static const qc::Dimension _dim = ConfiguratorType::Dim;
protected:
  MultiDimMultilevelArrayType _transformation;
  const bool _disableSaving;
public:
  RegistrationMultilevelDescentInterface ( const aol::ParameterParser &Parser ) :
      qc::RegistrationMultilevelDescentInterfaceBase<ConfiguratorType, ImageDOFType> ( Parser ),
      _transformation ( this->_grid, _dim ),
      _disableSaving ( Parser.hasVariable ( "saveDirectory" ) == false ) {
    _transformation.clear();

    if ( this->getParserReference().hasVariable ( "initial-deformation" ) ) {
      qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( _transformation, aol::FLAT_COPY );
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        phi[i].load ( this->getParserReference().getString ( "initial-deformation", i ).c_str() );
      _transformation.levRestrict ( 0, this->_curLevel );
    }
  }

  virtual ~RegistrationMultilevelDescentInterface( ) {}

  void setLevel ( const int Level ) {
    qc::RegistrationMultilevelDescentInterfaceBase<ConfiguratorType, ImageDOFType>::setLevel ( Level );
    this->_transformation.setCurLevel ( Level );
  }

  void prolongate( ) {
    if ( this->_curLevel < this->getMaxGridDepth() ) {
      if ( this->getParserReference().checkAndGetBool ( "resampleInsteadOfProlongateDeformation" ) ) {
        for ( int i = 0; i < this->_transformation.numComponents(); ++i )
          this->_transformation.getArray ( i, this->_curLevel + 1 ).resampleFrom ( this->_transformation.getArray ( i, this->_curLevel ) );
        this->_transformation.setCurLevel ( this->_curLevel + 1 );
      }
      else
        this->_transformation.levProlongate( );
      setLevel ( this->_curLevel + 1 );
    }
  }

  RealType getRegularizationWeight ( const int Level ) const {
    const RealType lambdaFactor = this->getParserReference().hasVariable ( "lambdaFactor" ) ? static_cast<RealType> ( pow ( this->getParserReference().getDouble ( "lambdaFactor" ), this->getParserReference().getInt ( "stopLevel" ) - Level ) ) : 1;
    return ( lambdaFactor * static_cast<RealType> ( this->getParserReference().getDouble ( "lambda" ) ) );
  }

  void logLevelStart ( const RealType Lambda ) const {
    cerr << "\n--------------------------------------------------------\n";
    cerr << "Registration on level " << this->_curLevel << " started";
    if ( aol::appeqAbsolute ( Lambda, static_cast<RealType> ( this->getParserReference().getDouble ( "lambda" ) ) ) == false )
      cerr << " (lambda " << Lambda << ")";
    cerr << "\n";
    cerr << "--------------------------------------------------------\n\n";
  }

  const MultiDimMultilevelArrayType& getTransformationReference ( ) const {
    return _transformation;
  }

  const InitType& getTransformationDOFInitializer ( ) const {
    return this->getInitializerRef();
  }

  string getDeformationFileNameSuffix ( ) const {
    return "_%d.dat.bz2";
  }

  void setTransformation ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) {
    for ( int comp = 0; comp < ConfiguratorType::Dim; ++comp )
      _transformation.getArray( comp, this->_maxDepth ) = Transformation[comp];
    _transformation.levRestrict ( 0, this->_maxDepth );
  }

  void getTransformation ( qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) const {
    for ( int comp = 0; comp < ConfiguratorType::Dim; ++comp )
      Transformation[comp] = _transformation.getArray( comp, this->_maxDepth );
  }

  RealType getTransformationNorm ( ) const {
    RealType norm = 0;
    for ( int comp = 0; comp < ConfiguratorType::Dim; ++comp )
      norm += _transformation.getArray( comp, this->_maxDepth ).norm();
    return norm;
  }

  void addTransformationTo ( qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) const {
    for ( int comp = 0; comp < ConfiguratorType::Dim; ++comp )
      Transformation[comp] += _transformation.getArray( comp, this->_maxDepth );
  }

  void setTransformationToZero ( ) {
    for ( int comp = 0; comp < ConfiguratorType::Dim; ++comp )
      _transformation.getArray( comp, this->_maxDepth ).setZero();
    _transformation.levRestrict ( 0, this->_maxDepth );
  }

  void setTransformationToTranslation ( const aol::Vec<_dim, RealType> &Translation ) {
    for ( int comp = 0; comp < ConfiguratorType::Dim; ++comp )
      _transformation.getArray( comp, this->_maxDepth ).setAll( Translation[comp] );
    _transformation.levRestrict ( 0, this->_maxDepth );
  }

  void composeDeformations ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation1,  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation2, qc::MultiArray<RealType, ConfiguratorType::Dim> &Composition ) {
    qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType> ( Transformation1, Transformation2, this->_grid, Composition );
  }

  void setTransformationToComposition ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation1,  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation2 ) {
    _transformation.setCurLevel ( _transformation.getDepth() );
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( _transformation, aol::FLAT_COPY );
    qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType> ( Transformation1, Transformation2, this->_grid, phi );
    _transformation.levRestrict ( 0, this->_maxDepth );
  }

  void applyTransformation ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false ) const {
    qc::DeformImage<ConfiguratorType> ( InputImage, this->_grid, DeformedImage, Transformation, true, ExtensionConstant, NearestNeighborInterpolation );
  }

  void applyCurrentTransformation ( const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage ) const {
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( _transformation, aol::FLAT_COPY );
    applyTransformation ( phi, InputImage, DeformedImage );
  }

  void applySavedTransformation ( const char *DefBaseName, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false, const RealType xDerivNormThreshold = 0 ) const {
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_grid );
    loadTransformationTo ( DefBaseName, phi );
    applyTransformation ( phi, InputImage, DeformedImage, ExtensionConstant, NearestNeighborInterpolation );
    if ( xDerivNormThreshold > 0 )
      throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  void loadTransformationTo ( const char *DefBaseName, qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) const {
    Transformation.load ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str() );
  }

  void saveTransformationTo ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation, const char *DefBaseName ) const {
    Transformation.save ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str(), qc::PGM_DOUBLE_BINARY );
  }

  void saveTransformation ( const char *FileName ) const {
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( _transformation, aol::FLAT_COPY );
    phi.save ( FileName, qc::PGM_DOUBLE_BINARY );
  }

  void loadTransformation ( const char *FileName ) {
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_grid );
    phi.load ( FileName );
    setTransformation ( phi );
  }

  void writeInitialRegistration ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    if ( this->getParserReference().checkAndGetBool ( "onlySaveDisplacement" ) )
      return;

    qc::RegistrationStepSaver<ConfiguratorType, ImageDOFType, true>
      stepSaver( this->_grid, *this->_curGrid, this->_org_reference[ this->_maxDepth ], this->_org_reference[ this->_curLevel ], this->_org_template[ this->_maxDepth ], this->getParserReference().getInt ( "checkboxWidth" ) );
    stepSaver.setSaveName ( aol::strprintf ( "_before_%02d", this->_curLevel ).c_str() );
    stepSaver.setSaveDirectory ( this->getSaveDirectory() );
    stepSaver.setSaveImagesAsDouble ( this->getParserReference().checkAndGetBool ( "saveAsDouble" ) );
    stepSaver.saveStep ( Phi, -1 );
  }

  void writeCurrentRegistration ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi, const int Iteration = -1 ) const {
    qc::RegistrationStepSaver<ConfiguratorType, ImageDOFType>
      stepSaver( this->_grid, *this->_curGrid, this->_org_reference[ this->_maxDepth ], this->_org_reference[ this->_curLevel ], this->_org_template[ this->_maxDepth ], this->getParserReference().getInt ( "checkboxWidth" ) );
    string numString = aol::strprintf ( "_%02d", this->_curLevel );
    stepSaver.setSaveName ( numString.c_str() );
    stepSaver.setSaveDirectory ( this->getSaveDirectory() );
    stepSaver.setSaveTimestepOffset ( 1 );
    stepSaver.setOnlySaveDisplacement ( this->getParserReference().checkAndGetBool ( "onlySaveDisplacement" ) );
    stepSaver.setSaveImagesAsDouble ( this->getParserReference().checkAndGetBool ( "saveAsDouble" ) );
    stepSaver.saveStep ( Phi, Iteration );
  }

  void writeDeformedMatch ( const aol::Vector<RealType> &Image, const char *Filename, const aol::MultiVector<RealType> &Phidofs ) const {
    ArrayType match( this->_grid );
    qc::deformImageWithCoarseDeformation<ConfiguratorType> ( Image, this->_grid, *this->_curGrid, match, Phidofs );
    qc::writeImage<RealType> ( this->_grid, match, Filename );
  }

  //! Write a checkbox image of \f$ A \f$ and \f$ B\circ\phi \f$, \f$ A \f$ and \f$ B \f$ are images on the finest level, and \f$ \phi \f$ a transformation on the current level
  void writeCheckBox ( const char *filename, const qc::Array< RealType > &ImageA, const qc::Array< RealType > &ImageB, const aol::MultiVector<RealType> &Phidofs, const int BoxWidth, const bool SliceView = false ) const {
    ArrayType checkBoxImage ( this->_grid );
    ArrayType match ( this->_grid );

    qc::deformImageWithCoarseDeformation<ConfiguratorType> ( ImageB, this->_grid, *this->_curGrid, match, Phidofs );
    qc::DataGenerator<ConfiguratorType> generator ( this->_grid );
    generator.generateCheckView ( checkBoxImage, ImageA, match, BoxWidth, SliceView );

    qc::writeImage<RealType> ( this->_grid, checkBoxImage, filename );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class DirichletRegularizationConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::StiffOp<ConfiguratorType> _stiff;
  const qc::DisplacementLengthEnergy<ConfiguratorType> _regE;
  const aol::DiagonalBlockOp<RealType> _regDE;
public:

  DirichletRegularizationConfigurator ( const typename ConfiguratorType::InitType &Grid,
                                        // Dummy argument, don't give it a name!
                                        const aol::ParameterParser &/*Parser*/ = *static_cast<aol::ParameterParser*> ( NULL ),
                                        // Dummy argument, don't give it a name!
                                        const aol::Vector<RealType> & = *static_cast<aol::Vector<RealType>*> ( NULL ) )
   : _grid ( Grid ),
     _stiff ( _grid,  aol::ASSEMBLED ),
     _regE ( _stiff ),
     _regDE ( _stiff ) {}

  template <typename RegistrationMultilevelDescentType>
  DirichletRegularizationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
   : _grid ( RegisMLD.getCurrentGrid() ),
     _stiff ( _grid,  aol::ASSEMBLED ),
     _regE ( _stiff ),
     _regDE ( _stiff ) {}

  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &getRegERef ( ) const {
    return _regE;
  }

  const aol::Op<aol::MultiVector<RealType> > &getRegDERef ( ) const {
    return _regDE;
  }
};

/**
 * \author Berkels
 */
template < typename ConfiguratorType >
class LaplaceConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const qc::LaplaceEnergy < ConfiguratorType > _laplaceEnergy;
  const aol::DerivativeWrapper < RealType, qc::LaplaceEnergy < ConfiguratorType >, aol::MultiVector < RealType > > _laplaceEnergyDerivative;
public:

  LaplaceConfigurator ( const typename ConfiguratorType::InitType &Grid,
                       const aol::ParameterParser &/*Parser*/ )
    : _grid ( Grid ),
      _laplaceEnergy ( _grid ),
      _laplaceEnergyDerivative ( _laplaceEnergy ) {}

  template < typename RegistrationMultilevelDescentType >
  LaplaceConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
    : _grid ( RegisMLD.getCurrentGrid () ),
      _laplaceEnergy ( _grid ),
      _laplaceEnergyDerivative ( _laplaceEnergy ) {}

  const aol::Op < aol::MultiVector < RealType >, aol::Scalar < RealType > > &getRegERef () const {
    return _laplaceEnergy;
  }

  const aol::Op < aol::MultiVector < RealType > > &getRegDERef () const {
    return _laplaceEnergyDerivative;
  }
};

/**
 * linear combination of Dirichlet and Laplace energy (Laplace energy term weighted with _laplaceWeight)
 *
 * \author Effland
 */
template < typename ConfiguratorType >
class DirichletLaplaceConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const RealType _laplaceWeight; //coefficient of Laplace term
  const qc::DirichletRegularizationConfigurator < ConfiguratorType > _dirichletRegConf;
  const qc::LaplaceEnergy < ConfiguratorType > _laplaceEnergy;
  const aol::DerivativeWrapper < RealType, qc::LaplaceEnergy < ConfiguratorType >, aol::MultiVector < RealType > > _laplaceEnergyDerivative;
  aol::LinCombOp < aol::MultiVector < RealType >, aol::Scalar < RealType > > _E;
  aol::LinCombOp < aol::MultiVector < RealType > > _DE;

  void init ( ) {
    _E.appendReference ( _dirichletRegConf.getRegERef () );
    _E.appendReference ( _laplaceEnergy, _laplaceWeight );
    _DE.appendReference ( _dirichletRegConf.getRegDERef () );
    _DE.appendReference ( _laplaceEnergyDerivative, _laplaceWeight );
  }
public:

  DirichletLaplaceConfigurator ( const typename ConfiguratorType::InitType &Grid,
                                 const aol::ParameterParser &Parser )
  : _grid ( Grid ),
    _laplaceWeight ( Parser.getDouble ( "laplaceWeight" ) ),
    _dirichletRegConf ( _grid ),
    _laplaceEnergy ( _grid ),
    _laplaceEnergyDerivative ( _laplaceEnergy ) {
      init();
  }

  template < typename RegistrationMultilevelDescentType >
  DirichletLaplaceConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
  : _grid ( RegisMLD.getCurrentGrid () ),
    _laplaceWeight ( RegisMLD.getParserReference ().getDouble ( "laplaceWeight" ) ),
    _dirichletRegConf ( RegisMLD ),
    _laplaceEnergy ( _grid ),
    _laplaceEnergyDerivative ( _laplaceEnergy ) {
      init();
  }

  const aol::Op < aol::MultiVector < RealType >, aol::Scalar < RealType > > &getRegERef () const {
    return _E;
  }

  const aol::Op < aol::MultiVector < RealType > > &getRegDERef () const {
    return _DE;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class ROFRegularizationConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  const RealType _delta;
  const aol::Vec<ConfiguratorType::Dim, RealType> _compScaling;
  const aol::IsoEnergyOp<ConfiguratorType> _tvE;
  const aol::VariationOfIsoEnergyOp<ConfiguratorType> _tvDE;
  const aol::ScaledScalarVecToScalarMVecOp<RealType, ConfiguratorType::Dim> _regE;
  const aol::ScaledVecToMVecOp<RealType, ConfiguratorType::Dim> _regDE;
public:

  template <typename RegistrationMultilevelDescentType>
  ROFRegularizationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
   : _grid ( RegisMLD.getCurrentGrid() ),
     _delta ( RegisMLD.getParserReference().getDouble ( "delta" ) ),
     _compScaling ( RegisMLD.getParserReference().hasVariable ( "regComponentScaling" )
                    // This cast does nothing except for working around a GCC bug.
                    ? static_cast<const aol::ParameterParser>(RegisMLD.getParserReference()).getRealVec<ConfiguratorType::Dim, RealType> ( "regComponentScaling" )
                    : aol::Vec<ConfiguratorType::Dim, RealType> ( aol::ZOTrait<RealType>::one ) ),
     _tvE ( _grid, _delta ),
     _tvDE ( _grid, _delta ),
     _regE ( _tvE, _compScaling ),
     _regDE ( _tvDE, _compScaling ) {}

  ROFRegularizationConfigurator ( const typename ConfiguratorType::InitType &Grid, const aol::ParameterParser &Parser )
    : _grid ( Grid ),
      _delta ( Parser.getDouble ( "delta" ) ),
      _compScaling ( Parser.hasVariable ( "regComponentScaling" )
                     ? Parser.getRealVec<ConfiguratorType::Dim, RealType> ( "regComponentScaling" )
                     : aol::Vec<ConfiguratorType::Dim, RealType> ( aol::ZOTrait<RealType>::one ) ),
      _tvE ( _grid, _delta ),
      _tvDE ( _grid, _delta ),
      _regE ( _tvE, _compScaling ),
      _regDE ( _tvDE, _compScaling ) {}

  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &getRegERef ( ) const {
    return _regE;
  }

  const aol::Op<aol::MultiVector<RealType> > &getRegDERef ( ) const {
    return _regDE;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class ImageDrivenAnisoTVVecNormRegularizationConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  typedef qc::PeronaMalikWeightingFunction<RealType> PeronaMalikFunctionType;
  PeronaMalikFunctionType _peronaMalikFunction;
  qc::ImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> _aTVVecNorm;
  qc::VariationOfImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> _aTVVecNormVar;
public:

  ImageDrivenAnisoTVVecNormRegularizationConfigurator ( const typename ConfiguratorType::InitType &Grid, const aol::ParameterParser &Parser, const aol::Vector<RealType> &ReferenceImage )
   : _grid ( Grid ),
     _peronaMalikFunction( 1, aol::ZOTrait<RealType>::one/aol::Sqr ( Parser.getDouble ( "mu" ) ) ),
     _aTVVecNorm( Grid, ReferenceImage, _peronaMalikFunction,  Parser.getDouble ( "nu" ), Parser.getDouble ( "delta" ) ),
     _aTVVecNormVar( Grid, ReferenceImage, _peronaMalikFunction,  Parser.getDouble ( "nu" ), Parser.getDouble ( "delta" ) ) { }

  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &getRegERef ( ) const {
    return _aTVVecNorm;
  }

  const aol::Op<aol::MultiVector<RealType> > &getRegDERef ( ) const {
    return _aTVVecNormVar;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class HyperelasticEnergyConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  typedef qc::HyperelasticEnergyDensityDefault<ConfiguratorType> HyperelasticEnergyDensityType;
  const HyperelasticEnergyDensityType _density;
  const qc::HyperelasticEnergy<ConfiguratorType, HyperelasticEnergyDensityType> _regE;
  const qc::HyperelasticGradient<ConfiguratorType, HyperelasticEnergyDensityType> _regDE;
public:

  template <typename RegistrationMultilevelDescentType>
  HyperelasticEnergyConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
   : _grid ( RegisMLD.getCurrentGrid() ),
     _density ( RegisMLD.getParserReference().getDouble ( "weightLengthEnergy" ),
                RegisMLD.getParserReference().getDouble ( "weightSurfaceEnergy" ),
                RegisMLD.getParserReference().getDouble ( "weightVolumeEnergy" ) ),
     _regE ( _grid, _density ),
     _regDE ( _grid, _density ) {}

  HyperelasticEnergyConfigurator ( const typename ConfiguratorType::InitType &Grid,
                                   const RealType WeightLengthEnergy,
                                   const RealType WeightSurfaceEnergy,
                                   const RealType WeightVolumeEnergy )
   : _grid ( Grid ),
     _density ( WeightLengthEnergy, WeightSurfaceEnergy, WeightVolumeEnergy ),
     _regE ( _grid, _density ),
     _regDE ( _grid, _density ) {}

  HyperelasticEnergyConfigurator ( const typename ConfiguratorType::InitType &Grid, const aol::ParameterParser &Parser )
   : _grid ( Grid ),
     _density ( Parser.getDouble ( "weightLengthEnergy" ),
                Parser.getDouble ( "weightSurfaceEnergy" ),
                Parser.getDouble ( "weightVolumeEnergy" ) ),
     _regE ( _grid, _density ),
     _regDE ( _grid, _density ) {}

  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &getRegERef ( ) const {
    return _regE;
  }

  const aol::Op<aol::MultiVector<RealType> > &getRegDERef ( ) const {
    return _regDE;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ImageDOFType, typename RegistrationConfiguratorType, typename RegularizationConfiguratorType = DirichletRegularizationConfigurator<ConfiguratorType>, typename GradientDescentType = aol::H1GradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, typename qc::MultilevelArrayTrait<typename ConfiguratorType::RealType, typename ConfiguratorType::InitType>::LinSmoothType> >
class StandardRegistration : public qc::RegistrationInterface<ConfiguratorType, ImageDOFType> {
protected:
  enum MINIMIZATION_ALGORITHM {
    GRADIENT_DESCENT,
    QUASI_NEWTON
  };
  typedef typename ConfiguratorType::RealType RealType;

  const RegistrationConfiguratorType &_regisConfig;
  const RegularizationConfiguratorType &_regulConfig;

  RealType _stopEpsilon;
  int _maxGDIterations;
  bool _validateDerivative;
  MINIMIZATION_ALGORITHM _minimizationAlgo;
public:
  StandardRegistration ( const ImageDOFType &Reference0,
                         const ImageDOFType &Template0,
                         const RegistrationConfiguratorType &RegisConfig,
                         const RegularizationConfiguratorType &RegulConfig,
                         const RealType Lambda )
    : qc::RegistrationInterface<ConfiguratorType, ImageDOFType> ( Reference0, Template0, Lambda ),
      _regisConfig ( RegisConfig ),
      _regulConfig ( RegulConfig ),
      _stopEpsilon ( aol::ZOTrait<RealType>::zero ),
      _maxGDIterations ( 1000 ),
      _validateDerivative ( false ),
      _minimizationAlgo ( GRADIENT_DESCENT ) {}

  StandardRegistration ( const qc::GridSize<ConfiguratorType::Dim> &Size,
                         const RegistrationConfiguratorType &RegisConfig,
                         const RegularizationConfiguratorType &RegulConfig,
                         const RealType Lambda )
    : qc::RegistrationInterface<ConfiguratorType, ImageDOFType> ( Size, Lambda ),
      _regisConfig ( RegisConfig ),
      _regulConfig ( RegulConfig ),
      _stopEpsilon ( aol::ZOTrait<RealType>::zero ),
      _maxGDIterations ( 1000 ),
      _validateDerivative ( false ),
      _minimizationAlgo ( GRADIENT_DESCENT ) {}

  RealType findTransformation ( aol::MultiVector<RealType> &Phi, const bool NoConsoleOutput = false, const char *EnergyPlotFile = NULL ) {
    // Check if the input data fulfills the requirements imposed by the chosen registration approach.
    _regisConfig.checkInput ( this->getRefImageReference(), this->getTemplImageReference() );

    aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;

    typename RegistrationConfiguratorType::Energy registrationEnergy ( this->_grid, this->getRefImageReference(), this->getTemplImageReference(), _regisConfig );

    E.appendReference ( registrationEnergy );
    E.appendReference ( _regulConfig.getRegERef(), this->_lambda );

    aol::LinCombOp<aol::MultiVector<RealType> > DE;
    typename RegistrationConfiguratorType::EnergyVariation variationOfRegisEnergy ( this->_grid, this->getRefImageReference(), this->getTemplImageReference(), _regisConfig );

    DE.appendReference ( variationOfRegisEnergy );
    DE.appendReference ( _regulConfig.getRegDERef(), this->_lambda );

    if ( _validateDerivative ) {
      aol::MultiVector<RealType> mtmp ( this->_dim, this->_grid.getNumberOfNodes() );
      aol::FirstDerivativeValidator<aol::MultiVector<RealType> > tester ( E, DE, this->_grid.H(), aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.0001 );
      tester.testDirection ( mtmp, aol::strprintf ( "test/test%06d", this->_grid.getNumberOfNodes() ).c_str() );
    }

    // Return the energy of our solution. May be useful to make a decision on what to do with the solution.
    switch ( _minimizationAlgo ) {
      case GRADIENT_DESCENT:
        return updateTransformationWithGD ( E, DE, Phi, NoConsoleOutput, EnergyPlotFile );
      case QUASI_NEWTON:
        return updateTransformationWithQuasiNewton ( E, DE, Phi, NoConsoleOutput, EnergyPlotFile );
      default:
        throw aol::Exception ( "Unknown minimization mode!", __FILE__, __LINE__ );
        return 0;
    }
  }

private:
  RealType updateTransformationWithGD ( const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                        const aol::Op<aol::MultiVector<RealType> > &DE,
                                        aol::MultiVector<RealType> &Phi, const bool NoConsoleOutput, const char *EnergyPlotFile ) {
    //aol::SimpleGradientDescent<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( this->_grid, DE, 100);
    typedef aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::MultiVector<RealType>, GradientDescentType> GDType;
    GDType gradientDescent_solver ( this->_grid, E, DE, _maxGDIterations, this->_tau, _stopEpsilon );
    gradientDescent_solver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG | ( NoConsoleOutput ? GDType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );

    std::ofstream out;
    if ( EnergyPlotFile ) {
      out.open ( EnergyPlotFile );
      gradientDescent_solver.setOutStream ( out );
    }

    gradientDescent_solver.applySingle ( Phi );
    this->_tau = gradientDescent_solver.getStartTau();

    return gradientDescent_solver.getEnergyAtLastPosition();
  }

  RealType updateTransformationWithQuasiNewton ( const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                const aol::Op<aol::MultiVector<RealType> > &DE,
                                                aol::MultiVector<RealType> &Phi, const bool NoConsoleOutput, const char * /*EnergyPlotFile*/ ) {

    typedef aol::QuasiNewtonBFGS<RealType, aol::MultiVector<RealType>, aol::MultiVector<RealType> > DescentType;
    DescentType descentSolver ( E, DE, _maxGDIterations, _stopEpsilon );
    if ( NoConsoleOutput )
      descentSolver.getNewtonInfo().setMegaQuietMode();
    descentSolver.setEnergyDecayBasedStopping ( true );
    descentSolver.applySingle ( Phi );
    return descentSolver.getNewtonInfo().getFinalResidual();
  }

public:

  void setValidateDerivative ( const bool ValidateDerivative ) {
    _validateDerivative = ValidateDerivative;
  }

  void setStopEpsilon ( const RealType StopEpsilon ) {
    _stopEpsilon = StopEpsilon;
  }

  void setMaxGDIterations ( const int MaxGDIterations ) {
    _maxGDIterations = MaxGDIterations;
  }

  void setMinimizationAlgo ( const string &MinimizationAlgo ) {
    if ( MinimizationAlgo == "GradientDescent" )
      _minimizationAlgo = GRADIENT_DESCENT;
    else if ( MinimizationAlgo == "QuasiNewton" )
      _minimizationAlgo = QUASI_NEWTON;
    else
      throw aol::Exception ( "Unknown minimization mode!", __FILE__, __LINE__ );
  }

  void setSettingsFromParser ( const aol::ParameterParser &Parser ) {
    if ( Parser.hasVariable ( "stopEpsilon" ) )
      setStopEpsilon( Parser.getDouble ( "stopEpsilon" ) );
    if ( Parser.hasVariable ( "maxGDIterations" ) )
      setMaxGDIterations( Parser.getInt ( "maxGDIterations" ) );
    if ( Parser.hasVariable ( "minimizationAlgo" ) )
      setMinimizationAlgo ( Parser.getString ( "minimizationAlgo" ) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename RegistrationConfiguratorType, typename _RegularizationConfiguratorType = DirichletRegularizationConfigurator<ConfiguratorType>, typename GradientDescentType = aol::H1GradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, typename qc::MultilevelArrayTrait<typename ConfiguratorType::RealType, typename ConfiguratorType::InitType>::LinSmoothType > >
class StandardRegistrationMultilevelDescent : public qc::RegistrationMultilevelDescentInterface<ConfiguratorType, typename RegistrationConfiguratorType::ImageDOFType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef _RegularizationConfiguratorType RegularizationConfiguratorType;
  typedef typename RegistrationConfiguratorType::ImageDOFType ImageDOFType;
  static const bool IsParametric = false;

  RealType _energyOfLastSolution;
  const bool _validateDerivative;
public:
  StandardRegistrationMultilevelDescent ( const aol::ParameterParser &Parser, const bool ValidateDerivative = false ) :
      qc::RegistrationMultilevelDescentInterface<ConfiguratorType, ImageDOFType> ( Parser ),
      _energyOfLastSolution ( aol::NumberTrait<RealType>::NaN ),
      _validateDerivative ( ValidateDerivative ) {}

  virtual ~StandardRegistrationMultilevelDescent( ) {}

  void descentOnCurrentGrid() {
    // make a multivector, referencing on the array in the multidimmultilevelarray _sol, which contains the deformation
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_transformation, aol::FLAT_COPY );

    const RegistrationConfiguratorType regisConfig ( *this );
    const RegularizationConfiguratorType regulConfig ( *this );
    const RealType lambda = this->getRegularizationWeight ( this->_curLevel );

    this->logLevelStart ( lambda );

    if ( this->_disableSaving == false )
      regisConfig.writeInitialRegistration( *this, phi );

    StandardRegistration<ConfiguratorType, ImageDOFType, RegistrationConfiguratorType, RegularizationConfiguratorType, GradientDescentType> stdRegistration ( this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ], regisConfig, regulConfig, lambda );
    stdRegistration.setValidateDerivative ( _validateDerivative );
    stdRegistration.setSettingsFromParser ( this->getParserReference() );
    _energyOfLastSolution = stdRegistration.findTransformation ( phi, false, this->_disableSaving ? NULL : aol::strprintf ( "%senergy_%02d.txt", this->getSaveDirectory(), this->_curLevel ).c_str() );

    if ( this->_disableSaving == false )
      regisConfig.writeCurrentRegistration( *this, phi );
  }

  RealType getEnergyOfLastSolution ( ) const {
    return _energyOfLastSolution;
  }

  aol::Vec2<RealType> computeEnergyComponents ( const aol::MultiVector<RealType> &Transformation ) const {
    aol::Vec2<RealType> result;
    const RegistrationConfiguratorType regisConfig ( *this );
    const RegularizationConfiguratorType regulConfig ( *this );
    typename aol::Scalar < typename ConfiguratorType::RealType > tempScalar;

    typename RegistrationConfiguratorType::Energy registrationEnergy ( this->getInitializerRef(), this->getRefImageReference(), this->getTemplImageReference(), regisConfig );
    registrationEnergy.apply ( Transformation, tempScalar );
    result[0] = tempScalar[0];

    regulConfig.getRegERef ().apply ( Transformation, tempScalar );
    result[1] = this->getRegularizationWeight ( this->getMaxGridDepth() ) * tempScalar[0];
    return result;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class BaseRegistrationConfigurator {
public:
  typedef typename ConfiguratorType::RealType RealType;

  BaseRegistrationConfigurator ( ) {}

  BaseRegistrationConfigurator ( const aol::ParameterParser &/*Parser*/ ) {}

  template <typename RegistrationMultilevelDescentType>
  void writeInitialRegistration ( const RegistrationMultilevelDescentType &RegisMLD,
                                  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    RegisMLD.writeInitialRegistration( Phi );
  }

  template <typename RegistrationMultilevelDescentType>
  void writeCurrentRegistration ( const RegistrationMultilevelDescentType &RegisMLD,
                                  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    RegisMLD.writeCurrentRegistration ( Phi );
  }

  template <typename ImageDOFType>
  void checkInput ( const ImageDOFType &/*ReferenceImage*/, const ImageDOFType &/*TemplateImage*/ ) const {}
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename GNFunctionType>
class StandardRegistrationMultilevelDescentGN : public qc::RegistrationMultilevelDescentInterface<ConfiguratorType, typename ConfiguratorType::ArrayType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef ArrayType ImageDOFType;
  static const bool IsParametric = false;

  RealType _energyOfLastSolution;
public:
  StandardRegistrationMultilevelDescentGN ( const aol::ParameterParser &Parser )
    : qc::RegistrationMultilevelDescentInterface<ConfiguratorType, ImageDOFType> ( Parser ),
      _energyOfLastSolution ( aol::NumberTrait<RealType>::NaN ) {}

  virtual ~StandardRegistrationMultilevelDescentGN( ) {}

  void descentOnCurrentGrid() {
    // make a multivector, referencing on the array in the multidimmultilevelarray _sol, which contains the deformation
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_transformation, aol::FLAT_COPY );

    const RealType lambda = this->getRegularizationWeight ( this->_curLevel );
    qc::BaseRegistrationConfigurator<ConfiguratorType> regisConfig;

    this->logLevelStart ( lambda );

    if ( this->_disableSaving == false )
      regisConfig.writeInitialRegistration( *this, phi );

    const GNFunctionType F ( this->getCurrentGrid(), this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ], lambda );

    typedef aol::SparseMatrix<RealType> MatType;
    typedef aol::SparseBlockMatrix<MatType> BlockMatType;
    aol::DerivativeWrapper<RealType, GNFunctionType, aol::MultiVector<RealType>, BlockMatType> DF ( F );

    aol::Vector<int> dimDomF;
    phi.getSizes ( dimDomF );
    aol::GaussNewtonAlgorithm<aol::MultiVector<RealType>, BlockMatType,       aol::GaussNewtonSparseNormalEquationsBlockSolver<RealType> > gnSolver ( F.getDimRangeF(), dimDomF, F, DF, this->getParserReference().getInt ( "maxGDIterations" ), this->getParserReference().getDouble ( "stopEpsilon" ) );
    gnSolver.applySingle ( phi );
    _energyOfLastSolution = gnSolver.getFNormSqrAtLastPosition();

    if ( false ) {
      aol::SquaredNormOp<aol::MultiVector<RealType> > E ( F, F.getDimRangeF() );
      aol::SquaredNormOpGradient<aol::MultiVector<RealType>, BlockMatType> DE ( F, DF, F.getDimRangeF(), F.getDimRangeF().size(), dimDomF.size() );
      aol::FirstDerivativeValidator<aol::MultiVector<RealType> >
      tester(E, DE, this->getCurrentGrid().H(), aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.0001);
      const string testBasename = aol::strprintf ( "test/test%06d", this->getCurrentGrid().getNumberOfNodes() );
      if ( aol::directoryExists ( "test" ) == false )
        aol::makeDirectory ( "test" );
      tester.testDirection ( phi, testBasename.c_str() );
      tester.setSkippingThreshold ( 1e-2 );
      tester.testAllDirections ( phi, testBasename.c_str() );
    }

    if ( this->_disableSaving == false )
      regisConfig.writeCurrentRegistration( *this, phi );
  }

  RealType getEnergyOfLastSolution ( ) const {
    return _energyOfLastSolution;
  }
};

/**
 * Converts a registration operator for gray scale images to color images by applying and summing the
 * gray scale op for each color channel.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename RegisEnergyType, typename OpDestType>
class RegistrationOpWrapper : public aol::LinCombOp<aol::MultiVector<typename ConfiguratorType::RealType>, OpDestType> {
  typedef typename ConfiguratorType::RealType RealType;
  aol::RandomAccessContainer<const RegisEnergyType> _regisOps;

public:
  template <typename RegistrationConfiguratorType>
  RegistrationOpWrapper ( const typename ConfiguratorType::InitType &Grid,
                          const aol::MultiVector<RealType> &ImR,
                          const aol::MultiVector<RealType> &ImT,
                          const RegistrationConfiguratorType &RegisConfig )  {
    for ( int i = 0; i < ImR.numComponents(); ++i ) {
      _regisOps.constructDatumAndPushBack ( Grid, ImR[i], ImT[i], RegisConfig );
      this->appendReference ( _regisOps[i] );
    }
  }
};

template <typename ConfiguratorType> class SSDEnergy;
template <typename ConfiguratorType> class SSDForce;

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class SSDRegistrationConfigurator : public BaseRegistrationConfigurator<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ImageDOFType;

  SSDRegistrationConfigurator ( ) : BaseRegistrationConfigurator<ConfiguratorType> ( ) {}

  explicit SSDRegistrationConfigurator ( const aol::ParameterParser &Parser ) : BaseRegistrationConfigurator<ConfiguratorType> ( Parser ) {}

  template <typename RegistrationMultilevelDescentType>
  explicit SSDRegistrationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
    : BaseRegistrationConfigurator<ConfiguratorType> ( RegisMLD.getParserReference() ) {}

  typedef SSDEnergy<ConfiguratorType> Energy;
  typedef SSDForce<ConfiguratorType> EnergyVariation;
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class RGBSSDRegistrationConfigurator : public BaseRegistrationConfigurator<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::MultiArray<RealType, ConfiguratorType::Dim, 3> ImageDOFType;

  explicit RGBSSDRegistrationConfigurator ( const aol::ParameterParser &Parser ) : BaseRegistrationConfigurator<ConfiguratorType> ( Parser ) {}

  template <typename RegistrationMultilevelDescentType>
  explicit RGBSSDRegistrationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
    : BaseRegistrationConfigurator<ConfiguratorType> ( RegisMLD.getParserReference() ) {}

  typedef RegistrationOpWrapper<ConfiguratorType, SSDEnergy<ConfiguratorType>, aol::Scalar<RealType> > Energy;
  typedef RegistrationOpWrapper<ConfiguratorType, SSDForce<ConfiguratorType>, aol::MultiVector<RealType> > EnergyVariation;
};

/**
 * \author Effland
 */
template < typename ConfiguratorType, typename ImageStorageType >
class RGBSSDAdditionalInformationRegistrationConfigurator : public qc::BaseRegistrationConfigurator < ConfiguratorType > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef ImageStorageType ImageDOFType;

  explicit RGBSSDAdditionalInformationRegistrationConfigurator ( const aol::ParameterParser &Parser )
  : qc::BaseRegistrationConfigurator < ConfiguratorType > ( Parser ) {
  }

  template < typename RegistrationMultilevelDescentType >
  explicit RGBSSDAdditionalInformationRegistrationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
  : qc::BaseRegistrationConfigurator < ConfiguratorType > ( RegisMLD.getParserReference () ) {
  }

  typedef qc::RegistrationOpWrapper < ConfiguratorType, qc::SSDEnergy < ConfiguratorType >, aol::Scalar < RealType > > Energy;
  typedef qc::RegistrationOpWrapper < ConfiguratorType, qc::SSDForce < ConfiguratorType >, aol::MultiVector < RealType > > EnergyVariation;
};

/**
  * \brief Implements the function \f$ f(u,v)=\frac12(u-v)^2 \f$ and its derivatives for use as fidelity term in image registration.
  *
  * This fidelity term is appropriate for Gaussian noise, where \f$ u \f$ is the noisefree and \f$ v \f$ the noisy image.
  *
  * \author Wirth
  * \ingroup AnalyticFunctions
  */
template <typename RealType>
class FidelitySSD{
public:
  FidelitySSD(){}

  //! returns \f$(u-v)^2/2 \f$
  inline RealType evaluate ( const RealType u, const RealType v ) const {
    return aol::Sqr( u - v ) / 2.0;
  }

  //! returns \f$u-v\f$
  inline RealType evaluateUDerivative ( const RealType u, const RealType v ) const {
    return ( u - v );
  }

  //! returns \f$(u-v,v-u)\f$
  inline aol::Vec2<RealType> evaluateDerivative ( const RealType u, const RealType v ) const {
    return aol::Vec2<RealType>( evaluateUDerivative ( u, v ), v - u );
  }

  //! returns I
  inline aol::Matrix22<RealType> evaluateSecondDerivative ( const RealType /*u*/, const RealType /*v*/ ) const {
    return aol::Matrix22<RealType>( 1, 0, 0, 1 );
  }

  //! returns 1
  inline RealType evaluateUSecondDerivative ( const RealType /*u*/, const RealType /*v*/ ) const {
    return 1;
  }

  static RealType maxUFactor ( const RealType UMax ) {
    return aol::Sqr ( UMax );
  }
};

/**
  * \brief Implements the function \f$ f(u,v)=u-v\log u \f$ and its derivatives for use as fidelity term in image registration.
  *
  * This fidelity term is appropriate for Poisson noise, where \f$ u \f$ is the noisefree and \f$ v \f$ the noisy image.
  *
  * \author Wirth
  * \ingroup AnalyticFunctions
  */
template <typename RealType>
class FidelityPoissonNoise{
public:
  FidelityPoissonNoise(){}

  //! returns \f$ u-v\log u \f$
  inline RealType evaluate ( const RealType u, const RealType v ) const {
    return ( u - v * log ( u ));
  }

  //! returns \f$1-v/u\f$
  inline RealType evaluateUDerivative ( const RealType u, const RealType v ) const {
    return ( 1 - v / u );
  }

  //! returns \f$(1-v/u,-\log u)\f$
  inline aol::Vec2<RealType> evaluateDerivative ( const RealType u, const RealType v ) const {
    return aol::Vec2<RealType>( evaluateUDerivative ( u, v ), - log ( u ) );
  }

  //! returns \f$ (v/u^2,-1/u;-1/u,0) \f$
  inline aol::Matrix22<RealType> evaluateSecondDerivative ( const RealType u, const RealType v ) const {
    return aol::Matrix22<RealType>( v/aol::Sqr(u), -1/u, -1/u, 0 );
  }

  //! returns \f$ v/u^2 \f$
  inline RealType evaluateUSecondDerivative ( const RealType u, const RealType v ) const {
    return v/aol::Sqr(u);
  }

  static RealType maxUFactor ( const RealType UMax ) {
    return UMax;
  }
};

/**
 * \brief Implements the function \f$ f(u,v)= u\log u+v\log v -(u+v)\log \left(\frac{u+v}{2}\right) \f$ and its derivatives.
 *
 * This function is appropriate as fidelity term for Poisson noise.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class FidelityPoissonMaximumLikelihood {
  inline RealType evalXLogX ( const RealType X ) const {
    return ( X > 0 ) ? X * log ( X ) : 0;
  }
public:
  FidelityPoissonMaximumLikelihood(){}

  //! returns \f$ u\log u+v\log v -(u+v)\log \left(\frac{u+v}{2}\right) \f$
  inline RealType evaluate ( const RealType u, const RealType v ) const {
    return evalXLogX ( u ) + evalXLogX ( v ) - 2 * evalXLogX ( 0.5 * ( u + v ) );
  }

  //! returns \f$log \left( \frac{2u}{u+v} \right)\f$
  inline RealType evaluateUDerivative ( const RealType u, const RealType v ) const {
    return log ( 2 * u / ( u + v ) );
  }

  //! returns \f$ \frac{v/u}{u+v} \f$
  inline RealType evaluateUSecondDerivative ( const RealType u, const RealType v ) const {
    return v/u/(u+v);
  }

  static RealType maxUFactor ( const RealType UMax ) {
    return UMax;
  }
};

/**
 * \brief Implements the function \f$ f(u,v)=-uv \f$ and its derivatives for use as fidelity term in image registration.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class FidelityCC{
public:
  FidelityCC(){}

  //! returns \f$-uv\f$
  inline RealType evaluate ( const RealType u, const RealType v ) const {
    return (- u * v );
  }

  //! returns \f$-v\f$
  inline RealType evaluateUDerivative ( const RealType u, const RealType v ) const {
    return ( - v );
  }

  //! returns \f$(-v,-u)\f$
  inline aol::Vec2<RealType> evaluateDerivative ( const RealType u, const RealType v ) const {
    return aol::Vec2<RealType>( evaluateUDerivative ( u, v ), -u );
  }

  inline aol::Matrix22<RealType> evaluateSecondDerivative ( const RealType /*u*/, const RealType /*v*/ ) const {
    return aol::Matrix22<RealType>( 0, -1, -1, 0 );
  }

  inline RealType evaluateUSecondDerivative ( const RealType /*u*/, const RealType /*v*/ ) const {
    return 0;
  }
};

/**
 * Data term for registration of a template image \f$ T \f$ to a reference image \f$ R \f$ via a deformation \f$ \Phi \f$.
 * Implements \f$ \int f((T\circ\Phi)(x),R(x)) dx \f$,
 * where \f$ \Phi=id+u \f$ is a deformation, the displacement \f$ u \f$ is given in the constructor,
 * \f$ R \f$ is given in the constructor, and \f$ T \f$ is given as argument in apply(Add).
 * \f$ f \f$ is passed as template argument (e.g. FidelitySSD or FidelityPoissonNoise).
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename FidelityType>
class RegDataTermOfImage : public qc::FENonlinDeformIntegrationInterface<ConfiguratorType, RegDataTermOfImage<ConfiguratorType, FidelityType> > {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _reference;
  const FidelityType &_fidelity;
public:

  RegDataTermOfImage( const typename ConfiguratorType::InitType &Grid,
                      const aol::Vector<RealType> &Reference,
                      const aol::MultiVector<RealType> &Displacement,
                      const FidelityType &Fidelity ) :
  qc::FENonlinDeformIntegrationInterface<ConfiguratorType, RegDataTermOfImage<ConfiguratorType, FidelityType> > ( Grid, Displacement ),
  _reference( Grid, Reference ),
  _fidelity ( Fidelity ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    return _fidelity.evaluate( DiscFunc.evaluate( TransformedEl, TransformedLocalCoord ), _reference.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

/**
 * Implements variation of \f$ \int f((T\circ\Phi)(x),R(x)) dx \f$ with respect to \f$ T \f$,
 * where the functional is explained in the class RegDataTermOfImage.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename FidelityType>
class RegDataTermOfImageVariation : public qc::FENonlinDeformOpInterface < ConfiguratorType, RegDataTermOfImageVariation<ConfiguratorType, FidelityType> > {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _reference;
  const FidelityType &_fidelity;
public:

  RegDataTermOfImageVariation ( const typename ConfiguratorType::InitType &Grid,
                                const aol::Vector<RealType> &Reference,
                                const aol::MultiVector<RealType> &Displacement,
                                const FidelityType &Fidelity ) :
  qc::FENonlinDeformOpInterface < ConfiguratorType, RegDataTermOfImageVariation<ConfiguratorType, FidelityType> > ( Grid, Displacement ),
    _reference( Grid, Reference ),
    _fidelity ( Fidelity ) {}

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         typename ConfiguratorType::RealType &NL ) const {
    NL = _fidelity.evaluateUDerivative( DiscFunc.evaluate ( TransformedEl, TransformedLocalCoord ), _reference.evaluateAtQuadPoint ( El, QuadPoint ) );
  }
};

/**
 * Implements variation of \f$ \int f((T\circ\Phi)(x),R(x)) dx \f$ with respect to \f$ \Phi \f$,
 * where the functional is explained in the class RegDataTermOfImage.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename FidelityType>
class RegDataTermOfDisplacementVariation : public aol::FENonlinVectorOpInterface < ConfiguratorType, 1, ConfiguratorType::Dim, RegDataTermOfDisplacementVariation<ConfiguratorType, FidelityType> > {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _reference;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
  const FidelityType &_fidelity;
public:

  RegDataTermOfDisplacementVariation( const typename ConfiguratorType::InitType &Grid,
                                      const aol::Vector<RealType> &Reference,
                                      const aol::MultiVector<RealType> &Displacement,
                                      const FidelityType &Fidelity ) :
  aol::FENonlinVectorOpInterface < ConfiguratorType, 1, ConfiguratorType::Dim, RegDataTermOfDisplacementVariation<ConfiguratorType, FidelityType> > ( Grid ),
    _reference( Grid, Reference ),
    _displacement( Grid, Displacement ),
    _fidelity ( Fidelity ) {}

  void getNonlinearity ( aol::auto_container<1,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::VecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType>( this->getConfigurator(), _displacement, El, QuadPoint, RefCoord, transformedEl, transformedCoord, coordinateWithinLimits );

    DiscFuncs[0].evaluateGradient( transformedEl, transformedCoord, NL );
    const RealType diff = _fidelity.evaluateUDerivative( DiscFuncs[0].evaluate ( transformedEl, transformedCoord ), _reference.evaluateAtQuadPoint ( El, QuadPoint ) );
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      if ( coordinateWithinLimits[i] )
        NL[i] *= diff;
      else
        NL[i] = 0;
  }
};

/**
 * \brief  \f$ \frac{1}{2}\int ((T\circ\Phi)(x)-R(x))^2 dx \f$
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class SSDEnergy : public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, SSDEnergy<ConfiguratorType>, ConfiguratorType::Dim > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;

  SSDEnergy ( const typename ConfiguratorType::InitType &Grid,
              const aol::Vector<RealType> &ImR,
              const aol::Vector<RealType> &ImT )
  : aol::FENonlinIntegrationVectorInterface<ConfiguratorType, SSDEnergy<ConfiguratorType>, ConfiguratorType::Dim > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  template <typename RegistrationConfiguratorType>
  SSDEnergy ( const typename ConfiguratorType::InitType &Grid,
              const aol::Vector<RealType> &ImR,
              const aol::Vector<RealType> &ImT,
              const RegistrationConfiguratorType &/*RegisConfig*/ )
  : aol::FENonlinIntegrationVectorInterface<ConfiguratorType, SSDEnergy<ConfiguratorType>, ConfiguratorType::Dim > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord );

    return 0.5 * aol::Sqr ( _t.evaluate(transformed_el, transformed_local_coord) - _r.evaluateAtQuadPoint(El, QuadPoint) );
  }
};

/**
 * \brief  \f$ \left( \int ((T\circ\Phi)(x)-R(x))\nabla T(x)\cdot \Psi_i(x) dx\right)_i \f$
 *
 * \author Berkels, Han
 */
template <typename ConfiguratorType>
class SSDForce : public aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, SSDForce<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;

  SSDForce ( const typename ConfiguratorType::InitType &Grid,
             const aol::Vector<RealType> &ImR,
             const aol::Vector<RealType> &ImT )
  : aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim, SSDForce<ConfiguratorType> > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  template <typename RegistrationConfiguratorType>
  SSDForce ( const typename ConfiguratorType::InitType &Grid,
             const aol::Vector<RealType> &ImR,
             const aol::Vector<RealType> &ImT,
             const RegistrationConfiguratorType &/*RegisConfig*/ )
  : aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim, SSDForce<ConfiguratorType> > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                       const typename ConfiguratorType::ElementType &El,
                       int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                       aol::Vec<ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;

    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord, coordinateWithinLimits );

    const RealType Td= _t.evaluate(transformed_el, transformed_local_coord);
    const RealType R = _r.evaluateAtQuadPoint(El, QuadPoint);
    _t.evaluateGradient( transformed_el, transformed_local_coord, NL );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      if ( coordinateWithinLimits[i] )
        NL[i] *= (Td-R);
      else
        NL[i] = 0;
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class SSDGNFuncFromImage : public aol::FELeastSquaresFunctionalInterface<ConfiguratorType, SSDGNFuncFromImage<ConfiguratorType>, aol::SparseMatrix<typename ConfiguratorType::RealType>, true> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
public:
  SSDGNFuncFromImage ( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &ImR,
                       const aol::MultiVector<RealType> &Displacement )
    : aol::FELeastSquaresFunctionalInterface<ConfiguratorType, SSDGNFuncFromImage<ConfiguratorType>, aol::SparseMatrix<typename ConfiguratorType::RealType>, true> ( Grid ),
     _r( Grid, ImR ),
     _displacement ( Grid, Displacement ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), _displacement, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord );

    return ( DiscFuncs.evaluate(transformedEl, transformedLocalCoord) - _r.evaluateAtQuadPoint(El, QuadPoint) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class VarOfSSDGNFunc : public aol::FELeastSquaresFunctionalVectorDerivativeInterface<ConfiguratorType, VarOfSSDGNFunc<ConfiguratorType>, aol::SparseBlockMatrix<aol::SparseMatrix<typename ConfiguratorType::RealType> > > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;
public:
  VarOfSSDGNFunc ( const typename ConfiguratorType::InitType &Grid,
                   const aol::Vector<RealType> &ImR,
                   const aol::Vector<RealType> &ImT )
    : aol::FELeastSquaresFunctionalVectorDerivativeInterface<ConfiguratorType, VarOfSSDGNFunc<ConfiguratorType>, aol::SparseBlockMatrix<aol::SparseMatrix<typename ConfiguratorType::RealType> > > ( Grid ),
     _r( Grid, ImR ),
     _t ( Grid, ImT ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim, RealType> &NL ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    _t.evaluateGradient ( transformedEl, transformedLocalCoord, NL );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      if ( coordinateWithinLimits[i] == false )
        NL[i] = 0;
    }
  }
};

template <typename ConfiguratorType> class CrossCorrelationEnergy;
template <typename ConfiguratorType> class CrossCorrelationForce;

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class CCRegistrationConfigurator : public BaseRegistrationConfigurator<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ImageDOFType;

  CCRegistrationConfigurator ( ) : BaseRegistrationConfigurator<ConfiguratorType> ( ) {}

  explicit CCRegistrationConfigurator ( const aol::ParameterParser &Parser ) : BaseRegistrationConfigurator<ConfiguratorType> ( Parser ) {}

  template <typename RegistrationMultilevelDescentType>
  explicit CCRegistrationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
  : BaseRegistrationConfigurator<ConfiguratorType> ( RegisMLD.getParserReference() ) {}

  typedef CrossCorrelationEnergy<ConfiguratorType> Energy;
  typedef CrossCorrelationForce<ConfiguratorType> EnergyVariation;
};

template <typename ConfiguratorType> class NormalizedCrossCorrelationEnergy;
template <typename ConfiguratorType> class VariationOfNormalizedCrossCorrelationEnergy;

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class NCCRegistrationConfigurator : public BaseRegistrationConfigurator<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ImageDOFType;

  NCCRegistrationConfigurator ( ) : BaseRegistrationConfigurator<ConfiguratorType> ( ) {}

  explicit NCCRegistrationConfigurator ( const aol::ParameterParser &Parser ) : BaseRegistrationConfigurator<ConfiguratorType> ( Parser ) {}

  template <typename RegistrationMultilevelDescentType>
  explicit NCCRegistrationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
    : BaseRegistrationConfigurator<ConfiguratorType> ( RegisMLD.getParserReference() ) {}

  typedef NormalizedCrossCorrelationEnergy<ConfiguratorType> Energy;
  typedef VariationOfNormalizedCrossCorrelationEnergy<ConfiguratorType> EnergyVariation;
};

/**
 * \brief  \f$ -\frac{1}{2}\int (T\circ\Phi)(x)R(x) dx \f$
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class CrossCorrelationEnergy : public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, CrossCorrelationEnergy<ConfiguratorType>, ConfiguratorType::Dim > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;

  CrossCorrelationEnergy ( const typename ConfiguratorType::InitType &Grid,
                           const aol::Vector<RealType> &ImR,
                           const aol::Vector<RealType> &ImT )
  : aol::FENonlinIntegrationVectorInterface<ConfiguratorType, CrossCorrelationEnergy<ConfiguratorType>, ConfiguratorType::Dim > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  template <typename RegistrationConfiguratorType>
  CrossCorrelationEnergy ( const typename ConfiguratorType::InitType &Grid,
                           const aol::Vector<RealType> &ImR,
                           const aol::Vector<RealType> &ImT,
                           const RegistrationConfiguratorType &/*RegisConfig*/ )
  : aol::FENonlinIntegrationVectorInterface<ConfiguratorType, CrossCorrelationEnergy<ConfiguratorType>, ConfiguratorType::Dim > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                             const typename ConfiguratorType::ElementType &El,
                             int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord );

    return - _t.evaluate(transformed_el, transformed_local_coord) * _r.evaluateAtQuadPoint(El, QuadPoint);
  }
};

/**
 * \brief  \f$ \left( - \int R(x)\nabla T(\Phi(x))\cdot \Psi_i(x) dx\right)_i \f$
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class CrossCorrelationForce : public aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, CrossCorrelationForce<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;

  CrossCorrelationForce ( const typename ConfiguratorType::InitType &Grid,
                          const aol::Vector<RealType> &ImR,
                          const aol::Vector<RealType> &ImT )
  : aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim, CrossCorrelationForce<ConfiguratorType> > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  template <typename RegistrationConfiguratorType>
  CrossCorrelationForce ( const typename ConfiguratorType::InitType &Grid,
                          const aol::Vector<RealType> &ImR,
                          const aol::Vector<RealType> &ImT,
                          const RegistrationConfiguratorType &/*RegisConfig*/ )
  : aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim, CrossCorrelationForce<ConfiguratorType> > ( Grid ),
    _r( Grid, ImR ),
    _t( Grid, ImT ) {}

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                       const typename ConfiguratorType::ElementType &El,
                       int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                       aol::Vec<ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    if ( !qc::transformCoord<ConfiguratorType> ( this->_initializer, DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ) {
      NL.setZero();
      return;
    }

    _t.evaluateGradient( transformed_el, transformed_local_coord, NL );
    NL *= -_r.evaluateAtQuadPoint(El, QuadPoint);
  }
};

/**
 * \brief  \f$ \left( - \int ( R(x) + T(\Phi(x)) E) \nabla T(\Phi(x))\cdot \Psi_i(x) dx\right)_i \f$
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class NormalizedCrossCorrelationForce : public aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, NormalizedCrossCorrelationForce<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;
  const RealType _NCCEnergy;

  NormalizedCrossCorrelationForce ( const typename ConfiguratorType::InitType &Grid,
                                    const aol::Vector<RealType> &ImR,
                                    const aol::Vector<RealType> &ImT,
                                    const RealType NCCEnergy )
    : aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim, NormalizedCrossCorrelationForce<ConfiguratorType> > ( Grid ),
      _r( Grid, ImR ),
      _t( Grid, ImT ),
      _NCCEnergy ( NCCEnergy ) {}

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !qc::transformCoord<ConfiguratorType, true> ( this->getConfigurator(), DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) ) {
      NL.setZero();
      return;
    }

    _t.evaluateGradient( transformedEl, transformedLocalCoord, NL );
    NL *= - ( _r.evaluateAtQuadPoint(El, QuadPoint) + _t.evaluate ( transformedEl, transformedLocalCoord ) * _NCCEnergy );
  }
};

/**
 * \brief  \f$ - \int \frac{(T\circ\Phi-\textrm{mean}(T\circ\Phi))(x)(R(x)-\textrm{mean}(R))}{\textrm{var}(T\circ\Phi)\textrm{var}(R)} dx \f$
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class NormalizedCrossCorrelationEnergy : public aol::StandardGenEnergyOp<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  const aol::MassOp<ConfiguratorType> _massOp;
  typename ConfiguratorType::ArrayType _normalizedR;
  const typename ConfiguratorType::ArrayType _t;
  mutable RealType _varOfLastDeformedT;
  mutable typename ConfiguratorType::ArrayType _lastNormalizedDeformedT;
public:
  NormalizedCrossCorrelationEnergy ( const typename ConfiguratorType::InitType &Grid,
                                     const aol::Vector<RealType> &ImR,
                                     const aol::Vector<RealType> &ImT,
                                     // Dummy argument so that this class can be used in qc::StandardRegistration.
                                     // It's potentially filled with an invalid reference, so don't give it a name!
                                     const NCCRegistrationConfigurator<ConfiguratorType> & = *static_cast<NCCRegistrationConfigurator<ConfiguratorType>*> ( NULL ) )
    : _grid ( Grid ),
      _massOp ( Grid, aol::ASSEMBLED ),
      _normalizedR( Grid ),
      _t( ImT, Grid, aol::FLAT_COPY ),
      _varOfLastDeformedT ( 0 ),
      _lastNormalizedDeformedT ( Grid ) {
    normalizeImageForNCC ( _massOp, ImR, _normalizedR );
  }

  //! Subtracts mean and then devides by variance. Furthermore, returns the variance.
  static RealType normalizeImageForNCC ( const aol::Op<aol::Vector<RealType> > &MassOp, const aol::Vector<RealType> &Image, aol::Vector<RealType> &NormalizedImage ) {
    typedef typename ConfiguratorType::RealType RealType;
    aol::Vector<RealType> temp ( Image, aol::STRUCT_COPY );
    MassOp.apply ( Image, temp );
    const RealType imageMean = temp.sum();
    NormalizedImage = Image;
    NormalizedImage.addToAll ( -imageMean );
    MassOp.apply ( NormalizedImage, temp );
    const RealType imageVar = sqrt ( NormalizedImage * temp );
    if ( aol::appeqAbsolute ( imageVar, aol::ZOTrait<RealType>::zero ) == false )
      NormalizedImage /= imageVar;
    return imageVar;
  }

  virtual void deformImage ( const typename ConfiguratorType::ArrayType &Image, typename ConfiguratorType::ArrayType &DeformedImage , const aol::MultiVector<RealType> &Phi) const {
    qc::FastDeformImage<ConfiguratorType, ConfiguratorType::Dim> ( Image, _grid, DeformedImage, Phi );
  }

  RealType computeEnergy ( const aol::MultiVector<RealType> &MArg ) const {
    typename ConfiguratorType::ArrayType deformedT ( _grid );
    // For NCC using zero extension is not a good idea, this would make calculating the mean and
    // the variance more complicated. Probably even R would need to be renormalized based on the
    // altered domain mask.
    deformImage ( _t, deformedT, MArg );
    _varOfLastDeformedT = normalizeImageForNCC ( _massOp, deformedT, _lastNormalizedDeformedT );
    // deformedImage is not needed anymore, we can use it as temp vector.
    _massOp.apply ( _lastNormalizedDeformedT, deformedT );
    this->_lastEnergy = - ( deformedT * _normalizedR );
    return ( this->_lastEnergy );
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    Dest += computeEnergy ( MArg );
  }

  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    Dest = computeEnergy ( Arg );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    aol::Scalar<RealType> energy;
    apply ( MArg, energy );
    aol::Vector<RealType> temp ( _normalizedR );
    temp.addMultiple ( _lastNormalizedDeformedT, energy[0] );
    if ( aol::appeqAbsolute ( _varOfLastDeformedT, aol::ZOTrait<RealType>::zero ) == false )
      temp /= _varOfLastDeformedT;

    CrossCorrelationForce<ConfiguratorType> force ( _grid, temp, _t );
    force.apply ( MArg, MDest );
  }

  const typename ConfiguratorType::ArrayType &getNormalizedRRef ( ) const {
    return _normalizedR;
  }

  const aol::MassOp<ConfiguratorType> &getMassOpRef ( ) const {
    return _massOp;
  }
};

template <typename ConfiguratorType>
class VariationOfNormalizedCrossCorrelationEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

private:
  NormalizedCrossCorrelationEnergy<ConfiguratorType> _E;
  aol::DerivativeWrapper<RealType, qc::NormalizedCrossCorrelationEnergy<ConfiguratorType>, aol::MultiVector<RealType> > _DE;

public:
  VariationOfNormalizedCrossCorrelationEnergy ( const typename ConfiguratorType::InitType &Grid,
                                                const aol::Vector<RealType> &ImR,
                                                const aol::Vector<RealType> &ImT,
                                                // Dummy argument so that this class can be used in qc::StandardRegistration.
                                                // It's potentially filled with an invalid reference, so don't give it a name!
                                                const NCCRegistrationConfigurator<ConfiguratorType> & = *static_cast<NCCRegistrationConfigurator<ConfiguratorType>*> ( NULL ) )
    : _E ( Grid, ImR, ImT ), _DE ( _E ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    _DE.applyAdd ( MArg, MDest );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class SSDDirichletGNFunc : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef aol::SparseMatrix<RealType> MatrixType;
  typedef aol::SparseBlockMatrix<MatrixType> BlockMatrixType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_r;
  const aol::Vector<RealType> &_t;
  const qc::VarOfSSDGNFunc<ConfiguratorType> _DF;
  aol::RandomAccessContainer<qc::DirichtletEnergyComponentGNFunc<ConfiguratorType> > _regF;
  aol::Vector<int> _dimRangeF;
public:
  SSDDirichletGNFunc ( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &ImR,
                       const aol::Vector<RealType> &ImT,
                       const RealType Lambda )
    : _grid ( Grid ),
      _r( ImR ),
      _t( ImT ),
      _DF ( _grid, ImR, ImT ),
      _dimRangeF ( 1 + ConfiguratorType::Dim * ConfiguratorType::Dim ) {
    _dimRangeF[0] = _DF.getDimRangeF ();
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _regF.constructDatumAndPushBack ( _grid, i, Lambda );

    for ( int j = 0; j < ConfiguratorType::Dim; ++j )
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        _dimRangeF[i+j*ConfiguratorType::Dim+1] = _regF[i].getDimRangeF ( );
  }

  void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    qc::SSDGNFuncFromImage<ConfiguratorType> dataF ( _grid, _r, Arg );
    dataF.apply ( _t, Dest[0] );
    for ( int j = 0; j < ConfiguratorType::Dim; ++j )
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        _regF[i].apply ( Arg[j], Dest[i+j*ConfiguratorType::Dim+1] );
  }

  void applyAdd ( const aol::MultiVector<RealType> &, aol::MultiVector<RealType> & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &Arg, BlockMatrixType &MatDest ) const {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      // Allocate matrices for Jacobian of the data term for each component of the deformation.
      MatDest.template allocateMatrixOfType<MatrixType> ( 0, i, _dimRangeF[0], Arg[i].size() );
      // Just make flat copies of the regularizer matrices.
      for ( int j = 0; j < ConfiguratorType::Dim; ++j ) {
        MatDest.template allocateMatrixAndCopyOtherMatrix<MatrixType>( i+j*ConfiguratorType::Dim+1, j, _regF[i].getMatrixRef(), aol::FLAT_COPY );
      }
    }

    // Data term derivative
    _DF.apply ( Arg, MatDest );
  }

  const aol::Vector<int> &getDimRangeF () const {
    return _dimRangeF;
  }
};

/**
 * \author Berkels
 *
 * \note Since this is derived from aol::StandardGenEnergyOp, child classes need to set
 *       this->_lastEnergy when the energy is computed.
 */
template <typename ConfiguratorType, typename ConfiguratorType1D>
class LineWiseRegisEnergyBase : public aol::StandardGenEnergyOp<aol::Vector<typename ConfiguratorType::RealType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const typename ConfiguratorType1D::InitType _grid1D;
  const int _numLines;
  const typename ConfiguratorType::ArrayType &_u0;
  const typename ConfiguratorType::ArrayType &_u;
  aol::RandomAccessContainer<const aol::Vector<RealType> > _u0Lines;
  aol::RandomAccessContainer<const aol::Vector<RealType> > _uLines;
public:
  LineWiseRegisEnergyBase ( const typename ConfiguratorType::InitType &Initializer,
                          const typename ConfiguratorType::ArrayType &U0,
                          const typename ConfiguratorType::ArrayType &U ):
  _grid( Initializer ),
  _grid1D ( aol::Vec3<int> ( U0.getNumX(), 1, 1 ) ),
  _numLines ( U0.getNumY() ),
  _u0 ( U0 ),
  _u ( U ) {
    for ( int line = 0; line < _numLines; ++line ) {
      _u0Lines.constructDatumAndPushBack (_u0.getRowDataPointer ( line ), _grid1D.getNumberOfNodes(), aol::FLAT_COPY );
      _uLines.constructDatumAndPushBack (_u.getRowDataPointer ( line ), _grid1D.getNumberOfNodes(), aol::FLAT_COPY );
    }
  }

  const typename ConfiguratorType1D::InitType& getGrid1DReference ( ) const {
    return _grid1D;
  }

  const aol::Vector<RealType>& getU0LineReference ( const int Line ) const {
    return _u0Lines[Line];
  }

  const aol::Vector<RealType>& getULineReference ( const int Line ) const {
    return _uLines[Line];
  }
};

/**
 * \author Berkels
 *
 * \note The displacement is interpreted relative to NumX-1, not max(NumX,NumY)-1.
 *       This has to be taken into account when using this class on a non-quadratic grid.
 */
template <typename ConfiguratorType, typename ConfiguratorType1D, typename RegistrationConfiguratorType1D>
class FELineShiftRegisEnergy : public qc::LineWiseRegisEnergyBase<ConfiguratorType, ConfiguratorType1D> {
  typedef typename ConfiguratorType::RealType RealType;
public:
  FELineShiftRegisEnergy ( const typename ConfiguratorType::InitType &Initializer,
                           const typename ConfiguratorType::ArrayType &U0,
                           const typename ConfiguratorType::ArrayType &U )
  : qc::LineWiseRegisEnergyBase<ConfiguratorType, ConfiguratorType1D> ( Initializer, U0, U ) {}

  virtual void apply ( const aol::Vector<RealType> &Y, aol::Scalar<RealType> &Dest ) const {
    const typename ConfiguratorType::ArrayType y ( Y, this->_grid, aol::FLAT_COPY );
    Dest.setZero();
    for ( int line = 0; line < this->_numLines; ++line ) {
      if ( aol::isNaN ( this->_u0Lines[line][0] ) || aol::isNaN ( this->_uLines[line][0] ) )
        continue;

      aol::MultiVector<RealType> yLine ( y.getRowDataPointer ( line ), this->_grid1D.getNumberOfNodes(), aol::FLAT_COPY );
      typename RegistrationConfiguratorType1D::Energy lineEnergy ( this->_grid1D,  this->_u0Lines[line],  this->_uLines[line] );
      lineEnergy.applyAdd ( yLine, Dest );
    }
    Dest /= this->_numLines;
    this->_lastEnergy = Dest[0];
  }

  void applyDerivative ( const aol::Vector<RealType> &Y, aol::Vector<RealType> &DY ) const {
    const typename ConfiguratorType::ArrayType yArray ( Y, this->_grid, aol::FLAT_COPY );
    const typename ConfiguratorType::ArrayType dY ( DY, this->_grid, aol::FLAT_COPY );
    for ( int line = 0; line < this->_numLines; ++line ) {
      if ( aol::isNaN ( this->_u0Lines[line][0] ) || aol::isNaN ( this->_uLines[line][0] ) )
        continue;

      aol::MultiVector<RealType> yLine ( yArray.getRowDataPointer ( line ), this->_grid1D.getNumberOfNodes(), aol::FLAT_COPY );
      aol::MultiVector<RealType> dYLine ( dY.getRowDataPointer ( line ), this->_grid1D.getNumberOfNodes(), aol::FLAT_COPY );
      typename RegistrationConfiguratorType1D::EnergyVariation lineEnergyDeriv ( this->_grid1D,  this->_u0Lines[line],  this->_uLines[line] );
      lineEnergyDeriv.apply ( yLine, dYLine );
      dYLine /= this->_numLines;
    }
  }

  /* Only implemented for SSD
  void applyDerivativeWRTU ( const typename ConfiguratorType::ArrayType &Y, typename ConfiguratorType::ArrayType &DU ) const {
    for ( int line = 0; line < this->_numLines; ++line ) {
      aol::MultiVector<RealType> yLine ( Y.getRowDataPointer ( line ), this->_grid1D.getNumberOfNodes(), aol::FLAT_COPY );
      aol::Vector<RealType> dULine ( DU.getRowDataPointer ( line ), this->_grid1D.getNumberOfNodes(), aol::FLAT_COPY );
      qc::VariationOfSSDEnergyWRTTemplate<ConfiguratorType1D> lineSSDDeriv ( this->_grid1D, this->_u0Lines[line], yLine );
      lineSSDDeriv.apply ( this->_uLines[line], dULine );
      dULine /= this->_numLines;
    }
  }
  */

  virtual void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( Arg, tmp );
    Dest += tmp;
  }
};

} // end namespace qc

#endif // __REGISTRATION_H
