#ifndef __QUOCTIMESTEPSAVER_H
#define __QUOCTIMESTEPSAVER_H

#include <timestepSaver.h>
#include <deformations.h>
#include <ChanVese.h>

namespace qc {

/**
 * Class to extract slices from a MultiVector that contains 3D arrays.
 *
 * \todo Find a better place for this class.
 * \todo Make the slice direction selectable.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class SliceExtractor {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiguratorType2D;
private:
  const typename ConfiguratorType::InitType &_grid;
  const aol::MultiVector<RealType> &_mVec;
  const int _numOfSlices, _numX, _numY;
  const qc::RectangularGrid<qc::QC_2D> _grid2D;
public:

  SliceExtractor ( const typename ConfiguratorType::InitType &Initializer,
                   const aol::MultiVector<RealType> &MVec )
    : _grid ( Initializer ),
      _mVec ( MVec ),
      _numOfSlices ( _grid.getNumZ() ),
      _numX ( _grid.getNumX() ),
      _numY ( _grid.getNumY() ),
      _grid2D ( aol::Vec3<int>( _numX, _numY, 1 ) ) {}

  void getSlice ( const int SliceNumber, aol::MultiVector<RealType> &MVecSlice ) const {
    for ( int i = 0; i < _mVec.numComponents(); i++ ){
      qc::ScalarArray<RealType, qc::QC_2D> tmp2d ( MVecSlice[i], _grid2D, aol::FLAT_COPY );
      qc::ScalarArray<RealType, qc::QC_3D> tmp3d ( _mVec[i], _grid, aol::FLAT_COPY );
      tmp3d.getSlice ( qc::QC_Z, SliceNumber, tmp2d );
    }
  }

  int getNumOfSlices () const {
    return _numOfSlices;
  }

  int getNumX() const {
    return _numX;
  }

  int getNumY() const {
    return _numY;
  }

  const qc::RectangularGrid<qc::QC_2D>& get2DGridReference () const {
    return _grid2D;
  }
};

/**
 * Class for Quoc specific safe functions.
 *
 * ATTENTION: Since this class is derived with aol::TimestepSaver as
 * virtual base class, you can't rely that the constructor of
 * aol::TimestepSaver you chose in the initialization list of a
 * QuocTimestepSaver constructor is called in classes derived from
 * QuocTimestepSaver. Therefore we only call the standard constructor
 * and configure aol::TimestepSaver by using the access functions.
 *
 * \author Berkels
 *
 */

template <typename ConfiguratorType>
class QuocTimestepSaver : public virtual aol::TimestepSaver<typename ConfiguratorType::RealType> {
  typedef typename ConfiguratorType::RealType RealType;
public:

  QuocTimestepSaver( ) : aol::TimestepSaver<RealType>(){}

  QuocTimestepSaver( const int k, const char *saveName, const bool quiet = false )
    : aol::TimestepSaver<RealType>(){
    this->setSaveTimestepOffset( k );
    this->setSaveName( saveName );
    this->setQuietMode ( quiet );
  }

  QuocTimestepSaver ( const int k, const int numSaveFirst, const char *saveName, const bool quiet = false )
    : aol::TimestepSaver<RealType>(){
    this->setSaveTimestepOffset( k );
    this->setNumberSaveFirstPics( numSaveFirst );
    this->setSaveName( saveName );
    this->setQuietMode ( quiet );
  }

  void saveDeformedImage ( const aol::Vector<RealType> &Img,
                           const typename ConfiguratorType::InitType &Grid,
                           const aol::MultiVector<RealType> &Deformation,
                           const int step = -1,
                           const int outerStep = -1 ) const {
    aol::Vector<RealType> tmp( Img, aol::STRUCT_COPY );
    qc::DeformImage<ConfiguratorType>( Img, Grid, tmp, Deformation);
    string filename;
    this->initBaseSaveNameString ( outerStep, step, filename );
    qc::writeImage<RealType> ( Grid, tmp, filename.c_str() );
    // since qc::writeImage appends the correct suffix on its own, we have to had it by hand
    if ( Grid.getDimOfWorld() == 2 ) {
      filename += ".png";
    }
    if ( Grid.getDimOfWorld() == 3 ) {
      filename += ".raw";
    }
    this->printSaveMessageToConsole ( filename.c_str() );

  }

  void saveDeformedImage ( const aol::Vector<RealType> &Img,
                           const typename ConfiguratorType::InitType &Grid,
                           const aol::MultiVector<RealType> &Deformation1,
                           const aol::MultiVector<RealType> &Deformation2,
                           const int step = -1,
                           const int outerStep = -1 ) const {
    aol::Vector<RealType> tmp( Img, aol::STRUCT_COPY );
    qc::DeformImage<ConfiguratorType>( Img, Grid, tmp, Deformation1);
    saveDeformedImage ( tmp, Grid, Deformation2, step, outerStep );
  }

  /**
   * Plots vectorfield with gnuplot in different colors. The different
   * color regions are specified with MaskVector.
   *
   * If levelset functions are supplied, their zero isolines will be
   * added to the plot.
   *
   * \author Berkels
   */
  void plotVectorField ( const aol::MultiVector< RealType > &vectorfield,
                         const typename ConfiguratorType::InitType &grid,
                         const aol::PlotOutFileType OutType,
                         const std::vector<const qc::BitArray<qc::QC_2D>*> &MaskVector,
                         const int outerStep = -1,
                         const int step = -1,
                         const aol::MultiVector<RealType> *LevelsetFunctions = NULL,
                         const RealType Isovalue = aol::NumberTrait<RealType>::zero,
                         const RealType Spacing = 0.05 ) const {
    if  ( this->checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        aol::PlotDataFileHandler<RealType> plotDataFileHandler;
        plotDataFileHandler.generateColoredVectorfieldData ( vectorfield, grid, MaskVector, Spacing );
        if ( LevelsetFunctions ){
          ConfiguratorType confDummy( grid );
          plotDataFileHandler.generateIsolineData( *LevelsetFunctions, grid, confDummy, Isovalue );
        }
        aol::Plotter<RealType> plotter;
        plotter.addPlotCommandsFromHandler( plotDataFileHandler );
        string filename;
        this->initBaseSaveNameString ( outerStep, step, filename );
        plotter.set_outfile_base_name ( filename.c_str() );
        plotter.genPlot( OutType, &filename );
        this->printSaveMessageToConsole ( filename.c_str() );
      } else {
        cerr << aol::color::red << "WARNING: plotVectorField not implemented for Dim != 2.\n" << aol::color::reset;
      }
    }
  }

  /**
   * Plots vectorfield with gnuplot in different colors. The different
   * color regions are specified with MaskVector.
   *
   * If levelset functions are supplied, their zero isolines will be
   * added to the plot.
   *
   * This function is needed because there is no standard conversion from
   * std::vector<qc::BitArray<qc::QC_2D>*> to std::vector<const qc::BitArray<qc::QC_2D>*>.
   *
   * \author Berkels
   */
  void plotVectorField ( const aol::MultiVector< RealType > &vectorfield,
                         const typename ConfiguratorType::InitType &grid,
                         const aol::PlotOutFileType OutType,
                         const std::vector<qc::BitArray<qc::QC_2D>*> &MaskVector,
                         const int outerStep = -1,
                         const int step = -1,
                         const aol::MultiVector<RealType> *LevelsetFunctions = NULL,
                         const RealType Isovalue = aol::NumberTrait<RealType>::zero ) const {
    std::vector<const qc::BitArray<qc::QC_2D>*> maskVectorConst( MaskVector.size() );
    for ( unsigned int i = 0; i<MaskVector.size(); i++ )
      maskVectorConst[i] = MaskVector[i];
    plotVectorField ( vectorfield, grid, OutType, maskVectorConst, outerStep, step, LevelsetFunctions, Isovalue );
  }

  /**
   * Plots vectorfield with gnuplot. A Mask can be given to hide parts
   * of vectorfield in the plot.
   *
   * \author Berkels
   */
  void plotVectorField ( const aol::MultiVector< RealType > &vectorfield,
                         const typename ConfiguratorType::InitType &grid,
                         const aol::PlotOutFileType OutType,
                         const int outerStep = -1,
                         const int step = -1,
                         const qc::BitArray<qc::QC_2D> *Mask = NULL ) const {
    std::vector<const qc::BitArray<qc::QC_2D>*> maskVector(1);
    maskVector[0] = Mask;
    plotVectorField ( vectorfield, grid, OutType, maskVector, outerStep, step );
  }

  /**
   * Plots Vectorfield with gnuplot in different colors. The different
   * color regions are based on the levelset functions.
   *
   * The plot is written in all formats selected in the OutTypeVec
   * vector. If DrawLevelLines == true the zero isolines of the levelset
   * functions will be added to the plot.
   *
   * \author Berkels
   */
  void plotVectorFieldSegmented ( const aol::MultiVector<RealType> &Vectorfield,
                                  const aol::MultiVector<RealType> &LevelsetFunctions,
                                  const typename ConfiguratorType::InitType &Grid,
                                  const std::vector<aol::PlotOutFileType> &OutTypeVec,
                                  const int OuterStep = -1,
                                  const int Step = -1,
                                  const bool DrawLevelLines = false,
                                  const RealType Isovalue = aol::NumberTrait<RealType>::zero ) const {
    if ( Grid.getDimOfWorld() == 2 ) {
      const int numberOfLevelsetFunctions = LevelsetFunctions.numComponents();
      std::vector<qc::BitArray<qc::QC_2D>*> maskVector( 1<<numberOfLevelsetFunctions );
      for ( unsigned int i = 0; i< maskVector.size(); i++ )
        maskVector[i] = new qc::BitArray<qc::QC_2D> ( qc::GridSize<qc::QC_2D> ( Grid ) );

      qc::DataGenerator<ConfiguratorType> generator( Grid );
      generator.generateMaskFromLevelsetfunctions ( maskVector, LevelsetFunctions, Isovalue );

      for ( unsigned int i = 0; i < OutTypeVec.size(); i++ )
        this->plotVectorField ( Vectorfield, Grid, OutTypeVec[i], maskVector, OuterStep, Step, (DrawLevelLines ? &LevelsetFunctions : NULL), Isovalue );

      for ( unsigned int i = 0; i< maskVector.size(); i++ )
        delete maskVector[i];
    }
    // Assume a 3D grid
    else {
      cerr << aol::color::red << "WARNING: plotVectorFieldSegmented not implemented for Dim != 2.\n" << aol::color::reset;
    }
  }
  /**
   * Plots Vectorfield with gnuplot in different colors. The different
   * color regions are based on the levelset functions.
   *
   * The plot is written in the format specified with OutType. If
   * DrawLevelLines == true the zero isolines of the levelset
   * functions will be added to the plot.
   *
   * \author Berkels
   */
  void plotVectorFieldSegmented ( const aol::MultiVector< RealType > &Vectorfield,
                                  const aol::MultiVector<RealType> &LevelsetFunctions,
                                  const typename ConfiguratorType::InitType &Grid,
                                  const aol::PlotOutFileType OutType,
                                  const int OuterStep = -1,
                                  const int Step = -1,
                                  const bool DrawLevelLines = false ) const {
    std::vector<aol::PlotOutFileType> outTypeVec(1);
    outTypeVec[0] = OutType;
    plotVectorFieldSegmented ( Vectorfield, LevelsetFunctions, Grid, outTypeVec, OuterStep, Step, DrawLevelLines );
  }

  void saveTimestepPNGDrawRedIsoline ( const int step,
                                       const double isovalue,
                                       const aol::Vector<RealType> &levelsetFunction,
                                       const aol::Vector<RealType> &img,
                                       const typename ConfiguratorType::InitType &grid ) const {
    saveTimestepPNGDrawRedIsoline ( -1, step, isovalue, levelsetFunction, img, grid );
  }
  void saveTimestepPNGDrawRedIsoline ( const int outerStep,
                                       const int step,
                                       const double isovalue,
                                       const aol::Vector<RealType> &levelsetFunction,
                                       const aol::Vector<RealType> &img,
                                       const typename ConfiguratorType::InitType &grid ) const {
    // To prevent code duplication, just call the function here, which draws
    // two isolines and deactivate the drawing of the second one.
    saveTimestepPNGDrawRedAndBlueIsoline ( outerStep, step, isovalue, 666, levelsetFunction, img, grid, false );
  }

  void saveTimestepPNGDrawRedAndBlueIsoline ( const int outerStep,
                                              const int step,
                                              const double isovalue1,
                                              const double isovalue2,
                                              const aol::Vector<RealType> &levelsetFunction,
                                              const aol::Vector<RealType> &img,
                                              const typename ConfiguratorType::InitType &grid,
                                              const bool DrawSecondIsoline = true ) const {
    if  ( this->checkSaveConditions ( step ) ) {
      qc::LevelSetDrawer<ConfiguratorType> drawer ( grid );
      qc::ScalarArray<RealType, qc::QC_2D> levelsetFunctionArray ( levelsetFunction, grid );

      qc::MultiArray<RealType, 2, 3 > colorImg ( grid );
      colorImg.setAll( img );

      // Draw the isovalue1 levelline in red.
      drawer.draw ( levelsetFunctionArray, colorImg, isovalue1, 0 );
      // Draw the isovalue2 levelline in blue.
      if( DrawSecondIsoline )
        drawer.draw ( levelsetFunctionArray, colorImg, isovalue2, 2 );

      string filename;
      this->initSaveNameStringAndSendMessageToConsole ( outerStep, step, ".png", filename );

      colorImg.setQuietMode ( true );
      colorImg.setOverflowHandling(aol::CLIP_THEN_SCALE, 0., 1.);
      colorImg.savePNG ( filename.c_str() );
    }
  }

  void saveIsolineWithGnuplot ( const int outerStep,
                                const int step,
                                const double isovalue,
                                const aol::PlotOutFileType OutType,
                                const aol::Vector<RealType> &levelsetFunction,
                                const typename ConfiguratorType::InitType &grid,
                                const bool HideBorder = false,
                                const aol::Vector<RealType> *PSecondLevelsetFunction = NULL ) const {
    if  ( this->checkSaveConditions ( step ) ) {
      ConfiguratorType confDummy( grid );

      aol::PlotDataFileHandler<RealType> plotDataFileHandler;
      plotDataFileHandler.generateIsolineData( levelsetFunction, grid, confDummy, isovalue);
      aol::Plotter<RealType> plotter;
      plotter.addPlotCommand( plotDataFileHandler.getDataFileNames()[0].c_str(), aol::GNUPLOT_LINE, 1 );
      if ( PSecondLevelsetFunction ) {
        plotDataFileHandler.generateIsolineData( *PSecondLevelsetFunction, grid, confDummy, isovalue);
        plotter.addPlotCommand( plotDataFileHandler.getDataFileNames()[1].c_str(), aol::GNUPLOT_LINE, 2 );
      }
      string special = "set xrange [0:1]\n"
                       "set yrange [0:1]\n"
                       "unset xtics\n"
                       "unset ytics\n"
                       "set size square 1, 1.4285715\n"
                       "set style line 1 linewidth 10\n"
                       "set style line 2 linewidth 10 linecolor -1\n";
      if ( HideBorder )
        special += "unset border\n";
      plotter.setSpecial( special.c_str() );

      string filename;
      this->initBaseSaveNameString ( outerStep, step, filename );
      plotter.set_outfile_base_name ( filename.c_str() );
      plotter.genPlot( OutType, &filename );
      this->printSaveMessageToConsole ( filename.c_str() );
    }
  }

  void savePiecewiseConstantImageSegmented ( const aol::MultiVector<RealType> &MeanGrayValues,
                                             const aol::MultiVector<RealType> &LevelsetFunctions,
                                             const typename ConfiguratorType::InitType &Grid,
                                             const RealType Isovalue,
                                             const int Step = -1,
                                             const int OuterStep = -1 ) const {
    if  ( this->checkSaveConditions ( Step ) ) {
      // Generate and write a grayscale image showing the current segmentation.
      typename ConfiguratorType::InitType::OldAllNodeIterator fnit;
      qc::FastILexMapper<ConfiguratorType::Dim> mapper ( Grid );
      aol::Vector<int> intVec( LevelsetFunctions.numComponents() );
      qc::ScalarArray<RealType, qc::QC_2D> imageArray ( Grid );
      for ( fnit = Grid._nBeginIt; fnit != Grid._nEndIt; ++fnit ) {
        for ( int i = 0; i < LevelsetFunctions.numComponents(); i++ )
          intVec[i] = (LevelsetFunctions[i].get ( mapper.getGlobalIndex ( *fnit ) )>= Isovalue) ? 1 : 0;
        imageArray.set ( *fnit, MeanGrayValues[aol::PowerSetIterator::getPositionNumberFromVector(intVec)][0] );
      }
      this->savePNG( imageArray, Grid, Step, OuterStep );
    }
  }

  void savePiecewiseConstantColorImageSegmented ( const aol::MultiVector<RealType> &MeanColorValues,
                                                  const aol::MultiVector<RealType> &LevelsetFunctions,
                                                  const typename ConfiguratorType::InitType &Grid,
                                                  const RealType Isovalue,
                                                  const int Step = -1,
                                                  const int OuterStep = -1 ) const {
    if  ( this->checkSaveConditions ( Step ) ) {
      // Generate and write a color image showing the current segmentation.
      typename ConfiguratorType::InitType::OldAllNodeIterator fnit;
      qc::FastILexMapper<ConfiguratorType::Dim> mapper ( Grid );
      aol::Vector<int> intVec( LevelsetFunctions.numComponents() );
      qc::MultiArray<RealType, 2, 3> imageMArray ( Grid );
      for ( fnit = Grid._nBeginIt; fnit != Grid._nEndIt; ++fnit ) {
        for ( int i = 0; i < LevelsetFunctions.numComponents(); i++ )
          intVec[i] = (LevelsetFunctions[i].get ( mapper.getGlobalIndex ( *fnit ) )>= Isovalue) ? 1 : 0;
        const int globalIndex = aol::PowerSetIterator::getPositionNumberFromVector(intVec);
        for ( int i = 0; i < imageMArray.numComponents(); i++ )
          imageMArray[i].set ( *fnit, MeanColorValues[globalIndex][i] );
      }
      this->savePNG ( imageMArray, Step, OuterStep );
    }
  }
};

/**
 * \author Berkels
 */
template <typename RealType, int _Dim, typename ArrayType = typename qc::ScalarArrayTrait<RealType, _Dim>::ArrayType>
class DefaultArraySaver : public aol::StepSaverBase<RealType, ArrayType> {
  const bool _scalePNGOutputIn2DTo01;
  const bool _onlySavePNGIn2D;
  const bool _dontSavePNGIn2D;
  const bool _saveRangeInPNGNameIn2D;
  RealType _enhanceContrastSaturationPercentage;
  static const qc::Dimension Dim = static_cast<qc::Dimension> ( _Dim );
  aol::Vec2<RealType> _clipMinMax;
protected:
  void doSaveStep ( const ArrayType &SaveInput, const int Iteration, const char *OverrideBaseSaveName ) const {
    if ( ( Dim != qc::QC_2D ) || ( _onlySavePNGIn2D == false ) )
      SaveInput.save( this->createSaveName ( "", qc::getDefaultArraySuffix( Dim ), Iteration, OverrideBaseSaveName ).c_str(), qc::PGM_DOUBLE_BINARY );
    if ( ( Dim == qc::QC_2D ) && ( _dontSavePNGIn2D == false ) ) {
      // Since ArrayType may differ from qc::ScalarArray<RealType, Dim>, we convert this by accessing the data pointer.
      qc::ScalarArray<RealType, qc::QC_2D> temp ( SaveInput.getNumX(), SaveInput.getNumY(), SaveInput.getData(), aol::FLAT_COPY );
      if ( _scalePNGOutputIn2DTo01 == false ) {
        temp.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
        temp.save( this->createSaveName ( "", ".png", Iteration, OverrideBaseSaveName ).c_str(), qc::PNG_2D );
      }
      else {
        if ( _enhanceContrastSaturationPercentage > 0 )
          temp.setOverflowHandlingToEnhanceContrast ( _enhanceContrastSaturationPercentage );
        else if ( _clipMinMax.norm() > 0 )
          temp.setOverflowHandling ( aol::CLIP_THEN_SCALE, _clipMinMax[0], _clipMinMax[1] );
        else
          temp.setOverflowHandlingToCurrentValueRange ( );
        if ( _saveRangeInPNGNameIn2D )
          temp.save( this->createSaveName ( "", aol::strprintf ( "-%f_%f.png", temp.getOverflowMin(), temp.getOverflowMax() ).c_str(), Iteration, OverrideBaseSaveName ).c_str(), qc::PNG_2D );
        else
          temp.save( this->createSaveName ( "", ".png", Iteration, OverrideBaseSaveName ).c_str(), qc::PNG_2D );
      }
    }
  };

public:
  DefaultArraySaver ( const bool ScalePNGOutputIn2DTo01 = false, const bool OnlySavePNGIn2D = false, const bool DontSavePNGIn2D = false,
                      const bool SaveRangeInPNGNameIn2D = true )
    : _scalePNGOutputIn2DTo01 ( ScalePNGOutputIn2DTo01 ),
      _onlySavePNGIn2D ( OnlySavePNGIn2D ) ,
      _dontSavePNGIn2D ( DontSavePNGIn2D ),
      _saveRangeInPNGNameIn2D ( SaveRangeInPNGNameIn2D ),
      _enhanceContrastSaturationPercentage ( 0 ),
      _clipMinMax ( 0, 0 ) {}

  void setEnhanceContrastSaturationPercentage ( const RealType EnhanceContrastSaturationPercentage ) {
    _enhanceContrastSaturationPercentage = EnhanceContrastSaturationPercentage;
  }

  void setClipMinMax ( const RealType Min, const RealType Max ) {
    _clipMinMax[0] = Min;
    _clipMinMax[1] = Max;
  }
};

} // end namespace qc

#endif // __QUOCTIMESTEPSAVER_H
