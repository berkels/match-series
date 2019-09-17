#ifndef __TIMESTEPSAVER_H
#define __TIMESTEPSAVER_H

#include <levelSetDrawer.h>
#include <colorWheel.h>
#include <aol.h>
#include <gnuplotter.h>
#include <configurators.h>
#include <imageTools.h>
#include <parameterParser.h>

namespace aol {

/**
 * A simple little class for saving every k'th timestep.
 * Just to have a clearer code in the main loop.
 * There are two possibilities to use a TimestepSaver:
   * 1) create an object TimestepSaver (good, if your saving happens in the main program)
   * 2) derive from this class and inherit everything (good, if it happens in an operator-class)
 * It is possible to configure the shell-output (switch off, or configure its colors).
 * Author: Nemitz
 */
template <typename RealType>
class TimestepSaver {

protected:
  int _k;                     // save every k'th timestep
  int _numSaveFirst;          // save the whole first #_saveFirst timesteps
  char _saveNameTrunc[1024];  // the name of the file to be saved without number and suffix
  bool _writeTimeSteps;       // switch the whole class off or on
  bool _writeAllTimeSteps;    // if true, all steps are written, no matter how _k and _numSaveFirst are set
  int _stepDigits;            // the number of digits, which are used for displaying steps in filename generation
  int _outerStepDigits;

  char _saveDirectory[1024];  // the name of the directory, where the output is saved
  bool _quiet;                // make outputs or not

  string _saveNameTruncColor; // color of the savename
  string _stdColor;           // color of the text around the savename
  string _saveOutputColor;    // color of the output of the save-method

public:

  // constructors
  TimestepSaver( ) {
    setDefaults();
  }

  TimestepSaver ( const int k, const char *saveName, const bool quiet = false ) {
    setDefaults();
    _k = k;
    if ( _k <= 0 ) throw aol::Exception ( "TimestepSaver: saveOffset k has to be positive.", __FILE__, __LINE__ );
    strncpy ( _saveNameTrunc, saveName, 1023 );
    _quiet = quiet;
  }

  TimestepSaver ( const int k, const int numSaveFirst, const char *saveName, const bool quiet = false ) {
    setDefaults();
    _k = k;
    if ( _k <= 0 ) throw aol::Exception ( "TimestepSaver: saveOffset k has to be positive.", __FILE__, __LINE__ );
    _numSaveFirst = numSaveFirst;
    strncpy ( _saveNameTrunc, saveName, 1023 );
    _quiet = quiet;
  }

  virtual ~TimestepSaver ( ) {}

  void setDefaults() {
    _k = 10;
    _numSaveFirst = 1;
    strcpy ( _saveNameTrunc, "PleaseSetSaveName" );
    strcpy ( _saveDirectory, "" );
    _saveNameTruncColor   = aol::color::blue;       // default values
    _stdColor        = aol::color::reset;
    _saveOutputColor = aol::color::green;
    _quiet = false;
    _writeTimeSteps = true;
    _writeAllTimeSteps = false;
    _stepDigits = 3;
    _outerStepDigits = 2;
  }

  // ------------------------- methods for accessing class-variables --------------------

  void setSaveTimestepOffset ( const int k ) {
    if ( k > 0 )
      _k = k;
    else
      if ( _k == 0 ) throw aol::Exception ( "TimestepSaver: saveOffset k has to be positive.", __FILE__, __LINE__ );
  }
  void setSaveName ( const char *saveName, const int outerSteps = -1 ) {
    if ( outerSteps != -1 )
      sprintf ( _saveNameTrunc, "%s_%02d", saveName, outerSteps );
    else {
      strncpy ( _saveNameTrunc, saveName, 1023 );
      _saveNameTrunc[1023] = 0;
    }
  }
  void setSaveName ( const char *saveName, const char *saveDir, const int outerSteps ) {
    sprintf ( _saveNameTrunc, "%s_%02d", saveName, outerSteps );
    strncpy ( _saveDirectory, saveDir, 1023 );
  }
  const char* getSaveName () const                                {
    return _saveNameTrunc;
  }
  //! \note The save directory name is assumed to be terminated with '/'.
  void setSaveDirectory ( const char *saveDir )                   {
    strncpy ( _saveDirectory, saveDir, 1023 );
    _saveDirectory[1023] = 0;
  }
  const char* getSaveDirectory () const                           {
    return _saveDirectory;
  }
  //! function will be overridden in class NewtonIterationBase.
  virtual void setQuietMode ( const bool quiet )                  {
    _quiet    = quiet;
  }
  void setSaveNameColor ( const string color )                    {
    _saveNameTruncColor   = color;
  }
  void setStdColor ( const string color )                         {
    _stdColor        = color;
  }
  void setSaveOutputColor ( const string color )                  {
    _saveOutputColor = color;
  }
  void setNumberSaveFirstPics ( const int first )                 {
    _numSaveFirst = first;
  }
  void setWriteTimeSteps ( const bool writeTS )                   {
    _writeTimeSteps = writeTS;
  }
  void setWriteAllTimeSteps ( const bool WriteAllTimeSteps )      {
    _writeAllTimeSteps = WriteAllTimeSteps;
  }
  void activateSaving()                                           {
    _writeTimeSteps = true;
  }
  void deactivateSaving()                                         {
    _writeTimeSteps = false;
  }
  void setStepDigits ( const int StepDigits )                     {
    _stepDigits = StepDigits;
  }
  void setStepDigitsFromMaxIter ( const int MaxIterations )       {
    _stepDigits = aol::countDigitsOfNumber ( MaxIterations );
  }
  void setOuterStepDigits ( const int OuterStepDigits )           {
    _outerStepDigits = OuterStepDigits;
  }
  void setOuterStepDigitsFromMaxIter ( const int MaxIterations )  {
    setOuterStepDigits ( aol::countDigitsOfNumber ( MaxIterations ) );
  }
  void activateSavingAndConfigure( const int NumSaveFirst,
                                   const int SaveOffset,
                                   const char *SaveName,
                                   const char *SaveDirectory,
                                   const int OuterSteps = -1 )    {
    activateSaving();
    setNumberSaveFirstPics( NumSaveFirst );
    setSaveTimestepOffset( SaveOffset );
    setSaveName( SaveName, OuterSteps );
    setSaveDirectory( SaveDirectory );
  }


  // ------------------------- methods for predefined color-schemes ---------------------

  void setAllBlack() {
    _saveNameTruncColor   = aol::color::black;
    _stdColor        = aol::color::black;
    _saveOutputColor = aol::color::black;
  }

protected:
  // --------------------- auxiliary methods for the saving ------------------------------


  void printSaveMessageToConsole ( const char *filename ) const {
    if ( !_quiet ) {
      cerr << _stdColor << "\n>>>>>> saving to file " << _saveNameTruncColor << filename;
      cerr << _stdColor << " <<<<<<\n\n" << _saveOutputColor;
    }
  }
  void initBaseSaveNameString ( const int outerStep, const int step, string &dest, const char* baseName ) const {
    aol::SimpleFormat outerStepsFormat ( _outerStepDigits, 0, ios::fixed | ios::right, '0' );
    aol::SimpleFormat stepsFormat ( _stepDigits, 0, ios::fixed | ios::right, '0' );
    dest = _saveDirectory;
    dest += baseName; //_saveNameTrunc;
    if ( outerStep != -1 ) {
      dest += "_";
      dest += outerStepsFormat ( outerStep );
    }
    if ( step != -1 ) {
      dest += "_";
      dest += stepsFormat ( step );
    }
  }
  void initSaveNameStringAndSendMessageToConsole ( const int outerStep,
                                                   const int step,
                                                   const char *suffix,
                                                   string &dest ) const {
    initBaseSaveNameString ( outerStep, step, dest );
    dest += suffix;
    printSaveMessageToConsole ( dest.c_str() );
  }
  void initBaseSaveNameString ( const int outerStep, const int step, string &dest ) const {
    initBaseSaveNameString ( outerStep, step, dest, _saveNameTrunc );
  }

  void initSaveNameStringMinMaxAndSendMessageToConsole ( const string baseName,
                                                         const aol::Vector<RealType> &img,
                                                         const char *suffix,
                                                         string &dest ) const {
        dest = baseName;
        char temp[1024];
        sprintf ( temp, "_%.5f_%.5f", img.getMinValue(), img.getMaxValue() );
        dest += temp;
        dest += suffix;
//         cerr << "("<<img.getMinValue() << ", "<< img.getMaxValue()<<") ";
        printSaveMessageToConsole ( dest.c_str() );
      }

  void initSaveNameStringMinMaxAndSendMessageToConsole ( const int outerStep,
                                                         const int step,
                                                         const aol::Vector<RealType> &img,
                                                         const char *suffix,
                                                         string &dest ) const {
    string baseName;
    initBaseSaveNameString ( outerStep, step, baseName );
    initSaveNameStringMinMaxAndSendMessageToConsole ( baseName, img, suffix, dest );
  }

//   initSaveNameStringMinMaxAndSendMessageToConsole ( baseSaveName, outerStep, step, img, ".pgm", "" );

  void initSaveNameStringMinMaxAndSendMessageToConsole ( const char* baseName,
                                                         const int outerStep,
                                                         const int step,
                                                         const aol::Vector<RealType> &img,
                                                         const char *suffix,
                                                         string &dest ) const {
    string baseNameString;
    initBaseSaveNameString ( outerStep, step, baseNameString, baseName );
    initSaveNameStringMinMaxAndSendMessageToConsole ( baseNameString, img, suffix, dest );
  }

public:

  bool checkSaveConditions ( const int step ) const {
    if  ( _writeTimeSteps && ( step % _k == 0 || step < _numSaveFirst || _writeAllTimeSteps ) ) return true;
    else return false;
  }

  // -------------------------------------------------------------------------------------
  // the main methods: Save a ScalarArray with the defined name
  // and makes some useful outputs.

  //! save 2d-pic as pgm (short values = 0..255)
  void saveTimestepPGM ( const int step, const qc::ScalarArray<RealType, qc::QC_2D> &img ) const {
    if  ( checkSaveConditions ( step ) ) {
      string filename;
      initSaveNameStringAndSendMessageToConsole ( -1, step, ".pgm", filename );
      img.save ( filename.c_str(), qc::PGM_UNSIGNED_CHAR_BINARY );
      if ( !_quiet ) cerr << endl << aol::color::reset;
    }
  }

  void saveTimestepPGM ( const int step,
                         const aol::Vector<RealType> &img,
                         const qc::GridStructure &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( grid );
        tmp = img;
        string filename;
        initSaveNameStringAndSendMessageToConsole ( -1, step, ".pgm", filename );

        tmp.save ( filename.c_str() );
        if ( !_quiet ) cerr << endl << aol::color::reset;
      }  else {
        cerr << aol::color::red << "\n\nWARNING: saveTimestepPGM not implemented yet in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  void saveTimestepPNG ( const int step,
                         const aol::Vector<RealType> &img,
                         const qc::GridStructure &grid,
                         const RealType Center,
                         const RealType Radius ) const {
    saveTimestepPPM ( -1, step, img, grid, Center, Radius );
  }

  void saveTimestepPNG ( const int outerStep,
                         const int step,
                         const aol::Vector<RealType> &img,
                         const qc::GridStructure &grid,
                         const RealType Center,
                         const RealType Radius ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( grid );
        tmp = img;
        string filename;
        initSaveNameStringMinMaxAndSendMessageToConsole ( outerStep, step, img, ".png", filename );

        const RealType minVal = Center - Radius;
        const RealType maxVal = Center + Radius;
        tmp.clamp ( minVal, maxVal );
        qc::MultiArray<unsigned char, qc::QC_2D, 3> bufArray ( grid );
        qc::convertScalarToColorUsingcolorTrans<RealType, qc::QC_2D> ( tmp, bufArray, aol::HSV_BLUE_TO_RED, minVal, maxVal );
        bufArray.savePNG ( filename.c_str() );

        if ( !_quiet ) cerr << endl << aol::color::reset;
      }  else {
        cerr << aol::color::red<< "\n\nWARNING: saveTimestepPPM not implemented yet in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  //! this method transforms the values of img, such that the full range of [0,1] is filled
  //! and saves these values multiplied by 255 to a pgm files.
  void saveTimestepPGMTransformTo01 ( const int step,
                                      const aol::Vector<RealType> &img,
                                      const qc::GridStructure &grid ) const {
    saveTimestepPGMTransformTo01 ( -1, step, img, grid );
  }

  void saveTimestepPGMTransformTo01 ( const int outerStep,
                                      const int step,
                                      const aol::Vector<RealType> &img,
                                      const qc::GridStructure &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( grid );
        tmp = img;
        string filename;
        initSaveNameStringMinMaxAndSendMessageToConsole ( outerStep, step, img, ".pgm", filename );

        tmp.addToAll ( - tmp.getMinValue() );
        tmp /= tmp.getMaxValue();
        tmp.setQuietMode ( true );
        tmp.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );

        tmp.save ( filename.c_str(), qc::PGM_UNSIGNED_CHAR_BINARY );
        if ( !_quiet ) cerr << endl << aol::color::reset;
      }  else {
        cerr << aol::color::red<< "\n\nWARNING: saveTimestepPGMTransformTo01 not implemented yet in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  void saveTimestepPNGTransformTo01 ( const int step,
                                      const aol::Vector<RealType> &img,
                                      const qc::GridStructure &grid ) const {
    saveTimestepPNGTransformTo01 ( -1, step, img, grid );
  }

  void saveTimestepPNGTransformTo01 ( const int outerStep,
                                      const int step,
                                      const aol::Vector<RealType> &img,
                                      const qc::GridStructure &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( grid );
        tmp = img;
        string filename;
        initSaveNameStringMinMaxAndSendMessageToConsole ( outerStep, step, img, ".png", filename );

        tmp.scaleValuesToAB ( -0.5, 0.5 );
        qc::MultiArray<unsigned char, qc::QC_2D, 3> bufArray ( grid );
        qc::convertScalarToColorUsingcolorTrans<RealType, qc::QC_2D> ( tmp, bufArray, aol::HSV_BLUE_TO_RED );
        bufArray.savePNG ( filename.c_str() );

        if ( !_quiet )
          cerr << endl << aol::color::reset;
      } else {
        cerr << aol::color::red<< "\n\nWARNING: saveTimestepPNGTransformTo01 not implemented yet in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  //! this method expects, that the data is contained in [0,1]. So it will be multiplied by
  //! 255 to receive reasonable values for a pgm-file.
  void saveTimestepPGMMult255 ( const int step,
                                const aol::Vector<RealType> &img,
                                const qc::GridStructure &grid ) const {
    saveTimestepPGMMult255 ( -1, step, img, grid );
  }

  void saveTimestepPGMMult255 ( const int outerStep,
                                const int step,
                                const aol::Vector<RealType> &img,
                                const qc::GridStructure &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( img, grid, aol::DEEP_COPY );
        string filename;
        initSaveNameStringAndSendMessageToConsole ( outerStep, step, ".pgm", filename );

        tmp *= 255.;
        tmp.setQuietMode ( true );
        tmp.save ( filename.c_str(), qc::PGM_UNSIGNED_CHAR_BINARY );
        if ( !_quiet ) cerr << endl << aol::color::reset;
      } else {
        cerr << aol::color::red<< "\n\nWARNING: saveTimestepPGMMult255 not implemented yet in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  //! save 2d- or 3d-image as bz2 (double values) - no outer iteration
  template <typename GridType>
  void saveTimestepBZ2 ( const int step,
                         const aol::Vector<RealType> &img,
                         const GridType &grid ) const {
    saveTimestepBZ2 ( -1, step, img, grid );
  }

  //! save 2d- or 3d-image as bz2 (double values) - with outer iteration
  template <typename GridType>
  void saveTimestepBZ2 ( const int outerStep,
                         const int step,
                         const aol::Vector<RealType> &img,
                         const GridType &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      string filename;
      initSaveNameStringAndSendMessageToConsole ( outerStep, step, ".bz2", filename );
      if ( grid.getDimOfWorld() == 1 ) {
        qc::ScalarArray<RealType, qc::QC_1D> tmp ( img, grid.getNumX() );
        tmp.setQuietMode ( _quiet );
        tmp.save ( filename.c_str(), qc::PGM_DOUBLE_BINARY );
      }
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( img, grid.getNumX(), grid.getNumY() );
        tmp.setQuietMode ( _quiet );
        tmp.save ( filename.c_str(), qc::PGM_DOUBLE_BINARY );
      }
      if ( grid.getDimOfWorld() == 3 ) {
        qc::ScalarArray<RealType, qc::QC_3D> tmp ( img, grid.getNumX(), grid.getNumY(), grid.getNumZ() );
        tmp.setQuietMode ( _quiet );
        tmp.save ( filename.c_str(), qc::PGM_DOUBLE_BINARY );
      }
      if ( !_quiet ) cerr << endl << aol::color::reset;
    }
  }

  //! save 2d- or 3d-image as bz2 (double values) - no outer iteration
  void saveTimestepMinMaxBZ2 ( const int step,
                         const aol::Vector<RealType> &img,
                         const qc::GridStructure &grid,
                         const char *baseSaveName ) const {
    saveTimestepBZ2 ( -1, step, img, grid, baseSaveName );
  }

  //! save 2d- or 3d-image as bz2 (double values) - with outer iteration
  void saveTimestepMinMaxBZ2 ( const int outerStep,
                         const int step,
                         const aol::Vector<RealType> &img,
                         const qc::GridStructure &grid,
                         const char *baseSaveName ) const {
    if  ( checkSaveConditions ( step ) ) {
      string filename;
      initSaveNameStringMinMaxAndSendMessageToConsole ( baseSaveName, outerStep, step, img, ".bz2", filename );
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( img, grid );
        tmp.setQuietMode ( _quiet );
        tmp.save ( filename.c_str(), qc::PGM_DOUBLE_BINARY );
      }
      if ( grid.getDimOfWorld() == 3 ) {
        qc::ScalarArray<RealType, qc::QC_3D> tmp ( img, grid );
        tmp.setQuietMode ( _quiet );
        tmp.save ( filename.c_str(), qc::PGM_DOUBLE_BINARY );
      }
      if ( !_quiet ) cerr << endl << aol::color::reset;
    }
  }

  template<typename GridType>
  void saveTimestepBZ2 ( const int step,
                         const aol::MultiVector<RealType> &img,
                         const GridType &grid ) const {
    for( int i = 0; i < img.numComponents(); i++ )
      saveTimestepBZ2 ( step, i, img[i], grid );
  }

  void saveTimestepColorWheel ( const int step,
                                const aol::Vector<RealType> &img,
                                const qc::GridStructure &grid ) const {
    saveTimestepColorWheel ( -1, step, img, grid );
  }

  void saveTimestepColorWheel ( const int outerStep,
                                const int step,
                                const aol::Vector<RealType> &img,
                                const qc::GridStructure &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {

        qc::ScalarArray<RealType, qc::QC_2D> tmp ( img, grid );

        qc::ColorWheel<RealType> colorWheel ( tmp, 0., 1., true );
        string filename;
        initSaveNameStringMinMaxAndSendMessageToConsole ( outerStep, step, img, ".ppm", filename );
        colorWheel.saveColoredWheel ( filename.c_str() );
        initSaveNameStringMinMaxAndSendMessageToConsole ( outerStep, step, img, ".pgm", filename );
        colorWheel.saveGrayWheel ( filename.c_str() );

      }
      if ( grid.getDimOfWorld() == 3 ) {
        cerr << aol::color::red<< "\n\nWARNING: saveTimestepColorWheel doesn't make sense in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  void saveTimestepColorWheel ( const int step,
                                const aol::Vector<RealType> &img,
                                const qc::GridStructure &grid,
                                const char *baseSaveName ) const {
                                  saveTimestepColorWheel ( -1, step, img, grid, baseSaveName );
                                }

  void saveTimestepColorWheel ( const int outerStep,
                                const int step,
                                const aol::Vector<RealType> &img,
                                const qc::GridStructure &grid,
                                const char *baseSaveName ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( img, grid );
        aol::SimpleFormat valueFormat ( 1, 4, ios::fixed | ios::right, '0' );

        string filename;
        initSaveNameStringMinMaxAndSendMessageToConsole ( baseSaveName, outerStep, step, img, ".pgm", filename );

        qc::ColorWheel<RealType> colorWheel ( tmp, 0., 1., true );
        colorWheel.saveGrayWheel ( filename.c_str() );

      }
      if ( grid.getDimOfWorld() == 3 ) {
        cerr << aol::color::red<< "\n\nWARNING: saveTimestepColorWheel doesn't make sense in 3D!!!\n\n" << aol::color::reset;
      }
    }
  }

  //! save 2d-pic as pgm (short values = 0..255) with marked isoline
  //! only for pics with dimension (2^n + 1)^2
  void saveTimestepPGMDrawIsoline ( const int step,
                                    const double isovalue,
                                    const qc::ScalarArray<RealType, qc::QC_2D> &img ) const {
    if  ( checkSaveConditions ( step ) ) {
      // get the grid belonging to the pic and draw the level line
      qc::GridDefinition grid ( qc::logBaseTwo ( img.getNumX() ) );
      //! \todo Check if this is the appropriate configurator.
      typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiType;
      qc::LevelSetDrawer<ConfiType> drawer ( grid );
      qc::ScalarArray<RealType, qc::QC_2D> tmp ( img.getNumX(), img.getNumY() );
      tmp = img;
      drawer.draw ( img, tmp, isovalue );

      string filename;
      initSaveNameStringAndSendMessageToConsole ( -1, step, ".pgm", filename );

      img.save ( filename.c_str() );
      if ( !_quiet ) cerr << endl << aol::color::reset;
    }
  }

  //! save 2d-pic as bz2 (double values) with marked isoline
  void saveTimestepBZ2DrawIsoline ( const int step,
                                    const double isovalue,
                                    const qc::ScalarArray<RealType, qc::QC_2D> &img ) const {
    if  ( checkSaveConditions ( step ) ) {
      // get the grid belonging to the pic and draw the level line
      qc::GridDefinition grid ( qc::logBaseTwo ( img.getNumX() ) );
      //! \todo Check if this is the appropriate configurator.
      typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiType;
      qc::LevelSetDrawer<ConfiType> drawer ( grid );
      qc::ScalarArray<RealType, qc::QC_2D> tmp ( img.getNumX(), img.getNumY() );
      tmp = img;
      drawer.draw ( img, tmp, isovalue );

      string filename;
      initSaveNameStringAndSendMessageToConsole ( -1, step, ".bz2", filename );

      tmp.save ( filename.c_str() , qc::PGM_DOUBLE_BINARY );
      if ( !_quiet ) cerr << aol::color::reset << endl;
    }
  }

  //! save Vector or MultiVector data on arbitrary Quoc grids
  template<typename VectorType, typename GridType>
  void saveTimestepVectorOrMultiVector ( const int outerStep,
                                         const int step,
                                         const VectorType &img,
                                         const GridType &grid ) const {
    if  ( checkSaveConditions ( step ) ) {
      string filename;
      initBaseSaveNameString ( outerStep, step, filename );
      qc::writeImage<RealType> ( grid, img, filename.c_str() );

      if ( !_quiet ) cerr << aol::color::reset << endl;
    }
  }

  template<typename VectorType, typename GridType>
  void saveTimestepVectorOrMultiVector ( const int step,
                                         const VectorType &img,
                                         const GridType &grid ) const {
    saveTimestepVectorOrMultiVector ( -1, step, img, grid );
  }

  void savePNG ( const aol::Vector<RealType> &img,
                 const qc::GridStructure &grid,
                 const int step = -1,
                 const int outerStep = -1 ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> tmp ( grid );
        tmp = img;
        string filename;
        initSaveNameStringAndSendMessageToConsole ( outerStep, step, ".png", filename );
        tmp.setQuietMode ( true );
        tmp.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
        tmp.savePNG ( filename.c_str() );
        if ( !_quiet ) cerr << endl << aol::color::reset;
      }  else {
        cerr << aol::color::red<< "\n\nWARNING: You can't write a 3D array to PNG!\n\n" << aol::color::reset;
      }
    }
  }

  void savePNG ( const qc::MultiArray<RealType, 2, 3> &VImage,
                 const int step = -1,
                 const int outerStep = -1,
                 const int suffixIndex = -1 ) const {
    if  ( checkSaveConditions ( step ) ) {
      string filename;
      if ( suffixIndex == -1 ) {
        initSaveNameStringAndSendMessageToConsole ( outerStep, step, ".png", filename );
      }
      else {
        string temp = strprintf( "_%03d.png", suffixIndex );
        initSaveNameStringAndSendMessageToConsole ( outerStep, step, temp.c_str(), filename );
      }

      qc::MultiArray<RealType, 2, 3> vImageTemp ( VImage, aol::FLAT_COPY );
      vImageTemp.setQuietMode ( true );
      vImageTemp.setOverflowHandling(aol::CLIP_THEN_SCALE, 0., 1.);
      vImageTemp.savePNG ( filename.c_str() );

      if ( !_quiet ) cerr << endl << aol::color::reset;
    }
  }

  void writeColorField ( const aol::MultiVector< RealType > &vectorfield,
                         const qc::GridStructure &grid,
                         const int outerStep = -1,
                         const int step = -1 ) const {
    if  ( checkSaveConditions ( step ) ) {
      if ( grid.getDimOfWorld() == 2 ) {
        string filename;
        initSaveNameStringAndSendMessageToConsole ( outerStep, step, ".png", filename );
        qc::Array<RealType> vectorfieldX ( vectorfield[0], grid );
        qc::Array<RealType> vectorfieldY ( vectorfield[1], grid );
        qc::writeColorField<RealType>( vectorfieldX, vectorfieldY, filename, -1. );
      }

      else if ( grid.getDimOfWorld() == qc::QC_3D ) {
        char time_and_suffix[1024];
        qc::Array<RealType> vectorfield_0 ( vectorfield[0], grid );
        qc::Array<RealType> vectorfield_1 ( vectorfield[1], grid );
        qc::Array<RealType> vectorfieldX2d ( grid.getNumX(), grid.getNumY() );
        qc::Array<RealType> vectorfieldY2d ( grid.getNumX(), grid.getNumY() );

        for( int i = 0; i < grid.getNumZ(); i++ ){
          sprintf( time_and_suffix, "_%03d.png", i );
          cerr << endl << time_and_suffix <<endl;
          string filename;
          initSaveNameStringAndSendMessageToConsole ( outerStep, step, time_and_suffix, filename );
          vectorfield_0.getSlice ( qc::QC_Z, i, vectorfieldX2d );
          vectorfield_1.getSlice ( qc::QC_Z, i, vectorfieldY2d );
          qc::writeColorField<RealType>( vectorfieldX2d, vectorfieldY2d, filename, -1. );
        }
      }
      else {
        throw aol::UnimplementedCodeException( "writeColorField not implemented for Dim != 2 and Dim != 3", __FILE__, __LINE__ );
      }
    }
  }

  void writeMultiVectorToTextFile ( const aol::MultiVector<RealType> &MVec,
                                    const char *FileName,
                                    const int Iterations,
                                    const int MaxIterations,
                                    const int StartComponent = 0,
                                    const int EndComponent = -1 ) const{
    char fn[1024];
    sprintf( fn, "%s%s", this->_saveDirectory, FileName);
    ofstream outfile(fn, ofstream::out | ofstream::app);
    aol::SimpleFormat iterationsFormat(aol::countDigitsOfNumber(MaxIterations), 0, ios::fixed | ios::right);
    outfile << iterationsFormat(Iterations);

    int endComponent = MVec.numComponents();
    if ( EndComponent != -1 )
      endComponent = EndComponent;
    for ( int i = StartComponent; i < endComponent; i++ )
      outfile << " " << MVec[i];
    outfile << endl;
    outfile.close();
  }
};

/**
 * \author Berkels
 * \ingroup Utilities
 */
template <typename RealType, typename SaveInputType>
class StepSaverBase : protected aol::TimestepSaver<RealType> {
  int _defaultOuterIteration;
public:
  StepSaverBase ( ) : aol::TimestepSaver<RealType> ( ), _defaultOuterIteration ( -1 ) {
    // On default write all time steps. The default behavior of aol::TimestepSaver is just too confusing.
    setSaveTimestepOffset ( 1 );
  }

  virtual ~StepSaverBase ( ) {}

  using aol::TimestepSaver<RealType>::activateSavingAndConfigure;
  using aol::TimestepSaver<RealType>::setSaveName;
  using aol::TimestepSaver<RealType>::getSaveName;
  using aol::TimestepSaver<RealType>::setStepDigitsFromMaxIter;
  using aol::TimestepSaver<RealType>::setOuterStepDigitsFromMaxIter;
  using aol::TimestepSaver<RealType>::getSaveDirectory;
  using aol::TimestepSaver<RealType>::setSaveTimestepOffset;

  void saveStep ( const SaveInputType &SaveInput, const int Iteration, const char *OverrideBaseSaveName = NULL ) const {
    if ( this->checkSaveConditions ( Iteration ) )
      doSaveStep ( SaveInput, Iteration, OverrideBaseSaveName );
  };

  void setAndCreateSaveDirectory ( const char *SaveDirectory ) {
    string saveDirectory = SaveDirectory;
    aol::makeDirectory( saveDirectory.c_str() );
    this->setSaveDirectory( saveDirectory.c_str() );
  }

  void setSaveDirectory ( const char *SaveDirectory ) {
    string saveDirectory = SaveDirectory;
    // aol::TimestepSaver expects the directory name to be terminated with '/', so make sure it is.
    if ( saveDirectory[saveDirectory.length()-1] != '/' )
      saveDirectory += "/";
    aol::TimestepSaver<RealType>::setSaveDirectory ( saveDirectory.c_str() );
  }

  void initFromParser ( const aol::ParameterParser &Parser, const bool OnlyGetDirectory = false, const char *SaveDirectorySuffix = NULL ) {
    if ( Parser.hasVariable ( "saveDirectory" ) ) {
      string saveDir = Parser.getStringExpandTilde( "saveDirectory" );
      if ( SaveDirectorySuffix ) {
        // The suffix creates a subdirectory, so we have to create the top directory first.
        if ( saveDir.at ( saveDir.length() -1 ) == '/' )
          aol::makeDirectory( saveDir.c_str() );
        saveDir += SaveDirectorySuffix;
      }
      setAndCreateSaveDirectory ( saveDir.c_str() );
      Parser.dumpToFile( "parameter-dump.txt", this->getSaveDirectory() );
      if ( OnlyGetDirectory == false ) {
        this->setNumberSaveFirstPics( Parser.getInt( "numSaveFirst" ) );
        this->setSaveTimestepOffset( Parser.getInt( "saveOffset" ) );
        this->setStepDigitsFromMaxIter( Parser.getInt( "maxiterations" ) );
      }
      this->activateSaving();
    }
    else
      this->deactivateSaving();
  }

  void copyFrom ( const StepSaverBase<RealType, SaveInputType> &Other ) {
    static_cast<aol::TimestepSaver<RealType>&> ( *this ) = Other;
    _defaultOuterIteration = Other._defaultOuterIteration;
  }

  string createSaveName ( const char* Prefix, const char* Suffix, const int Iteration = -1, const char *OverrideBaseSaveName = NULL, const int OuterIteration = -1 ) const {
    string out;
    string saveName = Prefix;
    saveName += OverrideBaseSaveName ? OverrideBaseSaveName : this->_saveNameTrunc;
    this->initBaseSaveNameString ( ( ( OuterIteration != -1 ) ? OuterIteration : _defaultOuterIteration), Iteration, out, saveName.c_str() );
    out += Suffix;
    return out;
  }

  void setDefaultOuterIteration ( const int DefaultOuterIteration ) {
    _defaultOuterIteration = DefaultOuterIteration;
  }

protected:
  virtual void doSaveStep ( const SaveInputType &SaveInput, const int Iteration, const char *OverrideBaseSaveName ) const = 0;
};

}   // end namespace

#endif

