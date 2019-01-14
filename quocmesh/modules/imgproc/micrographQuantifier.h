#ifndef __MICROGRAPHQUANTIFIER_H
#define __MICROGRAPHQUANTIFIER_H

#include <aol.h>
#include <atomFinder.h>
#include <gradientDescent.h>
#include <stats.h>

namespace im {

template <typename RealType>
RealType convertElectronsPerPixelToCoulombPerSquareCentimeter ( const int IncidentElectronsPerPixel, const RealType PixelSizeInPicometres ) {
  return IncidentElectronsPerPixel * 1.6021766208 * 10.0 / ( aol::Sqr<RealType> ( PixelSizeInPicometres ) ); // canceled 1e-19 from elementary charge against 1e-20 conversion factor from pm^2 to cm^2
}

template <typename RealType>
RealType convertAngstromToNanoMeter ( const RealType Angstrom ) {
  return 0.1 * Angstrom;
}

template <typename RealType>
RealType convertAngstromToPicoMeter ( const RealType Angstrom ) {
  return 100 * Angstrom;
}

template <typename RealType>
RealType convertPicoMeterToNanoMeter ( const RealType PicoMeter ) {
  return 1e-3 * PicoMeter;
}
 
template <typename RealType>
void plotNanoMeterScaleOverlay ( const char* BasePath, const qc::GridSize<qc::QC_2D> &SizeInPixels, const RealType PixelSizeInNanoMeters, const int BarLengthInNanoMeters = 1,
                                 const RealType RelativeXOffset = 0.05, const RealType RelativeYOffset = 0.05,
                                 const std::string &Font = "Verdana", const int FontSize = 14, const int LineWidth = 10 ) {
  
  const RealType barLengthInPixels = static_cast<RealType> ( BarLengthInNanoMeters ) / PixelSizeInNanoMeters;
  
  aol::RandomAccessContainer<aol::Vec<2, RealType> > points;
  aol::Vec<2, RealType> point;
  point[0] = RelativeXOffset * SizeInPixels.getNumX ( );
  point[1] = RelativeYOffset * SizeInPixels.getNumY ( );
  points.pushBack ( point );
  point[0] += barLengthInPixels;
  points.pushBack ( point );
  aol::Plotter<RealType> plotter;
  aol::PlotDataFileHandler<RealType> plotHandler;
  plotHandler.generateCurvePlot ( points, false, "", aol::strprintf ( "lc \"white\" lw %d", LineWidth ) );
  plotter.addPlotCommandsFromHandler ( plotHandler );
  plotter.setBackgroundPNGSettings ( SizeInPixels.getNumX ( ), SizeInPixels.getNumY ( ) );
  plotter.setXRange ( 0, SizeInPixels.getNumX ( )-1 );
  plotter.setYRange ( 0, SizeInPixels.getNumY ( )-1 );
  plotter.addText ( aol::strprintf ( "%d nm", BarLengthInNanoMeters ), 0.5 * ( points[0][0] + points[1][0] ), points[0][1] + 5 + FontSize, aol::GNUPLOT_TEXT_CENTER, Font, FontSize, "white" );
  plotter.set_outfile_base_name ( BasePath );
  plotter.genPDF ( );
}

  
/**
 * Quantifies electron micrographs (e.g. STEM images) based on Precision, Fidelity, and (Mis-) Detection fraction.
 *
 * \author Mevenkamp
 * \ingroup ElectronMicroscopy
 */
template <typename _RealType,
          typename _MatrixType = aol::FullMatrix<_RealType>,
          typename _LinearRegressionType = aol::LinearRegressionQR<_RealType>,
          typename _ScalarPictureType = qc::ScalarArray<_RealType, qc::QC_2D>,
          typename _ColoredPictureType = qc::MultiArray<_RealType, qc::QC_2D, 3> >
class MicrographQuantifier {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef _ScalarPictureType PictureType;
  typedef _ColoredPictureType ColoredPictureType;
protected:
  std::string _outputDir, _outputDirGT, _outputDirEstimate;
  bool _quietMode, _diskOutput;
public:
  MicrographQuantifier ( const std::string &OutputDir = "", const bool Quiet = true )
    : _outputDir ( OutputDir ), _quietMode ( Quiet ), _diskOutput ( OutputDir != "" ) {
    setOutputDir ( OutputDir );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Centers,
                                            const aol::MultiVector<RealType> &LatticeVectors,
                                            const RealType PeriodDelta, const RealType AngleDelta,
                                            const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const {
    const RealType periodX = LatticeVectors[0].norm ( );
    const RealType periodY = LatticeVectors[1].norm ( );
    const RealType angleX = atan2 ( LatticeVectors[0][1], LatticeVectors[0][0] ) * 180.0 / aol::NumberTrait<RealType>::pi;
    const RealType angleY = atan2 ( LatticeVectors[1][1], LatticeVectors[1][0] ) * 180.0 / aol::NumberTrait<RealType>::pi;
    return getPrecisions ( Centers, periodX, periodY, PeriodDelta, angleX, angleY, AngleDelta, GnuplotXMatches, GnuplotYMatches );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Centers,
                                            const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                            const RealType AngleX, const RealType AngleY, const RealType AngleDelta,
                                            const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const {
    aol::MultiVector<RealType> distances;
    getInterAtomicDistances ( distances, Centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta, GnuplotXMatches, GnuplotYMatches );
    return getPrecisions ( distances );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const std::string &CentersCSVSrcPath,
                                            const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                            const RealType AngleX, const RealType AngleY, const RealType AngleDelta ) const {
    aol::MultiVector<RealType> centers;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.readCentersFromCSV ( centers, CentersCSVSrcPath );
    return getPrecisions ( centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser,
                                            const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const {
    RealType periodX, periodY, periodDelta, angleX, angleY, angleDelta;
    readPrecisionAnalysisParameters ( Parser, periodX, periodY, periodDelta, angleX, angleY, angleDelta );
    return getPrecisions ( Centers, periodX, periodY, periodDelta, angleX, angleY, angleDelta, GnuplotXMatches, GnuplotYMatches );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const std::string &CentersCSVSrcPath, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.readCentersFromCSV ( centers, CentersCSVSrcPath );
    return getPrecisions ( centers, Parser );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.getAtomPositions ( centers, Data, Parser );
    return getPrecisions ( centers, Parser );
  }
  
  RealType getFidelity ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const;
  
  RealType getFidelity ( const std::string &CentersRefCSVSrcPath, const std::string &CentersEstimateCSVSrcPath ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.readCentersFromCSV ( centersRef, CentersRefCSVSrcPath );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.readCentersFromCSV ( centersEstimate, CentersEstimateCSVSrcPath );
    return getFidelity ( centersRef, centersEstimate );
  }
  
  RealType getFidelity ( const PictureType &Ref, const PictureType &Estimate, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Ref, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getFidelity ( centersRef, centersEstimate );
  }
  
  int getNumCorrespondences ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const;
  
  RealType getDetectionFraction ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
    return static_cast<RealType> ( getNumCorrespondences ( CentersRef, CentersEstimate ) ) / CentersRef.numComponents ( );
  }
  
  RealType getDetectionFraction ( const std::string &CentersRefCSVSrcPath, const std::string &CentersEstimateCSVSrcPath ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.readCentersFromCSV ( centersRef, CentersRefCSVSrcPath );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.readCentersFromCSV ( centersEstimate, CentersEstimateCSVSrcPath );
    return getDetectionFraction ( centersRef, centersEstimate );
  }
  
  RealType getDetectionFraction ( const PictureType &Ref, const PictureType &Estimate, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Ref, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getDetectionFraction ( centersRef, centersEstimate );
  }
  
  RealType getMisdetectionFraction ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
    return 1 - static_cast<RealType> ( getNumCorrespondences ( CentersRef, CentersEstimate ) ) / CentersEstimate.numComponents ( );
  }
  
  RealType getMisdetectionFraction ( const std::string &CentersRefCSVSrcPath, const std::string &CentersEstimateCSVSrcPath ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.readCentersFromCSV ( centersRef, CentersRefCSVSrcPath );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.readCentersFromCSV ( centersEstimate, CentersEstimateCSVSrcPath );
    return getMisdetectionFraction ( centersRef, centersEstimate );
  }
  
  RealType getMisdetectionFraction ( const PictureType &Reference, const PictureType &Estimate, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Reference, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getMisdetectionFraction ( centersRef, centersEstimate );
  }
  
  void saveStatistics ( const std::string &Path, const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers, centersInitial;
    aol::Vector<RealType> intensities, intensitiesInitial;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.analyzeAtoms ( centers, intensities, centersInitial, intensitiesInitial, Data, Parser );
    saveStatistics ( aol::strprintf ( "%s%s_initial.txt", aol::getPath ( Path.c_str ( ) ).c_str ( ), aol::getBaseFileName ( Path.c_str ( ) ).c_str ( ) ), aol::MultiVector<RealType> ( 0, 0 ), centersInitial, Parser );
    saveHistogram ( aol::strprintf ( "%sintensitiesHist_initial.csv", aol::getPath ( Path.c_str ( ) ).c_str ( ) ).c_str ( ), intensitiesInitial, 100, Data.getMinValue ( ), Data.getMaxValue ( ) );
    saveStatistics ( Path, aol::MultiVector<RealType> ( 0, 0 ), centers, Parser );
    saveHistogram ( aol::strprintf ( "%sintensitiesHist.csv", aol::getPath ( Path.c_str ( ) ).c_str ( ) ).c_str ( ), intensities, 100, Data.getMinValue ( ), Data.getMaxValue ( ) );
  }
  
  const std::string getStatistics ( const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.getAtomPositions ( centers, Data, Parser );
    return getStatistics ( aol::MultiVector<RealType> ( 0, 0 ), centers, Parser );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser ) const {
    saveStatistics ( Path, aol::MultiVector<RealType> ( 0 , 0 ), Centers, Parser );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser ) const {
    return getStatistics ( aol::MultiVector<RealType> ( 0 , 0 ), Centers, Parser );
  }
  
  void saveStatistics ( const std::string &Path, const PictureType &Reference, const PictureType &Estimate, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    aol::ParameterParser parser ( Parser );
    if ( Parser.hasVariable ( "gammaGT" ) )
      parser.changeVariableValue ( "gamma", Parser.getDouble ( "gammaGT" ) );
    atomFinder.getAtomPositions ( centersRef, Reference, parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    const RealType psnr = aol::PSNR<RealType> ( Reference, Estimate );
    saveStatistics ( Path, centersRef, centersEstimate, Parser, psnr );
  }
  
  const std::string getStatistics ( const PictureType &Reference, const PictureType &Estimate, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    aol::ParameterParser parser ( Parser );
    if ( Parser.hasVariable ( "gammaGT" ) )
      parser.changeVariableValue ( "gamma", Parser.getDouble ( "gammaGT" ) );
    atomFinder.getAtomPositions ( centersRef, Reference, parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    const RealType psnr = aol::PSNR<RealType> ( Reference, Estimate );
    return getStatistics ( centersRef, centersEstimate, Parser, psnr );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate, const aol::ParameterParser &Parser,
                        const RealType PSNR = aol::NumberTrait<RealType>::NaN ) const {
    ofstream txtFile;
    txtFile.open ( Path.c_str ( ) );
    txtFile << getStatistics ( CentersRef, CentersEstimate, Parser, PSNR );
    txtFile.close ( );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate, const aol::ParameterParser &Parser,
                                    const RealType PSNR = aol::NumberTrait<RealType>::NaN ) const {
    RealType periodX, periodY, periodDelta, angleX, angleY, angleDelta, pixelSizeInPicometres, incidentElectronsPerPixel;
    readPrecisionAnalysisParameters ( Parser, periodX, periodY, periodDelta, angleX, angleY, angleDelta, pixelSizeInPicometres, incidentElectronsPerPixel );
    return getStatistics ( CentersRef, CentersEstimate, periodX, periodY, periodDelta, angleX, angleY, angleDelta, pixelSizeInPicometres, incidentElectronsPerPixel, PSNR );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &Centers,
                        const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0,
                        const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0,
                        const RealType PixelSizeInPicometres = 0, const RealType IncidentElectronsPerPixel = 0,
                        const RealType PSNR = aol::NumberTrait<RealType>::NaN ) const {
    saveStatistics ( Path, aol::MultiVector<RealType> ( 0, 0 ), Centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta, PixelSizeInPicometres, IncidentElectronsPerPixel, PSNR );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &Centers,
                                    const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0,
                                    const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0,
                                    const RealType PixelSizeInPicometres = 0, const RealType IncidentElectronsPerPixel = 0,
                                    const RealType PSNR = aol::NumberTrait<RealType>::NaN ) const {
    return getStatistics ( aol::MultiVector<RealType> ( 0, 0 ), Centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta, PixelSizeInPicometres, IncidentElectronsPerPixel, PSNR );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate,
                        const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0,
                        const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0,
                        const RealType PixelSizeInPicometres = 0, const RealType IncidentElectronsPerPixel = 0,
                        const RealType PSNR = aol::NumberTrait<RealType>::NaN ) const {
    ofstream txtFile;
    txtFile.open ( Path.c_str ( ) );
    txtFile << getStatistics ( CentersRef, CentersEstimate, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta, PixelSizeInPicometres, IncidentElectronsPerPixel, PSNR );
    txtFile.close ( );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate,
                                    const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0,
                                    const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0,
                                    const RealType PixelSizeInPicometres = 0, const RealType IncidentElectronsPerPixel = 0,
                                    const RealType PSNR = aol::NumberTrait<RealType>::NaN ) const {
    aol::MultiVector<RealType> distances;
    aol::Vec2<RealType> precisions, meanDistances;
    aol::Vec2<int> numNeighbors;
    
    if ( PeriodX > 0 && PeriodY > 0 ) {
      getInterAtomicDistances ( distances, CentersEstimate, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
      precisions.set ( getPrecisions ( distances ) );
      meanDistances.set ( distances[0].getMeanValue ( ), distances[1].getMeanValue ( ) );
      numNeighbors.set ( distances[0].size ( ), distances[1].size ( ) );
    }
    
    std::stringstream ss;
    ss << "---------- General ------------" << std::endl;
    if ( aol::isFinite<RealType> ( PSNR ) )
      ss << "PSNR:                              " << PSNR << " dB" << std::endl;
    if ( IncidentElectronsPerPixel > 0 && PixelSizeInPicometres > 0 )
      ss << "Dose:                              " << convertElectronsPerPixelToCoulombPerSquareCentimeter<RealType> ( IncidentElectronsPerPixel, PixelSizeInPicometres ) << " C/cm^2" << std::endl;
    ss << "Number of atoms:                   " << CentersEstimate.numComponents ( ) << std::endl;
    if ( CentersRef.numComponents ( ) > 0 )
      ss << "Number of atoms (ref.):            " << CentersRef.numComponents ( ) << std::endl;
    
    if ( PeriodX > 0 && PeriodY > 0 ) {
      ss << std::endl << std::endl;
      ss << "---------- Precision ----------" << std::endl;
      ss << "Number of horizontal neighbors:    " << numNeighbors[0] << std::endl;
      ss << "Number of vertical neighbors:      " << numNeighbors[1] << std::endl;
      ss << "Mean horizontal atom separation:   " << meanDistances[0] << " px"
         << ( PixelSizeInPicometres > 0 ? aol::strprintf ( "; %f pm", meanDistances[0] * PixelSizeInPicometres ) : "" ) << std::endl;
      ss << "Mean vertical atom separation:     " << meanDistances[1] << " px"
         << ( PixelSizeInPicometres > 0 ? aol::strprintf ( "; %f pm", meanDistances[1] * PixelSizeInPicometres ) : "" ) << std::endl;
      ss << "Horizontal precision:              " << precisions[0] << " px"
         << ( PixelSizeInPicometres > 0 ? aol::strprintf ( "; %f pm", precisions[0] * PixelSizeInPicometres ) : "" ) << std::endl;
      ss << "Vertical precision:                " << precisions[1] << " px"
         << ( PixelSizeInPicometres > 0 ? aol::strprintf ( "; %f pm", precisions[1] * PixelSizeInPicometres ) : "" ) << std::endl;
      ss << "Total precision:                   " << precisions.norm ( ) << " px"
         << ( PixelSizeInPicometres > 0 ? aol::strprintf ( "; %f pm", precisions.norm ( ) * PixelSizeInPicometres ) : "" ) << std::endl;
    }
    
    if ( CentersRef.numComponents ( ) > 0 ) {
      const RealType detectionFraction = getDetectionFraction ( CentersRef, CentersEstimate );
      const RealType misdetectionFraction = getMisdetectionFraction ( CentersRef, CentersEstimate );
      const RealType fidelity = getFidelity ( CentersRef, CentersEstimate );
      
      ss << std::endl << std::endl;
      ss << "---------- Fidelity  ----------" << std::endl;
      ss << "Detection fraction:                " << detectionFraction << std::endl;
      ss << "Misdetection fraction:             " << misdetectionFraction << std::endl;
      ss << "Fidelity:                          " << fidelity << " px"
         << ( PixelSizeInPicometres > 0 ? aol::strprintf ( "; %f pm", fidelity * PixelSizeInPicometres ) : "" ) << std::endl;
    }
    
    return ss.str ( );
  }
  
  void saveQuantificationBinary ( const std::string &Path, const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers, centersInitial;
    aol::Vector<RealType> intensities, intensitiesInitial;
    im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirEstimate, _quietMode );
    atomFinder.analyzeAtoms ( centers, intensities, centersInitial, intensitiesInitial, Data, Parser );
    
    RealType periodX, periodY, periodDelta, angleX, angleY, angleDelta, pixelSizeInPicometres, incidentElectronsPerPixel;
    readPrecisionAnalysisParameters ( Parser, periodX, periodY, periodDelta, angleX, angleY, angleDelta, pixelSizeInPicometres, incidentElectronsPerPixel );
    aol::MultiVector<RealType> distances;
    getInterAtomicDistances ( distances, centers, periodX, periodY, periodDelta, angleX, angleY, angleDelta );
    aol::Vec2<RealType> precisions = getPrecisions ( distances );
    if ( pixelSizeInPicometres > 0 ) precisions *= pixelSizeInPicometres;
    
    aol::Vector<RealType> quantificationVec ( 4 );
    quantificationVec[0] = precisions[0];
    quantificationVec[1] = precisions[1];
    
    if ( Parser.hasVariable ( "groundTruthPath" ) ) {
      PictureType gt ( Parser.getString ( "groundTruthPath" ).c_str ( ) );
      
      aol::MultiVector<RealType> centersRef, centersRefInitial;
      aol::Vector<RealType> intensitiesRef, intensitiesRefInitial;
      im::AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
      atomFinder.analyzeAtoms ( centersRef, intensitiesRef, centersRefInitial, intensitiesRefInitial, gt, Parser );
      
      quantificationVec[2] = getDetectionFraction ( centersRef, centers );
      quantificationVec[3] = getMisdetectionFraction ( centersRef, centers );
    }
    
    quantificationVec.saveToFile ( Path.c_str ( ) );
  }
  
  void saveHistogram ( const std::string &Path, const aol::Vector<RealType> &Intensities,
                       const int NumBins,
                       const RealType Min, const RealType Max ) const {
    std::vector<std::pair<RealType, int> > hist;
    Intensities.createHistogramOfValues ( hist, NumBins, Min, Max );
    std::ofstream txtFile ( Path.c_str ( ) );
    for ( unsigned int i=0; i<hist.size ( ) ; ++i )
      txtFile << hist[i].first << ", " << hist[i].second << std::endl;
    txtFile.close ( );
  }
  
  void setOutputDir ( const std::string &OutputDir ) {
    _outputDir = OutputDir;
    _diskOutput = OutputDir != "";
    
    if ( _diskOutput ) {
      std::stringstream ss;
      ss << OutputDir << "/gt";
      _outputDirGT = ss.str ( );
      aol::makeDirectory ( _outputDirGT.c_str ( ) );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << OutputDir << "/estimate";
      _outputDirEstimate = ss.str ( );
      aol::makeDirectory ( _outputDirEstimate.c_str ( ) );
    }
  }
  
  void setQuietMode ( const bool Quiet = true ) {
    _quietMode = Quiet;
  }
private:
  void readPrecisionAnalysisParameters ( const aol::ParameterParser &Parser,
                                         RealType &PeriodX, RealType &PeriodY, RealType &PeriodDelta,
                                         RealType &AngleX, RealType &AngleY, RealType &AngleDelta ) const {
    RealType pixelSizeInPicometres, incidentElectronsPerPixel;
    readPrecisionAnalysisParameters ( Parser, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta, pixelSizeInPicometres, incidentElectronsPerPixel );
  }
  
  void readPrecisionAnalysisParameters ( const aol::ParameterParser &Parser,
                                         RealType &PeriodX, RealType &PeriodY, RealType &PeriodDelta,
                                         RealType &AngleX, RealType &AngleY, RealType &AngleDelta,
                                         RealType &PixelSizeInPicometres, RealType &IncidentElectronsPerPixel ) const {
    PeriodX = Parser.getDoubleOrDefault ( "periodX", 0 );
    PeriodY = Parser.getDoubleOrDefault ( "periodY", 0 );
    PeriodDelta = Parser.getDoubleOrDefault ( "periodDelta", 0 );
    AngleX = Parser.getDoubleOrDefault ( "angleX", 0 );
    AngleY = Parser.getDoubleOrDefault ( "angleY", 0 );
    AngleDelta = Parser.getDoubleOrDefault ( "angleDelta", 0 );
    PixelSizeInPicometres = Parser.getDoubleOrDefault ( "pixelSizeInPicometres", 0 );
    IncidentElectronsPerPixel = Parser.getDoubleOrDefault ( "incidentElectronsPerPixel", 0 );
  }

public:
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Distances ) const {
    if ( Distances.numComponents() < 2 ) throw aol::Exception ( "Distances must have at least two components (x/y distances)", __FILE__, __LINE__ );
    return aol::Vec2<RealType> ( Distances[0].getStdDev ( ), Distances[1].getStdDev ( ) );
  }

  const aol::Vec2<RealType> getPrecisionsComponentsAveraged ( const aol::MultiVector<RealType> &Distances ) const {
    if ( Distances.numComponents() != 6 )
      throw aol::Exception ( "Unexpected number of components", __FILE__, __LINE__ );

    return aol::Vec2<RealType> ( sqrt( aol::Sqr( Distances[2].getStdDev() ) + aol::Sqr( Distances[4].getStdDev() ) ), sqrt( aol::Sqr( Distances[3].getStdDev() ) + aol::Sqr( Distances[5].getStdDev() ) ) );
  }


  RealType getTotalPrecision ( const aol::MultiVector<RealType> &Distances ) const {
    return getPrecisions ( Distances ).norm ( );
  }
  
  void getInterAtomicDistances ( aol::MultiVector<RealType> &Distances, const aol::MultiVector<RealType> &Centers,
                                 const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                 const RealType AngleX, const RealType AngleY, const RealType AngleDelta,
                                 const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const;
  
  void getInterAtomicDistances ( aol::MultiVector<RealType> &Distances, const aol::MultiVector<RealType> &Centers,
                                const aol::ParameterParser &Parser,
                                const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const {
    RealType periodX, periodY, periodDelta, angleX, angleY, angleDelta;
    readPrecisionAnalysisParameters ( Parser, periodX, periodY, periodDelta, angleX, angleY, angleDelta );
    return getInterAtomicDistances ( Distances, Centers, periodX, periodY, periodDelta, angleX, angleY, angleDelta, GnuplotXMatches, GnuplotYMatches );
  }

private:
  void getCorrespondences ( aol::Vector<short> &Correspondences, const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const;
};
  
  
template <typename RealType>
bool testMicrographQuantifier ( const bool QuietMode = false ) {
  std::cerr << setprecision ( 16 ) << "Testing MicrographQuantifier.." << std::endl;
  const RealType threshold = 1e-12;
  bool allTestsOK = true;
  
  // Simple grid errors (exclusively in either x or y direction)
  const int nx = 11, ny = 11, n = nx * ny;
  aol::MultiVector<RealType> atomCentersReferenceX ( n, 2 ), atomCentersReferenceY ( n, 2 );
  const RealType epsX = 0.1, epsY = 0.2;
  for ( int y=0; y<ny ; ++y ) {
    for ( int x=0; x<nx ; ++x ) {
      atomCentersReferenceX[y * nx + x][0] = x + ( ( x % 2 == 0 ) ? 1 : -1 ) * epsX;
      atomCentersReferenceX[y * nx + x][1] = y;
      atomCentersReferenceY[y * nx + x][0] = x;
      atomCentersReferenceY[y * nx + x][1] = y + ( ( y % 2 == 0 ) ? 1 : -1 ) * epsY;
    }
  }
  aol::MultiVector<RealType> atomCenters ( atomCentersReferenceX );
  for ( int i=0; i<n ; ++i ) atomCenters[i][0] += epsX;  // shift all atoms by epsX                       ( fidelity )
  atomCenters[0][0] = -10;                               // make the first ground truth atom misdetected  ( detection fraction & misdetection fraction )
  atomCenters.resize ( n+1, 2 );                         // add a wrongly detected atom                   ( misdetection fraction )
  atomCenters[n][0] = -10;
  atomCenters[n][1] = -10;
  
  // Random grid errors (for testing total precision)
  aol::RandomGenerator rng;
  aol::MultiVector<RealType> atomCentersReference ( n, 2 );
  const RealType eps = 0.1;
  for ( int y=0; y<ny ; ++y ) {
    for ( int x=0; x<nx ; ++x ) {
      atomCentersReference[y * nx + x][0] = x + rng.rReal<RealType> ( -eps, eps );
      atomCentersReference[y * nx + x][1] = y + rng.rReal<RealType> ( -eps, eps );
    }
  }
  aol::MultiVector<RealType> distances ( 2, 0 );
  for ( int y=0; y<ny ; ++y ) {
    for ( int x=0; x<nx-1 ; ++x ) {
      aol::Vec2<RealType> diff ( atomCentersReference[y * nx + x][0] - atomCentersReference[y * nx + (x+1)][0],
                                 atomCentersReference[y * nx + x][1] - atomCentersReference[y * nx + (x+1)][1] );
      distances[0].pushBack ( diff.norm ( ) );
    }
  }
  for ( int y=1; y<ny ; ++y ) {
    for ( int x=0; x<nx ; ++x ) {
      aol::Vec2<RealType> diff ( atomCentersReference[y * nx + x][0] - atomCentersReference[(y-1) * nx + x][0],
                                 atomCentersReference[y * nx + x][1] - atomCentersReference[(y-1) * nx + x][1] );
      distances[1].pushBack ( diff.norm ( ) );
    }
  }
  const aol::Vec2<RealType> precisionReference = aol::Vec2<RealType> ( distances[0].getStdDev ( ), distances[1].getStdDev ( ) );
  
  MicrographQuantifier<RealType> micrographQuantifier;
  std::cerr << "Testing fidelity.. ";
  const RealType fidelityReference = epsX;
  const RealType fidelity = micrographQuantifier.getFidelity ( atomCentersReferenceX, atomCenters );
  const bool fidelityPassed = aol::Abs<RealType> ( fidelity - fidelityReference ) <= threshold;
  std::cerr << ( fidelityPassed ? "OK." : "FAIL!" ) << std::endl;
  if ( !QuietMode && !fidelityPassed ) {
    std::cerr << "Fidelity (ref):                                   " << fidelityReference << std::endl;
    std::cerr << "Fidelity (MicrographQuantifier):                  " << fidelity << std::endl;
  }
  allTestsOK = allTestsOK && fidelityPassed;
  
  std::cerr << "Testing detection fraction.. ";
  const RealType detectionFractionReference = static_cast<RealType> ( n - 1 ) / static_cast<RealType> ( n );
  const RealType detectionFraction = micrographQuantifier.getDetectionFraction ( atomCentersReferenceX, atomCenters );
  const bool detectionFractionPassed = aol::Abs<RealType> ( detectionFraction - detectionFractionReference ) <= threshold;
  std::cerr << ( detectionFractionPassed ? "OK." : "FAIL!" ) << std::endl;
  if ( !QuietMode && !detectionFractionPassed ) {
    std::cerr << "Detection fraction (ref):                         " << detectionFractionReference << std::endl;
    std::cerr << "Detection fraction (MicrographQuantifier):        " << detectionFraction << std::endl;
  }
  allTestsOK = allTestsOK && detectionFractionPassed;
  
  std::cerr << "Testing misdetection fraction.. ";
  const RealType misdetectionFractionReference = 2.0 / static_cast<RealType> ( n+1 );
  const RealType misdetectionFraction = micrographQuantifier.getMisdetectionFraction ( atomCentersReferenceX, atomCenters );
  const bool misdetectionFractionPassed = aol::Abs<RealType> ( misdetectionFraction - misdetectionFractionReference ) <= threshold;
  std::cerr << ( misdetectionFractionPassed ? "OK." : "FAIL!" ) << std::endl;
  if ( !QuietMode && !misdetectionFractionPassed ) {
    std::cerr << "Misdetection fraction (ref):                      " << misdetectionFractionReference << std::endl;
    std::cerr << "Misdetection fraction (MicrographQuantifier):     " << misdetectionFraction << std::endl;
  }
  allTestsOK = allTestsOK && misdetectionFractionPassed;
  
  std::cerr << "Testing precision.. ";
  const int nxDistances = (nx-1) * ny, nyDistances = (ny-1) * nx;
  const RealType biasFactorX = sqrt ( static_cast<RealType> ( nxDistances ) / static_cast<RealType> ( nxDistances-1 ) );
  const RealType biasFactorY = sqrt ( static_cast<RealType> ( nyDistances ) / static_cast<RealType> ( nyDistances-1 ) );
  const RealType precisionXReference = biasFactorX * 2 * epsX;
  const RealType precisionYReference = biasFactorY * 2 * epsY;
  const aol::Vec2<RealType> precisionX = micrographQuantifier.getPrecisions ( atomCentersReferenceX, 1, 1, 2.1 * epsX, 0, 90, 3 );
  const aol::Vec2<RealType> precisionY = micrographQuantifier.getPrecisions ( atomCentersReferenceY, 1, 1, 2.1 * epsY, 0, 90, 3 );
  const aol::Vec2<RealType> precision = micrographQuantifier.getPrecisions ( atomCentersReference, 1, 1, 2.1 * eps, 0, 90, 15 );
  const bool precisionPassed = aol::Abs<RealType> ( precisionX[0] - precisionXReference ) <= threshold && aol::Abs<RealType> ( precisionX[1] ) <= threshold
                            && aol::Abs<RealType> ( precisionY[1] - precisionYReference ) <= threshold && aol::Abs<RealType> ( precisionY[0] ) <= threshold
                            && aol::Abs<RealType> ( precision[0] - precisionReference[0] ) <= threshold && aol::Abs<RealType> ( precision[1] - precisionReference[1] ) <= threshold;
  std::cerr << ( precisionPassed ? "OK." : "FAIL!" ) << std::endl;
  if ( !QuietMode && !precisionPassed ) {
    std::cerr << "Precision x-errs (ref):                           " << precisionXReference << "; " << 0 << std::endl;
    std::cerr << "Precision x-errs (MicrographQuantifier):          " << precisionX[0] << "; " << precisionX[1] << std::endl;
    std::cerr << "Precision y-errs (ref):                           " << precisionYReference << "; " << 0 << std::endl;
    std::cerr << "Precision y-errs (MicrographQuantifier):          " << precisionY[1] << "; " << precisionY[0] << std::endl;
    std::cerr << "Precision combined-errs (ref):                    " << precisionReference[0] << "; " << precisionReference[1] << std::endl;
    std::cerr << "Precision combined-errs (MicrographQuantifier):   " << precision[0] << "; " << precision[1] << std::endl;
  }
  allTestsOK = allTestsOK && precisionPassed;
  std::cerr << "Tests on MicrographQuantifier: " << ( allTestsOK ? "OK." : "FAIL!" ) << std::endl;
  
  return allTestsOK;
}
  
/**
 * \brief Objective functional to determine the crystal lattice in im::Lattice.
 *
 * \author Berkels
 */
template <typename RealType>
class LatticeEnergy : public aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > {
  const aol::RandomAccessContainer<aol::Vec<2, RealType> > &_points;
  const aol::RandomAccessContainer<aol::Vec3<int> > &_pointLatticeIndices;
public:
  LatticeEnergy ( const aol::RandomAccessContainer<aol::Vec<2, RealType> > &Points,
                  const aol::RandomAccessContainer<aol::Vec3<int> > &PointLatticeIndices )
    : _points ( Points ), _pointLatticeIndices ( PointLatticeIndices ) { }

  static void MVecToMotifAndDirection ( const aol::MultiVector<RealType> &MArg,
                                        aol::RandomAccessContainer<aol::Vec<2, RealType> > &Motif,
                                        aol::Matrix22<RealType> &LatticeDirections ) {
    const int motifSize = MArg.numComponents() - 2;
    Motif.reallocate ( motifSize );
    for ( int i = 0; i < 2; ++i ) {
      for ( int k = 0; k < motifSize; ++k )
        Motif[k][i] = MArg[k][i];
      for ( int j = 0; j < 2; ++j )
        LatticeDirections[i][j] = MArg[i+motifSize][j];
    }
  }

  static void MotifAndDirectionToMVec ( const aol::RandomAccessContainer<aol::Vec<2, RealType> > &Motif,
                                        const aol::Matrix22<RealType> &LatticeDirections,
                                        aol::MultiVector<RealType> &MArg ) {
    const int motifSize = Motif.size();
    MArg.reallocate ( 2 + motifSize, 2 );
    for ( int i = 0; i < 2; ++i ) {
      for ( int k = 0; k < motifSize; ++k )
        MArg[k][i] = Motif[k][i];
      for ( int j = 0; j < 2; ++j )
        MArg[i+motifSize][j] = LatticeDirections[i][j];
    }
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest  ) const {
    aol::RandomAccessContainer<aol::Vec<2, RealType> > motif;
    aol::Matrix22<RealType> latticeDirections;
    MVecToMotifAndDirection ( MArg, motif, latticeDirections );

    RealType energy = 0;
    for ( int i = 0; i < _points.size(); ++i ) {
      energy += ( motif[_pointLatticeIndices[i][2]]
                  + static_cast<RealType> ( _pointLatticeIndices[i][0] ) * latticeDirections[0]
                  + static_cast<RealType> ( _pointLatticeIndices[i][1] ) * latticeDirections[1]
                  - _points[i] ).normSqr();
    }
    Dest[0] += 0.5 * energy;
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    MDest.setZero();
    aol::RandomAccessContainer<aol::Vec<2, RealType> > motif;
    aol::Matrix22<RealType> latticeDirections;
    MVecToMotifAndDirection ( MArg, motif, latticeDirections );
    const int motifSize = motif.size();

    for ( int i = 0; i < _points.size(); ++i ) {
      aol::Vec2<RealType> diff ( motif[_pointLatticeIndices[i][2]]
                                 + static_cast<RealType> ( _pointLatticeIndices[i][0] ) * latticeDirections[0]
                                 + static_cast<RealType> ( _pointLatticeIndices[i][1] ) * latticeDirections[1]
                                 - _points[i] );
      for ( int j = 0; j < 2; ++j ) {
        MDest[_pointLatticeIndices[i][2]][j] += diff[j];
        for ( int k = 0; k < 2; ++k )
          MDest[motifSize+k][j] += _pointLatticeIndices[i][k] * diff[j];
      }
    }
  }
};

/**
 * \brief Estimates the underlying crytal lattice from the precison matches 
 *        output by im::MicrographQuantifier.
 *
 * \author Berkels
 */
template <typename RealType>
class Lattice {
  aol::RandomAccessContainer<aol::Vec<2, RealType> > _points;
  aol::RandomAccessContainer<aol::Vec3<int> > _pointLatticeIndices;
  aol::Matrix22<RealType> _latticeDirections;
  aol::RandomAccessContainer<aol::Vec<2, RealType> > _motif;
  aol::Matrix22<RealType> _bbox;

  void addPoint ( const aol::Vec<2, RealType> &Point, const RealType MinDistance = 0.01 ) {
    for ( int i = 0; i < _points.size(); ++i ) {
      if ( ( _points[i] - Point ).norm() < MinDistance )
        return;
    }
    _points.pushBack ( Point );
  }

  void addPointsFromPrecMatch ( const aol::MultiVector<RealType> &PrecMatches ) {
    for ( int i = 0; i < PrecMatches.numComponents(); ++i ) {
      aol::Vec2<RealType> point ( PrecMatches[i][0], PrecMatches[i][1] );
      addPoint ( point );
      for ( int j = 0; j < 2; ++j )
        point[j] += PrecMatches[i][j+2];
      addPoint ( point );
    }
  }

  aol::Matrix22<RealType> getBoundingBox ( ) const {
    aol::Matrix22<RealType> bbox;
    bbox[0] = bbox[1] = _points[0];

    for ( int i = 1; i < _points.size(); ++i ) {
      for ( int j = 0; j < 2; ++j ) {
        if ( _points[i][j] < bbox[0][j] )
          bbox[0][j] = _points[i][j];
        if ( _points[i][j] > bbox[1][j] )
          bbox[1][j] = _points[i][j];
      }
    }
    return bbox;
  }

  aol::Vec2<int> estimateLatticeLoopBounds ( ) const {
    RealType bboxInfWidth = aol::Max ( _bbox[1][0] - _bbox[0][0], _bbox[1][1] - _bbox[0][1]  );
    return aol::Vec2<int> ( std::ceil ( 0.5 * bboxInfWidth / _latticeDirections[0].norm() ),
                            std::ceil ( 0.5 * bboxInfWidth / _latticeDirections[1].norm() ) );
  }

  int getCentralPointIndex ( ) const {
    aol::Matrix22<RealType> bbox = getBoundingBox ( );
    aol::Vec2<RealType> center ( ( bbox[0] + bbox[1] ) / 2 );

    int centralPointIndex = 0;
    RealType distToCenter = ( _points[0] - center ).norm();

    for ( int i = 1; i < _points.size(); ++i ) {
      RealType curDistToCenter = ( _points[i] - center ).norm();
      if ( curDistToCenter < distToCenter ) {
        distToCenter = curDistToCenter;
        centralPointIndex = i;
      }
    }
    return centralPointIndex;
  }

  aol::Vec2<RealType> getLatticePoint ( const aol::Vec3<int> &Index ) const {
    aol::Vec2<RealType> latticePoint;
    for ( int c = 0; c < 2; ++c )
      latticePoint[c] = _motif[Index[2]][c] + Index[0] * _latticeDirections[0][c] + Index[1] * _latticeDirections[1][c];
    return latticePoint;
  }

  aol::Vec3<int> findLatticeIndex ( aol::Vec<2, RealType> &Point ) const {
    // Extremely naive and non-efficient way to find the closest point on the lattice.
    // Should still be sufficiently fast is the number of points is small.
    RealType distToClosesPoint = ( Point - _motif[0] ).norm();
    aol::Vec3<int> index ( 0, 0, 0 );

    aol::Vec2<int> loopBounds = estimateLatticeLoopBounds ( );
    for ( int i = -loopBounds[0];  i <=  loopBounds[0]; ++i ) {
      for ( int j = -loopBounds[1];  j <=  loopBounds[1]; ++j ) {
        for ( int k = 0; k < _motif.size(); ++k ) {
          aol::Vec2<RealType> latticePoint = getLatticePoint ( aol::Vec3<int> ( i, j, k ) );
          RealType curDist = ( Point - latticePoint ).norm();
          if ( curDist < distToClosesPoint ) {
            distToClosesPoint = curDist;
            index[0] = i;
            index[1] = j;
            index[2] = k;
          }
        }
      }
    }
    return index;
  }

  //! Very simple and error prone way to estimate the lattice direction by checking
  //! the closest points to the center. Only works if the crystal motif contians just
  //! one element.
  void estimateLatticeDirections ( ) {
    const aol::Vec<2, RealType> &center = _motif[0];
    std::vector<std::pair<RealType, int> > distances;
    for ( int i = 0; i < _points.size(); ++i ) {
      distances.push_back ( std::pair<RealType, int> ( ( _points[i] - center ).norm() , i) );
    }
    std::sort( distances.begin(), distances.end() );

    _latticeDirections[0] = _points[distances[1].second] - center;
    aol::Vec<2, RealType> dirOne = _latticeDirections[0];
    dirOne.normalize();

    for ( int i = 2; i < _points.size(); ++i ) {
      _latticeDirections[1] = _points[distances[i].second] - center;
      aol::Vec<2, RealType> dirTwo = _latticeDirections[1];
      dirTwo.normalize();
      cerr << i << " " << dirOne * dirTwo << endl;
      if ( aol::Abs ( dirOne * dirTwo ) < 0.5 )
        break;
    }
  }

  void estimateMotif ( ) {
    aol::Vec<2, RealType> center = _points[getCentralPointIndex()];

    const aol::Vec<2, RealType> paraMeanVec = 0.5 * ( _latticeDirections[0] + _latticeDirections[1] );
    aol::Parallelogram<RealType> para ( center - 0.1 * paraMeanVec , _latticeDirections[0], _latticeDirections[1] );
    for ( int i = 0; i < _points.size(); ++i ) {
      if ( para.contains ( _points[i] ) ) {
        _motif.pushBack ( _points[i] );
      }
    }
  }

  void findLatticeIndices ( ) {
    _pointLatticeIndices.reallocate ( _points.size() );
    for ( int i = 0; i < _points.size(); ++i )
      _pointLatticeIndices[i] = findLatticeIndex ( _points[i] );
  }

public:
  Lattice ( const char *PrecXMatchFilename, const char *PrecYMatchFilename ) {
    aol::MultiVector<RealType> precMatches;
    readPrecMatchFile ( PrecXMatchFilename, precMatches );
    _latticeDirections[0] = getMeanShiftFromPrecMatch ( precMatches );
    addPointsFromPrecMatch ( precMatches );
    readPrecMatchFile ( PrecYMatchFilename, precMatches );
    _latticeDirections[1] = getMeanShiftFromPrecMatch ( precMatches );
    addPointsFromPrecMatch ( precMatches );

    if ( _points.size() == 0 ) {
      cerr << "Lattice: No points found, no lattice can be detected.\n";
      return;
    }

    estimateMotif();
    if ( _motif.size() == 0 ) {
      cerr << "Lattice: No motif found, no lattice can be detected.\n";
      return;
    }
    _bbox = getBoundingBox ( );

    findLatticeIndices ( );

    aol::MultiVector<RealType> mtmp;
    LatticeEnergy<RealType>::MotifAndDirectionToMVec ( _motif, _latticeDirections, mtmp );
    LatticeEnergy<RealType> E ( _points, _pointLatticeIndices );
    aol::DerivativeWrapper<RealType, LatticeEnergy<RealType>, aol::MultiVector<RealType> > DE ( E );
    aol::GridlessGradientDescent<RealType, aol::MultiVector<RealType> > gradDescent ( E, DE, 1000 );
    gradDescent.applySingle ( mtmp );
    LatticeEnergy<RealType>::MVecToMotifAndDirection ( mtmp, _motif, _latticeDirections );

    /*
    aol::FirstDerivativeValidator<aol::MultiVector<RealType> > tester ( E, DE, 1., aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.001 );
    tester.testAllDirections ( mtmp, "test/test" );
    */

    cerr << _latticeDirections << endl;
  }

  aol::Vec2<RealType> getAvgDist ( const bool SaveDatFiles = false, const char *SaveDirectory = NULL ) {
    if ( _points.size() == 0 ) {
      cerr << "Lattice::getAvgDist: No points found.\n";
      return aol::Vec2<RealType> ( aol::NumberTrait<RealType>::NaN );
    }
    
    if ( _points.size() != _pointLatticeIndices.size() ) {
      cerr << "Lattice::getAvgDist: Invalid number of point lattice indices.\n";
      return aol::Vec2<RealType> ( aol::NumberTrait<RealType>::NaN );
    }

    aol::RandomAccessContainer<aol::Vec<2, RealType> > diffVecs ( _points.size() );
    aol::Vec2<RealType> avgDist;
    for ( int i = 0; i < _points.size(); ++i ) {
      diffVecs[i] = ( static_cast<aol::Vec<2, RealType> >( getLatticePoint ( _pointLatticeIndices[i] ) ) - _points[i] );
      for ( int c = 0; c < 2; ++c )
        avgDist[c] += aol::Sqr ( diffVecs[i][c] );
    }
    for ( int c = 0; c < 2; ++c )
      avgDist[c] = sqrt ( avgDist[c] ) / _points.size();

    if ( SaveDatFiles ) {
      const string saveDirectory = SaveDirectory ? SaveDirectory : "";
      aol::Vec2<RealType>::setFormat ( aol::detailedFormat );

      exportPointsToGnuplot ( _points, saveDirectory + "points.dat" );
      exportPointsToGnuplot ( _motif, saveDirectory + "motif.dat" );

      std::ofstream outLat ( ( saveDirectory + "lattice.dat" ).c_str()  );
      aol::Vec2<int> loopBounds = estimateLatticeLoopBounds ( );
      for ( int i = -loopBounds[0];  i <=  loopBounds[0]; ++i ) {
        for ( int j = -loopBounds[1];  j <=  loopBounds[1]; ++j ) {
          for ( int k = 0; k < _motif.size(); ++k ) {
            outLat << getLatticePoint ( aol::Vec3<int> ( i, j, k ) ) << endl;
          }
        }
      }

      std::ofstream outUnitCell ( ( saveDirectory + "unitCell.dat" ).c_str()  );
      outUnitCell << _motif[0] << " " << _latticeDirections[0] << endl;
      outUnitCell << _motif[0] << " " << _latticeDirections[1] << endl;

      std::ofstream outDiffVecs ( ( saveDirectory + "diffVecs.dat" ).c_str()  );
      for ( int i = 0; i < _points.size(); ++i )
        outDiffVecs << getLatticePoint ( _pointLatticeIndices[i] ) << " " << static_cast<RealType> ( 20 ) * diffVecs[i] << endl;
      std::ofstream gnuplotdat ( ( saveDirectory + "latticeDiff.gp" ).c_str()  );
      gnuplotdat << "set terminal postscript eps color\n";
      gnuplotdat << "set output \"" << saveDirectory << "latticeDiff.eps\"\n";
      gnuplotdat << "unset xtics\n";
      gnuplotdat << "unset ytics\n";
      gnuplotdat << "set size square 1, 1.4285715\n";
      gnuplotdat << "plot \"" << saveDirectory << "points.dat\" w p, \"" << saveDirectory << "lattice.dat\" w p, \"" << saveDirectory << "motif.dat\" w p, \"" << saveDirectory << "unitCell.dat\" w vec, \"" << saveDirectory << "diffVecs.dat\" w vec\n";
      gnuplotdat.close();
      aol::runGnuplot ( ( saveDirectory + "latticeDiff.gp" ).c_str() );
    }
    return avgDist;
  }

  static void readPrecMatchFile ( const char *Filename, aol::MultiVector<RealType> &PrecMatches ) {
    std::ifstream in ( Filename );
    if ( !in.good() )
      throw aol::FileException ( aol::strprintf ( "cannot open file %s for input.", Filename ).c_str(), __FILE__, __LINE__ );

    aol::Vector<RealType> precMatchesVec;
    while ( in.good() ) {
      RealType v;
      in >> v;
      if ( in.good() )
        precMatchesVec.pushBack( v );
    }
    PrecMatches.reallocate ( precMatchesVec.size() / 4, 4 );
    PrecMatches.copySplitFrom( precMatchesVec );
  }

  static aol::Vec2<RealType> getMeanShiftFromPrecMatch ( const aol::MultiVector<RealType> &PrecMatches ) {
    aol::Vec2<RealType> meanShift ( 0, 0 );
    for ( int i = 0; i < PrecMatches.numComponents(); ++i ) {
      for ( int j = 0; j < 2; ++j ) {
        meanShift[j] += PrecMatches[i][j+2];
      }
    }
    meanShift /= PrecMatches.numComponents();
    return meanShift;
  }

  static void exportPointsToGnuplot ( const aol::RandomAccessContainer<aol::Vec<2, RealType> > &Points, const string &Filename ) {
    std::ofstream out ( Filename.c_str() );
    for ( int i = 0; i < Points.size(); ++i )
      out << Points[i] << endl;
  }
};

/**
 * \author Berkels
 */
template <typename RealType>
void analyzePrecisionFromCenters ( const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser,
                                   const char* BackgroudnImageFilename, const char *BaseSaveName,
                                   const bool SaveMatchedLattice = false ) {
  const std::string precXMatchFile = string ( BaseSaveName ) + "XMatch.dat";
  const std::string precYMatchFile = string ( BaseSaveName ) + "YMatch.dat";
  const std::string gnuplotdatFile = string ( BaseSaveName ) + "Gnuplot.dat";


  im::MicrographQuantifier<RealType> micrographQuantifier ( "", false );
  aol::MultiVector<RealType> distances;
  micrographQuantifier.getInterAtomicDistances ( distances, Centers, Parser,
                                                precXMatchFile.c_str(),
                                                precYMatchFile.c_str() );

  const aol::Vec2<RealType> prec = micrographQuantifier.getPrecisionsComponentsAveraged ( distances );

  aol::plotPrecisionMatches
    ( gnuplotdatFile, string ( BaseSaveName ) + "Match.eps",
      BackgroudnImageFilename, precXMatchFile, precYMatchFile );

  cerr << "Precision in pixels: " << prec << endl;

  std::ofstream outPrec ( ( string ( BaseSaveName ) + "MatchNum.dat" ).c_str() );
  const aol::Vec2<RealType> precNew = micrographQuantifier.getPrecisions ( distances );
  outPrec << precNew[0] << " " << precNew[1] << " ";
  outPrec << prec[0] << " " << prec[1] << " ";
  im::Lattice<RealType> lattice ( precXMatchFile.c_str(), precYMatchFile.c_str() );
  const aol::Vec2<RealType> avgDist = lattice.getAvgDist ( SaveMatchedLattice, BaseSaveName );
  outPrec << avgDist[0] << " " << avgDist[1] << endl;
  outPrec.close();
}

} // namespace im

#endif
