#ifndef __COLORCONV_H
#define __COLORCONV_H

#include <array.h>
#include <matrix.h>
#include <kernel3d.h>
#include <auxiliary.h>
#include <bzipiostream.h>
#include <multiArray.h>
#include <rectangularGrid.h>
#include <gnuplotter.h>
#include <parameterParser.h>
#include <convolution.h>
#include <rgbColorMap.h>
#include <connectedComponents.h>

namespace im {

/**
 *\brief Convert an RGB image to an HSI image
 * \author Tatano
 */
template<typename RealType>
void convertRGBtoHSI(const qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &inputArray, qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &outputArray){
  RealType rgb[3] = {0, 0, 0}, hsi[3];
  for (int i=0; i<outputArray[0].size(); i++) {
    rgb[0] = inputArray[0][i];
    rgb[1] = inputArray[1][i];
    rgb[2] = inputArray[2][i];
    aol::RGBColorMap<RealType>::rgb2hsi(rgb, hsi);
    outputArray[0][i] = hsi[0];
    outputArray[1][i] = hsi[1];
    outputArray[2][i] = hsi[2];
  }
}

/**
 *\brief Convert an HSI image to an RGB image
 * \author Tatano
 */
template<typename RealType>
void convertHSItoRGB(const qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &inputArray, qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &outputArray){
  RealType rgb[3] = {0, 0, 0}, hsi[3];
  for (int i=0; i<outputArray[0].size(); i++) {
    hsi[0] = inputArray[0][i];
    hsi[1] = inputArray[1][i];
    hsi[2] = inputArray[2][i];
    aol::RGBColorMap<RealType>::hsi2rgb( hsi, rgb);
    outputArray[0][i] = rgb[0];
    outputArray[1][i] = rgb[1];
    outputArray[2][i] = rgb[2];
  }
}

/**
 *\brief Convert an RGB image to a gray scale image using the transformation in Sridhar, Digital Image Processing, OUP, 2011
 * \author Tatano
 */
template <typename RealType>
void convertRGBToGrayScale(const qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &rgbArray, qc::ScalarArray<RealType, qc::QC_2D> &grayArray){
  if (grayArray.size() != rgbArray[0].size()) {
    throw aol::Exception("dimensions mismatch!", __FILE__, __LINE__ );
  }
  RealType rgb[3] = {0, 0, 0};
  for (int i=0; i<grayArray.size(); i++) {
    rgb[0] = rgbArray[0][i];
    rgb[1] = rgbArray[1][i];
    rgb[2] = rgbArray[2][i];
    aol::RGBColorMap<RealType>::rgb2gray(rgb, grayArray[i]);
  }
}
  
/**
 *\brief Convert an RGB image to a YCbCr image using the transformation in Sridhar, Digital Image Processing, OUP, 2011
 * \author Tatano
 */
template<typename RealType>
void convertRGBToYCbCr(const qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &rgbArray, qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &YCbCrArray) {
  if (YCbCrArray[0].size() != rgbArray[0].size()) {
    throw aol::Exception("dimensions mismatch!", __FILE__, __LINE__ );
  }
  RealType rgb[3] = {0, 0, 0}, YCbCr[3];
  for (int i=0; i<rgbArray[0].size(); i++) {
    rgb[0] = rgbArray[0][i];
    rgb[1] = rgbArray[1][i];
    rgb[2] = rgbArray[2][i];
    aol::RGBColorMap<RealType>::rgb2YCbCr(rgb, YCbCr);
    YCbCrArray[0][i] = YCbCr[0];
    YCbCrArray[1][i] = YCbCr[1];
    YCbCrArray[2][i] = YCbCr[2];
  }
}
  
/**
 *\brief Desaturate an RGB image using the transformation in Sridhar, Digital Image Processing, OUP, 2011. Value must be between 0 and 1.
 * \author Tatano
*/
template<typename RealType>
void desaturateImageTo(const qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &inputArray, const RealType value, qc::MultiArray<RealType, qc::QC_2D, qc::QC_3D> &satArray) {
  qc::ScalarArray<RealType, qc::QC_2D> temp(inputArray[0], aol::STRUCT_COPY);
  convertRGBToGrayScale(inputArray, temp);
  for (int i=0; i<temp.size(); i++) {
    satArray[0][i] = temp[i] + value * (inputArray[0][i] - temp[i]);
    satArray[1][i] = temp[i] + value * (inputArray[1][i] - temp[i]);
    satArray[2][i] = temp[i] + value * (inputArray[2][i] - temp[i]);
  }
}

} // end of namespace im.

#endif
