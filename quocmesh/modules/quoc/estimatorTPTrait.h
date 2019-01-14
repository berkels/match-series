#ifndef __ESTIMATORTPTRAIT_H
#define __ESTIMATORTPTRAIT_H

#include "estimatorTP.h"

// forward declaration
template <typename T> class EstimatorTP2d;
template <typename T> class EstimatorTP3d;

/**
 * Trait to select the proper estimator op based on the Configurator template.
 *
 *  \author Paetz
 */
template <typename T, qc::Dimension Dim>
class EstimatorTPTrait {};

template <typename T>
class EstimatorTPTrait<T, qc::QC_2D> {
public:
  typedef EstimatorTP2d<T> EstimatorType;
};

template <typename T>
class EstimatorTPTrait<T, qc::QC_3D> {
public:
  typedef EstimatorTP3d<T> EstimatorType;
};

#endif
