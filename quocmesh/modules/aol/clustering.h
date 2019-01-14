#ifndef __CLUSTERING_H
#define __CLUSTERING_H

#include <aol.h>
#include <scalarArray.h>
#include <vec.h>
#include <vectorExtensions.h>
#include <randomGenerator.h>
#include <progressBar.h>
#include <eigenvectors.h>


namespace aol {


class KMEANS_INIT_METHOD {
public:
  static const int RNG  = 0;    // Initialize clusters with random pairwise different points from the given data
  static const int KPP  = 1;
  static const int PCA  = 2;    // Initialize clusters with first principal components
  static const int MEC  = 3;    // Initialize clusters with maximum edge-weighted clique within the complete graph spanned by the data with Euclidean distances as edge weights
  static const int MAN  = 4;    // Initialize clusters with manually specified centers
  
  static const int NUM  = 5;
  
  static std::string getIdentifier ( const int InitMethod ) {
    if ( InitMethod == RNG ) return         "rng";
    else if ( InitMethod == KPP ) return    "kpp";
    else if ( InitMethod == PCA ) return    "pca";
    else if ( InitMethod == MEC ) return    "mec";
    else if ( InitMethod == MAN ) return    "man";
    else throw aol::Exception ( "Did not recognize initialization method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int InitMethod ) {
    if ( InitMethod == RNG ) return         "Random";
    else if ( InitMethod == KPP ) return    "K-means++";
    else if ( InitMethod == PCA ) return    "Principal Component Analysis";
    else if ( InitMethod == MEC ) return    "Maximum edge-weighted clique";
    else if ( InitMethod == MAN ) return    "Manual";
    else throw aol::Exception ( "Did not recognize initialization method!", __FILE__, __LINE__ );
  }
};


/**
 * \brief Provides a class for k-means clustering
 * 
 * The class can be applied to a vector of scalars and find k scalar clusters
 * It can also be applied to a multivector and find k vector-valued clusters
 *
 * \author Mevenkamp
 */
template <typename _RealType>
class KMeansClusterer {
  typedef _RealType RealType;
protected:
  aol::RandomGenerator _randomGenerator;
  bool _quietMode;
  aol::MultiVector<RealType> _manuallySpecifiedCenters;
public:
  KMeansClusterer ( ) : _quietMode ( true ) {
//    _randomGenerator.randomize ( );    // if "true" randomness is required, uncomment this line
  }
  
  // Cluster multiple times each time initializing with RNG and return best result (without cluster labels)
  RealType applyMultipleRNG ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, const short NumClusters,
                              const int NumPasses = 100, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::Vector<int> clusterLabels;
    return applyMultipleRNG ( Input, Clusters, NumClusters, clusterLabels, NumPasses, MaxIt, Epsilon );
  }
  
  // Cluster multiple times each time initializing with RNG and return best result
  RealType applyMultipleRNG ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, const short NumClusters, aol::Vector<int> &ClusterLabels,
                              const int NumPasses = 100, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    RealType leastResidual = aol::NumberTrait<RealType>::Inf, residual;
    aol::MultiVector<RealType> bestClusters ( Input.numComponents ( ), NumClusters );
    aol::Vector<int> bestClusterLabels ( Input[0].size ( ) );
    aol::ProgressBar<> progressBar ( "Clustering", std::cerr );
    if ( !_quietMode ) progressBar.start ( NumPasses );
    
    for ( int i=0; i<NumPasses ; ++i ) {
      residual = apply ( Input, Clusters, NumClusters, ClusterLabels, KMEANS_INIT_METHOD::RNG, MaxIt, Epsilon );
      if ( residual < leastResidual ) {
        leastResidual = residual;
        bestClusters = Clusters;
        bestClusterLabels = ClusterLabels;
      }
      if ( !_quietMode ) progressBar++;
    }
    if ( !_quietMode ) progressBar.finish ( );
    
    Clusters = bestClusters;
    ClusterLabels = bestClusterLabels;
    
    return leastResidual;
  }

  // Clustering of vectors (without output of cluster labels)
  RealType apply ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, const short NumClusters,
                   const int InitMethod = KMEANS_INIT_METHOD::RNG, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::Vector<int> clusterLabels;
    return apply ( Input, Clusters, NumClusters, clusterLabels, InitMethod, MaxIt, Epsilon );
  }
  
  // Clustering of vectors (including output of cluster labels for each input vector)
  RealType apply ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, const short NumClusters, aol::Vector<int> &ClusterLabels,
                   const int InitMethod = KMEANS_INIT_METHOD::RNG, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    const int dim = Input.numComponents ( );
    if ( dim <= 0 ) throw aol::Exception ( "Received empty MultiVector, but expected at least one component!", __FILE__, __LINE__ );
    
    ClusterLabels.resize ( Input[0].size ( ) );
    aol::MultiVector<RealType> clustersNew ( dim, NumClusters );
    initializeClusters ( Input, clustersNew, InitMethod );
    
    // Refine clusters until difference between old and new clusters is less than Epsilon
    aol::MultiVector<RealType> clustersOld ( dim, NumClusters ), clustersDiff ( dim, NumClusters );
    aol::Vector<int> clusterSizes ( NumClusters );
    int numIt = 0;
    RealType residual;
    aol::Vector<RealType> val ( dim );
    do {
      clustersOld = clustersNew;
      
      // Update clusters
      clusterSizes.setZero ( );
      clustersNew.setZero ( );
      for ( int j=0; j<Input[0].size ( ) ; ++j ) {
        Input.getTo ( j, val );
        const int nearestClusterIdx = getNearestClusterIndex ( val, clustersOld );
        ClusterLabels[j] = nearestClusterIdx;
        clusterSizes[nearestClusterIdx]++;
        for ( int c=0; c<dim ; ++c )
          clustersNew[c][nearestClusterIdx] += Input[c][j];
      }
      for ( int c=0; c<dim ; ++c )
        for ( int i=0; i<clustersNew[0].size ( ) ; ++i )
          if ( clusterSizes[i] > 0 ) clustersNew[c][i] /= clusterSizes[i];
      
      // Calculate difference between old and new clusters
      clustersDiff = clustersOld;
      clustersDiff -= clustersNew;
      numIt++;
      
      // Calculate residual
      residual = 0.0;
      for ( int j=0; j<Input[0].size ( ) ; ++j )
        for ( int c=0; c<dim ; ++c )
          residual += aol::Sqr<RealType> ( clustersNew[c][ClusterLabels[j]] - Input[c][j] );
      residual = sqrt ( residual );
      
      if ( !_quietMode ) cerr << "KMeans: Iteration " << numIt << ": " << clustersNew << endl;
    } while ( clustersDiff.norm ( ) > Epsilon && numIt < MaxIt );
    
    Clusters.resize ( Input.numComponents ( ), NumClusters );
    Clusters = clustersNew;
    
    return residual;
  }
  
  // Clustering of scalars (without output of cluster labels)
  RealType apply ( const aol::Vector<RealType> &Input, aol::Vector<RealType> &Clusters, const short NumClusters,
               const int InitMethod = KMEANS_INIT_METHOD::RNG, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::Vector<int> clusterLabels;
    return apply ( Input, Clusters, NumClusters, clusterLabels, InitMethod, MaxIt, Epsilon );
  }
  
  // Clustering of scalars (including output of cluster labels for each input element)
  RealType apply ( const aol::Vector<RealType> &Input, aol::Vector<RealType> &Clusters, const short NumClusters, aol::Vector<int> &ClusterLabels,
               const int InitMethod = KMEANS_INIT_METHOD::RNG, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::MultiVector<RealType> input ( 1, Input.size ( ) );
    input[0] = Input;
    aol::MultiVector<RealType> clusters;
    RealType residual = apply ( input, clusters, NumClusters, ClusterLabels, InitMethod, MaxIt, Epsilon );
    Clusters.resize ( NumClusters );
    Clusters = clusters[0];
    return residual;
  }
  
  void initializeClusters ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, const int InitMethod ) {
    const int numClusters = Clusters[0].size ( );
    const int dim = Input.numComponents ( );
    
    if ( InitMethod == KMEANS_INIT_METHOD::RNG ) {
      // Initialize clusters randomly
      aol::Vector<int> randClusterIndices ( numClusters );
      _randomGenerator.rIntVecPairwiseDifferent ( randClusterIndices, 0, Input[0].size ( ) );
      for ( int c=0; c<dim ; ++c )
        for ( int i=0; i<numClusters; ++i )
          Clusters[c][i] = Input[c][randClusterIndices[i]];
    } else if ( InitMethod == KMEANS_INIT_METHOD::KPP ) {
      aol::Vector<RealType> dSqrs ( Input[0].size ( ), 1.0 );
      RealType dSqrsSum = dSqrs.sum ( );
      int clusterIdx = 0;
      while ( clusterIdx < numClusters ) {
        int i = 0;
        const RealType x = _randomGenerator.rReal<RealType> ( );
        RealType limit = dSqrs[0] / dSqrsSum;
        while ( x > limit ) {
          ++i;
          limit += dSqrs[i] / dSqrsSum;
        }
        for ( int c=0; c<dim ; ++c )
          Clusters[c][clusterIdx] = Input[c][i];
        
        ++clusterIdx;
        
        if ( clusterIdx < numClusters ) {
          // Update dSqrs
          RealType minDistSqr = aol::NumberTrait<RealType>::Inf, distSqr;
          for ( int i=0; i<Input[0].size ( ) ; ++i ) {
            for ( int j=0; j<clusterIdx ; ++j ) {
              distSqr = 0;
              for ( int c=0; c<dim ; ++c ) distSqr += aol::Sqr<RealType> ( Clusters[c][j] - Input[c][i] );
              if ( distSqr < minDistSqr ) minDistSqr = distSqr;
            }
            dSqrs[i] = minDistSqr;
          }
          dSqrsSum = dSqrs.sum ( );
        }
      }
    } else if ( InitMethod == KMEANS_INIT_METHOD::PCA ) {
      throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
    } else if ( InitMethod == KMEANS_INIT_METHOD::MEC ) {
      // Initialize clusters with maximum edge-weighted clique within graph spanned by the data points with Euclidean distances as edge weights
      // Exact solution is NP hard to obtain, thus approximate with O(N) greedy algorithm
      // Step 1: Add point with largest norm
      RealType maxNorm = 0;
      int maxIdx = 0;
      aol::Vector<RealType> vec ( dim );
      for ( int i=0; i<Input[0].size ( ) ; ++i ) {
        Input.getTo ( i, vec );
        if ( vec.norm ( ) > maxNorm ) {
          maxIdx = i;
          maxNorm = vec.norm ( );
        }
      }
      for ( int c=0; c<dim ; ++c ) Clusters[c][0] = Input[c][maxIdx];
      // Step 2: Iteratively add points that result in the MEC for the current size, i.e. the points with largest summed distance to all points already in the clique
      for ( int k=1 ; k<numClusters ; ++k ) {
        RealType maxMinDist = 0;
        int maxIdx = 0;
        aol::Vector<RealType> cluster ( dim );
        for ( int i=0; i<Input[0].size ( ) ; ++i ) {
          aol::Vector<RealType> dists ( k );
          for ( int j=0; j<k ; ++j ) {
            Input.getTo ( i, vec );
            Clusters.getTo ( j, cluster );
            vec -= cluster;
            dists[j] = vec.norm ( );
          }
          if ( dists.getMinValue ( ) > maxMinDist ) {
            maxIdx = i;
            maxMinDist = dists.getMinValue ( );
          }
        }
        for ( int c=0; c<dim ; ++c ) Clusters[c][k] = Input[c][maxIdx];
      }
    } else if ( InitMethod == KMEANS_INIT_METHOD::MAN ) {
      if ( _manuallySpecifiedCenters.numComponents ( ) != dim || _manuallySpecifiedCenters[0].size ( ) != numClusters )
        throw aol::Exception ( "Centers were not specified manually or have the wrong dimensions", __FILE__, __LINE__ );
      
      for ( int c=0; c<dim ; ++c )
        for ( int i=0; i<numClusters; ++i )
          Clusters[c][i] = _manuallySpecifiedCenters[c][i];
    } else throw aol::Exception ( "Did not recognize k-means initialization method!", __FILE__, __LINE__ );
  }
  
  int getNearestClusterIndex ( const aol::Vector<RealType> &Value, const aol::MultiVector<RealType> &Clusters ) {
    RealType minDistance = aol::NumberTrait<RealType>::Inf;
    int dim = Clusters.numComponents ( ), nearestClusterIdx = -1;
    for ( int i=0; i<Clusters[0].size ( ); ++i ) {
      RealType distance = 0;
      for ( int c=0; c<dim ; ++c ) distance += aol::Sqr<RealType> ( Clusters[c][i] - Value[c] );
      if ( distance < minDistance ) {
        minDistance = distance;
        nearestClusterIdx = i;
      }
    }
    
    if ( nearestClusterIdx < 0 ) throw aol::Exception ( "Distance of specified Value to each Cluster is Infinity", __FILE__, __LINE__ );
    
    return nearestClusterIdx;
  }
  
  void setQuietMode ( bool Quiet ) {
    _quietMode = Quiet;
  }
  
  void setInitialCenters ( const aol::MultiVector<RealType> &Centers ) {
    _manuallySpecifiedCenters.resize ( Centers.numComponents ( ), Centers[0].size ( ) );
    for ( int c=0; c<Centers.numComponents ( ) ; ++c )
      for ( int i=0; i<Centers[0].size ( ); ++i )
        _manuallySpecifiedCenters[c][i] = Centers[c][i];
  }
};



enum PMEANS_CRITERION_TYPE {
  AIC,
  AICc,
  BIC
};

enum PMEANS_DECISION_TYPE {
  MIN_CRITERION,
  MIN_NUMCLUSTERS_SIGNIFICANT_CRITERIA
};

/**
 * \brief Class for k-means clustering with automatic (probabilistic) choice of k (number of clusters)
 *
 * This class runs k-means several times using different k.
 * The clustering that minimizes a given probabilistic criterion is chosen.
 * Available criteria: Akaike's Information Criterion (AIC), AIC with correction for finite sample sizes (AICc) and Bayesian Information Criterion (BIC)
 * Available decision types: (1) minimize criterion, (2) = select smallest k from those that yield a criterion not significantly different from the minimizer
 *
 * \author Mevenkamp
 */
template <typename _RealType>
class PMeansClusterer {
  typedef _RealType RealType;
protected:
  int _maxNumClusters;
  const RealType _significance;
  bool _quietMode;
public:
  PMeansClusterer ( const int MaxNumClusters = 10, const RealType Significance = 0.1, const bool Quiet = true )
    : _maxNumClusters ( MaxNumClusters ), _significance ( Significance ), _quietMode ( Quiet ) { }
  
  // Clustering of vectors (without output of cluster labels)
  RealType apply ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters,
                   const PMEANS_CRITERION_TYPE CriterionType = AIC, const PMEANS_DECISION_TYPE DecisionType = MIN_CRITERION,
                   const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::Vector<int> clusterLabels;
    return apply ( Input, Clusters, clusterLabels, CriterionType, DecisionType, MaxIt, Epsilon );
  }
  
  // Clustering of vectors (including output of cluster labels for each input vector)
  RealType apply ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, aol::Vector<int> &ClusterLabels,
                   const PMEANS_CRITERION_TYPE CriterionType = AIC, const PMEANS_DECISION_TYPE DecisionType = MIN_CRITERION,
                   const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    KMeansClusterer<RealType> kMeansClusterer;
    aol::Vector<RealType> Cs, residuals;
    aol::VectorContainer<aol::MultiVector<RealType> > clustersVec;
    aol::VectorContainer<aol::Vector<int> > clusterLabelsVec;
    for ( int k=1; k<=_maxNumClusters ; ++k ) {
      aol::MultiVector<RealType> clusters;
      aol::Vector<int> clusterLabels;
      RealType residual = kMeansClusterer.applyMultipleRNG ( Input, clusters, k, clusterLabels, MaxIt, Epsilon );
      
      if ( getMinClusterSize ( clusterLabels ) <= 1 ) {
        _maxNumClusters = k-1;
        break;
      }
      
      residuals.pushBack ( residual );
      clustersVec.pushBack ( clusters );
      clusterLabelsVec.pushBack ( clusterLabels );
      Cs.pushBack ( getCriterion ( CriterionType, Input, clusters, clusterLabels ) );
      if ( !_quietMode )
        std::cerr << "#clusters=" << k << "; residual=" << residuals[k-1] << "; C=" << Cs[k-1] << std::endl;
    }
    
    const int kOpt = getOptimalNumClusters ( clustersVec, Cs, DecisionType );
    Clusters.reallocate ( clustersVec[kOpt-1].numComponents ( ), clustersVec[kOpt-1][0].size ( ) );
    Clusters = clustersVec[kOpt-1];
    ClusterLabels.resize ( clusterLabelsVec[kOpt-1].size ( ) );
    ClusterLabels = clusterLabelsVec[kOpt-1];

    if ( !_quietMode )
      std::cerr << "Optimal clustering: #clusters=" << Clusters[0].size ( ) << "; residual=" << residuals[kOpt-1] << "; C=" << Cs[kOpt-1] << std::endl;
    
    return residuals[kOpt-1];
  }
  
  // Clustering of scalars (without output of cluster labels)
  RealType apply ( const aol::Vector<RealType> &Input, aol::Vector<RealType> &Clusters,
                   const PMEANS_CRITERION_TYPE CriterionType = AIC, const PMEANS_DECISION_TYPE DecisionType = MIN_CRITERION,
                   const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::Vector<int> clusterLabels;
    return apply ( Input, Clusters, clusterLabels, CriterionType, DecisionType, MaxIt, Epsilon );
  }
  
  // Clustering of scalars (including output of cluster labels for each input element)
  RealType apply ( const aol::Vector<RealType> &Input, aol::Vector<RealType> &Clusters, aol::Vector<int> &ClusterLabels,
                   const PMEANS_CRITERION_TYPE CriterionType = AIC, const PMEANS_DECISION_TYPE DecisionType = MIN_CRITERION,
                   const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    aol::MultiVector<RealType> input ( 1, Input.size ( ) );
    input[0] = Input;
    aol::MultiVector<RealType> clusters;
    RealType residual = apply ( input, clusters, ClusterLabels, CriterionType, DecisionType, MaxIt, Epsilon );
    Clusters.resize ( clusters[0].size ( ) );
    Clusters = clusters[0];
    return residual;
  }
  
  void setQuietMode ( bool Quiet ) {
    _quietMode = Quiet;
  }
protected:
  int getOptimalNumClusters ( aol::VectorContainer<aol::MultiVector<RealType> > &PClusters, aol::Vector<RealType> &Cs,
                              const PMEANS_DECISION_TYPE DecisionType = MIN_CRITERION ) {
    if ( DecisionType == MIN_CRITERION ) {
      return Cs.getMinIndexAndValue ( ).first + 1;
    } else if ( DecisionType == MIN_NUMCLUSTERS_SIGNIFICANT_CRITERIA ) {
      aol::Vector<int> significantNumClusters;
      const RealType Cmin = Cs.getMinValue ( );
      for ( int k=1; k<=PClusters.size ( ) ; ++k ) {
        const RealType relativeLikelihood = exp ( 0.5 * ( Cmin - Cs[k-1] ) );
        if ( relativeLikelihood > _significance )
          significantNumClusters.pushBack ( k );
        if ( !_quietMode )
          std::cerr << "Relative likelihood (k=" << k << ") = " << relativeLikelihood << std::endl;
      }
      return significantNumClusters.getMinValue ( );
    } else throw aol::Exception ( "Could not recognize specified decision type!", __FILE__, __LINE__ );
    return 1;
  }
  
  RealType getCriterion ( const PMEANS_CRITERION_TYPE CriterionType,
                          const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, aol::Vector<int> &ClusterLabels ) {
    if ( CriterionType == AIC ) return getAIC ( Input, Clusters, ClusterLabels );
    else if ( CriterionType == AICc ) return getAICc ( Input, Clusters, ClusterLabels );
    else if ( CriterionType == BIC ) return getBIC ( Input, Clusters, ClusterLabels );
    else throw aol::Exception ( "Could not recognize specified criterion type!", __FILE__, __LINE__ );
    return 0;
  }
  
  RealType getAIC ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, aol::Vector<int> &ClusterLabels ) {
    return -2 * getMaximumLogLikelihood ( Input, Clusters, ClusterLabels ) + 2 * getNumFreeParameters ( Clusters );
  }
  
  RealType getAICc ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, aol::Vector<int> &ClusterLabels ) {
    const int N = Input[0].size ( );
    const int kFree = getNumFreeParameters ( Clusters );
    return getAIC ( Input, Clusters, ClusterLabels ) + 2 * kFree * ( kFree + 1 ) / ( N - kFree - 1 );
  }
  
  RealType getBIC ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, aol::Vector<int> &ClusterLabels ) {
    const int N = Input[0].size ( );
    return -2 * getMaximumLogLikelihood ( Input, Clusters, ClusterLabels ) + getNumFreeParameters ( Clusters ) * log ( N );
  }
  
  RealType getMaximumLogLikelihood ( const aol::MultiVector<RealType> &Input, aol::MultiVector<RealType> &Clusters, aol::Vector<int> &ClusterLabels ) {
    const int N = Input[0].size ( );
    const int K = Clusters[0].size ( );
    const int d = Clusters.numComponents ( );
    
    aol::Vector<int> clusterSizes ( K );
    aol::Vector<RealType> clusterVariances ( K );
    for ( int i=0; i<N ; ++i ) {
      ++clusterSizes[ClusterLabels[i]];
      for ( int c=0; c<d ; ++c )
        clusterVariances[ClusterLabels[i]] += aol::Sqr<RealType> ( Clusters[c][ClusterLabels[i]] - Input[c][i] );
    }
    for ( int j=0; j<K ; ++j )
      clusterVariances[j] /= static_cast<RealType> ( clusterSizes[j] );
    
    RealType maximumLogLikelihood = 0;
    for ( int j=0; j<K ; ++j ) {
      const int Nj = clusterSizes[j];
      maximumLogLikelihood +=  Nj * log ( Nj ) - Nj * log ( N ) - 0.5 * Nj * ( log ( 2 * aol::NumberTrait<RealType>::pi ) + d * log ( clusterVariances[j] ) + 1 );
    }
    
    return maximumLogLikelihood;
  }
  
  int getNumFreeParameters ( aol::MultiVector<RealType> &Clusters ) {
    const int K = Clusters[0].size ( );
    const int d = Clusters.numComponents ( );
    return K - 1 + d * K + K;
  }
  
  int getMinClusterSize ( const aol::Vector<int> &ClusterLabels ) {
    const int N = ClusterLabels.size ( );
    aol::Vector<int> clusterSizes ( ClusterLabels.getMaxValue ( ) + 1 );
    for ( int i=0; i<N ; ++i )
      ++clusterSizes[ClusterLabels[i]];
  
    return clusterSizes.getMinValue ( );
  }
};


/**
 * \author Mevenkamp
 */
template <typename _RealType, typename _PictureType>
class KMeans2DSpatialIntensityClusterer {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  aol::RandomGenerator _randomGenerator;
  bool _quietMode;
public:
  KMeans2DSpatialIntensityClusterer ( ) : _quietMode ( true ) { }

  void apply ( const PictureType &Arg, aol::RandomAccessContainer<aol::Vec2<RealType> > &Dest, const short NumClusters, const int MaxIt = 100, const RealType Epsilon = 0.001 ) {
    // Initialize class member variables
    aol::RandomAccessContainer<aol::Vec2<RealType> > clustersOld ( NumClusters ), clustersNew ( NumClusters );
    aol::Vec2<RealType> clustersPosDiff;
    RealType clustersDiff;
    aol::Vector<RealType> clusterSizes ( NumClusters );

    // Initialize clusters
    aol::MultiVector<int> randomPixels ( 2, NumClusters );
    
//    _randomGenerator.randomize ( );    // if "true" randomness is required, uncomment this line
    aol::Vector<int> min ( 2 ), max ( 2 );
    max[0] = Arg.getNumX ( ), max[1] = Arg.getNumY ( );
    _randomGenerator.rIntMultiVecPairwiseDifferent ( randomPixels, min, max );
    for ( int i=0; i<NumClusters; ++i )
      clustersNew[i].set ( randomPixels[0][i], randomPixels[1][i] );

    int numIt = 0;
    aol::Vec2<short> pos;
    do {
      clustersOld = clustersNew;
      clusterSizes.setZero ( );
      for ( int i=0; i<clustersNew.size ( ) ; ++i ) clustersNew[i].setZero ( );
      for ( short i=0; i<Arg.getNumX ( ) ; ++i ) {
        for ( short j=0; j<Arg.getNumY ( ) ; ++j ) {
          pos.set ( i, j );
          const int nearestClusterIdx = getNearestClusterIndex ( pos, clustersOld, Arg );
          clusterSizes[nearestClusterIdx] += exp ( Arg.get ( pos ) );
          clustersNew[nearestClusterIdx][0] += pos[0] * exp ( Arg.get ( pos ) );
          clustersNew[nearestClusterIdx][1] += pos[1] * exp ( Arg.get ( pos ) );
        }
      }
      clustersDiff = 0;
      for ( int i=0; i<clustersNew.size ( ) ; ++i ) {
        clustersNew[i] /= clusterSizes[i];
        clustersPosDiff = clustersNew[i];
        clustersPosDiff -= clustersOld[i];
        clustersDiff += clustersPosDiff.normSqr ( );
      }
      clustersDiff = sqrt ( clustersDiff );

      numIt++;

      if ( !_quietMode ) std::cerr << "KMeans: Iteration " << numIt << ": " << std::endl;
      for ( int i=0; i<clustersNew.size ( ) ; ++i )
        cerr << clustersNew[i] << endl;
    } while ( clustersDiff > Epsilon && numIt < MaxIt );

    Dest.reallocate ( NumClusters );
    Dest = clustersNew;
  }

private:
  int getNearestClusterIndex ( const aol::Vec2<short> &Pos, const aol::RandomAccessContainer<aol::Vec2<RealType> > &Clusters, const PictureType &Data ) {
    RealType minDistance = aol::NumberTrait<RealType>::Inf;
    int nearestClusterIdx = -1;
    aol::Vec2<RealType> posDiff;
    for ( int i=0; i < Clusters.size ( ); ++i ) {
      posDiff.set ( Pos[0], Pos[1] );
      posDiff -= Clusters[i];
      RealType distance = aol::Abs ( exp ( Data.get ( round ( Clusters[i][0] ), round ( Clusters[i][0] ) ) ) - exp ( Data.get ( Pos ) ) ) * posDiff.norm ( );
      if ( distance < minDistance ) {
        minDistance = distance;
        nearestClusterIdx = i;
      }
    }
    return nearestClusterIdx;
  }
};


class KSUBSPACE_INIT_METHOD {
public:
  static const int RNG  = 0;    // Initialize sub-spaces with random pairwise different means from the given data and random bases
  
  static const int NUM  = 1;
  
  static std::string getIdentifier ( const int InitMethod ) {
    if ( InitMethod == RNG ) return         "rng";
    else throw aol::Exception ( "Did not recognize initialization method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int InitMethod ) {
    if ( InitMethod == RNG ) return         "Random";
    else throw aol::Exception ( "Did not recognize initialization method!", __FILE__, __LINE__ );
  }
};
  
/**
 * \brief Provides a class for k-subspace clustering (with known sub-space dimensions)
 *
 * \author Mevenkamp
 */
template <typename _RealType>
class KSubSpaceClusterer {
  typedef _RealType RealType;
protected:
  aol::RandomGenerator _randomGenerator;
  bool _quietMode;
  std::string _outputDir;
  qc::GridSize<qc::QC_2D> _labelDimensions2D;
public:
  KSubSpaceClusterer ( ) : _quietMode ( true ), _outputDir ( "" ), _labelDimensions2D ( static_cast<short> ( 0 ) ) { }
  
  // Cluster multiple times each time initializing with RNG and return best result
  RealType applyMultipleRNG ( const aol::MultiVector<RealType> &Input, const short NumClusters, const short &SubSpaceDim,
                              aol::Vector<int> &ClusterLabels,
                              const int NumPasses = 100, const int MaxIt = 100 ) {
    aol::Vector<short> subSpaceDims ( NumClusters );
    subSpaceDims.setAll ( SubSpaceDim );
    return applyMultipleRNG ( Input, NumClusters, subSpaceDims, ClusterLabels, NumPasses, MaxIt );
  }
  
  // Cluster multiple times each time initializing with RNG and return best result
  RealType applyMultipleRNG ( const aol::MultiVector<RealType> &Input, const short NumClusters, const aol::Vector<short> &SubSpaceDims,
                              aol::Vector<int> &ClusterLabels,
                              const int NumPasses = 100, const int MaxIt = 100 ) {
    RealType leastResidual = aol::NumberTrait<RealType>::Inf, residual;
    aol::Vector<int> bestClusterLabels ( Input[0].size ( ) );
    aol::ProgressBar<> progressBar ( "Clustering", std::cerr );
    if ( !_quietMode ) progressBar.start ( NumPasses );
    
    const std::string baseOutputDir = _outputDir;
    
    for ( int i=0; i<NumPasses ; ++i ) {
      if ( baseOutputDir != "" ) {
        _outputDir = aol::strprintf ( "%s/pass%d", baseOutputDir.c_str ( ), i );
        aol::makeDirectory ( _outputDir.c_str ( ) );
      }
      
      residual = apply ( Input, NumClusters, SubSpaceDims, ClusterLabels, KMEANS_INIT_METHOD::RNG, MaxIt );
      if ( residual < leastResidual ) {
        leastResidual = residual;
        bestClusterLabels = ClusterLabels;
      }
      if ( !_quietMode ) progressBar++;
    }
    if ( !_quietMode ) progressBar.finish ( );
    
    _outputDir = baseOutputDir;
    
    ClusterLabels = bestClusterLabels;
    
    return leastResidual;
  }
  
  RealType apply ( const aol::MultiVector<RealType> &Input, const short NumClusters, const short &SubSpaceDim,
                   aol::Vector<int> &ClusterLabels,
                   const int InitMethod = KSUBSPACE_INIT_METHOD::RNG, const int MaxIt = 100 ) {
    aol::Vector<short> subSpaceDims ( NumClusters );
    subSpaceDims.setAll ( SubSpaceDim );
    return apply ( Input, NumClusters, subSpaceDims, ClusterLabels, InitMethod, MaxIt );
  }
  
  RealType apply ( const aol::MultiVector<RealType> &Input, const short NumClusters, const aol::Vector<short> &SubSpaceDims,
                   aol::Vector<int> &ClusterLabels,
                   const int InitMethod = KSUBSPACE_INIT_METHOD::RNG, const int MaxIt = 100 ) {
    const int dim = Input.numComponents ( );
    if ( dim <= 0 ) throw aol::Exception ( "Received empty MultiVector, but expected at least one component!", __FILE__, __LINE__ );
    
    aol::MultiVector<RealType> meanValues ( NumClusters, dim );
    aol::RandomAccessContainer<aol::MultiVector<RealType> > subSpaces ( NumClusters );
    for ( int l=0; l<NumClusters ; ++l ) subSpaces[l].reallocate ( SubSpaceDims[l], dim );
    initializeSubSpaces ( Input, meanValues, subSpaces, InitMethod );
    
    return refineSubSpaces ( meanValues, subSpaces, Input, NumClusters, SubSpaceDims, ClusterLabels, MaxIt );
  }
  
  void initializeSubSpaces ( const aol::MultiVector<RealType> &Input,
                             aol::MultiVector<RealType> &MeanValues, aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces,
                             const int InitMethod ) {
    const int numClusters = SubSpaces.size ( );
    const int dim = Input.numComponents ( );
    
    if ( InitMethod == KSUBSPACE_INIT_METHOD::RNG ) {
      // Initialize sub-spaces randomly
      aol::Vector<int> randClusterIndices ( numClusters );
      _randomGenerator.rIntVecPairwiseDifferent ( randClusterIndices, 0, Input[0].size ( ) );
      aol::Vector<RealType> v ( dim );
      for ( int l=0; l<numClusters; ++l ) {
        Input.getTo ( randClusterIndices[l], v );
        MeanValues[l] = v;
        const int subSpaceDim = SubSpaces[l].numComponents ( );
        aol::FullMatrix<RealType> A ( dim, subSpaceDim );
        for ( int k=0; k<subSpaceDim ; ++k )
          for ( int j=0; j<dim ; ++j )
            A.set ( j, k, _randomGenerator.rReal<RealType> ( ) );
        aol::FullMatrix<RealType> Q ( dim, subSpaceDim ), R ( subSpaceDim, subSpaceDim );
        aol::QRDecomposeModifiedGramSchmidt<RealType> qrGramSchmidt;
        qrGramSchmidt.transform ( A, R, Q );
        for ( int k=0; k<subSpaceDim ; ++k )
          for ( int j=0; j<dim ; ++j )
            SubSpaces[l][k][j] = Q.get ( j, k );
      }
    } else throw aol::Exception ( "Did not recognize k-subspaces initialization method!", __FILE__, __LINE__ );
  }
  
  RealType refineSubSpaces ( aol::MultiVector<RealType> &MeanValues, aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces,
                             const aol::MultiVector<RealType> &Input, const short NumClusters, const aol::Vector<short> &SubSpaceDims,
                             aol::Vector<int> &ClusterLabels, const int MaxIt ) {
    const int dim = Input.numComponents ( );
    const int numPoints = Input[0].size ( );
    ClusterLabels.reallocate ( numPoints );
    
    aol::Vector<int> clusterLabelsDiff ( numPoints );
    int numIt = 0;
    RealType residual;
    int numClusterSwitches = 0;
    aol::Vector<RealType> v ( dim );
    do {
      clusterLabelsDiff = ClusterLabels;
      
      // Update cluster labels
      for ( int i=0; i<numPoints ; ++i ) {
        Input.getTo ( i, v );
        ClusterLabels[i] = getNearestClusterIndex ( v, MeanValues, SubSpaces );
      }
      
      updateSubSpaces ( MeanValues, SubSpaces, Input, NumClusters, SubSpaceDims, ClusterLabels );
      
      // Calculate difference between old and new clusters
      clusterLabelsDiff -= ClusterLabels;
      numClusterSwitches = clusterLabelsDiff.numNonZeroes ( );
      numIt++;
      
      // Calculate residual
      residual = 0;
      for ( int i=0; i<numPoints ; ++i ) {
        Input.getTo ( i, v );
        residual += distanceToSubSpace ( v, MeanValues[ClusterLabels[i]], SubSpaces[ClusterLabels[i]] );
      }
      
      if ( !_quietMode ) {
        cerr << "KSubSpaces: Iteration " << numIt << ": #switches = " << numClusterSwitches << endl;
        
        if ( _outputDir != "" && _labelDimensions2D.getNumberOfNodes ( ) == ClusterLabels.size ( ) ) {
          qc::ScalarArray<int, qc::QC_2D> clusterLabels2D ( _labelDimensions2D );
          for ( int i=0; i<ClusterLabels.size ( ) ; ++i ) clusterLabels2D[i] = ClusterLabels[i];
          clusterLabels2D.savePNG ( aol::strprintf ( "%s/labels2D_%d.png", _outputDir.c_str ( ), numIt ).c_str ( )  );
        }
      }
    } while ( numClusterSwitches > 0 && numIt < MaxIt );
    
    return residual;
  }
  
  RealType refineSubSpaces ( aol::MultiVector<RealType> &MeanValues, aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces,
                             const aol::MultiVector<RealType> &Input, const short NumClusters, const short SubSpaceDim,
                             aol::Vector<int> &ClusterLabels, const int MaxIt ) {
    aol::Vector<short> subSpaceDims ( NumClusters );
    subSpaceDims.setAll ( SubSpaceDim );
    
    return refineSubSpaces ( MeanValues, SubSpaces, Input, NumClusters, subSpaceDims, ClusterLabels, MaxIt );
  }
  
  void updateSubSpaces ( aol::MultiVector<RealType> &MeanValues, aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces,
                         const aol::MultiVector<RealType> &Input, const short NumClusters, const aol::Vector<short> &SubSpaceDims,
                         const aol::Vector<int> &ClusterLabels ) {
    const int dim = Input.numComponents ( );
    const int numPoints = Input[0].size ( );
    
    aol::Vector<RealType> v ( dim );
    for ( int l=0; l<NumClusters ; ++l ) {
      // Assemble vectors of l-th cluster
      aol::MultiVector<RealType> segVecs ( dim, numPoints );
      int n = 0;
      for ( int i=0; i<numPoints ; ++i ) {
        if ( ClusterLabels[i] == l ) {
          Input.getTo ( i, v );
          segVecs.set ( n, v );
          ++n;
        }
      }
      segVecs.resize ( dim, n );
      segVecs.getMeanComponents ( MeanValues[l] );
      aol::Vector<RealType> eigenVals;
      aol::getPCAEigenVecs<RealType> ( SubSpaces[l], eigenVals, segVecs, SubSpaceDims[l], 0.0, true );
    }
  }

  void updateSubSpaces ( aol::MultiVector<RealType> &MeanValues, aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces,
                         const aol::MultiVector<RealType> &Input, const short NumClusters, const short SubSpaceDim,
                         const aol::Vector<int> &ClusterLabels ) {
    aol::Vector<short> subSpaceDims ( NumClusters );
    subSpaceDims.setAll ( SubSpaceDim );
    updateSubSpaces ( MeanValues, SubSpaces, Input, NumClusters, subSpaceDims, ClusterLabels );
  }
  
  int getNearestClusterIndex ( const aol::Vector<RealType> &V,
                               const aol::MultiVector<RealType> &MeanValues,
                               const aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces ) {
    RealType minDistance = aol::NumberTrait<RealType>::Inf;
    int nearestClusterIdx = -1;
    for ( int i=0; i<SubSpaces.size ( ); ++i ) {
      const RealType distance = distanceToSubSpace ( V, MeanValues[i], SubSpaces[i] );
      if ( distance < minDistance ) {
        minDistance = distance;
        nearestClusterIdx = i;
      }
    }
    
    if ( nearestClusterIdx < 0 ) throw aol::Exception ( "Distance of specified Value to each Cluster is Infinity", __FILE__, __LINE__ );
    
    return nearestClusterIdx;
  }
  
  RealType distanceToSubSpace ( const aol::Vector<RealType> &V, const aol::Vector<RealType> &MeanValue, const aol::MultiVector<RealType> &SubSpace ) {
    aol::Vector<RealType> v ( V ), proj ( V.size ( ) );
    v -= MeanValue;
    for ( int k = 0; k < SubSpace.numComponents ( ) ; ++k ) proj.addMultiple ( SubSpace[k], v.dotProduct ( SubSpace[k] ) );
    v -= proj;
    
    return v.normSqr ( );
  }
  
  void setLabelDimensions2D ( const int NumX, const int NumY ) {
    _labelDimensions2D[0] = NumX;
    _labelDimensions2D[1] = NumY;
  }
  
  void setOutputDirectory ( const std::string &OutputDir ) {
    _outputDir = OutputDir;
  }
  
  void setQuietMode ( bool Quiet ) {
    _quietMode = Quiet;
  }
  
  void randomize ( ) {
    _randomGenerator.randomize ( );
  }
};
  
} // end namespace


#endif /* CLUSTERING_H_ */
