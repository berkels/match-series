#ifndef __CONNECTEDCOMPONENTS_H
#define __CONNECTEDCOMPONENTS_H

#include <aol.h>
#include <array.h>
#include <geom.h>
#include <multiArray.h>
#include <smallVec.h>
#include <vectorExtensions.h>


namespace qc {


/**
 * ConnectedComponentsLabeler is an algorithm that applies the Connected Component Labeling
 * alogrithm to an input qc::BitArray<qc::QC_2D>.
 *
 * \author Neil Brown, DAI
 * \author Judy Robertson, SELLIC OnLine
 *
 * Converted from Java by Berkels, cf.
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/label.htm
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/flatjavasrc/ImageLabel.java
 */
class ConnectedComponentsLabeler {
private:
  /**
   * Associate(equivalence) A with B.
   * A should be less than B to give some ordering (sorting)
   * if B is already associated with some other value, then propagate
   * down the list.
   */
  static void associate( aol::Vector<unsigned int> &Labels, unsigned int A, unsigned int B ) {
    if( A > B ) {
      associate( Labels, B, A );
      return;
    }
    if( ( A == B ) || ( Labels[ B ] == A ) ) return;
    if( Labels[ B ] == B ) {
      Labels[ B ] = A;
    } else {
      associate( Labels, Labels[ B ], A );
      if (Labels[ B ] > A) {             //***rbf new
        Labels[ B ] = A;
      }
    }
  }
  /**
   * Reduces the number of Labels.
   */
  static unsigned int reduce( aol::Vector<unsigned int> &Labels, unsigned int A ){
    if( Labels[ A ] == A ){
      return A;
    } else {
      return reduce( Labels, Labels[ A ] );
    }
  }
public:
  
  /**
   * doLabel applies the Labeling alogrithm
   *
   * \param Mask The input pixel array
   * \param LabelArray A pixel array containing the labelled image
   *
   * Returns the number of labels used.
   *
   * NB For images  0,0 is the top left corner.
   */
  static int doLabel ( const qc::BitArray<qc::QC_2D> &Mask, qc::ScalarArray<int, qc::QC_2D> &LabelArray ) {
    
    unsigned int nextlabel = 1;
    aol::Vec<4, bool> nbs;
    aol::Vec<4, unsigned int> nbls;
    unsigned int result = 0;
    const int numX = Mask.getNumX();
    const int numY = Mask.getNumY();
    // the most labels there can be is 1/2 of the points in checkerboard
    // part of the code assumes that there are at least three possible labels.
    aol::Vector<unsigned int> labels ( aol::Max ( numX * numY / 2, 3 ) );
    //initialise labels
    for (int i=0; i<labels.size(); i++) labels[ i ] = i;
    
    //now Label the image
    for ( int y = 0; y < numY; ++y ) {
      for ( int x = 0; x < numX; ++x ) {
        if ( Mask.get( x, y ) == false ) {
          result = 0;  //nothing here
        } else {
          
          //The 4 visited neighbours
          nbs[ 0 ] = ( x > 0 ) && Mask.get ( x-1, y );
          nbs[ 1 ] = ( y > 0 ) && Mask.get ( x, y-1 );
          nbs[ 2 ] = ( x > 0 ) && ( y > 0 ) && Mask.get ( x-1, y-1 );
          nbs[ 3 ] = ( x < numX - 1 ) && ( y > 0 ) && Mask.get ( x+1, y-1 );
          
          //Their corresponding labels
          nbls[ 0 ] = ( x > 0 ) ? LabelArray.get ( x-1, y ) : 0;
          nbls[ 1 ] = ( y > 0 ) ? LabelArray.get ( x, y-1 ) : 0;
          nbls[ 2 ] = ( ( x > 0 ) && ( y > 0 ) ) ? LabelArray.get ( x-1, y-1 ) : 0;
          nbls[ 3 ] = ( ( x < numX - 1 ) && ( y > 0 ) ) ? LabelArray.get ( x+1, y-1 ) : 0;
          
          //label the point
          if( (nbs[0] == nbs[1]) && (nbs[1] == nbs[2]) && (nbs[2] == nbs[3])
             && (nbs[0] == 0 )) {
            // all neighbours are 0 so gives this point a new label
            result = nextlabel;
            nextlabel++;
          } else { //one or more neighbours have already got labels
            int count = 0;
            int found = -1;
            for( int j=0; j<4; j++){
              if( nbs[ j ] != 0 ){
                count +=1;
                found = j;
              }
            }
            if( count == 1 ) {
              // only one neighbour has a label, so assign the same label to this.
              result = nbls[ found ];
            } else {
              // more than 1 neighbour has a label
              result = nbls[ found ];
              // Equivalence the connected points
              for( int j=0; j<4; j++){
                if( ( nbls[ j ] != 0 ) && (nbls[ j ] != result ) ){
                  associate( labels, nbls[ j ], result );
                }
              }
            }
          }
        }
        LabelArray.set ( x, y, result );
      }
    }
    //reduce labels ie 76=23=22=3 -> 76=3
    //done in reverse order to preserve sorting
    for( unsigned int i= labels.size()-1; i > 0; i-- ){
      labels[ i ] = reduce( labels, i );
    }
    
    /*now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
     this needs to be condensed down again, so that there is no wasted
     space eg in the above, the labels 3 and 4 are not used instead it jumps
     to 5.
     */
    aol::Vector<unsigned int> condensed ( nextlabel ); // cant be more than nextlabel labels
    
    unsigned int count = 0;
    for (unsigned int i=0; i< nextlabel; i++){
      if( i == labels[ i ] ) condensed[ i ] = count++;
    }
    // Record the number of labels
    const int numberOfLabels = count - 1;
    
    // now run back through our preliminary results, replacing the raw label
    // with the reduced and condensed one, and do the scaling and offsets too
    for (int i=0; i< LabelArray.size(); i++)
      LabelArray[i] = condensed[ labels[ LabelArray[ i ] ] ];
    
    return numberOfLabels;
  }
  
  /**
   * Extension of doLabel to 3D. Can also handle periodic domains.
   * \author Paul Springer (AICES)
   * \param[out] labelSize Stores the number of pixels belonging to each connected component. Label-size of label i is stored at position i (i.e., the last label-size is stored at labelSize[numberOfLabels]).
   * \param connectivity 0: 6-way connectivity (i.e. left, right, top, bottom, front, back) or 1: full connectivity(all 26 neighbors (i.e. with diagonal elements))
   * \return numberOfLabels Number of labels found (excluding label 0)
   */
  template<int PBC_X, int PBC_Y, int PBC_Z, int connectivity>
  static int doLabel3D ( const qc::BitArray<qc::QC_3D> &Mask, qc::ScalarArray<int, qc::QC_3D> &LabelArray, aol::Vector<unsigned int> &labelSize ) {
    
    unsigned int nextlabel = 1;
    //unsigned int result = 0;
    const int numX = Mask.getNumX();
    const int numY = Mask.getNumY();
    const int numZ = Mask.getNumZ();
    // the most labels there can be is 1/2 of the points in checkerboard
    aol::Vector<unsigned int> labels ( numX * numY * numZ / 2 );
    //initialise labels
    for (int i=0; i<labels.size(); i++) labels[ i ] = i;
    
    //now Label the image
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          
          if ( Mask.get( x, y, z ) == false ) {
            //result = 0;  //nothing here
          } else if ( LabelArray.get( x, y, z ) == 0 ) {
            
            int found = 0;
            int myLabel = 0;
            //loop over neighbors
            //1) look for a neighbor which already has a label
            //2) associate myLable with all my neihgbors
            for(int xx = -1; xx <= 1; ++xx){
              int tmpX;
              if(PBC_X){
                tmpX = (x + xx+numX)%numX;
              }else {
                if( x+xx >= numX ||  x+xx < 0) continue;
                tmpX = x + xx;
              }
              
              for(int yy = -1; yy <= 1; ++yy){
                
                //if yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                if(connectivity == 0 && xx*xx + yy*yy > 1)
                  continue;
                
                int tmpY;
                if(PBC_Y){
                  tmpY = (y + yy+numY)%numY;
                }else{
                  if( y+yy >= numY ||  y+yy < 0) continue;
                  tmpY = y+yy;
                }
                
                for(int zz = -1; zz <= 1; ++zz){
                  int tmpZ;
                  if(PBC_Z){
                    tmpZ = (z + zz+numZ)%numZ;
                  }else{
                    if( z+zz >= numZ ||  z+zz < 0) continue;
                    tmpZ = z+zz;
                  }
                  
                  //if zz*zz + yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                  if(connectivity == 0 && zz*zz + yy*yy + xx*xx > 1)
                    continue;
                  
                  if(Mask.get(tmpX, tmpY, tmpZ)){
                    int nbrLabel = LabelArray.get ( tmpX, tmpY, tmpZ );
                    if(nbrLabel){
                      if(found == 0){
                        found = 1;
                        myLabel = nbrLabel;
                        LabelArray.set ( x, y, z, myLabel);
                      }else{
                        associate( labels, nbrLabel, myLabel);
                      }
                    }
                  }
                }
              }
            }
            
            //none of my neighbors has a label yet
            if ( found == 0 ){
              myLabel = nextlabel;
              LabelArray.set ( x, y, z, myLabel);
              nextlabel++;
            }
            
            //loop over neighbors
            //assign same label to all neighbors which don't have a label yet
            for(int xx = -1; xx <= 1; ++xx){
              int tmpX;
              if(PBC_X){
                tmpX = (x + xx+numX)%numX;
              }else {
                if( x+xx >= numX ||  x+xx < 0) continue;
                tmpX = x + xx;
              }
              
              for(int yy = -1; yy <= 1; ++yy){
                
                //if yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                if(connectivity == 0 && xx*xx + yy*yy > 1)
                  continue;
                
                int tmpY;
                if(PBC_Y){
                  tmpY = (y + yy+numY)%numY;
                }else{
                  if( y+yy >= numY ||  y+yy < 0) continue;
                  tmpY = y+yy;
                }
                
                for(int zz = -1; zz <= 1; ++zz){
                  int tmpZ;
                  if(PBC_Z){
                    tmpZ = (z + zz+numZ)%numZ;
                  }else{
                    if( z+zz >= numZ ||  z+zz < 0) continue;
                    tmpZ = z+zz;
                  }
                  
                  //if zz*zz + yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                  if(connectivity == 0 && zz*zz + yy*yy + xx*xx > 1)
                    continue;
                  
                  if(Mask.get(tmpX, tmpY, tmpZ) && LabelArray.get ( tmpX, tmpY, tmpZ ) == 0){
                    LabelArray.set ( tmpX, tmpY, tmpZ, myLabel );
                  }
                }
              }
            }
            
          }
        }
      }
    }
    
    //reduce labels ie 76=23=22=3 -> 76=3
    //done in reverse order to preserve sorting
    for( unsigned int i= labels.size()-1; i > 0; i-- ){
      labels[ i ] = reduce( labels, i );
    }
    
    /*now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
     this needs to be condensed down again, so that there is no wasted
     space eg in the above, the labels 3 and 4 are not used instead it jumps
     to 5.
     */
    aol::Vector<unsigned int> condensed ( nextlabel ); // cant be more than nextlabel labels
    
    unsigned int count = 0;
    for (unsigned int i=0; i< nextlabel; i++){
      if( i == labels[ i ] ) condensed[ i ] = count++;
    }
    // Record the number of labels
    const int numberOfLabels = count - 1;
    
    // now run back through our preliminary results, replacing the raw label
    // with the reduced and condensed one, and do the scaling and offsets too
    for (int i=0; i< LabelArray.size(); i++){
      LabelArray[i] = condensed[ labels[ LabelArray[ i ] ] ];
      labelSize[LabelArray[i]]++;
    }
    
    return numberOfLabels;
  }
};
  
  
template<typename _DataType>
class RandomAccessStencil : protected aol::RandomAccessContainer<_DataType> {
public:
  typedef _DataType DataType;
  
  RandomAccessStencil ( ) {
  }
  
  RandomAccessStencil ( const DataType &N, const DataType &E, const DataType &S, const DataType &W ) {
    set ( N, E, S, W );
  }
  
  void set ( const DataType &N, const DataType &E, const DataType &S, const DataType &W ) {
    this->clear ( );
    this->pushBack ( N );
    this->pushBack ( E );
    this->pushBack ( S );
    this->pushBack ( W );
  }
  
  const DataType& N ( ) const {
    return (*this)[0];
  }
  
  const DataType& E ( ) const {
    return (*this)[1];
  }
  
  const DataType& S ( ) const {
    return (*this)[2];
  }
  
  const DataType& W ( ) const {
    return (*this)[3];
  }
  
  bool empty ( ) const {
    return ( this->_data.size ( ) < 4 );
  }
};


template<typename _RealType>
class Component : public std::set<aol::Vec2<short> > {
public:
  typedef _RealType RealType;
  
protected:
  RandomAccessStencil<short> _boundaryBox;
  
public:
  Component ( ) { }
  
  const aol::Vec2<short> getGeometricCenter ( ) const {
    // short is too small to store the distance.
    aol::Vec2<RealType> center ( 0.0, 0.0 );
    for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it ) {
      center[0] += (*it)[0];
      center[1] += (*it)[1];
    }
    center /= this->size ( );
    center[0] = aol::Rint ( center[0] );
    center[1] = aol::Rint ( center[1] );
    const aol::Vec2<short> res ( center );
    return res;
  }
  
  const aol::Vec<4, RealType> getMinMaxDiameterAndAngle ( ) {
    const aol::Vec2<short> geomCenter = getGeometricCenter ( );
    const Component bndry = boundary ( );
    const RealType maxD = 2 * aol::Max ( width ( ), height ( ) );
    RealType maxDiameter = 0.0, minDiameter = maxD, diameter = 0.0, maxAngle = 0.0, minAngle = 0.0;
    for ( int deg=0; deg<360 ; deg+=10 ) {
      const RealType rad = deg * aol::NumberTrait<RealType>::pi / 180.0;
      diameter = 0;
      for ( short sign=-1; sign<=1 ; sign+=2 ) {
        Component boundaryIntersections;
        for ( std::set<aol::Vec2<short> >::const_iterator it=bndry.begin ( ); it != bndry.end ( ) ; ++it ) {
          aol::LineSegment<RealType, qc::QC_2D> line ( aol::Vec2<RealType> ( -maxD * cos ( rad ), -maxD * sin ( rad ) ),
                                                      aol::Vec2<RealType> ( maxD * cos ( rad ), maxD * sin ( rad ) ) );
          aol::Vec2<RealType> pt ( (*it)[0] - geomCenter[0], (*it)[1] - geomCenter[1] ), projPt;
          line.projectTo ( pt, projPt );
          pt -= projPt;
          if ( pt.norm ( ) < 1 )
            boundaryIntersections.insert ( *it );
        }
        RealType intersectionDist = 0.0, maxIntersectionDist = 0.0;
        for ( std::set<aol::Vec2<short> >::const_iterator it=boundaryIntersections.begin ( ); it != boundaryIntersections.end ( ) ; ++it ) {
          intersectionDist = aol::Vec2<RealType> ( (*it)[0] - geomCenter[0], (*it)[1] - geomCenter[1] ).norm ( );
          if ( intersectionDist > maxIntersectionDist )
            maxIntersectionDist = intersectionDist;
        }
        diameter += maxIntersectionDist;
      }
      if ( diameter > maxDiameter ) {
        maxDiameter = diameter;
        maxAngle = rad;
      }
      if ( diameter < minDiameter ) {
        minDiameter = diameter;
        minAngle = rad;
      }
    }
    
    aol::Vec<4, RealType> res;
    res[0] = maxDiameter;
    res[1] = maxAngle;
    res[2] = minDiameter;
    res[3] = minAngle;
    return res;
  }
  
  const Component boundary ( ) const {
    Component bndry;
    aol::Vec2<short> pos;
    for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it ) {
      bool isInBoundary = false;
      for ( int dx=-1; dx<=1 && !isInBoundary ; dx+=2 ) {
        pos.set ( *it );
        pos[0] += dx;
        if ( this->find ( pos ) == this->end ( ) ) {
          bndry.insert ( *it );
          isInBoundary = true;
        }
      }
      for ( int dy=-1; dy<=1 && !isInBoundary ; dy+=2 ) {
        pos.set ( *it );
        pos[1] += dy;
        if ( this->find ( pos ) == this->end ( ) ) {
          bndry.insert ( *it );
          isInBoundary = true;
        }
      }
    }
    
    return bndry;
  }
  
  RealType naiveNearestNeighborSearchDistance ( const aol::Vec<2, short> &Point ) const {
    RealType distance = aol::NumberTrait<RealType>::Inf;
    for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it ) {
      // short is too small to store the distance.
      const aol::Vec2<RealType> distanceVec ( (*it)[0] - Point[0], (*it)[1] - Point[1] );
      distance = aol::Min ( distance, distanceVec.normSqr ( ) );
    }
    return sqrt ( distance );
  }
  
  const RandomAccessStencil<short>& getBoundaryIndices ( ) {
    if ( _boundaryBox.empty ( ) ) {
      short N, E, S, W;
      // find left-most X coordinate
      W = 1000; // TODO: find more appropriate definition of upper bound
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[0] < W )
          W = (*it)[0];
      // find right-most X coordinate
      E = -1000;
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[0] > E )
          E = (*it)[0];
      // find bottom-most Y coordinate
      S = 1000;
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[1] < S )
          S = (*it)[1];
      // find top-most Y coordinate
      N = -1000;
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[1] > N )
          N = (*it)[1];
      _boundaryBox.set ( N, E, S, W );
    }
    return _boundaryBox;
  }
  
  short width ( ) {
    RandomAccessStencil<short> boundaryIndices = getBoundaryIndices ( );
    return boundaryIndices.E ( ) - boundaryIndices.W ( );
  }
  
  short height ( ) {
    RandomAccessStencil<short> boundaryIndices = getBoundaryIndices ( );
    return boundaryIndices.N ( ) - boundaryIndices.S ( );
  }
};

  
template<typename _RealType>
class ComponentsCollection {
public:
  typedef _RealType RealType;
  
protected:
  aol::RandomAccessContainer<Component<RealType> > _components;
  qc::ScalarArray<int, qc::QC_2D> _labelArr;
  aol::RandomAccessContainer<Component<RealType> > _neighborhoods;
  qc::ScalarArray<int, qc::QC_2D> _neighborhoodArr;
  
  
public:
  ComponentsCollection ( ) {
  }
  
  ComponentsCollection ( const qc::ScalarArray<int, qc::QC_2D> &Labels )
  : _labelArr ( Labels ), _neighborhoodArr ( Labels ) {
    init( );
  }
  
  ComponentsCollection ( const qc::BitArray<qc::QC_2D> &Mask )
  : _labelArr ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) ), _neighborhoodArr ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) ) {
    ConnectedComponentsLabeler::doLabel ( Mask, _labelArr );
    init ( );
  }
  
  void initializeFrom ( const qc::BitArray<qc::QC_2D> &Mask ) {
    _labelArr.reallocate ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) );
    _neighborhoodArr.reallocate ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) );
    ConnectedComponentsLabeler::doLabel ( Mask, _labelArr );
    init ( );
  }
  
private:
  void init ( ) {
    for ( int i=0; i<=_labelArr.getMaxValue ( ) ; ++i )
      _components.pushBack ( Component<RealType> ( ) );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _labelArr ); it.notAtEnd ( ) ; ++it )
      if ( _labelArr.get ( *it ) )
        _components[_labelArr.get ( *it )].insert ( *it );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _neighborhoodArr ); it.notAtEnd ( ) ; ++it )
      _neighborhoodArr.set ( *it, 0 );
  }
  
public:
  void clearBoundaryComponents ( ) {
    qc::GridStructure grid ( qc::GridSize2d::createFrom ( _labelArr ) );
    for ( qc::GridStructure::AllBoundaryNodeIterator itArr ( grid ); itArr.notAtEnd ( ) ; ++itArr ) {
      const int componentNr = _labelArr.get ( *itArr );
      if ( componentNr ) {
        for ( typename Component<RealType>::const_iterator itSet = _components[componentNr].begin ( ); itSet != _components[componentNr].end ( ) ; ++itSet )
          _labelArr.set ( *itSet, 0 );
        _components[componentNr].clear ( );
      }
    }
  }
  
  void cropComponentsByGeometricCenter ( const aol::Vec2<short> &Origin, const aol::Vec2<short> &Size ) {
    const std::set<aol::Vec2<short> > centers = getGeometricCenters ( );
    for ( std::set<aol::Vec2<short> >::const_iterator it = centers.begin ( ); it != centers.end ( ) ; ++it ) {
      if ( (*it)[0] < Origin[0] || (*it)[0] >= Origin[0] + Size[0] || (*it)[1] < Origin[1] || (*it)[1] >= Origin[1] + Size[1] ) {
        const int componentNr = getComponentNumberByPosition ( *it );
        for ( typename Component<RealType>::const_iterator itSet = _components[componentNr].begin ( ); itSet != _components[componentNr].end ( ) ; ++itSet )
          _labelArr.set ( *itSet, 0 );
        _components[componentNr].clear ( );
      }
    }
  }
  
  void createPictureBWComponents ( qc::ScalarArray<RealType, qc::QC_2D> &Picture ) const {
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _labelArr ); it.notAtEnd ( ) ; ++it )
      Picture.set ( *it, ( _labelArr.get ( *it ) > 0 ) );
    Picture.setOverflowHandlingToCurrentValueRange ( );
  }
  
  void createPictureBWComponentsRedGeometricCenters ( qc::MultiArray<RealType, qc::QC_2D, 3> &Picture ) {
    qc::GridStructure grid ( qc::GridSize2d::createFrom ( _labelArr ) );
    qc::ScalarArray<RealType, qc::QC_2D> bwComponents ( grid );
    createPictureBWComponents ( bwComponents );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( grid ); it.notAtEnd ( ) ; ++it )
      for ( int k=0; k<3 ; ++k )
        Picture[k].set ( *it, bwComponents.get ( *it ) );
    const std::set<aol::Vec2<short> > geometricCenters = getGeometricCenters ( );
    for ( std::set<aol::Vec2<short> >::const_iterator it = geometricCenters.begin ( ); it != geometricCenters.end ( ) ; ++it ) {
      Picture[0].set ( *it, 1 );
      Picture[1].set ( *it, 0 );
      Picture[2].set ( *it, 0 );
    }
  }
  
  void createPictureBWComponentsRGBNeighborhoods ( qc::MultiArray<RealType, qc::QC_2D, 3> &Picture ) {
    int color = 0;
    for ( ComponentsCollection<RealType>::NonEmptyComponentsIterator it ( *this ); it.notAtEnd ( ) ; ++it ) {
      const Component<RealType> neighborhood = getNeighborhoodByComponentNumber ( it.getCurrentComponentNumber ( ) );
      for ( std::set<aol::Vec2<short> >::const_iterator nhIt=neighborhood.begin ( ); nhIt != neighborhood.end ( ) ; ++nhIt ) {
        for ( int k=0; k<3 ; ++k )
          Picture[k].set ( *nhIt, 0 );
        Picture[color].set ( *nhIt, 1 );
      }
      color++;
      color %= 3;
    }
    qc::GridStructure grid ( qc::GridSize2d::createFrom ( _labelArr ) );
    qc::ScalarArray<RealType, qc::QC_2D> bwComponents ( grid );
    createPictureBWComponents ( bwComponents );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( grid ); it.notAtEnd ( ) ; ++it )
      if ( bwComponents.get( *it ) )
        for ( int k=0; k<3 ; ++k )
          Picture[k].set ( *it, bwComponents.get ( *it ) );
  }
  
  int getComponentNumberByPosition ( const aol::Vec2<short int> &Coord ) const {
    return _labelArr.get ( Coord );
  }
  
  int getComponentNumberByPosition ( const int X, const int Y ) const {
    return _labelArr.get ( X, Y );
  }
  
  Component<RealType>& getComponentByNumber ( const int I ) {
    return _components[I];
  }
  
  Component<RealType>& getComponentByPosition ( const aol::Vec2<short int> &Coord ) {
    return _components[getComponentNumberByPosition( Coord )];
  }
  
  Component<RealType>& getComponentByPosition ( const int X, const int Y ) {
    return _components[getComponentNumberByPosition( X , Y )];
  }
  
  std::set<aol::Vec2<short> > getGeometricCenters ( ) {
    std::set<aol::Vec2<short> > res;
    for ( NonEmptyComponentsIterator it ( *this ); it.notAtEnd ( ) ; ++it )
      res.insert ( (*it).getGeometricCenter( ) );
    return res;
  }
  
  int getNumNonEmptyComponents ( ) {
    int n = 0;
    for ( NonEmptyComponentsIterator it ( *this ); it.notAtEnd ( ) ; ++it )
      ++n;
    return n;
  }
  
  const Component<RealType>& getLargestComponent ( ) {
    int n = 0;
    for ( int i=1; i<_components.size ( ) ; ++i ) {
      if ( _components[i].size ( ) > _components[n].size ( ) )
        n = i;
    }
    return _components[n];
  }
  
  int getComponentNumberByNeighborhoodPosition ( const aol::Vec2<short int> &Coord ) {
    if ( _neighborhoodArr.size ( ) == 0 )
      naiveCalculateNeighborhoods ( );
    return _neighborhoodArr.get ( Coord );
  }
  
  int getComponentNumberByNeighborhoodPosition ( const int X, const int Y ) {
    if ( _neighborhoods.size ( ) == 0 )
      naiveCalculateNeighborhoods ( );
    return _neighborhoodArr.get ( X, Y );
  }
  
  const Component<RealType>& getNeighborhoodByComponentNumber ( const int I ) {
    if ( _neighborhoods.size ( ) == 0 )
      naiveCalculateNeighborhoods ( );
    return _neighborhoods[I];
  }
  
  void setSquareNeighborhoods ( const int Size ) {
    aol::Vec2<short> pos;
    const short neighborhoodOffset = ( Size - 1 ) / 2;
    for ( int i=0; i<=_labelArr.getMaxValue ( ) ; ++i )
      _neighborhoods.pushBack ( Component<RealType> ( ) );
    const std::set<aol::Vec2<short> > centers = getGeometricCenters ( );
    for ( std::set<aol::Vec2<short> >::const_iterator it = centers.begin ( ); it != centers.end ( ) ; ++it ) {
      const int componentNr = getComponentNumberByPosition ( *it );
      for ( int i=-neighborhoodOffset; i<=neighborhoodOffset ; ++i )
        for ( int j=-neighborhoodOffset; j<=neighborhoodOffset ; ++j )
          if ( (*it)[0] + i >= 0 && (*it)[0] + i < _labelArr.getNumXYZ ( ) && (*it)[1] + j >= 0 && (*it)[1] + j < _labelArr.getNumXYZ ( ) ) {
            pos.set ( i, j );
            pos += *it;
            _neighborhoods[componentNr].insert ( pos );
            _neighborhoodArr.set ( pos, componentNr );
          }
    }
  }
  
  void naiveCalculateNeighborhoods ( ) {
    RealType distance, minDist;
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > arrIt ( _neighborhoodArr ); arrIt.notAtEnd ( ) ; ++arrIt ) {
      minDist = aol::NumberTrait<RealType>::Inf;
      for ( NonEmptyComponentsIterator compIt ( *this ); compIt.notAtEnd ( ) ; ++compIt ) {
        distance = (*compIt).naiveNearestNeighborSearchDistance ( *arrIt );
        if ( distance < minDist ) {
          minDist = distance;
          _neighborhoodArr.set ( *arrIt, compIt.getCurrentComponentNumber ( ) );
        }
      }
    }
    for ( int i=0; i<=_labelArr.getMaxValue ( ) ; ++i )
      _neighborhoods.pushBack ( Component<RealType> ( ) );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _neighborhoodArr ); it.notAtEnd ( ) ; ++it )
      if ( _neighborhoodArr.get ( *it ) )
        _neighborhoods[_neighborhoodArr.get ( *it )].insert ( *it );
  }
  
  class NonEmptyComponentsIterator {
  private:
    ComponentsCollection &_compCollection;
    int _curCompNr, _lastNonEmptyCompNr;
  public:
    typedef NonEmptyComponentsIterator Self;
    
    NonEmptyComponentsIterator ( ComponentsCollection &CompCollection ) : _compCollection ( CompCollection ) {
      _lastNonEmptyCompNr = 0;
      for ( int i=1; i<_compCollection._components.size ( ) ; ++i ) {
        if ( !_compCollection.getComponentByNumber ( i ).empty ( ) )
          _lastNonEmptyCompNr = i;
      }
      
      _curCompNr = 0;
      ++(*this);
    }
    
    bool notAtEnd ( ) const {
      return ( _curCompNr <= _lastNonEmptyCompNr );
    }
    
    Self& operator++ ( ) {
      ++_curCompNr;
      while ( notAtEnd ( ) && _compCollection.getComponentByNumber ( _curCompNr ).empty ( ) )
        ++_curCompNr;
      return *this;
    }
    
    Component<RealType>& operator* ( ) const {
      return _compCollection.getComponentByNumber ( _curCompNr );
    }
    
    int getCurrentComponentNumber ( ) const {
      return _curCompNr;
    }
  };
};
  
  
} // end namespace
  

#endif
