#include<adaptiveTriangMesh.h>

#if defined(__MINGW32_VERSION) || defined(_MSC_VER) || defined(__BCPLUSPLUS__)
#include <process.h>
#define GNUPLOT_PATH string("\"C:\\Program Files\\gnuplot\\bin\\wgnuplot.exe\" ")
#else
#include <unistd.h>
#define GNUPLOT_PATH string("gnuplot ")
#endif

// switch off warnings induced by system() calls in AdaptiveFETriangMesh<RealType, TriangleType>::plot2DProjection
WARNING_OFF( unused-result )

/**********************************************************************************************/
/*                                  class DartIterator                                */
/**********************************************************************************************/

//! moves to neighbouring triangle across _edge
template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::DartIterator::flipTriangle(){
  
    if( !canFlipTriangle() )
      throw aol::Exception ( "DartIterator::flipTriangle(): cannot flip this triangle!", __FILE__, __LINE__ ); 
  
    GlobalIndex newTriangleIndex = getNextTriangleIndex();
    GlobalIndex globalNodeIndex = _grid.getTriangNodeIdx( _triangle, _node);    
    GlobalIndex globalOtherNodesIndex = _grid.getTriangNodeIdx( _triangle, 3 - (_node + _edge) );

    // find out which node of the new triangle is ours
    LocalIndex newNodeIndex = IndexNotSet;
    for (LocalIndex i = 0; i < 3; ++i){
      if ( _grid.getTriangNodeIdx( newTriangleIndex, i) == globalNodeIndex ){
        newNodeIndex = i;
        break;
      }
    }
      
    // find out which node of the new triangle lies on our edge
    LocalIndex newOtherNodesIndex = IndexNotSet;
    for (LocalIndex i = 1; i < 3; ++i){
      if ( _grid.getTriangNodeIdx( newTriangleIndex,(newNodeIndex + i) % 3 ) == globalOtherNodesIndex ){
        newOtherNodesIndex = (newNodeIndex + i) % 3;
        break;
      }
    }


    if ( min( newNodeIndex, newOtherNodesIndex ) < 0 ) {
        stringstream err_msg;
        err_msg << "DartIterator::flipTriangle(): found inconsistency in grid: triangle " << newTriangleIndex << " should be neighbour of " << _triangle
                << ", but does not contain node " << globalNodeIndex << " or " << globalOtherNodesIndex << ".";
        throw aol::Exception ( err_msg.str(), __FILE__, __LINE__ ); 
    }

    LocalIndex newEdgeIndex = (-(newNodeIndex + newOtherNodesIndex)) % 3;
    if (newEdgeIndex < 0) newEdgeIndex += 3;

    // all information found, now set:
    _triangle = newTriangleIndex;
    _node = newNodeIndex;
    _edge = newEdgeIndex;
}



/**********************************************************************************************/
/*                              class BoundaryIterator                                */
/**********************************************************************************************/

//! Kostruktor aus einem grid
template< typename RealType, typename TriangleType>
AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator::BoundaryIterator( const MeshType & grid )
: _grid(grid)
, _dart(grid,0,0,0)
, _startIndex(0)
, _initialized(true)
{
  // Start looking for the boundary
  int found = 0;
  typename aol::FETriangMesh<RealType, TriangleType >::TriangleIterator tIter (_grid);
  while ( ( tIter.notAtEnd() ) && (!found) )  {
    for( int i = 0; i < 3; i++ )
      if ( _grid.getNeighbour( tIter.getIndex() , i ) == -1 ) found = i+1;
    if ( !found ) ++tIter;
  }
  if (!tIter.notAtEnd()) throw  aol::Exception ( "BoundaryIterator: Could not find boudary node, iterator will not work!", __FILE__, __LINE__ );

  _startIndex = tIter.getNodeIndex(found%3);
  _dart.set( tIter.getIndex(), found-1, found%3 );
}
// ------------------------------------------------------------------------------------------------

//! Iterator am Ende der Dreiecksliste?
template< typename RealType, typename TriangleType>
bool AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator::notAtEnd() const
{
  return (( _initialized ) || ( _dart.getGlobalNodeIndex() != _startIndex ));
}
// ------------------------------------------------------------------------------------------------

//! Inkrement-Operator
template< typename RealType, typename TriangleType>
typename AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator & AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator::operator++()
{
  _initialized = false;
  _dart.flipNode();
  _dart.flipEdge();
  while ( _dart.canFlipTriangle() )  {
    _dart.flipTriangle();
    _dart.flipEdge();
  }
  return *this;
}
// ------------------------------------------------------------------------------------------------

//! Dekrement-Operator
template< typename RealType, typename TriangleType>
typename AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator & AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator::operator--()
{
  _initialized = false;
  _dart.flipEdge();
  while ( _dart.canFlipTriangle() )  {
    _dart.flipTriangle();
    _dart.flipEdge();
  }
  _dart.flipNode();
  return *this;
}
// ------------------------------------------------------------------------------------------------

//! Knotenindex am Anfang des Segments
template< typename RealType, typename TriangleType>
GlobalIndex AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator::getStartIndex() const
{
  return _dart.getGlobalNodeIndex();
}
// ------------------------------------------------------------------------------------------------

//! Knotenindex am Ende des Segments
template< typename RealType, typename TriangleType>
GlobalIndex AdaptiveFETriangMesh<RealType, TriangleType>::BoundaryIterator::getEndIndex() const
{
  return _dart.getNextNodeGlobalIndex();
}
// ------------------------------------------------------------------------------------------------



//! plot projection of *this to (x,y)-plane in desired format, possibly with encoding neighbouring relations and refinement markers
//! For FormatType JPG there is the option to have a rectangular output with width/height = relWidthHeight and min(width, height) = 800 pixel.
template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::plot2DProjection( string filename, FormatType format, bool neighbours, bool marked, double relWidthHeight ) const {
  
    fstream file;

    switch ( format )
    {
        case MATLAB:
        {
            // print triangulation into three files used by MatLab to display the triangulation:
            // filename + "_x.dat" and filename + "_y.dat" for x- and y-coordinates of the
            // grid nodes, and filename + "_tri.dat" for a 3-by-m-matrix containing the global
            // node indices. The triangulation can be displayed in MatLab loading all three files
            // and using "triplot(tri,x,y)".
            //filename = filename.substr( 0, filename.rfind(".") );

            ofstream x_file((filename + "_x.dat").c_str()),
                     y_file((filename + "_y.dat").c_str());

            for (NodeIteratorType nIter( *this ); nIter.notAtEnd(); ++nIter)
            {
                x_file << (*nIter)[0] << endl;
                y_file << (*nIter)[1] << endl;
            }
            x_file.close();
            y_file.close();

            ofstream tri_file((filename + "_tri.dat").c_str());
            for (ElementIteratorType tIter = *this; tIter.notAtEnd(); ++tIter)
                tri_file << tIter.getNodeIndex(0) + 1 << "\t"
                         << tIter.getNodeIndex(1) + 1 << "\t"
                         << tIter.getNodeIndex(2) + 1 << endl;
            tri_file.close();
            break;
        }

        case GNUPLOT:
        {
            //filename = filename.substr( 0, filename.rfind(".") );
            file.open( (filename+".txt").c_str(), fstream::out );
            file << "set title 'Planar mesh'" << endl;
            file << "set size square" << endl;
            file << "unset border" << endl;
            file << "set xtics nomirror" << endl;
            file << "set ytics nomirror" << endl;
            file << "plot ";
            if ( marked )     file << "'" << filename << "_m.txt' notitle with filledcurves lw 1 lc rgb \"red\",";
            if ( neighbours ) file << "'" << filename << "_n.txt' notitle with lines lw 1 lc rgb \"green\",";
            file << "'" << filename << "_t.txt' notitle with lines lw 2 lc rgb \"black\"" << endl;
            file.close();

            file.open( (filename+"_t.txt").c_str(), fstream::out );
            if ( file.is_open() )
            {
                for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                {
                    file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl;
                    file << iter->getNode(1)[0] << " " << iter->getNode(1)[1] << endl;
                    file << iter->getNode(2)[0] << " " << iter->getNode(2)[1] << endl;
                    file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl << endl;
                }
                file.close();
                cout << "To see the plot use 'gnuplot " << filename << ".txt -'" << endl;
            }
            else
            {
                cerr << "Error opening file " << filename << endl;
            }

            if ( neighbours )
            {
                file.open( (filename+"_n.txt").c_str(), fstream::out );
                if ( file.is_open() )
                {
                    for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    {
		        aol::Vec3<RealType> midpoint = iter->centerOfMass ( );
                         for (LocalIndex i = 0; i < 3; ++i)
                            if (this->getNeighbour(iter->getIndex(), i) != IndexNotSet)
                            {
                                aol::Vec3<RealType> otherMidpoint = TriangleType( *this, this->getNeighbour(iter->getIndex(), i) ).centerOfMass ( );
                                file << midpoint[0] << " " << midpoint[1] << endl;
                                file << otherMidpoint[0] << " " << otherMidpoint[1] << endl << endl;
                            }
                    }
                    file.close();
                }
                else
                {
                    cerr << "Error opening file " << filename << endl;
                }
            }

            if ( marked )
            {
                file.open( (filename+"_m.txt").c_str(), fstream::out );
                if ( file.is_open() )
                {
                    for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    {
                        if ( isMarkedForRefinement( iter->getIndex() ) )
                        {
                            file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl;
                            file << iter->getNode(1)[0] << " " << iter->getNode(1)[1] << endl;
                            file << iter->getNode(2)[0] << " " << iter->getNode(2)[1] << endl;
                            file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl << endl;
                        }
                    }
                    file.close();
                }
                else
                {
                    cerr << "Error opening file " << filename << endl;
                }
            }

            break;
        }

        case POSTSCRIPT:
        {
            //filename = filename.substr( 0, filename.rfind(".") ) + ".eps";
            stringstream ss; ss <<
#ifdef _MSC_VER
                _getpid();
#else
                getpid();
#endif
            string tmpfile = ".gp_" + ss.str() + ".tmp";
            file.open( tmpfile.c_str(), fstream::out );
            if ( file.is_open() )
            {
                file << "set terminal postscript enhanced color eps" << endl;
                if ( filename == "output" ) file << "set output \"output.eps\"" << endl;
                else file << "set output \"" << filename << "\"" << endl;
                file << "set title '{/*2 Planar mesh}'" << endl;
		file << "set size ratio -1" << endl;
                file << "unset border" << endl;
                file << "set xtics nomirror" << endl;
                file << "set ytics nomirror" << endl;
                file << "plot ";
                if ( marked )     file << "'-' notitle with filledcurves lw 1 lc rgb \"red\",";
                if ( neighbours ) file << "'-' notitle with lines lw 1 lc rgb \"green\",";
                file << "'-' notitle with lines lt 1 lw 2 lc rgb \"black\"" << endl;

                if ( marked ) {
                    for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    {
                        if ( isMarkedForRefinement( iter->getIndex() ) )
                        {
                            file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl;
                            file << iter->getNode(1)[0] << " " << iter->getNode(1)[1] << endl;
                            file << iter->getNode(2)[0] << " " << iter->getNode(2)[1] << endl;
                            file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl << endl;
                        }
                    }
                    file << "e" << endl;
                }

                if ( neighbours ) {
                    for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    {
                        aol::Vec3<RealType> midpoint = iter->centerOfMass ( );
                        for (LocalIndex i = 0; i < 3; ++i)
                            if (this->getNeighbour(iter->getIndex(), i) != IndexNotSet)
                            {
                                aol::Vec3<RealType> otherMidpoint = TriangleType( *this, this->getNeighbour(iter->getIndex(), i) ).centerOfMass ( );
                                file << midpoint[0] << " " << midpoint[1] << endl;
                                file << otherMidpoint[0] << " " << otherMidpoint[1] << endl << endl;
                            }
                    }
                    file << "e" << endl;
                }

                for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                {
                    file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl;
                    file << iter->getNode(1)[0] << " " << iter->getNode(1)[1] << endl;
                    file << iter->getNode(2)[0] << " " << iter->getNode(2)[1] << endl;
                    file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl << endl;
                }
                file << "e" << endl;

                file.close();
                system( (GNUPLOT_PATH + tmpfile).c_str() );
                remove( tmpfile.c_str() );
            }
            else
            {
                cerr << "Error opening file " << filename << endl;
            }
            break;
        }

        case JPG:
        {
	    //For FormatType JPG there is the option to have a rectangular output with width/height = relWidthHeight and min(width, height) = 800 pixel.

            //filename = filename.substr( 0, filename.rfind(".") ) + ".jpg";
            stringstream ss; ss <<
#ifdef _MSC_VER
                _getpid();
#else
                getpid();
#endif
            string tmpfile = ".gp_" + ss.str() + ".tmp";
            file.open( tmpfile.c_str(), fstream::out );
            if ( file.is_open() )
            {
		int pixelsX, pixelsY;
	        if(relWidthHeight > 1.0){
		  pixelsX = 800*relWidthHeight;
		  pixelsY = 800;
		} else {
		  pixelsX = 800;
		  pixelsY = aol::Rint(800/relWidthHeight);
		}
                file << "set terminal jpeg large size " << pixelsX << "," << pixelsY << " enhanced" << endl;
                if ( filename == "output" ) file << "set output \"output.jpg\"" << endl;
                else file << "set output \"" << filename << "\"" << endl;
                file << "set title 'Planar mesh'" << endl;
                file << "unset border" << endl;
		file << "set size ratio -1" << endl;
                file << "set xtics nomirror" << endl;
                file << "set ytics nomirror" << endl;
                file << "plot ";
                if ( marked )     file << "'-' notitle with filledcurves lw 1 lc rgb \"red\",";
                if ( neighbours ) file << "'-' notitle with lines lw 1 lc rgb \"green\",";
                file << "'-' notitle with lines lw 2 lc rgb \"black\"" << endl;

                if ( marked ) {
                    for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    {
                        if ( isMarkedForRefinement( iter->getIndex() ) )
                        {
                            file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl;
                            file << iter->getNode(1)[0] << " " << iter->getNode(1)[1] << endl;
                            file << iter->getNode(2)[0] << " " << iter->getNode(2)[1] << endl;
                            file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl << endl;
                        }
                    }
                    file << "e" << endl;
                }

                if ( neighbours ) {
                    for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    {
                        aol::Vec3<RealType> midpoint = iter->centerOfMass ( );
                        for (LocalIndex i = 0; i < 3; ++i)
                            if (this->getNeighbour(iter->getIndex(), i) != IndexNotSet)
                            {
                                aol::Vec3<RealType> otherMidpoint = TriangleType( *this, this->getNeighbour(iter->getIndex(), i) ).centerOfMass ( );
                                file << midpoint[0] << " " << midpoint[1] << endl;
                                file << otherMidpoint[0] << " " << otherMidpoint[1] << endl << endl;
                            }
                    }
                    file << "e" << endl;
                }

                for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                {
                    file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl;
                    file << iter->getNode(1)[0] << " " << iter->getNode(1)[1] << endl;
                    file << iter->getNode(2)[0] << " " << iter->getNode(2)[1] << endl;
                    file << iter->getNode(0)[0] << " " << iter->getNode(0)[1] << endl << endl;
                }
                file << "e" << endl;

                file.close();
                system( (GNUPLOT_PATH + tmpfile).c_str() );
                remove( tmpfile.c_str() );
            }
            else
            {
                cerr << "Error opening file " << filename << endl;
            }
            break;
        }

        case PLY:
        {
            //filename = filename.substr( 0, filename.rfind(".") ) + ".ply";
            file.open( filename.c_str(), fstream::out );
            if ( file.is_open() )
            {
                file << "ply\n"
                     << "format ascii 1.0\n"
                     << "element vertex " << this->getNumVertices () << endl
                     << "property float x\n"
                     << "property float y\n"
                     << "property float z\n"
                     << "element face "   << this->getNumFaces() << endl
                     << "property list uchar int vertex_index\n"
                     << "end_header\n";
                for ( NodeIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    file << (*iter)[0] << " " << (*iter)[1] << " " << 0.0 << endl;

                for ( ElementIteratorType iter ( (*this) ); iter.notAtEnd(); ++iter )
                    file << "3 " << iter.getNodeIndex(0) << " " << iter.getNodeIndex(1) << " " << iter.getNodeIndex(2) << endl;
            }
            else
            {
                cerr << "Error opening file " << filename << endl;
            }
            break;
        }
      
	default: throw aol::Exception ( "AdaptiveFETriangMesh<>::plot(): unknown format!", __FILE__, __LINE__ ); 

    }

    return;
}

//!
template< typename RealType, typename TriangleType>
LocalIndex AdaptiveFETriangMesh<RealType, TriangleType>::getLongestEdgeIndex( GlobalIndex triangle ) const {
  // get global node indices of triangle
  aol::Vec3<int> triangIdx = this->getTriang( triangle );

  // assume that edge 0 is longest edge
  LocalIndex localIdx = 0;
  aol::Vec3<RealType> edge = this->getVertex( triangIdx[1] );
  edge -= this->getVertex( triangIdx[2] );
  RealType maxLength = edge.normSqr();  
 
  // now check if edge 1 is longer
  edge = this->getVertex( triangIdx[2] );
  edge -= this->getVertex( triangIdx[0] );
  RealType candidate = edge.normSqr();
  
  if( candidate > maxLength ){
    localIdx = 1;
    maxLength = candidate;
  }

  // now check if edge 2 is longer
  edge = this->getVertex( triangIdx[0] );
  edge -= this->getVertex( triangIdx[1] );
  candidate = edge.normSqr();
  
  if( candidate > maxLength )
    localIdx = 2;
  
  return localIdx;  
}
 
//! refine all triangles via subdivision (into four congruent sub triangles)
template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::subdivide()
{
  // adaptive meshes need neighbor information, but this is ensured in the constructor

  // Merke die Anzahl der Dreiecke vor der Verfeinerung
  int triangleCount = this->getNumTriangs();

  // Fuer jedes Dreieck: Kein TriangleIterator, da sich das Gitter aendern wird
  for ( int tIndex = 0; tIndex < triangleCount; tIndex++ )
  {
    // Sammle Daten fuer neue Knoten (Mittelpunkte), alte Knoten (Ecken) und Nachbarindices
    int newNodes[3];
    int oldNodes[3];
    int neighbours[6];
    // Fuer jede Seite
    for ( int i = 0; i < 3; i++ )
    {
      // Laufe von Knoten 0 in Richtung Knoten 1 usw.
      DartIterator d ( *this, tIndex, (2+i)%3, i );
      oldNodes[i] = d.getGlobalNodeIndex();
      // Falls diese Seite am Rand liegt kann ein neuer Knoten hinzugefuegt werden
      if ( !d.canFlipTriangle() )
      {
        newNodes[i]       = addEdgeMidpoint( d );
        neighbours[2*i  ] = IndexNotSet;
        neighbours[2*i+1] = IndexNotSet;
      }
      else
      {
        // Falls das andere Dreieck noch nicht verfeinert wurde, kann ein neuer Knoten hinzugefuegt werden
        if ( d.getNextTriangleIndex() > d.getGlobalTriangleIndex() )
        {
          // Errechne die Indices der Nachbardreiecke (evtl. noch gar nicht vorhanden!)
          int nindex = triangleCount + 3 * d.getNextTriangleIndex();
          newNodes[i] = addEdgeMidpoint( d );
          // Beachte Reihenfolge in der Dreiecke hinzugefuegt werden
          d.flipTriangle();
          neighbours[2*i  ] = nindex + d.getLocalNodeIndex();
          neighbours[2*i+1] = nindex + d.getNextNodeLocalIndex();
        }
        // Ansonsten existiert das andere Dreieck bereits
        else
        {
          // Suche nach dem Eckpunkt der auf aktueller Kante liegt
          aol::Vec3<RealType> midpoint =  this->getVertex ( d.getGlobalNodeIndex() );
          midpoint += this->getVertex ( d.getNextNodeGlobalIndex() );
          midpoint /= 2.;
          for ( int j = 0; j < 3; j++ )
          {
            if ( midpoint == this->getVertex( this->getTriangNodeIdx( d.getNextTriangleIndex(), j ) ) )
            {
              newNodes[i] = this->getTriangNodeIdx( d.getNextTriangleIndex(), j );
              // Bewege den Iterator ins naechste Dreieck und dort zur aktuellen Seite
              DartIterator dOther1( *this, d.getNextTriangleIndex(), (j+1)%3, j );
              dOther1.flipTriangle();
              dOther1.flipEdge();
              dOther1.flipNode();
              // Bewege einen zweiten in die andere Richtung
              DartIterator dOther2( *this, d.getNextTriangleIndex(), (j+1)%3, j );
              dOther2.flipEdge();
              dOther2.flipTriangle();
              // Entweder die Eckpunkte stimmen ueberein...
              if ( dOther1.getGlobalNodeIndex() == d.getGlobalNodeIndex() )
              {
                neighbours[2*i  ] = dOther1.getGlobalTriangleIndex();
                neighbours[2*i+1] = dOther2.getGlobalTriangleIndex();
              }
              // ...oder nicht
              else
              {
                neighbours[2*i  ] =  dOther2.getGlobalTriangleIndex();
                neighbours[2*i+1] =  dOther1.getGlobalTriangleIndex();
              }
            }
          }
        }
      }
    }
    // Fuege 3 neue Dreiecke in den Ecken hinzu und mache aus dem alten das zentrale Dreieck
    for ( int i = 0; i < 3; i++ )
    {
      // Knoten 0 ist ehemalige Ecke
      GlobalIndex T_neuIndex = this->pushBackTriang ( aol::Vec3<int>( oldNodes[i], newNodes[i], newNodes[(2+i)%3] ) );
      // Setze Nachbarn
      this->_neighbour_.push_back( aol::Vec3<int> ( tIndex, neighbours[(2*i+5)%6], neighbours[2*i] ) ); // arg2: das ist 2*i-1 im geeigneten Sinne
      this->getTriang( tIndex )[ i ] = newNodes[i];
      this->setNeighbour ( tIndex, (i+1)%3, T_neuIndex );
    }
  }

  this->makeOrientationConsistent();
}
// ------------------------------------------------------------------------------------------------

//! refine all triangles via bisection
template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::refineAll()
{
    for ( ElementIteratorType tIter = *this; tIter.notAtEnd(); ++tIter )
      mark( tIter.getIndex() );
    refineMarkedTriangles();
}
// ------------------------------------------------------------------------------------------------

//! refine all MARKED triangles via bisection
template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::refineMarkedTriangles()
{
    for ( ElementIteratorType tIter = *this; tIter.notAtEnd(); ++tIter )
      if ( isMarkedForRefinement( tIter.getIndex() ) ) 
	refine( tIter.getIndex() );
    this->makeOrientationConsistent();
}
// ------------------------------------------------------------------------------------------------

//! refinement
//! \param  triangleToRefine    Globaler Index des zu verfeinernden Dreiecks
template< typename RealType, typename TriangleType>
GlobalIndex AdaptiveFETriangMesh<RealType, TriangleType>::refine( GlobalIndex triangleToRefine )
{
    DartIterator d( *this, triangleToRefine, getLongestEdgeIndex( triangleToRefine ) );
    return refine( d );

}
// ------------------------------------------------------------------------------------------------

//! \brief results in at least one bisection of d.triangle along d.edge.
//!
//! modifies d.triangle such that a new node in the middle of d.edge
//! is inserted and d.triangle has d.node, the new point and the point
//! opposite to d.edge after the function's end.
//! The function returns the index of a new triangle that lies on the other half
//! of d.edge and is correctly connected with all triangles that were neighboured to
//! d.triangle before calling bisect().
//! If d.edge was not the longest edge of d.triangle, more bisections along other edges are made
//! before d.edge is cut in two. Indices of these inbetween-created triangles are not
//! returned.
template< typename RealType, typename TriangleType>
GlobalIndex AdaptiveFETriangMesh<RealType, TriangleType>::refine( const DartIterator& d )
{
#ifdef DEBUGMODE
cerr << "refine " << d.getGlobalTriangleIndex() << " along edge " << d.getLocalEdgeIndex() << "...";
#endif

    // first cut the neighbour on edge l in two triangles:
    GlobalIndex P_neu = addEdgeMidpoint( d ); 
    DartIterator d_prime = d;
    bool hasNeighbourOn_d_prime = d_prime.canFlipTriangle();
    GlobalIndex T_l_neu = IndexNotSet;

    if ( hasNeighbourOn_d_prime )
    {
        d_prime.flipTriangle();
#ifdef DEBUGMODE
cerr << "has first to refine " << d_prime.getGlobalTriangleIndex() << "...";
#endif
        T_l_neu = refineOnlyThis(d_prime, P_neu);
    }
#ifdef DEBUGMODE
cerr << "Now refine " << d.getGlobalTriangleIndex() << "...";
#endif
    
    // create new triangle T_neu and connect with neighbours
    GlobalIndex T_neu = refineOnlyThis(d, P_neu);

    // take care of neighbours
    if (hasNeighbourOn_d_prime)
    {
#ifdef DEBUGMODE    
cerr << "Now take care of neighbors of " << d.getGlobalTriangleIndex() << "...";
#endif
        // set d_prime to edge which lies in T_l_neu, opposite to edge d2.edge of T_neu
        d_prime.flipNode();
        d_prime.flipEdge();
        // flipTriangle() here will always work, as we have just created this
        // triangle. Thus, no canFlipTriangle() call is needed.
        d_prime.flipTriangle();
        d_prime.flipEdge();
        // set neighbour GlobalIndex in T_l_neu
#ifdef DEBUGMODE
cerr << "set neighbor(a) " << d_prime.getGlobalTriangleIndex() << " over edge " << d_prime.getLocalEdgeIndex() << " meets " << T_neu << "...";
#endif
	this->setNeighbour ( d_prime.getGlobalTriangleIndex(), d_prime.getLocalEdgeIndex(), T_neu );

        // move d1 to T_neu
        DartIterator d1 = d; d1.flipNode(); d1.flipEdge();
        d1.flipTriangle();
        d1.flipEdge();
        // set neighbour in T_neu
#ifdef DEBUGMODE
cerr << "set neighbor(b) " << d1.getGlobalTriangleIndex() << " over edge " << d1.getLocalEdgeIndex() << " meets " << T_l_neu << "...";
#endif
	this->setNeighbour ( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_l_neu );
    }
#ifdef DEBUGMODE
cerr << "done (" <<d.getGlobalTriangleIndex() << " to "<<T_neu<<").\n\n";
#endif

    return T_neu;
}
// ------------------------------------------------------------------------------------------------

template< typename RealType, typename TriangleType>
GlobalIndex AdaptiveFETriangMesh<RealType, TriangleType>::refineOnlyThis( const DartIterator& d, GlobalIndex midpoint)
{
#ifdef DEBUGMODE
cerr << "refineOnlyThis "<<d.getGlobalTriangleIndex() <<"...";
#endif

    // NOTE: in comments for this function, we write T for *this.
    LocalIndex l = getLongestEdgeIndex( d.getGlobalTriangleIndex() );
    if (l != d.getLocalEdgeIndex())
    {
#ifdef DEBUGMODE
cerr << " not longest edge...";
#endif

        DartIterator d_l(*this, d.getGlobalTriangleIndex(), l);
        LocalIndex commonNode = d.getCommonNodeLocalIndex(d_l);
        DartIterator d_prime( *this, d.getGlobalTriangleIndex(), l, commonNode);
        // bisect T along the longest edge
        refine(d_prime);
        // because d_prime.node also lies on d.edge, T is now modifies in such a way
        // that d.edge has remained unchanged. Proceed by recursion until d.edge is
        // the longest edge in T.
        return refineOnlyThis(d, midpoint);
    }

    // create DartIterators to all edges, starting from l
    DartIterator d1 = d;  d1.flipNode(); d1.flipEdge();
    DartIterator d2 = d1; d2.flipNode(); d2.flipEdge();

    // create new triangle T_neu and connect with neighbours
    GlobalIndex T_neuIndex = this->pushBackTriang ( aol::Vec3<int>( midpoint, d1.getGlobalNodeIndex(), d2.getGlobalNodeIndex() ) );
#ifdef DEBUGMODE
    cerr << " add new face (" << midpoint << ", " << d1.getGlobalNodeIndex() << ", " << d2.getGlobalNodeIndex() << ")...";
    cerr << " set neighbours (" << d1.getNextTriangleIndex()<< ", " <<d.getGlobalTriangleIndex()<< ", " << d.getNextTriangleIndex() << ") ... ";
#endif

    this->_neighbour_.push_back( aol::Vec3<int> ( d1.getNextTriangleIndex(), d.getGlobalTriangleIndex(), d.getNextTriangleIndex() ) );


    // set neighbour GlobalIndex in triangle lying opposite to d1.edge to T_neu
    DartIterator d1_prime = d1;
    if (d1_prime.canFlipTriangle())
    {
        d1_prime.flipTriangle();
	this->setNeighbour ( d1_prime.getGlobalTriangleIndex(), d1_prime.getLocalEdgeIndex(), T_neuIndex );
#ifdef DEBUGMODE
cerr << "set neighbor(c) " << d1_prime.getGlobalTriangleIndex() << " over edge " << d1_prime.getLocalEdgeIndex() << " meets " << T_neuIndex << "...";
#endif
    }

    // connect ourselves (i. e., connect T) to T_neu:
    this->setNeighbour( d1.getGlobalTriangleIndex(), d1.getLocalEdgeIndex(), T_neuIndex );
#ifdef DEBUGMODE    
cerr << " unmark " << d1.getGlobalTriangleIndex() << "...";    
#endif
    // do not have to refine this triangle again
    unmark( d1.getGlobalTriangleIndex() );
    
    // old triangle now has a new vertex, ie the midpoint
    this->getTriang( d1.getGlobalTriangleIndex() )[ d1.getLocalNodeIndex() ] = midpoint;

#ifdef DEBUGMODE
cerr << "done(" <<d.getGlobalTriangleIndex() << " to "<<T_neuIndex<<").\n";
#endif
    return T_neuIndex;
}

// ------------------------------------------------------------------------------------------------

template< typename RealType, typename TriangleType>
GlobalIndex AdaptiveFETriangMesh<RealType, TriangleType>::addEdgeMidpoint( const DartIterator& d ) {    
  DartIterator d1 = d;
  d1.flipNode();
  // get edge midpoint between nodes n1 and n2
  aol::Vec3<RealType> newVertex =  this->getVertex ( d.getGlobalNodeIndex() );
  newVertex += this->getVertex ( d1.getGlobalNodeIndex() );
  newVertex /= 2.;
  // add new vertex n = (n1 + n2)/2
  int newNodeIdx = this->pushBackVertex( newVertex ); 
  // update interpolation map by element ( n -> (n1,n2) )
  _interpolationMap.insert( std::pair< int, aol::Vec2<int> >( newNodeIdx, aol::Vec2<int>( d.getGlobalNodeIndex(), d1.getGlobalNodeIndex() ) ) );
  return newNodeIdx;
}

// ------------------------------------------------------------------------------------------------

template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::prolongateLinearly( aol::Vector<RealType>& function ) const {
  int oldSize = function.size();
  function.resize( this->getNumVertices() );
  
  std::map< int, aol::Vec2<int> >::const_iterator iter;
  for( int i = oldSize; i < function.size(); i++ ){
    // find new node index in map
    iter = _interpolationMap.find( i );
    if( iter == _interpolationMap.end() )
      throw aol::Exception ( "AdaptiveFETriangMesh<>::prolongateLinearly(): unknown vertex!", __FILE__, __LINE__ ); 
    // interpolate linearly
    function[i] = 0.5 * function[iter->second[0]] + 0.5 * function[iter->second[1]];
  }  
}

template< typename RealType, typename TriangleType>
void AdaptiveFETriangMesh<RealType, TriangleType>::prolongateConst( aol::Vector<RealType>& function, const RealType constant ) const {
  int oldSize = function.size();
  function.resize( this->getNumVertices() );
  
  std::map< int, aol::Vec2<int> >::const_iterator iter;
  for( int i = oldSize; i < function.size(); i++ ){
    // find new node index in map
    iter = _interpolationMap.find( i );
    if( iter == _interpolationMap.end() )
      throw aol::Exception ( "AdaptiveFETriangMesh<>::prolongateConst(): unknown vertex!", __FILE__, __LINE__ ); 
    // interpolate linearly
    function[i] = constant;
  }  
}

WARNING_ON( unused-result )

/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
template class AdaptiveFETriangMesh< float, aol::TriangBaseElement<float> >;
template class AdaptiveFETriangMesh< double, aol::TriangBaseElement<double> >;

template class AdaptiveFETriangMesh< float, aol::DKTPlateElement<float> >;
template class AdaptiveFETriangMesh< double, aol::DKTPlateElement<double> >;
