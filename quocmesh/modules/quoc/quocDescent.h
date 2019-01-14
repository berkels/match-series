#ifndef __QUOCDESCENT_H
#define __QUOCDESCENT_H

#include <quoc.h>
#include <gridBase.h>

namespace qc {

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class MultilevelDescentInterface {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  typedef typename ConfiguratorType::InitType InitType;
  const int _maxDepth;
  const InitType _grid;
  const InitType *_curGrid;
  int _curLevel;
  bool _logLevelStart;
  vector<const InitType*> _grids;
public:
  MultilevelDescentInterface ( const int MaxDepth )
   : _maxDepth ( MaxDepth ),
     _grid ( MaxDepth, ConfiguratorType::Dim ),
     _curGrid ( NULL ),
     _curLevel ( MaxDepth ),
     _logLevelStart ( false ) {

    for ( int level = 0; level <= _grid.getGridDepth(); level++ ) {
      _grids.push_back( new InitType( level, ConfiguratorType::Dim ) );
    }

    _curGrid = _grids[ _curLevel ];
  }

  virtual ~MultilevelDescentInterface( ) {
    for ( int level = 0; level <= _grid.getGridDepth(); level++ ) {
     delete _grids[ level ];
   }
  }

  RealType H() const {
    return _curGrid->H();
  }

  int getLevel( ) const {
    return _curLevel;
  }

  int getMaxGridDepth( ) const {
    return _grid.getGridDepth();
  }

  const InitType& getInitializerRef( ) const {
    return _grid;
  }

  const InitType& getCurrentGrid () const {
    return *this->_curGrid;
  }

  virtual void setLevel( const int Level ) {
    QUOC_ASSERT ( (Level <= _grid.getGridDepth()) && (Level >= 0) );
    _curLevel = Level;
    _curGrid = _grids[ _curLevel ];
  }

  void setLogLevelStart ( const bool LogLevelStart ) {
    _logLevelStart = LogLevelStart;
  }

  virtual void prolongate( ) {
      throw aol::UnimplementedCodeException ( "prolongate must be implemented in the derived class first", __FILE__, __LINE__ );
  };

  virtual void descentOnCurrentGrid( ) {
      throw aol::UnimplementedCodeException ( "descentOnCurrentGrid must be implemented in the derived class first", __FILE__, __LINE__ );
  };
  
  void solve( const int StartLevel, const int StopLevel, const int LevelIncrement = 1 ) {
    if ( StartLevel > StopLevel ) {
      cerr << "StartLevel > StopLevel. Nothing to do.\n";
      return;
    }

    setLevel( StartLevel );
    for ( int level = StartLevel; level <= StopLevel; ) {
      if ( _logLevelStart ) {
        cerr << "\n--------------------------------------------------------------------------------\n";
        cerr << "Descent on level " << getLevel() << " started";
        cerr << "\n";
        cerr << "--------------------------------------------------------------------------------\n\n";
      }

      descentOnCurrentGrid();
      for ( int i = 0; i < LevelIncrement; i++ ) {
        if ( level < StopLevel )
          prolongate( );
        ++level;
      }
    }
  }
};

} // end namespace qc

#endif // __QUOCDESCENT_H
