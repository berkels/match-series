#ifndef __OCTREE_H
#define __OCTREE_H

#include <quadTree.h>

#include <quoc.h>
#include <aol.h>

namespace qc {

class OcTree;

// necessary?:
typedef unsigned char byte;

const int OC_NORMAL   = 1;
const int OC_COARSE   = 2;
const int OC_FINE     = 4;
const int OC_ALL_FINE =  8;
const int OC_COARSE_NORMAL = OC_COARSE | OC_NORMAL;
const int OC_ALL           = OC_FINE | OC_COARSE | OC_NORMAL;

const int TN_DIR_RIGHT = 1;
const int TN_DIR_UP    = 2;
const int TN_DIR_BACK  = 4;

class OcGridField {
  friend class OcTree;
protected:
  inline OcGridField ( unsigned int Step, unsigned int NumX, unsigned int NumY, unsigned int NumZ );

  inline ~OcGridField();

  void setElement ( unsigned int X, unsigned int Y, unsigned int Z ) {
    /*if ( X % step != 0 ) {
      cerr << "ERROR in OcGridField::setElement X = " << X << " not in this grid...\n";
      }
      if ( Y % step != 0 ) {
      cerr << "ERROR in OcGridField::setElement Y = " << Y << " not in this grid...\n";
      }
      if ( Z % step != 0 ) {
      cerr << "ERROR in OcGridField::setElement Z = " << Z << " not in this grid...\n";
      } */

    byte x = X >> base;
    byte y = Y >> base;
    byte z = Z >> base;
    byte byteMask = 1;
    if ( X & checkBit ) byteMask <<= 1;
    if ( Y & checkBit ) byteMask <<= 2;
    if ( Z & checkBit ) byteMask <<= 4;

    bitField[ ( ( ( z << Ypot ) + y ) << Xpot ) + x ] |= byteMask;
  }

  void clearElement ( unsigned int X, unsigned int Y, unsigned int Z ) {
    /*
      if ( X % step != 0 ) {
      cerr << "ERROR in OcGridField::clearElement X = " << X << " not in this grid...\n";
      }
      if ( Y % step != 0 ) {
      cerr << "ERROR in OcGridField::clearElement Y = " << Y << " not in this grid...\n";
      }
      if ( Z % step != 0 ) {
      cerr << "ERROR in OcGridField::clearElement Z = " << Z << " not in this grid...\n";
      }*/

    byte x = X >> base;
    byte y = Y >> base;
    byte z = Z >> base;
    byte byteMask = 1;
    if ( X & checkBit ) byteMask <<= 1;
    if ( Y & checkBit ) byteMask <<= 2;
    if ( Z & checkBit ) byteMask <<= 4;

    bitField[ ( ( ( z << Ypot ) + y ) << Xpot ) + x ] &= ~byteMask;
  }


  bool checkElement ( unsigned int X, unsigned int Y, unsigned int Z ) {
    /*if ( X % step != 0 ) {
      cerr << "ERROR in OcGridField::checkElement X = " << X << " not in this grid...\n";
      }
      if ( Y % step != 0 ) {
      cerr << "ERROR in OcGridField::checkElement Y = " << Y << " not in this grid...\n";
      }
      if ( Z % step != 0 ) {
      cerr << "ERROR in OcGridField::checkElement Z = " << Z << " not in this grid...\n";
      } */


    byte x = X >> base;
    byte y = Y >> base;
    byte z = Z >> base;
    byte byteMask = 1;
    if ( X & checkBit ) byteMask <<= 1;
    if ( Y & checkBit ) byteMask <<= 2;
    if ( Z & checkBit ) byteMask <<= 4;
    return ( ( bitField[ ( ( ( z << Ypot ) + y ) << Xpot ) + x ] & byteMask ) != 0 );
  }

  byte getBitMask ( unsigned int X, unsigned int Y, unsigned int Z ) {
    byte x = X >> base;
    byte y = Y >> base;
    byte z = Z >> base;
    unsigned int byteIndex = ( ( ( z << Ypot ) + y ) << Xpot ) + x;
    return bitField[ byteIndex ];
  }

  inline void save ( ostream &out );

  inline void read ( istream &in );

  inline void clearAll();

private:

  unsigned checkBit;
  unsigned base;
  unsigned Xpot, Ypot;
  unsigned step, numX, numY, numZ;
  byte *bitField;
};

OcGridField::OcGridField ( unsigned Step, unsigned NumX, unsigned NumY, unsigned NumZ ) {
  step = Step;
  numX = NumX;
  numY = NumY;
  numZ = NumZ;

  base = qc::logBaseTwo ( Step ) + 1;
  Xpot = qc::logBaseTwo ( numX ) - 1;
  Ypot = qc::logBaseTwo ( numY ) - 1;
  checkBit = 1 << ( base - 1 );

  bitField = new byte[ ( NumX * NumY * NumZ ) >> 3 ];
  if ( bitField == NULL ) {
    cerr << "ERROR in OcGridField::OcGridField: could not allocate the Bitfield!\n";
  }
}

OcGridField::~OcGridField() {
  if ( bitField )
    delete[] bitField;
}

void OcGridField::clearAll() {
  memset ( bitField, 0, ( numX * numY * numZ ) >> 3 );
}

void OcGridField::save ( ostream &out ) {
  out.write ( reinterpret_cast< char* >( &step ), sizeof ( unsigned int ) );
  out.write ( reinterpret_cast< char* >( &numX ), sizeof ( unsigned int ) );
  out.write ( reinterpret_cast< char* >( &numY ), sizeof ( unsigned int ) );
  out.write ( reinterpret_cast< char* >( &numZ ), sizeof ( unsigned int ) );
  out.write ( reinterpret_cast< char* >( bitField ), ( numX * numY * numZ ) >> 3 );
}

void OcGridField::read ( istream &in ) {
  in.read ( reinterpret_cast< char* >( &step ), sizeof ( unsigned int ) );
  in.read ( reinterpret_cast< char* >( &numX ), sizeof ( unsigned int ) );
  in.read ( reinterpret_cast< char* >( &numY ), sizeof ( unsigned int ) );
  in.read ( reinterpret_cast< char* >( &numZ ), sizeof ( unsigned int ) );
  in.read ( reinterpret_cast< char* >( bitField ), ( numX * numY * numZ ) >> 3 );
}

class OcTree {
public:
  inline OcTree ( unsigned NumX, unsigned NumY, unsigned NumZ, unsigned MaxLevel );

  inline ~OcTree();

  inline void setElement ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  inline void setElement ( int Dir, unsigned X, unsigned Y, unsigned Z, unsigned Level );

  int nodeInGrid ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
    int Step = 1 << ( maxLevel - Level );
    return ( ( X % Step ) == 0 && ( Y % Step ) == 0 && ( Z % Step ) == 0 );
  }

  void dumpElement ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
    cerr << " - Cell = ( " << X << ", " << Y << ", " << Z << " ) Level = " << Level << endl;
  }

  inline void setElementAndUnite ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  inline void setElementAndClearChildren ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  inline bool checkAllChildren ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  inline bool checkChildren ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  inline bool checkParent ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  inline void clearElement ( unsigned X, unsigned Y, unsigned Z );

  inline void clearElement ( unsigned X, unsigned Y, unsigned Z, unsigned Level );

  bool getElement ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
    if ( Level > maxLevel ) {
      cerr << "ERROR in getElement: Level = " << Level << " out of range... \n";
      return 0;
    }
    return gridFields[ Level ]->checkElement ( X, Y, Z );
  }

  byte   getChildMask ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
    return gridFields[ Level + 1 ]->getBitMask ( X, Y, Z );
  }

  inline bool upNeighbour ( unsigned X, unsigned Y, unsigned Z, unsigned Level, char Mode = OC_ALL );

  inline bool lowNeighbour ( unsigned X, unsigned Y, unsigned Z, unsigned Level, char Mode = OC_ALL );

  inline bool leftNeighbour ( unsigned X, unsigned Y, unsigned Z, unsigned Level, char Mode = OC_ALL );

  inline bool rightNeighbour ( unsigned X, unsigned Y, unsigned Z, unsigned Level, char Mode = OC_ALL );

  inline bool frontNeighbour ( unsigned X, unsigned Y, unsigned Z, unsigned Level, char Mode = OC_ALL );

  inline bool backNeighbour ( unsigned X, unsigned Y, unsigned Z, unsigned Level, char Mode = OC_ALL );

  inline void   uniteElements ( unsigned X, unsigned Y, unsigned Z, unsigned int Level );

  inline void clearAll();

  unsigned int getStep ( unsigned int Level ) {
    if ( Level > maxLevel ) {
      cerr << "ERROR in OcTree::getStep! Level = " << Level << " out of range..\n";
    }
    return ( 1 << ( maxLevel - Level ) );
  }

  // returns the coarsest(lowest) number of level in which the node (X,Y) appears
  inline unsigned int getCoarsestLevel ( unsigned X, unsigned Y, unsigned Z );

  unsigned int getMaxLevel() {
    return maxLevel;
  }

  unsigned int getNumX() {
    return numX;
  }

  unsigned int getNumY() {
    return numY;
  }

  unsigned int getNumZ() {
    return numZ;
  }

  inline void statistics ( void );

  inline void save ( char *FileName );

  inline void read ( char *FileName );
protected:
  byte leftNbMask;
  byte rightNbMask;
  byte backNbMask;
  byte frontNbMask;
  byte lowNbMask;
  byte upNbMask;

  OcGridField **gridFields;

  unsigned int numX, numY, numZ, maxLevel;
};


typedef OcGridField* POcGridField;

bool OcTree::upNeighbour ( unsigned int X, unsigned int Y, unsigned int Z, unsigned int Level, char Mode ) {
  unsigned int step = getStep ( Level );
  unsigned coarseStep = step << 1;
  unsigned coarseMask = ~ ( 1 << ( maxLevel - Level ) );

  Y += step;

  if ( Y >= numY ) return true;

  // coarser level check
  if ( ( Mode & OC_COARSE ) && ( Level > 0 ) && ( Y % ( coarseStep ) == 0 ) &&
       getElement ( X & coarseMask, Y, Z & coarseMask, Level - 1 ) ) return true;

  // this level
  if ( ( Mode & OC_NORMAL ) && ( getElement ( X, Y, Z, Level ) ) ) return true;

  // next finer level
  if ( Level < maxLevel ) {
    byte childMask = getChildMask ( X, Y, Z, Level );
    if ( Mode & OC_FINE ) {
      if ( childMask & lowNbMask )
        return true;
    } else if ( Mode & OC_ALL_FINE ) {
      if ( (childMask & lowNbMask) == lowNbMask )
        return true;
    }
  }

  return false;

}

bool OcTree::lowNeighbour ( unsigned int X, unsigned int Y, unsigned int Z, unsigned int Level, char Mode ) {
  unsigned int step = getStep ( Level );
  unsigned coarseStep = step << 1;
  unsigned coarseMask = ~ ( 1 << ( maxLevel - Level ) );

  if ( Y <= 0 ) return true;

  // coarser level check
  if ( ( Mode & OC_COARSE ) && ( Level > 0 ) && ( Y % ( coarseStep ) == 0 ) && ( Y >= coarseStep ) &&
       getElement ( X & coarseMask, Y - coarseStep, Z & coarseMask, Level - 1 ) ) return true;

  // this level
  if ( ( Mode & OC_NORMAL ) && ( Y >= step ) && ( getElement ( X, Y - step, Z,  Level ) ) ) return true;

  // next finer level
  if ( ( Level < maxLevel ) && ( Y >= step ) ) {
    byte childMask = getChildMask ( X, Y - step , Z, Level );
    if ( Mode & OC_FINE ) {
      if ( childMask & upNbMask )
        return true;
    } else if ( Mode & OC_ALL_FINE ) {
      if ( (childMask & upNbMask) == upNbMask )
        return true;
    }
  }

  return false;
}

bool OcTree::leftNeighbour ( unsigned int X, unsigned int Y, unsigned int Z, unsigned int Level, char Mode ) {
  unsigned int step = getStep ( Level );
  unsigned coarseStep = step << 1;
  unsigned coarseMask = ~ ( 1 << ( maxLevel - Level ) );


  if ( X <= 0 ) return true;

  // coarser level check
  if ( ( Mode & OC_COARSE ) && ( Level > 0 ) && ( X % ( coarseStep ) == 0 ) && ( X >= coarseStep ) &&
       getElement ( X - coarseStep, Y & coarseMask, Z & coarseMask, Level - 1 ) ) return true;

  // this level
  if ( ( Mode & OC_NORMAL ) && ( X >= step ) && ( getElement ( X - step, Y, Z, Level ) ) ) return true;

  // next finer level
  if ( ( Level < maxLevel ) && ( X >= step ) ) {
    byte childMask = getChildMask ( X - step, Y, Z, Level );
    if ( Mode & OC_FINE ) {
      if ( childMask & rightNbMask )
        return true;
    } else if ( Mode & OC_ALL_FINE ) {
      if ( (childMask & rightNbMask) == rightNbMask )
        return true;
    }
  }

  return false;
}

bool OcTree::rightNeighbour ( unsigned int X, unsigned int Y, unsigned int Z, unsigned int Level, char Mode ) {
  unsigned int step = getStep ( Level );
  unsigned coarseStep = step << 1;
  unsigned coarseMask = ~ ( 1 << ( maxLevel - Level ) );

  X += step;

  if ( X >= numX ) return true;

  // coarser level check
  if ( ( Mode & OC_COARSE ) && ( Level > 0 ) && ( X % ( coarseStep ) == 0 ) &&
       getElement ( X, Y & coarseMask, Z & coarseMask, Level - 1 ) ) return true;

  // this level
  if ( ( Mode & OC_NORMAL ) && ( getElement ( X, Y, Z, Level ) ) ) return true;

  // next finer level
  if ( Level < maxLevel ) {
    byte childMask = getChildMask ( X, Y, Z, Level );
    if ( Mode & OC_FINE ) {
      if ( childMask & leftNbMask )
        return true;
    } else if ( Mode & OC_ALL_FINE ) {
      if ( (childMask & leftNbMask) == leftNbMask )
        return true;
    }
  }

  return false;
}

bool OcTree::frontNeighbour ( unsigned int X, unsigned int Y, unsigned int Z, unsigned int Level, char Mode ) {
  unsigned int step = getStep ( Level );
  unsigned coarseStep = step << 1;
  unsigned coarseMask = ~ ( 1 << ( maxLevel - Level ) );


  if ( Z == 0 ) return true;

  // coarser level check
  if ( ( Mode & OC_COARSE ) && ( Level > 0 ) && ( Z % ( coarseStep ) == 0 ) && ( Z >= coarseStep ) &&
       getElement ( X & coarseMask, Y & coarseMask, Z - coarseStep, Level - 1 ) ) return true;

  // this level
  if ( ( Mode & OC_NORMAL ) && ( Z >= step ) && ( getElement ( X, Y, Z - step, Level ) ) ) return true;

  // next finer level
  if ( ( Level < maxLevel ) && ( Z >= step ) ) {
    byte childMask = getChildMask ( X, Y, Z - step, Level );
    if ( Mode & OC_FINE ) {
      if ( childMask & backNbMask )
        return true;
    } else if ( Mode & OC_ALL_FINE ) {
      if ( (childMask & backNbMask) == backNbMask )
        return true;
    }
  }
  return false;
}

bool OcTree::backNeighbour ( unsigned int X, unsigned int Y, unsigned int Z, unsigned int Level, char Mode ) {
  unsigned int step = getStep ( Level );
  unsigned coarseStep = step << 1;
  unsigned coarseMask = ~ ( 1 << ( maxLevel - Level ) );

  Z += step;

  if ( Z >= numZ ) return true;

  // coarser level check
  if ( ( Mode & OC_COARSE ) && ( Level > 0 ) && ( Z % ( coarseStep ) == 0 ) &&
       getElement ( X & coarseMask, Y & coarseMask, Z, Level - 1 ) ) return true;

  // this level
  if ( ( Mode & OC_NORMAL ) && ( getElement ( X, Y, Z, Level ) ) ) return true;

  // next finer level
  if ( Level < maxLevel ) {
    byte childMask = getChildMask ( X, Y, Z, Level );
    if ( Mode & OC_FINE ) {
      if ( childMask & frontNbMask )
        return true;
    } else if ( Mode & OC_ALL_FINE ) {
      if ( (childMask & frontNbMask) == frontNbMask )
        return true;
    }
  }
  return false;
}


void OcTree::setElement ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  /*if ( checkChildren( X, Y, Z, Level ) ) {
    cerr << "WARNING in OcTree::setElement: a child of this element exists already!\n";
    dumpElement( X, Y, Z, Level );
  }

  if ( checkParent( X, Y, Z, Level ) ) {
    cerr << "WARNING in OcTree::setElement: parent of this element exists already!\n";
    dumpElement( X, Y, Z, Level );
  }

  if ( Level > maxLevel ) {
    cerr << "ERROR in OcTree::setElement( " << X << ", " << Y << ", " << Z << ", "
  << Level << " ): Level out of range...\n";
    return;
  } */
  gridFields[ Level ]->setElement ( X, Y, Z );
}

void OcTree::setElement ( int Dir, unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  int Step = 1 << ( maxLevel - Level );
  int NewX = X, NewY = Y, NewZ = Z;
  if ( ! ( Dir & TN_DIR_RIGHT ) ) {
    NewX -= Step;
  }
  if ( ! ( Dir & TN_DIR_UP ) ) {
    NewY -= Step;
  }
  if ( ! ( Dir & TN_DIR_BACK ) ) {
    NewZ -= Step;
  }
  setElement ( NewX, NewY, NewZ, Level );
}

bool OcTree::checkAllChildren ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  if ( Level >= maxLevel ) return false;

  unsigned Step = 1 << ( maxLevel - ( Level + 1 ) );

  return ( getElement ( X, Y, Z, Level + 1 )
           && getElement ( X + Step, Y, Z, Level + 1 )
           && getElement ( X, Y + Step, Z, Level + 1 )
           && getElement ( X + Step, Y + Step, Z, Level + 1 )
           && getElement ( X, Y, Z + Step, Level + 1 )
           && getElement ( X + Step, Y, Z + Step, Level + 1 )
           && getElement ( X, Y + Step, Z + Step, Level + 1 )
           && getElement ( X + Step, Y + Step, Z + Step, Level + 1 ) );
}

bool OcTree::checkChildren ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  if ( Level >= maxLevel ) return false;

  unsigned Step = 1 << ( maxLevel - ( Level + 1 ) );

  return ( getElement ( X, Y, Z, Level + 1 )
           || getElement ( X + Step, Y, Z, Level + 1 )
           || getElement ( X, Y + Step, Z, Level + 1 )
           || getElement ( X + Step, Y + Step, Z, Level + 1 )
           || getElement ( X, Y, Z + Step, Level + 1 )
           || getElement ( X + Step, Y, Z + Step, Level + 1 )
           || getElement ( X, Y + Step, Z + Step, Level + 1 )
           || getElement ( X + Step, Y + Step, Z + Step, Level + 1 ) );
}

bool OcTree::checkParent ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  if ( Level <= 0 ) return false;

  unsigned Step = 1 << ( maxLevel - ( Level - 1 ) );

  return getElement ( X - X % Step, Y - Y % Step, Z - Z % Step, Level - 1 );
}

void OcTree::clearElement ( unsigned X, unsigned Y, unsigned Z ) {
  // search for grid with the Element
  for ( unsigned level = 0 ; level <= maxLevel; level++ ) {

    if ( ( X % getStep ( level ) == 0 )
         && ( Y % getStep ( level ) == 0 )
         && ( Z % getStep ( level ) == 0 ) ) {

      if ( gridFields[ level ]->checkElement ( X, Y, Z ) == true ) {
        gridFields[ level ]->clearElement ( X, Y, Z );
        level = maxLevel; // exit loop
      }
    }
  }
}

void OcTree::clearElement ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  if ( Level <= maxLevel ) {
    gridFields[ Level ]->clearElement ( X, Y, Z );
  } else {
    cerr << "OcTree::clearElement Level not in range...\n";
  }
}

OcTree::OcTree ( unsigned NumX, unsigned NumY, unsigned NumZ,
                 unsigned MaxLevel ) {
  numX = NumX;
  numY = NumY;
  numZ = NumZ;

  if ( MaxLevel == 0 ) {
    cerr << "ERROR in OcTree::qc::OcTree: GridNumZ < 1 \n";
  }

  maxLevel = MaxLevel;

  // create all the gridfields now
  gridFields = new POcGridField[ MaxLevel + 1 ];

  for ( unsigned level = 0; level <= MaxLevel ; level++ ) {
    gridFields[ level ] = new OcGridField ( getStep ( level ),
                                            NumX / getStep ( level ),
                                            NumY / getStep ( level ),
                                            NumZ / getStep ( level ) );
  }

  leftNbMask = ( 1 << 0 ) | ( 1 << 2 ) | ( 1 << 4 ) | ( 1 << 6 );
  rightNbMask =  leftNbMask << 1;
  lowNbMask = ( 1 << 0 ) | ( 1 << 1 ) | ( 1 << 4 ) | ( 1 << 5 );
  upNbMask = lowNbMask << 2;
  frontNbMask = ( 1 << 0 ) | ( 1 << 1 ) | ( 1 << 2 ) | ( 1 << 3 );
  backNbMask = frontNbMask << 4;
}

OcTree::~OcTree() {
  for ( unsigned level = 0; level <= maxLevel; level++ ) {
    delete gridFields[ level ];
  }
  delete[] gridFields;
}

void OcTree::setElementAndUnite ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  setElement ( X, Y, Z, Level );
  if ( Level > maxLevel ) return;

  unsigned coarseX, coarseY, coarseZ;
  unsigned coarseStep = getStep ( Level - 1 );

  coarseX = X - X % coarseStep;
  coarseY = Y - Y % coarseStep;
  coarseZ = Z - Z % coarseStep;

  if ( checkAllChildren ( coarseX, coarseY, coarseZ, Level - 1 ) ) {
    uniteElements ( coarseX, coarseY, coarseZ, Level - 1 );
  }
}

void OcTree::setElementAndClearChildren ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  setElement ( X, Y, Z, Level );
  if ( Level >= maxLevel ) return;


  unsigned x, y, z, CoarseStep;

  CoarseStep = getStep ( Level );
  for ( unsigned level = Level + 1; level <= maxLevel; level++ ) {
    unsigned Step = getStep ( level );

    for ( x = X; x < X + CoarseStep; x += Step )
      for ( y = Y; y < Y + CoarseStep; y += Step )
        for ( z = Z; z < Z + CoarseStep; z += Step )
          clearElement ( x, y, z, level ) ;
  }
}

void OcTree::clearAll() {
  for ( unsigned level = 0; level <= maxLevel; level++ ) {
    gridFields[ level ]->clearAll();
  }
}

void OcTree::uniteElements ( unsigned X, unsigned Y, unsigned Z, unsigned Level ) {
  if ( Level > maxLevel ) return;

  unsigned HalfStep = + getStep ( Level + 1 );

  clearElement ( X, Y, Z, Level + 1 );
  clearElement ( X + HalfStep, Y, Z, Level + 1 );
  clearElement ( X, Y + HalfStep, Z, Level + 1 );
  clearElement ( X + HalfStep, Y + HalfStep, Z,  Level + 1 );

  clearElement ( X, Y, Z + HalfStep, Level + 1 );
  clearElement ( X + HalfStep, Y, Z + HalfStep, Level + 1 );
  clearElement ( X, Y + HalfStep, Z + HalfStep, Level + 1 );
  clearElement ( X + HalfStep, Y + HalfStep, Z + HalfStep,  Level + 1 );

  setElement ( X, Y, Z, Level );

}

unsigned OcTree::getCoarsestLevel ( unsigned X, unsigned Y, unsigned /*Z*/ ) {
  unsigned Level, CoarseLevel = 0;

  for ( Level = 0; Level <= maxLevel; Level ++ ) {
    if ( ( X % getStep ( Level ) ) == 0 ) {
      CoarseLevel = Level;
      break;
    }
  }

  if ( CoarseLevel < maxLevel ) {
    for ( Level = CoarseLevel; Level <= maxLevel; Level ++ ) {
      if ( ( Y % getStep ( Level ) ) == 0 ) {
        CoarseLevel = Level;
        break;
      }
    }
  }

  return CoarseLevel;
}

void OcTree::statistics ( void ) {
  unsigned X, Y, Z, Step, TotCells = 0, TotHierCells = 0, Cells, HierCells;
  cout << "L   HierCells  Cells" << endl;
  cout << "=== ========== ==========" << endl;
  for ( unsigned Level = 0; Level <= maxLevel; Level ++ ) {

    Step = getStep ( Level );

    Cells = 0;
    HierCells = 0;

    for ( X = 0; X + Step < numX; X += Step ) {
      for ( Y = 0; Y + Step < numY; Y += Step ) {
        for ( Z = 0; Z + Step < numZ; Z += Step ) {
          if ( getElement ( X, Y, Z, Level ) ) {
            Cells += Step * Step * Step;
            HierCells++;
          }
        }
      }
    }
    TotCells += Cells;
    TotHierCells += HierCells;

    cout << setiosflags ( ios::right ) << setw ( 3 ) << Level
    << " " << setw ( 10 ) << HierCells << " " << setiosflags ( ios::right ) << setw ( 10 ) << Cells << " " << endl;
  }
  cout << setiosflags ( ios::right ) << "TOT:" << setw ( 10 ) << TotHierCells
  << " " << setiosflags ( ios::right ) << setw ( 10 ) << TotCells << endl
  << " eff = " << float ( TotHierCells ) / float ( TotCells ) << endl;
}

void OcTree::save ( char *FileName ) {
  ofstream out ( FileName, ios::binary );

  out.write ( reinterpret_cast< char* >( &numX ), sizeof ( unsigned ) );
  out.write ( reinterpret_cast< char* >( &numY ), sizeof ( unsigned ) );
  out.write ( reinterpret_cast< char* >( &numZ ), sizeof ( unsigned ) );
  out.write ( reinterpret_cast< char* >( &maxLevel ), sizeof ( unsigned ) );

  for ( unsigned int l = 0; l <= maxLevel; l++ ) {
    gridFields[ l ]->save ( out );
  }
}

void OcTree::read ( char *FileName ) {
  ifstream in ( FileName, ios::binary );

  if ( in.fail() ) {
    cerr << "could not open filename " << FileName << endl;
    return;
  }

  in.read ( reinterpret_cast< char* >( &numX ), sizeof ( unsigned ) );
  in.read ( reinterpret_cast< char* >( &numY ), sizeof ( unsigned ) );
  in.read ( reinterpret_cast< char* >( &numZ ), sizeof ( unsigned ) );
  in.read ( reinterpret_cast< char* >( &maxLevel ), sizeof ( unsigned ) );

  for ( unsigned int l = 0; l <= maxLevel; l++ ) {
    gridFields[ l ]->read ( in );
  }
}

}

#endif

