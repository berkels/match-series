#include <quadTree.h>

qc::QuadGridField::QuadGridField ( int Step, int NumX, int NumY ) {
  step = Step;
  numX = NumX;
  numY = NumY;

  bitField = new aol::BitVector ( NumX * NumY );
  if ( bitField == NULL ) {
    cerr << "ERROR in qc::QuadGridField::qc::QuadGridField: could not allocate the Bitfield!\n";
  }
}

qc::QuadGridField::~QuadGridField() {
  delete bitField;
}

void qc::QuadGridField::setElement ( int X, int Y ) {
  if ( X % step != 0 ) {
    cerr << "ERROR in qc::QuadGridField::setElement X = " << X << " not in this grid...\n";
  }
  if ( Y % step != 0 ) {
    cerr << "ERROR in qc::QuadGridField::setElement Y = " << Y << " not in this grid...\n";
  }

  int x = X / step;
  int y = Y / step;

  bitField->set ( y * numX + x , true );
}

void qc::QuadGridField::clearElement ( int X, int Y ) {
  if ( X % step != 0 ) {
    cerr << "ERROR in qc::QuadGridField::clearElement X = " << X << " not in this grid...\n";
  }
  if ( Y % step != 0 ) {
    cerr << "ERROR in qc::QuadGridField::clearElement Y = " << Y << " not in this grid...\n";
  }

  int x = X / step;
  int y = Y / step;

  bitField->set ( y * numX + x , true );
}

void qc::QuadGridField::clearAll() {
  bitField->setZero();
}

bool qc::QuadGridField::checkElement ( int X, int Y ) {
  if ( X % step != 0 ) {
    cerr << "ERROR in qc::QuadGridField::checkElement X = " << X << " not in this grid...\n";
  }
  if ( Y % step != 0 ) {
    cerr << "ERROR in qc::QuadGridField::checkElement Y = " << Y << " not in this grid...\n";
  }

  int x = X / step;
  int y = Y / step;

  return bitField->get ( y * numX + x );
}




qc::QuadTree::QuadTree ( int GridDepth ) {
  width = height = 1 << ( GridDepth );

  if ( GridDepth < 1 ) {
    cerr << "ERROR in qc::QuadTree::qc::qc::QuadTree: Depth < 1 \n";
  }

  maxLevel = GridDepth;

  // create all the gridfields now

  gridFields = new PQuadGridField[ GridDepth + 1];

  for ( int level = 0; level <= GridDepth; level++ ) {
    gridFields[ level ] = new qc::QuadGridField ( getStep ( level ),
                                              width / getStep ( level ),
                                              height / getStep ( level ) );
  }
}

qc::QuadTree::~QuadTree() {
  for ( int level = 0; level <= maxLevel; level++ ) {
    delete gridFields[ level ];
  }
  delete[] gridFields;
}

void qc::QuadTree::setElement ( int X, int Y, int Level ) {
  if ( Level > maxLevel ) {
    cerr << "ERROR in qc::QuadTree::setElement(): Level out of range...\n";
    return;
  }
  gridFields[ Level ]->setElement ( X, Y );
}

void qc::QuadTree::setElementAndUnite ( int X, int Y, int Level ) {
  setElement ( X, Y, Level );
  if ( Level > maxLevel ) return;

  int coarseX, coarseY;
  int coarseStep = getStep ( Level - 1 );

  coarseX = X - X % coarseStep;
  coarseY = Y - Y % coarseStep;

  if ( checkForUnite ( coarseX, coarseY, Level - 1 ) ) {
    uniteElements ( coarseX, coarseY, Level - 1 );
  }
}

void qc::QuadTree::setElementAndClearChildren ( int X, int Y, int Level ) {
  setElement ( X, Y, Level );
  if ( Level >= maxLevel ) return;

  int x, y, CoarseStep;

  CoarseStep = getStep ( Level );
  for ( int level = Level + 1; level <= maxLevel; level++ ) {
    int Step = getStep ( level );

    for ( x = X; x < X + CoarseStep; x += Step )
      for ( y = Y; y < Y + CoarseStep; y += Step )
        clearElement ( x, y, level ) ;
  }
}

bool qc::QuadTree::checkForUnite ( int X, int Y, int Level ) {
  if ( Level >= maxLevel ) return false;

  int Step = getStep ( Level + 1 );

  return ( getElement ( X, Y, Level + 1 )
           && getElement ( X + Step, Y, Level + 1 )
           && getElement ( X, Y + Step, Level + 1 )
           && getElement ( X + Step, Y + Step, Level + 1 ) );
}

void qc::QuadTree::clearElement ( int X, int Y ) {
  // search for grid with the qc::Element
  for ( int level = 0 ; level <= maxLevel; level++ ) {

    if ( ( X % getStep ( level ) == 0 ) && ( Y % getStep ( level ) == 0 ) ) {

      if ( gridFields[ level ]->checkElement ( X, Y ) == true ) {
        gridFields[ level ]->clearElement ( X, Y );
        level = maxLevel; // exit loop
      }

    }

  }
}

void qc::QuadTree::clearElement ( int X, int Y, int Level ) {
  if ( Level <= maxLevel ) {
    gridFields[ Level ]->clearElement ( X, Y );
  } else {
    cerr << "qc::QuadTree::clearElement Level not in range...\n";
  }
}

bool qc::QuadTree::getElement ( int X, int Y, int Level ) {
  if ( Level > maxLevel ) {
    cerr << "ERROR in getElement: Level = " << Level << " out of range... \n";
    return 0;
  }
  return gridFields[ Level ]->checkElement ( X, Y );
}

void qc::QuadTree::clearAll() {
  for ( int level = 0; level <= maxLevel; level++ ) {
    gridFields[ level ]->clearAll();
  }
}

bool qc::QuadTree::upNeighbour ( int X, int Y, int Level ) {
  int step = getStep ( Level );
  Y += step;

  if ( Y >= height - 1 ) return false;

  if ( ( Level > 0 ) && ( Y % ( 2*step ) == 0 ) && ( Y + 2*step < height ) &&
       getElement ( X - X % ( 2*step ), Y, Level - 1 ) ) return true;

  if ( ( Y + step < height ) && ( getElement ( X, Y, Level ) ) ) return true;

  if ( ( Level < maxLevel ) && ( Y + step / 2 < height ) &&
       ( ( getElement ( X, Y, Level + 1 ) || getElement ( X + step / 2, Y, Level + 1 ) ) ) ) return true;

  return false;

}

bool qc::QuadTree::lowNeighbour ( int X, int Y, int Level ) {
  int step = getStep ( Level );

  if ( Y == 0 ) return false;

  if ( ( Level > 0 ) && ( Y % ( 2*step ) == 0 ) && ( Y >= 2*step ) &&
       getElement ( X - X % ( 2*step ), Y - 2*step, Level - 1 ) ) return true;

  if ( ( Y >= step ) && ( getElement ( X, Y - step, Level ) ) ) return true;

  if ( ( Level < maxLevel ) && ( Y >= step / 2 ) &&
       ( ( getElement ( X, Y - step / 2, Level + 1 ) || getElement ( X + step / 2, Y - step / 2, Level + 1 ) ) ) ) return true;

  return false;
}

bool qc::QuadTree::leftNeighbour ( int X, int Y, int Level ) {
  int step = getStep ( Level );

  if ( X == 0 ) return false;

  if ( ( Level > 0 ) && ( X % ( step << 1 ) == 0 ) && ( X >= ( step << 1 ) ) &&
       getElement ( X - ( step << 1 ), Y - Y % ( step << 1 ), Level - 1 ) ) return true;

  if ( ( X >= step ) && ( getElement ( X - step, Y, Level ) ) ) return true;

  if ( ( Level < maxLevel ) && ( X >= step / 2 ) &&
       ( ( getElement ( X - step / 2, Y, Level + 1 ) || getElement ( X - step / 2, Y + step / 2, Level + 1 ) ) ) ) return true;

  return false;
}

bool qc::QuadTree::rightNeighbour ( int X, int Y, int Level ) {
  int step = getStep ( Level );
  X += step;

  if ( X >= width - 1 ) return false;

  if ( ( Level > 0 ) && ( X % ( 2*step ) == 0 ) && ( X + 2*step < width ) &&
       getElement ( X, Y - Y % ( 2*step ), Level - 1 ) ) return true;

  if ( ( X + step < width ) && ( getElement ( X, Y, Level ) ) ) return true;

  if ( ( Level < maxLevel ) && ( X + step / 2 < width ) &&
       ( ( getElement ( X, Y, Level + 1 ) || getElement ( X, Y + step / 2, Level + 1 ) ) ) ) return true;

  return false;
}

void qc::QuadTree::uniteElements ( int X, int Y, int Level ) {
  if ( Level > maxLevel ) return;

  clearElement ( X, Y, Level + 1 );
  clearElement ( X + getStep ( Level + 1 ), Y, Level + 1 );
  clearElement ( X, Y + getStep ( Level + 1 ), Level + 1 );
  clearElement ( X + getStep ( Level + 1 ), Y + getStep ( Level + 1 ), Level + 1 );
  setElementAndUnite ( X, Y, Level );

}


int qc::QuadTree::getCoarsestLevel ( int X, int Y ) {
  int Level, CoarseLevel = 0;

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

void qc::QuadTree::statistics ( void ) {
  int X, Y, Step, TotCells = 0, TotHierCells = 0, Cells, HierCells;
  cout << "L   HierCells  Cells" << endl;
  cout << "=== ========== ==========" << endl;
  for ( int Level = 0; Level <= maxLevel; Level ++ ) {

    Step = getStep ( Level );

    Cells = 0;
    HierCells = 0;

    for ( X = 0; X + Step < width; X += Step ) {
      for ( Y = 0; Y + Step < height; Y += Step ) {
        if ( getElement ( X, Y, Level ) ) {
          Cells += Step * Step;
          HierCells++;
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


