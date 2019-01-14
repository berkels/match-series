#include <imageTools.h>

namespace qc {

void cleanMask ( qc::BitArray<qc::QC_2D> &Mask, const int DropComponentsSmallerThan, const bool DropBoundaryComponents ) {
  qc::GridStructure grid ( qc::GridSize2d::createFrom ( Mask ) );
  qc::ScalarArray<int, qc::QC_2D> labelArray ( grid );
  const int numberOfLabels = qc::ConnectedComponentsLabeler::doLabel ( Mask, labelArray );

  std::set<int> boundaryComponents;
  if ( DropBoundaryComponents ) {
    for ( qc::GridStructure::AllBoundaryNodeIterator it = grid; it.notAtEnd(); ++it ) {
      if ( Mask.get ( *it ) )
        boundaryComponents.insert ( labelArray.get ( *it ) );
    }
  }

  for ( int i = 0; i < numberOfLabels; ++i ) {
    const int componentIndex = i + 1;
    const int componentSize = labelArray.numOccurence ( componentIndex );
    if ( ( componentSize < DropComponentsSmallerThan ) || ( boundaryComponents.find ( componentIndex ) != boundaryComponents.end() ) ) {
      for ( int j = 0; j < Mask.size(); ++j ) {
        if ( labelArray[j] == ( i + 1 ) )
          Mask.set ( j, false );
      }
    }
  }
}

void cleanMask3D ( qc::BitArray<qc::QC_3D> &Mask, const int DropComponentsSmallerThan, const bool DropBoundaryComponents ) {
  qc::GridStructure grid ( qc::GridSize3d::createFrom ( Mask ) );
  qc::ScalarArray<int, qc::QC_3D> labelArray ( grid );
  aol::Vector<unsigned int> labelSize(Mask.size());
  const int numberOfLabels = qc::ConnectedComponentsLabeler::doLabel3D<0, 0, 0, 1> ( Mask, labelArray, labelSize );
  cerr << "numberOfLabels = " << numberOfLabels << endl;
  std::set<int> boundaryComponents;
  if ( DropBoundaryComponents ) {
    for ( qc::GridStructure::AllBoundaryNodeIterator it = grid; it.notAtEnd(); ++it ) {
      if ( Mask.get ( *it ) )
        boundaryComponents.insert ( labelArray.get ( *it ) );
    }
  }
    
  for ( int i = 0; i < numberOfLabels; ++i ) {
    const int componentIndex = i + 1;
    const int componentSize = labelArray.numOccurence ( componentIndex );
    if ( ( componentSize < DropComponentsSmallerThan ) || ( boundaryComponents.find ( componentIndex ) != boundaryComponents.end() ) ) {
      for ( int j = 0; j < Mask.size(); ++j ) {
        if ( labelArray[j] == ( i + 1 ) )
          Mask.set ( j, false );
      }
    }
  }
}

void getBiggestComponent ( qc::BitArray<qc::QC_2D> &Mask ) {
  qc::GridStructure grid ( qc::GridSize2d::createFrom ( Mask ) );
  qc::ScalarArray<int, qc::QC_2D> labelArray ( grid );
  const int numberOfLabels = qc::ConnectedComponentsLabeler::doLabel ( Mask, labelArray );
  int biggestComponentSize = -1;
  aol::Vector<int> componentSize(numberOfLabels);
  for ( int i = 0; i < numberOfLabels; ++i ) {
    const int componentIndex = i + 1;
    componentSize[i] = labelArray.numOccurence ( componentIndex );
  }
  biggestComponentSize = componentSize.getMaxValue();
  for ( int i = 0; i < numberOfLabels; ++i ) {
    if ( ( componentSize[i] < biggestComponentSize ) ) {
      for ( int j = 0; j < Mask.size(); ++j ) {
        if ( labelArray[j] == ( i + 1 ) )
          Mask.set ( j, false );
      }
    }
  }
}

void convertMaskToSet ( const qc::BitArray<qc::QC_2D> &Mask, std::set<aol::Vec<2, short> > &MaskSet ) {
  MaskSet.clear();
  for ( qc::RectangularIterator<qc::QC_2D, aol::Vec<2, short> > it ( Mask ); it.notAtEnd(); ++it ) {
    if ( Mask.get ( *it ) )
      MaskSet.insert ( *it );
  }
}

int computeMean ( const aol::Vector<int> &Histo, int StartIndex, int EndIndex ) {
  int firstMoment = 0;
  for ( int i = StartIndex; i <= EndIndex; ++i )
    firstMoment += Histo[i];
  return firstMoment;
}

int computeFirstMoment ( const aol::Vector<int> &Histo, int StartIndex, int EndIndex ) {
  int firstMoment = 0;
  for ( int i = StartIndex; i <= EndIndex; ++i )
    firstMoment += i * Histo[i];
  return firstMoment;
}

} // end of namespace qc.
