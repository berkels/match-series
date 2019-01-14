#include <multiDObject.h>

namespace qc {

template<>
void qc::MultiDStorageObject<qc::QC_3D>::reallocate ( const GridStructure &Grid ) {
  QUOC_ASSERT ( Grid.getDimOfWorld() == QC_3D );
  this->reallocate ( Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
}


template<>
void qc::MultiDStorageObject<qc::QC_3D>::resize ( const GridStructure &Grid ) {
  QUOC_ASSERT ( Grid.getDimOfWorld() == QC_3D );
  resize ( Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
}

}
