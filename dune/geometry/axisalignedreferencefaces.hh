#ifndef DUNE_GEOMETRY_AXISALIGNEDREFERENCEFACES_HH
#define DUNE_GEOMETRY_AXISALIGNEDREFERENCEFACES_HH

#include <cassert>

#include <memory>

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 5)
#include <dune/geometry/type.hh>
#else
#include <dune/geometry/genericgeometry/topologytypes.hh>
#endif

namespace Dune
{

  // axisAlignedReferenceFaces
  // -------------------------

  inline void axisAlignedReferenceFaces ( int dim, unsigned int *faceIndices, unsigned int *numFaces )
  {
    // initialize to dimension 0
    numFaces[ 0 ] = 0u;
    unsigned int numTopologies = 1u;

    for( int d = 0; d < dim; ++d )
    {
      for( unsigned int i = 0; i < numTopologies; ++i )
      {
        // setup the prism (index: numTopologies+i)
        for( int j = 0; j < d; ++j )
          faceIndices[ dim*(numTopologies + i) + j ] = faceIndices[ dim*i + j ];
        faceIndices[ dim*(numTopologies + i) + d ] = numFaces[ i ];
        numFaces[ numTopologies + i ] = numFaces[ i ]+2u;

        // setup the pyramid (index: numTopologies)
        for( int j = 0; j < d; ++j )
          ++faceIndices[ dim*i + j ];
        faceIndices[ dim*i + d ] = 0u;
        numFaces[ i ] += 1u + static_cast< unsigned int >( d == 0 );
      }

      // update numTopoloties
      numTopologies *= 2u;
    }
  }

  inline void axisAlignedReferenceFaces ( int dim, unsigned int *faceIndices )
  {
#if DUNE_VERSION_NEWER( DUNE_GEOMETRY, 2, 5 )
    const unsigned int numTopo = Impl::numTopologies( dim );
#else
    const unsigned int numTopo = GenericGeometry::numTopologies( dim );
#endif
    std::unique_ptr< unsigned int[] > numFaces( new unsigned int[ numTopo ] );
    axisAlignedReferenceFaces( dim, faceIndices, numFaces.get() );
  }

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_AXISALIGNEDREFERENCEFACES_HH
