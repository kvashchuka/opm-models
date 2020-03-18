#ifndef DUNE_OPTIM_COMMON_FIELDMATRIX_HH
#define DUNE_OPTIM_COMMON_FIELDMATRIX_HH

#include <dune/common/fmatrix.hh>

namespace Dune
{

  template< class X, class Field, int rows, int cols, class Y >
  Field multiply ( const X &x, const FieldMatrix< Field, rows, cols > &A, const Y &y )
  {
    Field z( 0 );
    for( int i = 0; i < rows; ++i )
    {
      Field zi( 0 );
      for( int j = 0; j < cols; ++j )
        zi += A[ i ][ j ] * x[ j ];
      z += y[ i ] * zi;
    }
    return z;
  }

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_COMMON_FIELDMATRIX_HH
