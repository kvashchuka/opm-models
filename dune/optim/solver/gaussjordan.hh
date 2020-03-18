#ifndef DUNE_OPTIM_SOLVER_DENSEMATRIX_HH
#define DUNE_OPTIM_SOLVER_DENSEMATRIX_HH

#include <cassert>

#include <type_traits>

#include <dune/common/densematrix.hh>
#include <dune/common/math.hh>

namespace Dune
{

  namespace Optim
  {

    // invertGaussJordan
    // -----------------

    template< class M >
    bool invertGaussJordan ( DenseMatrix< M > &matrix, const typename FieldTraits< DenseMatrix< M > >::real_type &tolerance )
    {
      typedef std::remove_const_t< typename FieldTraits< DenseMatrix< M > >::field_type > Field;
      typedef std::remove_const_t< typename FieldTraits< DenseMatrix< M > >::real_type > Real;
      typedef typename DenseMatrix< M >::size_type Size;
      typedef typename DenseMatrix< M >::row_type Row;

      using std::abs;
      using std::swap;

      assert( matrix.rows() == matrix.cols() );
      const Size size = matrix.rows();
      if( size == 0 )
        return true;

      std::vector< Size > p( size );
      for( Size j = 0; j < size; ++j )
        p[ j ] = j;
      for( Size j = 0; j < size; ++j )
      {
        // pivot search
        Size r = j;
        Real max = abs( matrix[ j ][ j ] );
        for( Size i = j+1; i < size; ++i )
        {
          Real absij = abs( matrix[ i ][ j ] );
          if( absij > max )
          {
            max = absij;
            r = i;
          }
        }
        if( max < tolerance )
          return false;

        // row swap
        if( r > j )
        {
          for( Size k = 0; k < size; ++k )
            swap( matrix[ j ][ k ], matrix[ r ][ k ] );
          swap( p[ j ], p[ r ] );
        }

        // transformation
        Field hr = Field( 1 ) / matrix[ j ][ j ];
        for( Size i = 0; i < size; ++i )
          matrix[ i ][ j ] *= hr;
        matrix[ j ][ j ] = hr;
        for( Size k = 0; k < size; ++k )
        {
          if( k == j )
            continue;
          for( Size i = 0; i < size; ++i )
          {
            if( i == j )
              continue;
            matrix[ i ][ k ] -= matrix[ i ][ j ]*matrix[ j ][ k ];
          }
          matrix[ j ][ k ] *= -hr;
        }
      }

      // column exchange
      Row hv( matrix[ 0 ] );
      for( Size i = 0; i < size; ++i )
      {
        for( Size k = 0; k < size; ++k )
          hv[ p[ k ] ] = matrix[ i ][ k ];
        for( Size k = 0; k < size; ++k )
          matrix[ i ][ k ] = hv[ k ];
      }
      return true;
    }



    // GaussJordanSolver
    // -----------------

    template< class Matrix >
    struct GaussJordanSolver
    {
      typedef typename FieldTraits< Matrix >::field_type Field;
      typedef typename FieldTraits< Matrix >::real_type Real;

      struct Inverse
      {
        template< class A >
        Inverse ( const A &a, Real epsilon )
          : matrix_( a ), epsilon_( epsilon )
        {
          invertGaussJordan( matrix_, epsilon_ );
        }

        template< class X, class Y >
        void mv ( const X &x, Y &y )
        {
          matrix_.mv( x, y );
        }

        template< class X, class Y >
        void mtv ( const X &x, Y &y )
        {
          matrix_.mtv( x, y );
        }

        template< class A >
        void updateRow ( int k, const A &a )
        {
          matrix_ = a;
          invertGaussJordan( matrix_, epsilon_ );
        }

      private:
        Matrix matrix_;
        Real epsilon_;
      };

      explicit GaussJordanSolver ( Real epsilon = 1e-8 ) : epsilon_( epsilon ) {}

      template< class A >
      Inverse operator() ( const A &a ) const
      {
        return Inverse( a, epsilon_ );
      }

    private:
      Real epsilon_;
    };

  } // namespace Optim

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_SOLVER_DENSEMATRIX_HH
