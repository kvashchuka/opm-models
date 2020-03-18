#ifndef DUNE_OPTIM_SOLVER_GAUSS_HH
#define DUNE_OPTIM_SOLVER_GAUSS_HH

#include <cstddef>

#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/dynvector.hh>

#include <dune/optim/solver/fullydynmatrix.hh>

namespace Dune
{

  namespace Optim
  {

    // PivotizedGaussSolver
    // --------------------

    template< class F, class A = std::allocator< F > >
    struct PivotizedGaussSolver
    {
      typedef F Field;
      typedef A Allocator;

      explicit PivotizedGaussSolver ( Allocator allocator = Allocator() )
        : epsilon_( 1e-10 ), allocator_( std::move( allocator ) )
      {}

      explicit PivotizedGaussSolver ( Field epsilon, Allocator allocator = Allocator() )
        : epsilon_( std::move( epsilon ) ), allocator_( std::move( allocator() ) )
      {}

      template< class Matrix, class RangeVector, class DomainVector >
      DomainVector &solve ( const Matrix &matrix, const RangeVector &b, DomainVector &x ) const
      {
        FullyDynamicMatrix< Field, Allocator > matrixCopy( matrix.N(), matrix.M(), allocator() );
        matrixCopy = matrix;

        DynamicVector< Field, Allocator > bCopy( b.size(), allocator() );
        for( std::size_t i = 0; i < b.size(); ++i )
          bCopy[ i ] = b[ i ];

        doSolve( matrixCopy, bCopy, x );
        return x;
      }

      const Field &epsilon () const { return epsilon_; }

      const Allocator &allocator () const { return allocator_; }

    private:
      template< class DomainVector >
      void doSolve ( FullyDynamicMatrix< Field, Allocator > &matrix, DynamicVector< Field, Allocator > &b, DomainVector &x ) const
      {
        using std::abs;

        const std::size_t size = b.size();
        assert( (matrix.rows() == size) && (matrix.cols() == size) );
        assert( x.size() == size );

        std::vector< std::size_t, typename Allocator::template rebind< std::size_t >::other > col( allocator() );
        col.reserve( size );
        for( std::size_t i = 0; i < size; ++i )
          col.emplace_back( i );

        for( std::size_t i = 0; i < size; ++i )
        {
          // Find pivot element
          Field pivot( 0 );
          int prow = i, pcol = i;
          for( std::size_t r = i; r < size; ++r )
          {
            for( std::size_t c = i; c < size; ++c )
            {
              Field a = abs( matrix[ r ][ col[ c ] ] );
              if( a > pivot )
              {
                pivot = a;
                prow = r;
                pcol = c;
              }
            }
          }
          if( pivot < epsilon() )
          {
            for( std::size_t r = i; r < size; ++r )
            {
              if( std::abs( b[ r ] ) >= epsilon() )
                DUNE_THROW( MathError, "PivotizedGaussSolver: No solution to system (b = " << b[ r ] << ")." );
              b[ r ] = 0;
              matrix[ r ][ col[ r ] ] = 1;
            }
            break;
          }
          for( std::size_t j = 0; j < size; ++j )
            std::swap( matrix[ i ][ j ], matrix[ prow ][ j ] );
          std::swap( b[ i ], b[ prow ] );
          std::swap( col[ i ], col[ pcol ] );

          // clean the i-th column (though we don't actually write 0 there)
          for( std::size_t j = i+1; j < size; ++j )
          {
            const Field factor = matrix[ j ][ col[ i ] ] / matrix[ i ][ col[ i ] ];
            matrix[ j ][ col[ i ] ] = 0;
            for( std::size_t k = i+1; k < size; ++k )
              matrix[ j ][ col[ k ] ] -= factor * matrix[ i ][ col[ k ] ];
            b[ j ] -= factor * b[ i ];
          }
        }

        for( std::size_t i = size; i > 0; --i )
        {
          Field sum = 0;
          for( std::size_t j = i; j < size; ++j )
            sum += matrix[ i-1 ][ col[ j ] ] * x[ col[ j ] ];
          x[ col[ i-1 ] ] = (b[ i-1 ] - sum) / matrix[ i-1 ][ col[ i-1 ] ];
        }
      }

      Field epsilon_;
      Allocator allocator_;
    };



    // RowPivotizedGaussSolver
    // -----------------------

    template< class F, class A = std::allocator< F > >
    struct RowPivotizedGaussSolver
    {
      typedef F Field;
      typedef A Allocator;

      explicit RowPivotizedGaussSolver ( Allocator allocator = Allocator() )
        : epsilon_( 1e-10 ), allocator_( std::move( allocator ) )
      {}

      explicit RowPivotizedGaussSolver ( Field epsilon, Allocator allocator = Allocator() )
        : epsilon_( std::move( epsilon ) ), allocator_( std::move( allocator() ) )
      {}

      template< class Matrix, class RangeVector, class DomainVector >
      DomainVector &solve ( const Matrix &matrix, const RangeVector &b, DomainVector &x ) const
      {
        FullyDynamicMatrix< Field, Allocator > matrixCopy( matrix.N(), matrix.M(), allocator() );
        matrixCopy = matrix;

        DynamicVector< Field, Allocator > bCopy( b.size(), allocator() );
        for( std::size_t i = 0; i < b.size(); ++i )
          bCopy[ i ] = b[ i ];

        doSolve( matrixCopy, bCopy, x );
        return x;
      }

      const Field &epsilon () const { return epsilon_; }

      const Allocator &allocator () const { return allocator_; }

    private:
      template< class DomainVector >
      void doSolve ( FullyDynamicMatrix< Field, Allocator > &matrix, DynamicVector< Field, Allocator > &b, DomainVector &x ) const
      {
        using std::abs;

        const std::size_t size = b.size();
        assert( (matrix.rows() == size) && (matrix.cols() == size) );
        assert( x.size() == size );

        for( std::size_t i = 0; i < size; ++i )
        {
          // Find pivot element
          Field pivot( 0 );
          int prow = i;
          for( std::size_t r = i; r < size; ++r )
          {
            Field a = abs( matrix[ r ][ i ] );
            if( a > pivot )
            {
              pivot = a;
              prow = r;
            }
          }
          if( pivot < epsilon() )
          {
            for( std::size_t r = i; r < size; ++r )
            {
              if( std::abs( b[ r ] ) >= epsilon() )
                DUNE_THROW( MathError, "PivotizedGaussSolver: No solution to system (b = " << b[ r ] << ")." );
              b[ r ] = 0;
              matrix[ r ][ r ] = 1;
            }
            break;
          }
          for( int j = 0; j < size; ++j )
            std::swap( matrix[ i ][ j ], matrix[ prow ][ j ] );
          std::swap( b[ i ], b[ prow ] );

          // clean the i-th column (though we don't actually write 0 there)
          for( std::size_t j = i+1; j < size; ++j )
          {
            const Field factor = matrix[ j ][ i ] / matrix[ i ][ i ];
            matrix[ j ][ i ] = 0;
            for( std::size_t k = i+1; k < size; ++k )
              matrix[ j ][ k ] -= factor * matrix[ i ][ k ];
            b[ j ] -= factor * b[ i ];
          }
        }

        for( std::size_t i = size; i > 0; --i )
        {
          Field sum = 0;
          for( std::size_t j = i; j < size; ++j )
            sum += matrix[ i-1 ][ j ] * x[ j ];
          x[ i ] = (b[ i-1 ] - sum) / matrix[ i-1 ][ i-1 ];
        }
      }

      Field epsilon_;
      Allocator allocator_;
    };



    // GaussSolver
    // -----------

    template< class F, class A = std::allocator< F > >
    struct GaussSolver
    {
      typedef F Field;
      typedef A Allocator;

      explicit GaussSolver ( Allocator allocator = Allocator() )
        : epsilon_(  1e-10 ), allocator_( std::move( allocator ) )
      {}

      explicit GaussSolver ( Field epsilon, Allocator allocator = Allocator() )
        : epsilon_( std::move( epsilon ) ), allocator_( std::move( allocator() ) )
      {}

      template< class Matrix, class RangeVector, class DomainVector >
      DomainVector &solve ( const Matrix &matrix, const RangeVector &b, DomainVector &x ) const
      {
        FullyDynamicMatrix< Field, Allocator > matrixCopy( matrix.N(), matrix.M(), allocator() );
        matrixCopy = matrix;

        DynamicVector< Field, Allocator > bCopy( b.size(), allocator() );
        for( std::size_t i = 0; i < b.size(); ++i )
          bCopy[ i ] = b[ i ];

        doSolve( matrixCopy, bCopy, x );
        return x;
      }

      const Field &epsilon () const { return epsilon_; }

      const Allocator &allocator () const { return allocator_; }

    private:
      template< class DomainVector >
      void doSolve ( FullyDynamicMatrix< Field, Allocator > &matrix, DynamicVector< Field, Allocator > &b, DomainVector &x ) const
      {
        const std::size_t size = b.size();
        assert( (matrix.rows() == size) && (matrix.cols() == size) );
        assert( x.size() == size );

        for( std::size_t i = 0; i < size; ++i )
        {
          const Field pivot = matrix[ i ][ i ];
          if( pivot < epsilon() )
          {
            for( std::size_t r = i; r < size; ++r )
            {
              if( std::abs( b[ r ] ) >= epsilon() )
                DUNE_THROW( MathError, "GaussSolver: No solution to system (b = " << b[ r ] << ")." );
              b[ r ] = 0;
              matrix[ r ][ r ] = 1;
            }
            break;
          }

          // clean the i-th column (though we don't actually write 0 there)
          Field (&ri)[ size ] = matrix[ i ];
          for( std::size_t j = i+1; j < size; ++j )
          {
            Field (&rj)[ size ] = matrix[ j ];
            const Field factor = rj[ i ] / pivot;
            rj[ i ] = 0;
            for( std::size_t k = i+1; k < size; ++k )
              rj[ k ] -= factor * ri[ k ];
            b[ j ] -= factor * b[ i ];
          }
        }

        for( std::size_t i = size; i > 0; --i )
        {
          Field sum = 0;
          for( std::size_t j = i; j < size; ++j )
            sum += matrix[ i-1 ][ j ] * x[ j ];
          x[ i-1 ] = (b[ i-1 ] - sum) / matrix[ i-1 ][ i-1 ];
        }
      }

      Field epsilon_;
      Allocator allocator_;
    };

  } // namespace Optim

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_SOLVER_GAUSS_HH
