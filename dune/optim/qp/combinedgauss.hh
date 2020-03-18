#ifndef DUNE_OPTIM_QP_COMBINEDGAUSS_HH
#define DUNE_OPTIM_QP_COMBINEDGAUSS_HH

#include <cassert>
#include <cstddef>

#include <dune/fem/storage/vector.hh>

#include <dune/optim/common/matrix.hh>
#include <dune/optim/solver/gauss.hh>

namespace Dune
{

  namespace Optim
  {

    // CombinedGaussQP
    // ---------------

    template< class F >
    struct CombinedGaussQP
    {
      typedef F Field;

      template< class Hessian, class Gradient >
      struct Problem;

    private:
      template< class ConstraintArray >
      struct ConstraintMatrix;

      template< class ConstraintArray >
      struct ConstraintVector;

    public:
      template< class Hessian, class Gradient, class ConstraintArray,
                class DomainVector, class MultiplierVector >
      void operator() ( const Problem< Hessian, Gradient > &problem,
                        const ConstraintArray &constraints,
                        DomainVector &x, MultiplierVector &lambda ) const;

      template< class Hessian, class Gradient >
      static Problem< Hessian, Gradient >
      problem ( const Hessian &hessian, const Gradient &negGradient )
      {
        return Problem< Hessian, Gradient >( hessian, negGradient );
      }

    private:
      PivotizedGaussSolver< Field > solver_;
    };



    // CombinedGaussQP::Problem
    // ------------------------

    template< class F >
    template< class Hessian, class Gradient >
    struct CombinedGaussQP< F >::Problem
    {
      Problem ( const Hessian &h, const Gradient &nGrad )
        : hessian( h ), negGradient( nGrad )
      {
        assert( hessian.rows == negGradient.size() );
        assert( hessian.rows == hessian.cols );
      }

      std::size_t dimension () const { return negGradient.size(); }

      const Hessian &hessian;
      const Gradient &negGradient;
    };



    // CombinedGaussQP::ConstraintMatrix
    // ---------------------------------

    template< class F >
    template< class ConstraintArray >
    struct CombinedGaussQP< F >::ConstraintMatrix
    {
      typedef typename CombinedGaussQP< F >::Field Field;

      ConstraintMatrix ( std::size_t columns, const ConstraintArray &array )
        : columns_( columns ),
          array_( array )
      {}

      const Field &operator() ( std::size_t row, std::size_t column ) const
      {
        assert( column < columns_ );
        return array_[ row ].normal()[ column ];
      }

      std::size_t columns () const { return columns_; }
      std::size_t rows () const { return array_.size(); }

    private:
      std::size_t columns_;
      const ConstraintArray &array_;
    };



    // CombinedGaussQP::ConstraintVector
    // ---------------------------------

    template< class F >
    template< class ConstraintArray >
    struct CombinedGaussQP< F >::ConstraintVector
      : public VectorDefaultReadOnly< F, ConstraintVector< ConstraintArray > >
    {
      typedef typename CombinedGaussQP< F >::Field Field;

      explicit ConstraintVector ( const ConstraintArray &array ) : array_( array ) {}

      const Field &operator[] ( int index ) const { return array_[ index ].rhs(); }

      std::size_t size () const { return array_.size(); }

    private:
      const ConstraintArray &array_;
    };



    // Implementation of CombinedGaussQP
    // ---------------------------------

    template< class F >
    template< class Hessian, class Gradient, class ConstraintArray,
              class DomainVector, class MultiplierVector >
    inline void
    CombinedGaussQP< F >::operator() ( const Problem< Hessian, Gradient > &problem,
                                       const ConstraintArray &constraints,
                                       DomainVector &x, MultiplierVector &lambda ) const
    {
      typedef ConstraintMatrix< ConstraintArray > CMatrix;
      typedef TransposedMatrix< CMatrix > CTMatrix;
      typedef Dune::NullMatrix< Field > NullMatrix;
      typedef Dune::BlockMatrix< Hessian, CTMatrix, CMatrix, NullMatrix > BlockMatrix;

      typedef ConstraintVector< ConstraintArray > CVector;
      typedef Dune::Fem::CombinedVector< Gradient, CVector > RHS;

      const std::size_t dimension = problem.dimension();
      assert( x.size() == dimension );

      const std::size_t nConstraints = constraints.size();
      assert( lambda.size() == nConstraints );

      CMatrix C( dimension, constraints );
      CTMatrix CT( C );
      NullMatrix N( nConstraints, nConstraints );
      BlockMatrix M( problem.hessian, CT, C, N );

      CVector c( constraints );
      const RHS b( const_cast< Gradient & >( problem.negGradient ), c );

      Dune::Fem::CombinedVector< DomainVector, MultiplierVector > y( x, lambda );
      solver_.solve( M, b, y );
    }

  } // namespace Optim

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_QP_COMBINEDGAUSS_HH
