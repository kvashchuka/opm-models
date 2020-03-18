#ifndef DUNE_OPTIM_QP_REDUCEDGAUSS_HH
#define DUNE_OPTIM_QP_REDUCEDGAUSS_HH

#include <cassert>
#include <cstddef>

#include <dune/common/densematrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/optim/solver/gauss.hh>

namespace Dune
{

  namespace Optim
  {

    namespace __ReducedGaussQP
    {

      // ContraintMatrix
      // ---------------

      template< class Hessian, class Constraints >
      struct ConstraintMatrix
      {
        ConstraintMatrix ( const Hessian &hessianInverse, const Constraints &constraints )
          : hessianInverse_( hessianInverse ), constraints_( constraints )
        {}

        const Hessian &hessianInverse () const noexcept { return hessianInverse_; }
        const Constraints &constraints () const noexcept { return constraints_; }

        const std::size_t N () const { return constraints().size(); }
        const std::size_t M () const { return constraints().size(); }

      private:
        const Hessian &hessianInverse_;
        const Constraints &constraints_;
      };

    } // namespace __ReducedGaussQP



    // ReducedGaussQP
    // --------------

    template< class F >
    struct ReducedGaussQP
    {
      typedef F Field;

      template< class Hessian, class Gradient >
      struct Problem;

    private:
      template< class Hessian, class Gradient, class ConstraintArray >
      struct ConstraintVector;

    public:
      template< class Hessian, class Gradient, class ConstraintArray, class DomainVector, class MultiplierVector >
      void operator() ( const Problem< Hessian, Gradient > &problem, const ConstraintArray &constraints, DomainVector &x, MultiplierVector &lambda ) const;

      template< class Hessian, class Gradient >
      static Problem< Hessian, Gradient >
      problem ( const Hessian &hessianInverse, const Gradient &negGradient )
      {
        return Problem< Hessian, Gradient >( hessianInverse, negGradient );
      }

    private:
      PivotizedGaussSolver< Field > solver_;
    };



    // ReducedGaussQP::Problem
    // -----------------------

    template< class F >
    template< class Hessian, class Gradient >
    struct ReducedGaussQP< F >::Problem
    {
      Problem ( const Hessian &hInv, const Gradient &nGrad )
        : hessianInverse( hInv ), negGradient( nGrad ), inverseGradient( negGradient )
      {
        assert( hessianInverse.rows == negGradient.size() );
        assert( hessianInverse.rows == hessianInverse.cols );
        hessianInverse.mtv( negGradient, inverseGradient );
      }

      std::size_t dimension () const { return negGradient.size(); }

      const Hessian &hessianInverse;
      const Gradient &negGradient;
      Gradient inverseGradient;
    };



    // ReducedGaussQP::ConstraintVector
    // --------------------------------

    template< class F >
    template< class Hessian, class Gradient, class ConstraintArray >
    struct ReducedGaussQP< F >::ConstraintVector
    {
      typedef typename ReducedGaussQP< F >::Field Field;
      typedef typename ReducedGaussQP< F >::template Problem< Hessian, Gradient > Problem;

      ConstraintVector ( const Problem &problem, const ConstraintArray &array )
        : problem_( problem ), array_( array )
      {}

      Field operator[] ( int index ) const
      {
        return problem_.inverseGradient *  array_[ index ].normal() - array_[ index ].rhs();
      }

      std::size_t size () const { return array_.size(); }

    private:
      const Problem &problem_;
      const ConstraintArray &array_;
    };



    // Implementation of ReducedGaussQP
    // --------------------------------

    template< class F >
    template< class Hessian, class Gradient, class ConstraintArray, class DomainVector, class MultiplierVector >
    inline void
    ReducedGaussQP< F >::operator() ( const Problem< Hessian, Gradient > &problem, const ConstraintArray &constraints, DomainVector &x, MultiplierVector &lambda ) const
    {
      typedef __ReducedGaussQP::ConstraintMatrix< Hessian, ConstraintArray > ConstraintMatrix;
      typedef typename ReducedGaussQP< F >::template ConstraintVector< Hessian, Gradient, ConstraintArray > ConstraintVector;

      assert( x.size() == problem.dimension() );

      const std::size_t nConstraints = constraints.size();
      assert( lambda.size() == nConstraints );

      if( nConstraints > 0 )
      {
        ConstraintMatrix A( problem.hessianInverse, constraints );
        ConstraintVector b( problem, constraints );
        solver_.solve( A, b, lambda );
      }

      FieldVector< Field, Hessian::rows > y( problem.negGradient );
      for( std::size_t i = 0; i < nConstraints; ++i )
        y.axpy( -lambda[ i ], constraints[ i ].normal() );
      problem.hessianInverse.mv( y, x );
    }

  } // namespace Optim



  // DenseMatrixAssigner for Optim::__ReducedGaussQP::ConstraintMatrix
  // -----------------------------------------------------------------

  template< class DenseMatrix, class Hessian, class Constraints >
  class DenseMatrixAssigner< DenseMatrix, Optim::__ReducedGaussQP::ConstraintMatrix< Hessian, Constraints > >
  {
    typedef typename Hessian::field_type Field;

  public:
    static void apply ( DenseMatrix &denseMatrix, const Optim::__ReducedGaussQP::ConstraintMatrix< Hessian, Constraints > &constraintMatrix )
    {
      const Hessian &hessianInverse = constraintMatrix.hessianInverse();
      const Constraints &constraints = constraintMatrix.constraints();

      const std::size_t size = constraints.size();
      assert( (denseMatrix.N() == size) && (denseMatrix.M() == size) );
      FieldVector< Field, Hessian::rows > hin;
      for( std::size_t i = 0; i < size; ++i )
      {
        hessianInverse.mv( constraints[ i ].normal(), hin );
        for( std::size_t j = 0; j < size; ++j )
          denseMatrix[ i ][ j ] = hin * constraints[ j ].normal();
      }
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_QP_REDUCEDGAUSS_HH
