#ifndef DUNE_OPTIM_QP_HH
#define DUNE_OPTIM_QP_HH

#include <dune/optim/activeindexmapper.hh>
#include <dune/optim/common/smallobject.hh>
#include <dune/optim/common/densesubvector.hh>
#include <dune/optim/std/subarray.hh>

namespace Dune
{

  namespace Optim
  {

    // ActiveSetStrategy
    // -----------------

    template< class QP, bool verbose = false >
    struct ActiveSetStrategy
    {
      typedef typename QP::Field Field;

      explicit ActiveSetStrategy ( const Field &epsilon, const QP &qp = QP() )
        : epsilon_( epsilon ), qp_( qp )
      {}

      /**
       * \brief solve inequality-constrained QP problem
       *
       * \note The constraints are such that \f$n * x \le c\f$.
       **/
      template< class Problem, class ConstraintArray, class DomainVector, class ActiveIndexMapper, class ActiveMultiplierVector >
      void operator() ( const Problem &problem, const ConstraintArray &constraints, DomainVector &x, ActiveIndexMapper &active, ActiveMultiplierVector &lambda ) const;

      template< class Problem, class ConstraintArray, class DomainVector, class MultiplierVector >
      void operator() ( const Problem &problem, const ConstraintArray &constraints, DomainVector &x, MultiplierVector &lambda ) const;

    private:
      Field epsilon_;
      QP qp_;
    };



    // Implementation of ActiveSetStrategy
    // -----------------------------------

    template< class QP, bool verbose >
    template< class Problem, class ConstraintArray, class DomainVector, class ActiveIndexMapper, class ActiveMultiplierVector >
    inline void
    ActiveSetStrategy< QP, verbose >::operator() ( const Problem &problem, const ConstraintArray &constraints, DomainVector &x, ActiveIndexMapper &active, ActiveMultiplierVector &lambda ) const
    {
      typedef typename ConstraintArray::value_type Constraint;
      typedef typename ActiveIndexMapper::InactiveIterator InactiveIterator;

      const auto activeConstraints = Std::make_subarray( constraints, std::ref( active ) );

      DomainVector dx( x );
      while( true )
      {
        qp_( problem, activeConstraints, dx, lambda );
        dx -= x;

        if( dx*dx < epsilon_*epsilon_ )
        {
          // find an active constraint with negative Lagrange multiplier
          int q = -1;
          Field lambda_min = -epsilon_;
          const unsigned int nActive = active.size();
          for( unsigned int i = 0; i < nActive; ++i )
          {
            const Field &lambda_i = lambda[ i ];
            if( lambda_i >= lambda_min )
              continue;
            lambda_min = lambda_i;
            q = i;
          }
          if( q == -1 )
            return;

          if( verbose )
            std::cout << "Releasing constraint " << active[ q ] << std::endl;
          lambda[ q ] = 0;
          active.release( q );
        }
        else
        {
          Field alpha = std::numeric_limits< Field >::infinity();
          int p = 0;
          const InactiveIterator end = active.endInactive();
          for( InactiveIterator it = active.beginInactive(); it != end; ++it )
          {
            const Constraint &constraint = constraints[ *it ];
            const Field ndx = constraint.normal() * dx;
            if( ndx < epsilon_ )
              continue;
            const Field beta = -constraint.evaluate( x ) / ndx;
            if( beta >= alpha )
              continue;

            alpha = beta;
            p = *it;
          }

          if( alpha <= Field( 1 ) )
          {
            x.axpy( alpha, dx );

            if( verbose )
              std::cout << "Adding Constraint: " << p << " (alpha = " << alpha << ")" << std::endl;
            active.activate( p );
          }
          else
            x += dx;
        }
      }
    }


    template< class QP, bool verbose >
    template< class Problem, class ConstraintArray, class DomainVector, class MultiplierVector >
    inline void
    ActiveSetStrategy< QP, verbose >::operator() ( const Problem &problem, const ConstraintArray &constraints, DomainVector &x, MultiplierVector &lambda ) const
    {
      typedef Optim::ActiveIndexMapper< SmallObjectAllocator< unsigned int > > ActiveIndexMapper;

      const unsigned int dimension = problem.dimension();
      const unsigned int nConstraints = constraints.size();
      ActiveIndexMapper active( dimension, nConstraints );

      lambda = 0;
      auto lambdaActive = denseSubVector( lambda, std::ref( active ) );
      (*this)( problem, constraints, x, active, lambdaActive );
    }

  } // namespace Optim

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_QP_HH
