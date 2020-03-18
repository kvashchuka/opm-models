#ifndef DUNE_OPTIM_QP_GRADIENT_HH
#define DUNE_OPTIM_QP_GRADIENT_HH

#include <functional>
#include <iostream>
#include <limits>

#include <dune/optim/activeindexmapper.hh>
#include <dune/optim/common/smallobject.hh>
#include <dune/optim/common/densesubvector.hh>

namespace Dune
{

  namespace Optim
  {

    // ActiveSetGradientQP
    // -------------------

    template< class F, bool verbose = false >
    struct ActiveSetGradientQP
    {
      typedef F Field;


      template< class Hessian, class Gradient >
      struct Problem;

      explicit ActiveSetGradientQP ( const Field &epsilon )
        : epsilon_( epsilon )
      {}

      template< class Hessian, class Gradient, class ConstraintArray, class DomainVector, class ActiveIndexMapper, class ActiveMultiplierVector >
      void operator() ( const Problem< Hessian, Gradient > &problem,
                        const ConstraintArray &constraints,
                        DomainVector &x,
                        ActiveIndexMapper &active,
                        ActiveMultiplierVector &lambda ) const;

      template< class Hessian, class Gradient, class ConstraintArray, class DomainVector, class MultiplierVector >
      void operator() ( const Problem< Hessian, Gradient > &problem,
                        const ConstraintArray &constraints,
                        DomainVector &x,
                        MultiplierVector &lambda ) const;

      template< class Hessian, class Gradient >
      static Problem< Hessian, Gradient >
      problem ( const Hessian &hessianInverse, const Gradient &negGradient )
      {
        return Problem< Hessian, Gradient >( hessianInverse, negGradient );
      }

    private:
      Field epsilon_;
    };



    // ActiveSetGradientQP::Problem
    // ----------------------------

    template< class F, bool verbose >
    template< class Hessian, class Gradient >
    struct ActiveSetGradientQP< F, verbose >::Problem
    {
      Problem( const Hessian &h, const Gradient &nGrad )
        : hessian( h ), negGradient( nGrad )
      {
        assert( hessian.rows == negGradient.size() );
        assert( hessian.rows == hessian.cols );
      }

      unsigned int dimension () const { return negGradient.size(); }

      const Hessian &hessian;
      const Gradient &negGradient;
    };



    // Implementation of ActiveSetStrategy
    // -----------------------------------

    template< class F, bool verbose >
    template< class Hessian, class Gradient, class ConstraintArray, class DomainVector, class ActiveIndexMapper, class ActiveMultiplierVector >
    inline void
    ActiveSetGradientQP< F, verbose >::operator() ( const Problem< Hessian, Gradient > &problem, const ConstraintArray &constraints, DomainVector &x, ActiveIndexMapper &active, ActiveMultiplierVector &lambda ) const
    {
      typedef typename ConstraintArray::value_type Constraint;
      typedef typename ActiveIndexMapper::InactiveIterator InactiveIterator;

      DomainVector g( x );
      problem.hessian.mv( x, g );
      g -= problem.negGradient;

      DomainVector dx( g );
#if 0
      for( unsigned int i = 0; i < active.size(); ++i )
        dx.addScaled( lambda[ i ], constraints[ active[ i ] ].normal() );
      dx *= Field( -1 );
#endif

      DomainVector Adx( x );
      while( true )
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

        // can we release a constraint?
        if( q >= 0 )
        {
          if( verbose )
            std::cout << "Releasing constraint " << active[ q ] << std::endl;
          //dx.addScaled( lambda[ q ], constraints[ active[ q ] ].normal() );
          lambda[ q ] = 0;
          active.release( q );
        }

        dx = g;
        for( unsigned int i = 0; i < active.size(); ++i )
        {
          const Constraint &constraint = constraints[ active[ i ] ];
          const Field s = (g*constraint.normal()) / (constraint.normal()*constraint.normal());
          dx.axpy( -s, constraints[ active[ i ] ].normal() );
        }
        dx *= Field( -1 );

        const Field dxdx = dx*dx;
        if( dxdx < epsilon_*epsilon_ )
          return;

        // find minimal step that will hit a constraint
        Field alpha = std::numeric_limits< Field >::infinity();
        Field pndx = 0;
        int p = -1;
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
          pndx = ndx;
          p = *it;
        }

        // perform the update step
        problem.hessian.mv( dx, Adx );
        const Field gamma = dx * Adx;
        if( gamma > epsilon_*dxdx )
        {
          const Field beta = -(g*dx) / gamma;
          if( beta < alpha )
          {
            alpha = beta;
            p = -1;
          }
        }
        x.axpy( alpha, dx );
        g.axpy( alpha, Adx );

        // do we need to activate a constraint?
        if( p >= 0 )
        {
          if( verbose )
            std::cout << "Adding Constraint: " << p << "(alpha = " << alpha << ")" << std::endl;
          const unsigned int j = active.activate( p );
          lambda[ j ] = (alpha*gamma) / pndx;
          //dx.addScaled( -lambda[ j ], constraints[ p ].normal() );
        }
      }
    }


    template< class F, bool verbose >
    template< class Hessian, class Gradient, class ConstraintArray, class DomainVector, class MultiplierVector >
    inline void
    ActiveSetGradientQP< F, verbose >::operator() ( const Problem< Hessian, Gradient > &problem, const ConstraintArray &constraints, DomainVector &x, MultiplierVector &lambda ) const
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

#endif // #ifndef DUNE_OPTIM_QP_GRADIENT_HH
