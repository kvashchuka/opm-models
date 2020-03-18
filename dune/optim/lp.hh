#ifndef DUNE_OPTIM_LP_HH
#define DUNE_OPTIM_LP_HH

#include <dune/common/densematrix.hh>

#include <dune/optim/activeindexmapper.hh>
#include <dune/optim/common/smallobject.hh>
#include <dune/optim/std/subarray.hh>

namespace Dune
{

  namespace Optim
  {

    namespace __LinearProgramming
    {

      // ContraintMatrix
      // ---------------

      template< class Constraints, class ActiveIndexMapper >
      struct ConstraintMatrix
      {
        ConstraintMatrix ( const Constraints &constraints, const ActiveIndexMapper &active )
          : constraints_( constraints ), active_( active )
        {}

        const Constraints &constraints () const noexcept { return constraints_; }
        const ActiveIndexMapper &active () const { return active_; }

        const std::size_t N () const { return active().size(); }
        const std::size_t M () const { return active().size(); }

      private:
        const Constraints &constraints_;
        const ActiveIndexMapper &active_;
      };

    } // namespace __LinearProgramming



    // LinearProgramming
    // -----------------

    template< class LinearSolver, bool verbose = false >
    struct LinearProgramming
    {
      typedef typename LinearSolver::Field Field;
      typedef Opm::MathToolbox<Field> Toolbox;


      explicit LinearProgramming ( const Field &epsilon, const LinearSolver &linearSolver = LinearSolver() ) : epsilon_( epsilon ), linearSolver_( linearSolver ) {}

      /**
       * \brief solve inequality-constrained LP problem
       *
       * \note The constraints are such that \f$n * x \le c\f$.
       **/
      template< class DomainVector, class ConstraintArray, class EvalVector, class ActiveIndexMapper >
      void operator() ( const DomainVector &descent, const ConstraintArray &constraints, EvalVector &x, ActiveIndexMapper &active ) const;

    private:
      Field epsilon_;
      LinearSolver linearSolver_;
    };



    // Implementation of LinearProgramming
    // -----------------------------------

    template< class LinearSolver, bool verbose >
    template< class DomainVector, class ConstraintArray, class Jacobian, class ActiveIndexMapper >
    inline void
    LinearProgramming< LinearSolver, verbose >::operator() ( const DomainVector &descent, const ConstraintArray &constraints, Jacobian &x, ActiveIndexMapper &active ) const
    {
      typedef typename ConstraintArray::value_type Constraint;
      typedef typename ActiveIndexMapper::InactiveIterator InactiveIterator;

      __LinearProgramming::ConstraintMatrix< ConstraintArray, ActiveIndexMapper > constraintMatrix( constraints, active );

      assert( active.size() == x.size() );
      auto linearInverse = linearSolver_( constraintMatrix );

      Jacobian dx( x );
      Jacobian lambda( x );
      int counter = 0;
      while( counter < 10 ) //max 10
      {
        // compute Langrange multipliers
        linearInverse.mtv( descent, lambda );

        // find a constraint with negative Lagrange multiplier
        int q = -1;
        Field lambda_min = -epsilon_;
        for( unsigned int i = 0; i < active.size(); ++i )
        {
          const Field &lambda_i = lambda[ i ];
          if( Toolbox::value(lambda_i) >= Toolbox::value(lambda_min) )
            continue;
          lambda_min = lambda_i;
          q = i;
        }
        if( q == -1 )
          return;
        if( verbose )
          std::cout << "Releasing constraint " << active[ q ] << std::endl;

        // compute search direction
        lambda = Field( 0 );
        lambda[ q ] = Field( -1 );
        linearInverse.mv( lambda, dx );

        // find the nearest constraint in the descent direction
        Field alpha = std::numeric_limits< double >::infinity();
        int p = 0;
        const InactiveIterator end = active.endInactive();
        for( InactiveIterator it = active.beginInactive(); it != end; ++it )
        {
          const Constraint &constraint = constraints[ *it ];
          const Field ndx = constraint.normal() * dx;
          if( Toolbox::value(ndx) < Toolbox::value(epsilon_) )
            continue;
          const Field beta = -constraint.evaluate( x ) / ndx;
          if( Toolbox::value(beta) >= Toolbox::value(alpha) )
            continue;

          alpha = beta;
          p = *it;
        }
        if( verbose )
          std::cout << "Adding Constraint: " << p << " (alpha = " << alpha << ")" << std::endl;

        //if( true )
        //  std::cout << "Adding Constraint: " << p << " (alpha = " << alpha << ")";

//        std::cout << " dx ";
//        for (int i = 0; i < 2; ++i) {
//            std::cout << dx[i] << " ";
//            for (int j = 0; j < 2; ++j) {
//                std::cout << dx[i].derivative(j) << " ";
//            }
//        }
//        std::cout << std::endl;

        // walk towards this constraint
        x.axpy( alpha, dx );

        // update active indices and linear solver
        active.update( q, p );
        linearInverse.updateRow( q, constraintMatrix );
        counter ++;
      }
    }

  } // namespace Optim



  // DenseMatrixAssigner for Optim::__LinearProgramming::ConstraintMatrix
  // ---------------------------------------------------------------------

  template< class DenseMatrix, class Constraints, class ActiveIndexMapper >
  struct DenseMatrixAssigner< DenseMatrix, Optim::__LinearProgramming::ConstraintMatrix< Constraints, ActiveIndexMapper > >
  {
    static void apply ( DenseMatrix &denseMatrix, const Optim::__LinearProgramming::ConstraintMatrix< Constraints, ActiveIndexMapper > &constraintMatrix )
    {
      const Constraints &constraints = constraintMatrix.constraints();
      const ActiveIndexMapper &active = constraintMatrix.active();
      typedef typename FieldTraits< DenseMatrix >::field_type Field;
      typedef Opm::MathToolbox<Field> Toolbox;
      const std::size_t size = active.size();
      assert( (denseMatrix.N() == size) && (denseMatrix.M() == size) );
      for( std::size_t i = 0; i < size; ++i ) {
          for (std::size_t j = 0; j < constraints[ active[ i ] ].normal().size(); ++j)
              denseMatrix[ i ][ j ] = Toolbox::createConstant(constraints[ active[ i ] ].normal()[j]);


      }
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_LP_HH
