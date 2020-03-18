#ifndef DUNE_FV_FUNCTION_PIECEWISELINEAR_HH
#define DUNE_FV_FUNCTION_PIECEWISELINEAR_HH

#include <type_traits>
#include <utility>

#include <dune/common/fvector.hh>

#include <dune/common/getreference.hh>

namespace Dune
{

  namespace FV
  {

    // PiecewiseLinearFunction
    // -----------------------

    template< class GridView, class Mapper, class AverageVector, class JacobianVector >
    class PiecewiseLinearFunction
    {
      typedef PiecewiseLinearFunction< GridView, Mapper, AverageVector, JacobianVector > This;

    public:
      typedef FieldVector< typename GridView::ctype, GridView::dimension > GlobalCoordinate;

      typedef typename GridView::template Codim< 0 >::Entity Entity;
      typedef typename GridView::template Codim< 0 >::Geometry Geometry;

      typedef typename getReferredType< const AverageVector & >::value_type Value;

    private:
      typedef typename getReferredType< const JacobianVector & >::value_type Jacobian;

      struct LocalFunction
      {
        LocalFunction ( const GlobalCoordinate &center, Value average, Jacobian jacobian )
          : center_( center ), average_( std::move( average ) ), jacobian_( std::move( jacobian ) )
        {}

        Value operator() ( const GlobalCoordinate &x ) const noexcept
        {
          Value value = average_;
          jacobian_.umv( x - center_, value );
          return std::move( value );
        }

      private:
        GlobalCoordinate center_;
        Value average_;
        Jacobian jacobian_;
      };

    public:
      PiecewiseLinearFunction ( const Mapper &mapper, AverageVector averages, JacobianVector jacobians ) noexcept
        : mapper_( mapper ), averages_( std::move( averages ) ), jacobians_( std::move( jacobians ) )
      {}

      friend LocalFunction localFunction ( const This &f, const Entity &entity )
      {
        return localFunction( f, entity, entity.geometry() );
      }

      friend LocalFunction localFunction ( const This &f, const Entity &entity, const Geometry &geometry )
      {
        const std::size_t index = f.mapper_.index( entity );
        return LocalFunction( geometry.center(), getReference( f.averages_ )[ index ], getReference( f.jacobians_ )[ index ] );
      }

    private:
      const Mapper &mapper_;
      AverageVector averages_;
      JacobianVector jacobians_;
    };



    // piecewiseLinearFunction
    // -----------------------

    template< class GridView, class Mapper, class AverageVector, class JacobianVector >
    inline static PiecewiseLinearFunction< GridView, Mapper, typename std::decay< AverageVector >::type, typename std::decay< JacobianVector >::type >
    piecewiseLinearFunction ( const GridView &gridView, const Mapper &mapper, AverageVector averages, JacobianVector jacobians ) noexcept
    {
      return PiecewiseLinearFunction< GridView, Mapper, typename std::decay< AverageVector >::type, typename std::decay< JacobianVector >::type >( mapper, std::move( averages ), std::move( jacobians ) );
    }

  } // namespace FV

} // namespace Dune

#endif // #ifndef DUNE_FV_FUNCTION_PIECEWISELINEAR_HH
