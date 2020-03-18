#ifndef DUNE_FV_LEASTSQUARESRECONSTRUCTION_HH
#define DUNE_FV_LEASTSQUARESRECONSTRUCTION_HH

#include <cassert>
#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/common/getreference.hh>

#include <dune/grid/utility/vectorcommdatahandle.hh>

#include <dune/fv/function/piecewiselinear.hh>

//#include <fem-fv/limiterutility.hh>
//#include <fem-fv/limitermodel.hh>

namespace Dune
{

  namespace FV
  {

    namespace Limiters
    {
      // base class for Limiters, mainly reading limitEps
      struct LimiterFunctionBase
      {
        const double limitEps_;
        LimiterFunctionBase()
          : limitEps_( getEpsilon() )
        {}

        //! get tolerance for shock detector
        static double getEpsilon()
        {
          // default value is 1e-8
          return Dune::Fem::Parameter::getValue("femdg.limiter.limiteps", double(1e-8) );
        }

        //! return epsilon for limiting
        double epsilon () const { return limitEps_; }

      protected:
        void printInfo(const std::string& name ) const
        {
          //if( Parameter::verbose() )
          {
            std::cout << "LimiterFunction: " << name << " with limitEps = " << limitEps_ << std::endl;
          }
        }
      };

      template <class Field>
      struct MinModLimiter : public LimiterFunctionBase
      {
        using LimiterFunctionBase :: limitEps_ ;
        typedef Field FieldType;

        MinModLimiter() : LimiterFunctionBase()
        {
          printInfo( "minmod" );
        }

        FieldType operator ()(const FieldType& g, const FieldType& d ) const
        {
          typedef Opm::MathToolbox<FieldType> Toolbox;
          const FieldType absG = Toolbox::abs( g );

          // avoid rounding errors
          // gradient of linear functions is nerly zero
          if ( absG < 1e-10 ) return 1;

          // if product smaller than zero set localFactor to zero
          // otherwise d/g until 1
          const FieldType localFactor =
                ( (d * g) < 0.0 ) ? 0.0 : Toolbox::min( Toolbox::createConstant(1.0) , ( d / g ) );

          return localFactor;
        }
      };

      template <class Field>
      struct SuperBeeLimiter : public LimiterFunctionBase
      {
        using LimiterFunctionBase :: limitEps_ ;
        typedef Field FieldType;

        SuperBeeLimiter() : LimiterFunctionBase()
        {
          printInfo( "superbee" );
        }

        FieldType operator ()(const FieldType& g, const FieldType& d ) const
        {
          // avoid rounding errors
          // gradient of linear functions is nerly zero
          if ( std::abs( g ) < limitEps_ ) return 1;

          const FieldType r = d / g ;

          // if product smaller than zero set localFactor to zero
          // otherwise d/g until 1
          const FieldType localFactor =
                ( (g * d) < 0.0 ) ? 0.0 :
                ( std::max( std::min( 2.0 * r , 1.0 ), std::min( r, 2.0 ) ) );

          return localFactor;
        }
      };

      template <class Field>
      struct VanLeerLimiter : public LimiterFunctionBase
      {
        using LimiterFunctionBase :: limitEps_ ;
        typedef Field FieldType;

        VanLeerLimiter() : LimiterFunctionBase()
        {
          printInfo( "vanLeer" );
        }

        FieldType operator ()(const FieldType& g, const FieldType& d ) const
        {
          const FieldType absG = std::abs( g );
          const FieldType absD = std::abs( d );

          // avoid rounding errors
          if ( (absG < limitEps_) && (absD < limitEps_) ) return 1;
          if ( absG < limitEps_ ) return 1;

          const FieldType r = d / g ;
          const FieldType absR = std::abs( r );

          return ( absR + r ) / ( 1.0 + absR );
        }
      };



    }


    // LeastSquaresReconstruction
    // --------------------------

    template< class GV, class SV, class BV >
    class LeastSquaresReconstruction
    {
      typedef LeastSquaresReconstruction< GV, SV, BV > This;

    public:
      typedef GV GridView;
      typedef SV StateVector;
      typedef BV BoundaryValue;

      typedef FieldVector< typename GridView::ctype, GridView::dimensionworld > GlobalCoordinate;

      typedef typename GridView::Intersection Intersection;

      typedef typename FieldTraits< StateVector >::field_type Field;
      typedef FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > Jacobian;

      template< class Mapper, class Vector >
      using Reconstruction = PiecewiseLinearFunction< GridView, Mapper, Vector, std::vector< Jacobian > >;

      LeastSquaresReconstruction ( const GridView &gridView, BoundaryValue boundaryValue )
        : gridView_( gridView ), boundaryValue_( std::move( boundaryValue ) )
      {
      }

      template< class Mapper, class Vector >
      Reconstruction< Mapper, Vector > operator() ( const Mapper &mapper, Vector u ) const
      {
        std::vector< Jacobian > du;
        (*this)( mapper, getReference( u ), du );
        return Reconstruction< Mapper, Vector >( mapper, std::move( u ), std::move( du ) );
      }

      template< class Mapper, class Vector >
      void operator () ( const Mapper &mapper, const Vector &u, std::vector< Jacobian > &du ) const
      {
        du.resize( mapper.size() );

        const auto end = gridView().template end<0, Dune::InteriorBorder_Partition>();
        for(auto it = gridView().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it)
        {
          const auto element = *it;

          const std::size_t elIndex = mapper.index( element );
          const GlobalCoordinate elCenter = element.geometry().center();

          Dune::FieldMatrix< Field, GlobalCoordinate::dimension, GlobalCoordinate::dimension > hessianInverse( 0 );
          Dune::FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > negGradient( 0 );
          const auto iend = gridView().iend( element );
          for( auto iit = gridView().ibegin( element ); iit != iend; ++iit )
          {
            const auto intersection = *iit ;
            const GlobalCoordinate iCenter = intersection.geometry().center();
            if( intersection.boundary() )
            {
              const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
              const StateVector du = boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] ) - u[ elIndex ];
              const GlobalCoordinate dx = iCenter - elCenter;
              for( int i = 0; i < GlobalCoordinate::dimension; ++i )
                hessianInverse[ i ].axpy( dx[ i ], dx );
              for( int i = 0; i < StateVector::dimension; ++i )
                negGradient[ i ].axpy( du[ i ], dx );
            }
            else if( intersection.neighbor() )
            {
              const auto neighbor = intersection.outside();
              const GlobalCoordinate dx = neighbor.geometry().center() - elCenter;
              const StateVector du = u[ mapper.index( neighbor ) ] - u[ elIndex ];
              for( int i = 0; i < GlobalCoordinate::dimension; ++i )
                hessianInverse[ i ].axpy( dx[ i ], dx );
              for( int i = 0; i < StateVector::dimension; ++i )
                negGradient[ i ].axpy( du[ i ], dx );
            }
          }

          hessianInverse.invert();
          for( int j = 0; j < StateVector::dimension; ++j )
            hessianInverse.mv( negGradient[ j ], du[ elIndex ][ j ] );
        }

        auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        gridView().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
      }

      const GridView &gridView () const { return gridView_; }

    private:
      GridView gridView_;
      BoundaryValue boundaryValue_;
    };



    // centralDifference
    // -----------------

    template< class SV, class GV, class BV >
    inline static LeastSquaresReconstruction< GV, SV, BV > leastSquaresReconstruction ( const GV &gridView, BV boundaryValue )
    {
      return LeastSquaresReconstruction< GV, SV, BV >( gridView, std::move( boundaryValue ) );
    }



    // LimitedLeastSquaresReconstruction
    // ---------------------------------

    template< class GV, class SV, class BV >
    class LimitedLeastSquaresReconstruction
    {
      typedef LimitedLeastSquaresReconstruction< GV, SV, BV > This;

    public:
      typedef GV GridView;
      typedef SV StateVector;
      typedef BV BoundaryValue;

      typedef FieldVector< typename GridView::ctype, GridView::dimensionworld > GlobalCoordinate;

      typedef typename GridView::Intersection Intersection;

      typedef typename FieldTraits< StateVector >::field_type Field;
      typedef FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > Jacobian;

      typedef Limiters::MinModLimiter< Field >  LimiterFunction;
      typedef Opm::MathToolbox<Field> Toolbox;

      typedef Dune::FieldMatrix< double, GlobalCoordinate::dimension, GlobalCoordinate::dimension > HessianInverseType;

      template< class Mapper, class Vector >
      using Reconstruction = PiecewiseLinearFunction< GridView, Mapper, Vector, std::vector< Jacobian > >;

      LimitedLeastSquaresReconstruction ( const GridView &gridView, BoundaryValue boundaryValue )
        : gridView_( gridView ), boundaryValue_( std::move( boundaryValue ) ), limiterFunction_(),
          hessianInverse_( gridView.indexSet().size( 0 ), HessianInverseType(0 ) )
      {
          const auto end = gridView.template end<0, Dune::InteriorBorder_Partition>();
          for(auto it = gridView.template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it)
          {
              const auto element = *it;
              const int elIdx = gridView.indexSet().index( element );
              const GlobalCoordinate elCenter = element.geometry().center();
              const auto iend = gridView.iend( element );
              auto& hessianInverse = hessianInverse_[ elIdx ];
              for( auto iit = gridView.ibegin( element ); iit != iend; ++iit )
              {
                  const auto intersection = *iit ;
                  const GlobalCoordinate iCenter = intersection.geometry().center();
                  if( intersection.boundary() )
                  {
                      const GlobalCoordinate dx = iCenter - elCenter;
                      for( int i = 0; i < GlobalCoordinate::dimension; ++i )
                          hessianInverse[ i ].axpy( dx[ i ], dx );

                  } else {
                      const auto neighbor = intersection.outside();
                      const GlobalCoordinate dx = neighbor.geometry().center() - elCenter;
                      for( int i = 0; i < GlobalCoordinate::dimension; ++i )
                          hessianInverse[ i ].axpy( dx[ i ], dx );
                  }
              }
              hessianInverse.invert();
          }

      }

      template< class Mapper, class Vector >
      Reconstruction< Mapper, Vector > operator() ( const Mapper &mapper, Vector u ) const
      {
        std::vector< Jacobian > du;
        (*this)( mapper, getReference( u ), du );
        return Reconstruction< Mapper, Vector >( mapper, std::move( u ), std::move( du ) );
      }

      template< class Mapper, class Vector >
      void operator () ( const Mapper &mapper, const Vector &u, std::vector< Jacobian > &du ) const
      {
        Vector factor;
        this->operator()( mapper, u, du, factor );
      }

      template< class Mapper, class ElementContext, class Vector >
      void operator () ( const Mapper &mapper, ElementContext& elemCtx,
                         const Vector &u, std::vector< Jacobian > &du, Vector& factor, const bool recompute = true ) const
      {
        //du.resize( u.size() );

        std::vector< std::pair< GlobalCoordinate, StateVector > > differences;

        size_t duCounter = 0;
        const auto end = gridView().template end<0, Dune::InteriorBorder_Partition>();
        for(auto it = gridView().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it)
        {
          const auto element = *it;
          const std::size_t elIndex = mapper.index( element );
          const GlobalCoordinate elCenter = element.geometry().center();

          //const unsigned int entityIndex = dofMapper_.index( entity );

          elemCtx.updateStencil(element);

          const auto& hessianInverse = hessianInverse_[ elIndex ];
          const auto& stencil = elemCtx.stencil(0);
          size_t numDof = stencil.numDof();
          std::vector< Dune::FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > > negGradient( 0 );
          for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx, ++duCounter)
          {
            differences.clear();
            Dune::FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > negGradient( 0 );

            unsigned int nbIndexLocal = 1;
            const auto iend = gridView().iend( element );
            for( auto iit = gridView().ibegin( element ); iit != iend; ++iit, ++nbIndexLocal )
            {
              const auto intersection = *iit ;
              const GlobalCoordinate iCenter = intersection.geometry().center();
              if( intersection.boundary() )
              /*{
                const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
                StateVector uBv =  boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] );
                for( int j = 0; j < StateVector::dimension; ++j )
                {
                    uBv[j] = Toolbox::createConstant(Toolbox::value(uBv[j]));
                }
                StateVector diff = u[ elIndex ] - u[ elIndex ];
                if( dofIdx > 0 )
                {
                  for( int j = 0; j < StateVector::dimension; ++j )
                  {
                    diff[ j ] = Toolbox::createConstant(Toolbox::value(diff[ j ]));
                  }
                }
                const GlobalCoordinate dx = iCenter - elCenter;
                for( int i = 0; i < StateVector::dimension; ++i )
                  negGradient[ i ].axpy( diff[ i ], dx );
                differences.emplace_back( dx, diff );
              }
              else*/ if( intersection.neighbor() )
              {
                const auto neighbor = intersection.outside();
                const GlobalCoordinate dx = neighbor.geometry().center() - elCenter;
                StateVector uNb = u[ mapper.index( neighbor ) ];
                StateVector uEl = u[ elIndex ];
                for( int j = 0; j < StateVector::dimension; ++j )
                {
                  if( dofIdx > 0 )
                    uEl[j] = Toolbox::createConstant(Toolbox::value(uEl[j]));

                  if( dofIdx != nbIndexLocal )
                    uNb[j] = Toolbox::createConstant(Toolbox::value(uNb[j]));
                }
                StateVector diff = uNb - uEl;
                for( int i = 0; i < StateVector::dimension; ++i )
                  negGradient[ i ].axpy( diff[ i ], dx );
                differences.emplace_back( dx, diff );
              }
            }

            for( int j = 0; j < StateVector::dimension; ++j )
              hessianInverse.mv( negGradient[ j ], du[ duCounter ][ j ] );

#if 1
            if( ! factor.empty() )
            {
              StateVector& localfactor = factor[ elIndex ];
              if( recompute )
              {
                localfactor = 1;
                du[ duCounter ] = limit( differences, du[ duCounter ], localfactor );
              }
              else
              {
                for( int j = 0; j < StateVector::dimension; ++j )
                  du[duCounter][ j ] *= localfactor[ j ];
              }

            }
            else

#endif
            {
              StateVector localfactor( 1 );
              du[ duCounter ] = limit( differences, du[ duCounter ], localfactor );
            }
          } // end for over dofIdx

        }
        //auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        //gridView().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
      }

      const GridView &gridView () const { return gridView_; }

    private:
      inline Jacobian limit ( const std::vector< std::pair< GlobalCoordinate, StateVector > > &differences, Jacobian du, StateVector& factor ) const
      {
        using std::abs;
        for( const auto &difference : differences )
        {
          StateVector slope;
          du.mv( difference.first, slope );
          for( int j = 0; j < StateVector::dimension; ++j )
          {
            const Field localFactor = limiterFunction_( slope[ j ], difference.second[ j ] );
            factor[ j ] = std::min( factor[ j ], localFactor );
            /*
            if( slope[ j ] * difference.second[ j ] <= 0 )
              return Jacobian( 0 );

            if( abs( difference.second[ j ] ) < abs( slope[ j ] ) )
              factor[ j ] = std::min( factor[ j ], difference.second[ j ] / slope[ j ] );
            */
          }
        }

        for( int j = 0; j < StateVector::dimension; ++j )
          du[ j ] *= factor[ j ];
        return du;
      }

      GridView gridView_;
      BoundaryValue boundaryValue_;
      LimiterFunction limiterFunction_;
      std::vector< HessianInverseType > hessianInverse_;

    };



    // limitedLeastSquaresReconstruction
    // ---------------------------------

    template< class SV, class GV, class BV >
    inline static LimitedLeastSquaresReconstruction< GV, SV, BV > limitedLeastSquaresReconstruction ( const GV &gridView, BV boundaryValue )
    {
      return LimitedLeastSquaresReconstruction< GV, SV, BV >( gridView, std::move( boundaryValue ) );
    }


    // LimitedLeastSquaresReconstruction
    // ---------------------------------

    template< class GV, class SV, class BV >
    class LimitedLeastSquaresReconstructionLocal
    {
      typedef LimitedLeastSquaresReconstructionLocal< GV, SV, BV > This;

    public:
      typedef GV GridView;
      typedef SV StateVector;
      typedef BV BoundaryValue;

      typedef FieldVector< typename GridView::ctype, GridView::dimensionworld > GlobalCoordinate;

      typedef typename GridView::Intersection Intersection;

      typedef typename FieldTraits< StateVector >::field_type Field;
      typedef FieldVector< Field, GlobalCoordinate::dimension > Jacobian;

      typedef Limiters::MinModLimiter< Field >  LimiterFunction;
      typedef Opm::MathToolbox<Field> Toolbox;

      typedef Dune::FieldMatrix< double, GlobalCoordinate::dimension, GlobalCoordinate::dimension > HessianInverseType;

      template< class Mapper, class Vector >
      using Reconstruction = PiecewiseLinearFunction< GridView, Mapper, Vector, Jacobian >;

      LimitedLeastSquaresReconstructionLocal ( const GridView &gridView, BoundaryValue boundaryValue )
        : gridView_( gridView ), boundaryValue_( std::move( boundaryValue ) ), limiterFunction_(),
          hessianInverse_( gridView.indexSet().size( 0 ), HessianInverseType(0 ) )
      {
          const auto end = gridView.template end<0, Dune::InteriorBorder_Partition>();
          for(auto it = gridView.template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it)
          {
              const auto element = *it;
              const int elIdx = gridView.indexSet().index( element );
              const GlobalCoordinate elCenter = element.geometry().center();
              const auto iend = gridView.iend( element );
              auto& hessianInverse = hessianInverse_[ elIdx ];
              for( auto iit = gridView.ibegin( element ); iit != iend; ++iit )
              {
                  const auto intersection = *iit ;
                  const GlobalCoordinate iCenter = intersection.geometry().center();
                  if( intersection.boundary() )
                  {
                      GlobalCoordinate dx = iCenter - elCenter;
                      dx *= boundaryFactor_;

                      for( int i = 0; i < GlobalCoordinate::dimension; ++i )
                          hessianInverse[ i ].axpy( dx[ i ], dx );

                  }
                  else {
                      const auto neighbor = intersection.outside();
                      const GlobalCoordinate dx = neighbor.geometry().center() - elCenter;
                      for( int i = 0; i < GlobalCoordinate::dimension; ++i )
                          hessianInverse[ i ].axpy( dx[ i ], dx );
                  }
              }
              hessianInverse.invert();
          }
      }

      template< class ElementContext, class Vector >
      Reconstruction< ElementContext, Vector > operator() ( const ElementContext& elemCtx, unsigned elIndex, unsigned timeIdx, Vector u ) const
      {
        Jacobian du;
        (*this)( elemCtx, elIndex, timeIdx, getReference( u ), du );
        return Reconstruction< ElementContext, Vector >( elemCtx, elIndex, timeIdx, std::move( u ), std::move( du ) );
      }

      template< class ElementContext, class Vector >
      void operator () ( const ElementContext& elemCtx,  unsigned elIndex, unsigned timeIdx, const Vector &u, Jacobian &du ) const
      {
        Field factor;
        this->operator()( elemCtx, elIndex, timeIdx , u, du, factor );
      }

      template< class ElementContext, class Vector >
      void operator () ( const ElementContext& elemCtx, unsigned elIndex, unsigned timeIdx, const Vector &u, Jacobian &du, Field& factor, const bool recompute = true ) const
      {
        std::vector< std::pair< GlobalCoordinate, StateVector > > differences;
        std::vector< std::pair< StateVector , GlobalCoordinate > > neighbours;

        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto element = stencil.element(elIndex);
        const size_t globalIdx = elemCtx.globalSpaceIndex( elIndex, 0 );
        const auto& hessianInverse = hessianInverse_[ globalIdx ];
        const size_t stencilSize = stencil.numDof();
        size_t neighboursCounter = 0;

        const GlobalCoordinate elCenter = element.geometry().center();
        Dune::FieldVector< Field, GlobalCoordinate::dimension > negGradient( 0 );
        differences.clear();
        neighbours.clear();
        const auto iend = gridView().iend( element );
        for( auto iit = gridView().ibegin( element ); iit != iend; ++iit ) {
            const auto intersection = *iit;
            const GlobalCoordinate iCenter = intersection.geometry().center();

            if (intersection.neighbor()) {
                const auto neighbor = intersection.outside();
                //const size_t outsideIdx = stencil.global
                //const size_t insideIdx = intersection.indexInInside();
                const size_t neighborIdx = stencil.globalToLocal(neighbor);
                const GlobalCoordinate dx = neighbor.geometry().center() - elCenter;
                const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
                const StateVector uBv = u[neighborIdx];

                const StateVector du1 = u[neighborIdx] - u[elIndex];
                negGradient.axpy(du1, dx);
                differences.emplace_back(dx, du1);
                neighbours.emplace_back(uBv, iNormal);
                neighboursCounter++;
            }
        }

        if ( (stencilSize - neighboursCounter) > 1.0 ) {
            const auto iend = gridView().iend( element );
            for( auto iit = gridView().ibegin( element ); iit != iend; ++iit ) {
                const auto intersection = *iit;
                const GlobalCoordinate iCenter = intersection.geometry().center();

                if( intersection.boundary() )
                {
                    const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
                    GlobalCoordinate iNormalOpposite;
                    for (auto ii=0; ii < iNormal.size(); ++ii)
                        iNormalOpposite [ii] = (-1.0) * iNormal[ii];
                    for( const auto &neighbour : neighbours ) {
                        double product = neighbour.second[0]*iNormalOpposite[0] + neighbour.second[1]*iNormalOpposite[1] + neighbour.second[2]*iNormalOpposite[2];
                        if (product < 0.0) {

                            const StateVector du1 = neighbour.first - u[elIndex];
                            const GlobalCoordinate dx = (iCenter - elCenter) + (iCenter - elCenter);

                            negGradient.axpy(du1, dx);
                            differences.emplace_back(dx, du1);
                        }
                    }
                    }
            }
        }


//            if( intersection.boundary() )
//            {
//                const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
//
//                StateVector uBv = u[elIndex];
//
//                const StateVector du1 = uBv - u[ elIndex ];
//                const GlobalCoordinate dx = boundaryFactor_ * (iCenter - elCenter);
//
//                negGradient.axpy( du1, dx );
//                differences.emplace_back( dx, du1 );
//
//            }
//            else if( intersection.neighbor() )
//            {
//                const auto neighbor = intersection.outside();
//                //const size_t outsideIdx = stencil.global
//                //const size_t insideIdx = intersection.indexInInside();
//                const size_t neighborIdx = stencil.globalToLocal(neighbor);
//                const GlobalCoordinate dx = neighbor.geometry().center() - elCenter;
//
//                const StateVector du1 = u[ neighborIdx ] - u[ elIndex ];
//                negGradient.axpy( du1, dx );
//                differences.emplace_back( dx, du1 );
//            }
//        }

        hessianInverse.mv( negGradient, du );

        if( recompute )
        {
            factor = 1;
            du = limit( differences, du, factor );
        }
        else
        {
            du *= factor;
        }


        //auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        //gridView().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
    }

      const GridView &gridView () const { return gridView_; }

    private:
      inline Jacobian limit ( const std::vector< std::pair< GlobalCoordinate, StateVector > > &differences, Jacobian du, StateVector& factor ) const
      {
        using std::abs;
        for( const auto &difference : differences )
        {
          StateVector slope = du * difference.first;
          //du.mv( difference.first, slope );

          const Field localFactor = limiterFunction_( slope, difference.second );
          factor = std::min( factor, localFactor );
          /*
            if( slope[ j ] * difference.second[ j ] <= 0 )
              return Jacobian( 0 );

            if( abs( difference.second[ j ] ) < abs( slope[ j ] ) )
              factor[ j ] = std::min( factor[ j ], difference.second[ j ] / slope[ j ] );
              */
        }

        du *= factor;
        return du;
      }

      GridView gridView_;
      BoundaryValue boundaryValue_;
      LimiterFunction limiterFunction_;
      std::vector< HessianInverseType > hessianInverse_;
      const double boundaryFactor_ = 2.0;

    };



    // limitedLeastSquaresReconstruction
    // ---------------------------------

    template< class SV, class GV, class BV >
    inline static LimitedLeastSquaresReconstructionLocal< GV, SV, BV > limitedLeastSquaresReconstructionLocal ( const GV &gridView, BV boundaryValue )
    {
      return LimitedLeastSquaresReconstructionLocal< GV, SV, BV >( gridView, std::move( boundaryValue ) );
    }

  } // namespace FV

} // namespace Dune

#endif // #ifndef DUNE_FV_LEASTSQUARESRECONSTRUCTION_HH
