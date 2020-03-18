#ifndef DUNE_FV_LPRECONSTRUCTION_HH
#define DUNE_FV_LPRECONSTRUCTION_HH

#include <cassert>
#include <cstddef>

#include <numeric>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#if ! DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 5 )
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>
#endif

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/common/getreference.hh>

#include <dune/geometry/axisalignedreferencefaces.hh>

#include <dune/grid/utility/vectorcommdatahandle.hh>

#include <dune/optim/activeindexmapper.hh>
#include <dune/optim/common/smallobject.hh>
#include <dune/optim/constraint/linear.hh>
#include <dune/optim/lp.hh>
#include <dune/optim/solver/gaussjordan.hh>

#include <dune/fv/function/piecewiselinear.hh>

namespace Dune
{

namespace FV
{

// LPReconstruction
// ----------------

/**
     * \class LPReconstruction
     * \brief Minmod-type reconstruction based on linear programming problem
     *
     * The LPReconstruction was suggested in the following paper:
     * \code
     * @article{Chen:IntegratedLinearReconstruction,
     *   author  = {Chen, L. and Li, R.},
     *   title   = {An Integrated Linear Reconstruction for Finite Volume scheme
     *              on Unstructured Grids},
     *   journal = {J. Sci. Comput.},
     *   year    = {2016},
     *   volume  = {68},
     *   pages   = {1172--1197},
     *   doi     = {10.1007/s10915-016-0173-1}
     * }
     * \endcode
     **/
template< class GV, class SV, class BV >
class LPReconstruction
{
    typedef LPReconstruction< GV, SV, BV > This;

public:
    typedef GV GridView;
    typedef SV StateVector;
    typedef BV BoundaryValue;

    typedef FieldVector< typename GridView::ctype, GridView::dimensionworld > GlobalCoordinate;

    typedef typename GridView::Intersection Intersection;

    typedef typename FieldTraits< StateVector >::field_type Field;
    typedef typename FieldTraits< StateVector >::real_type Real;
    typedef FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > Jacobian;

    typedef Opm::MathToolbox<Field> Toolbox;


    template< class Mapper, class Vector >
    using Reconstruction = PiecewiseLinearFunction< GridView, Mapper, Vector, std::vector< Jacobian > >;


    static const int dimension = GridView::dimension;

private:
    typedef Optim::LinearConstraintEval< GlobalCoordinate, Field > Constraint;
    typedef std::vector< Constraint > Constraints;

    typedef Optim::GaussJordanSolver< FieldMatrix< Field, GlobalCoordinate::dimension, GlobalCoordinate::dimension > > Solver;
    typedef Optim::LinearProgramming< Solver, false > LP;

public:
    LPReconstruction ( const GridView &gridView, BoundaryValue boundaryValue, Real tolerance )
        : gridView_( gridView ),
          boundaryValue_( std::move( boundaryValue ) ),
          tolerance_( std::move( tolerance ) ),
          lp_( tolerance_ ),
          faceAxes_( LocalGeometryTypeIndex::size( dimension ) )
    {
#if DUNE_VERSION_NEWER( DUNE_GEOMETRY, 2, 5 )
        const unsigned int numTopo = Impl::numTopologies( dimension );
#else
        const unsigned int numTopo = GenericGeometry::numTopologies( dimension );
#endif
        std::unique_ptr< unsigned int[] > faceIndices( new unsigned int[ dimension * numTopo ] );
        std::unique_ptr< unsigned int[] > numFaces( new unsigned int[ numTopo ] );
        axisAlignedReferenceFaces( dimension, faceIndices.get(), numFaces.get() );
        for( unsigned int topologyId = 0; topologyId < numTopo; ++topologyId )
        {
            const GeometryType type( topologyId, dimension );
            std::vector< unsigned int > &faceAxes = faceAxes_[ LocalGeometryTypeIndex::index( type ) ];
            faceAxes.resize( numFaces[ topologyId ], dimension );
            for( int i = 0; i < dimension; ++i )
                faceAxes[ faceIndices[ topologyId*dimension + i ] ] = i;
        }
    }

    template< class Mapper, class Vector >
    Reconstruction< Mapper, Vector > operator() ( const Mapper &mapper, Vector u ) const
    {
        std::vector< Jacobian > du;
        (*this)( mapper, getReference( u ), du );
        return Reconstruction< Mapper, Vector >( mapper, std::move( u ), std::move( du ) );
    }

    template< class Mapper, class ElementContext, class Vector >
    void operator () ( const Mapper &mapper, ElementContext& elemCtx, const Vector &u, std::vector< Jacobian > &du ) const
    {
        //du.resize( u.size() );

        size_t duCounter = 0;

        std::vector< std::pair< GlobalCoordinate, StateVector > > differences;
        Constraints constraints;

        const auto end = gridView().template end< 0, Dune::InteriorBorder_Partition >();
        for( auto it = gridView().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it )
            //for( const auto element : elements( gridView(), Partitions::interiorBorder ) )
        {
            const auto element = *it ;
            elemCtx.updateStencil(element);

            const std::size_t elIndex = mapper.index( element );
            const GlobalCoordinate elCenter = element.geometry().center();

            const auto& stencil = elemCtx.stencil(0);
            size_t numDof = stencil.numDof();
            for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx, ++duCounter)
            {

            std::array< unsigned int, dimension+1 > select;
            select.fill(0);
            const std::vector< unsigned int > &faceAxes = faceAxes_[ LocalGeometryTypeIndex::index( element.type() ) ];
            if( !faceAxes.empty() )
            {
                unsigned int nbIndexLocal = 1;
                differences.clear();
                const auto iend = gridView().iend( element );
                for( auto iit = gridView().ibegin( element ); iit != iend; ++iit, ++nbIndexLocal )
                    //for( const auto intersection : intersections( gridView(), element ) )
                {
                    const auto intersection = *iit;

                    select[ faceAxes[ intersection.indexInInside() ] ] = differences.size();

                    if( intersection.boundary() )
                    {
                        const GlobalCoordinate iCenter = intersection.geometry().center();
                        const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
                        StateVector uBnd = boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] );
                        for( int j = 0; j < StateVector::dimension; ++j )
                        {
                            uBnd[j] = Toolbox::createConstant(Toolbox::value(uBnd[j]));
                        }
                        StateVector diff = uBnd - u[ elIndex ];
                        if( dofIdx > 0 )
                        {
                            for( int j = 0; j < StateVector::dimension; ++j )
                            {
                                diff[ j ] = Toolbox::createConstant(Toolbox::value(diff[ j ]));
                            }
                        }
                        differences.emplace_back( iCenter - elCenter, diff );
                    }
                    else if( intersection.neighbor() )
                    {
                        const auto neighbor = intersection.outside();
                        const GlobalCoordinate nbCenter = neighbor.geometry().center();
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
                        differences.emplace_back( nbCenter - elCenter, diff );
                    }
                }
            }
            else
            {
                Dune::ReservedVector< GlobalCoordinate, dimension > onb;

                unsigned int nbIndexLocal = 1;

                differences.clear();
                const auto iend = gridView().iend( element );
                for( auto iit = gridView().ibegin( element ); iit != iend; ++iit, ++nbIndexLocal )
                    //for( const auto intersection : intersections( gridView(), element ) )
                {
                    const auto intersection = *iit;

                    select[ onb.size() ] = differences.size();

                    if( intersection.boundary() )
                    {
                        const GlobalCoordinate iCenter = intersection.geometry().center();
                        const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
                        StateVector uBnd = boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] );
                        StateVector diff = uBnd - u[ elIndex ];
                        if( dofIdx > 0 )
                        {
                            for( int j = 0; j < StateVector::dimension; ++j )
                            {
                                diff[ j ] = Toolbox::createConstant(Toolbox::value(diff[ j ]));
                            }
                        }
                        differences.emplace_back( iCenter - elCenter, diff );
                    }
                    else if( intersection.neighbor() )
                    {
                        const auto neighbor = intersection.outside();
                        const GlobalCoordinate nbCenter = neighbor.geometry().center();
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
                        differences.emplace_back( nbCenter - elCenter, diff );
                    }

                    if( onb.size() < dimension )
                    {
                        GlobalCoordinate dx = differences.back().first;
                        for( const GlobalCoordinate &v : onb )
                            dx.axpy( -(dx*v), v );

                        const auto dxNorm = dx.two_norm();
                        if( dxNorm >= tolerance_ )
                            onb.push_back( dx /= dxNorm );
                    }
                }
            }

            // reserve memory for constraints
            const std::size_t numConstraints = differences.size();
            constraints.resize( 2u*numConstraints );
            Optim::ActiveIndexMapper< SmallObjectAllocator< unsigned int > > active( GlobalCoordinate::dimension, constraints.size() );
            for( int j = 0; j < StateVector::dimension; ++j ) {
                GlobalCoordinate negGradient(0);
                for (std::size_t i = 0u; i < numConstraints; ++i) {
                    const Field sign = (differences[i].second[j] >= 0 ? 1 : -1);

                    negGradient.axpy(Toolbox::value(sign), differences[i].first);

                    constraints[2 * i].normal() = differences[i].first;
                    constraints[2 * i].normal() *= Toolbox::value(sign);
                    constraints[2 * i].rhs() = sign * differences[i].second[j];

                    constraints[2 * i + 1].normal() = constraints[2 * i].normal();
                    constraints[2 * i + 1].normal() *= Toolbox::value(Field(-1));
                    constraints[2 * i + 1].rhs() = 0;
                }

                // activate GlobalCoordinate::dimension constraints active in the origin
                active.clear();
                for (int i = 0; i < dimension; ++i)
                    active.activate(2 * select[i] + 1);

                // solve
                du[duCounter][j] = 0;
                lp_( negGradient, constraints, du[ duCounter ][ j ], active );
              //  std::cout << " du[ " << duCounter << " ][ " << j << " ] " << du[ duCounter ][ j ] << std::endl;
            }

            }


        }

 //       std::cout << " ducounter " << duCounter << std::endl;

        //auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        //gridView().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
    }

    const GridView &gridView () const { return gridView_; }

private:
    GridView gridView_;
    BoundaryValue boundaryValue_;
    Real tolerance_;
    LP lp_;
    std::vector< std::vector< unsigned int > > faceAxes_;
};



// lpReconstruction
// ----------------

template< class SV, class GV, class BV >
inline static LPReconstruction< GV, SV, BV > lpReconstruction ( const GV &gridView, BV boundaryValue, typename FieldTraits< SV >::real_type tolerance )
{
    return LPReconstruction< GV, SV, BV >( gridView, std::move( boundaryValue ), std::move( tolerance ) );
}

// LPReconstruction
// ----------------

/**
   * \class LPReconstruction
   * \brief Minmod-type reconstruction based on linear programming problem
   *
   * The LPReconstruction was suggested in the following paper:
   * \code
   * @article{Chen:IntegratedLinearReconstruction,
   *   author  = {Chen, L. and Li, R.},
   *   title   = {An Integrated Linear Reconstruction for Finite Volume scheme
   *              on Unstructured Grids},
   *   journal = {J. Sci. Comput.},
   *   year    = {2016},
   *   volume  = {68},
   *   pages   = {1172--1197},
   *   doi     = {10.1007/s10915-016-0173-1}
   * }
   * \endcode
   **/
template< class GV, class SV, class BV >
class LPReconstructionLocal
{
    typedef LPReconstructionLocal< GV, SV, BV > This;

public:
    typedef GV GridView;
    typedef SV StateVector;
    typedef BV BoundaryValue;

    typedef FieldVector< typename GridView::ctype, GridView::dimensionworld > GlobalCoordinate;

    typedef typename GridView::Intersection Intersection;

    typedef typename FieldTraits< StateVector >::field_type Field;
    typedef typename FieldTraits< StateVector >::real_type Real;
    typedef FieldVector< Field, GlobalCoordinate::dimension > Jacobian;

    typedef Opm::MathToolbox<Field> Toolbox;


    template< class Mapper, class Vector >
    using Reconstruction = PiecewiseLinearFunction< GridView, Mapper, Vector, std::vector< Jacobian > >;


    static const int dimension = GridView::dimension;

private:
    typedef Optim::LinearConstraintEval< GlobalCoordinate, Field > Constraint;
    typedef std::vector< Constraint > Constraints;

    typedef Optim::GaussJordanSolver< FieldMatrix< Field, GlobalCoordinate::dimension, GlobalCoordinate::dimension > > Solver;
    typedef Optim::LinearProgramming< Solver, false > LP;

public:
    LPReconstructionLocal ( const GridView &gridView, BoundaryValue boundaryValue, Real tolerance )
        : gridView_( gridView ),
          boundaryValue_( std::move( boundaryValue ) ),
          tolerance_( std::move( tolerance ) ),
          lp_( tolerance_ ),
          faceAxes_( LocalGeometryTypeIndex::size( dimension ) )
    {
#if DUNE_VERSION_NEWER( DUNE_GEOMETRY, 2, 5 )
        const unsigned int numTopo = Impl::numTopologies( dimension );
#else
        const unsigned int numTopo = GenericGeometry::numTopologies( dimension );
#endif
        std::unique_ptr< unsigned int[] > faceIndices( new unsigned int[ dimension * numTopo ] );
        std::unique_ptr< unsigned int[] > numFaces( new unsigned int[ numTopo ] );
        axisAlignedReferenceFaces( dimension, faceIndices.get(), numFaces.get() );
        for( unsigned int topologyId = 0; topologyId < numTopo; ++topologyId )
        {
            const GeometryType type( topologyId, dimension );
            std::vector< unsigned int > &faceAxes = faceAxes_[ LocalGeometryTypeIndex::index( type ) ];
            faceAxes.resize( numFaces[ topologyId ], dimension );
            for( int i = 0; i < dimension; ++i )
                faceAxes[ faceIndices[ topologyId*dimension + i ] ] = i;
        }

        differences_.reserve( 100 );
        constraints_.reserve( 100 );
    }

    template< class ElementContext, class Vector >
    Reconstruction< ElementContext, Vector > operator() ( const ElementContext& elemCtx, unsigned elIndex, unsigned timeIdx, Vector u ) const
    {
        Jacobian du;
        (*this)( elemCtx, elIndex, timeIdx, getReference( u ), du );
        return Reconstruction< ElementContext, Vector >( elemCtx, elIndex, timeIdx, std::move( u ), std::move( du ) );
    }

    template< class ElementContext, class Vector >
    void operator () ( const ElementContext& elemCtx, unsigned elIndex, unsigned timeIdx, const Vector &u, Jacobian &du ) const
    {
        //std::vector< std::pair< GlobalCoordinate, StateVector > > differences;
        //Constraints constraints;
        std::vector< std::pair< GlobalCoordinate, StateVector > >& differences = differences_;
        Constraints& constraints = constraints_;
        constraints.clear();

        const auto& stencil = elemCtx.stencil(timeIdx);

        const auto element = stencil.element(elIndex);

        const GlobalCoordinate elCenter = element.geometry().center();

        std::array< unsigned int, dimension+1 > select;
        select.fill(0);
        const std::vector< unsigned int > &faceAxes = faceAxes_[ LocalGeometryTypeIndex::index( element.type() ) ];
        if( !faceAxes.empty() )
        {
            differences.clear();
            const auto iend = gridView().iend( element );
            for( auto iit = gridView().ibegin( element ); iit != iend; ++iit )
            {
                const auto intersection = *iit;

                select[ faceAxes[ intersection.indexInInside() ] ] = differences.size();

                if( intersection.boundary() )
                {
                    const GlobalCoordinate iCenter = intersection.geometry().center();
                    //const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
                    //StateVector uBnd = u[ elIndex ];
                    // differences.emplace_back( iCenter - elCenter, uBnd - u[ elIndex ] );
                    differences.emplace_back( iCenter - elCenter, StateVector(0) );
                }
                else if( intersection.neighbor() )
                {
                    const auto neighbor = intersection.outside();
                    const GlobalCoordinate nbCenter = neighbor.geometry().center();
                    const size_t neighborIdx = stencil.globalToLocal( neighbor );
                    differences.emplace_back( nbCenter - elCenter, u[ neighborIdx ] - u[ elIndex ] );
                }
            }
        }
        else
        {
            Dune::ReservedVector< GlobalCoordinate, dimension > onb;

            differences.clear();
            const auto iend = gridView().iend( element );
            for( auto iit = gridView().ibegin( element ); iit != iend; ++iit )
            {
                const auto intersection = *iit;

                select[ onb.size() ] = differences.size();

                if( intersection.boundary() )
                {
                    const GlobalCoordinate iCenter = intersection.geometry().center();
                    //StateVector uBnd = u[ elIndex ];
                    //differences.emplace_back( iCenter - elCenter, uBnd - u[ elIndex ] );
                    differences.emplace_back( iCenter - elCenter, StateVector(0) );
                }
                else if( intersection.neighbor() )
                {
                    const auto neighbor = intersection.outside();
                    //const GlobalCoordinate nbCenter = neighbor.geometry().center();
                    const size_t neighborIdx = stencil.globalToLocal( neighbor );
                    const GlobalCoordinate& nbCenter = stencil.subControlVolume( neighborIdx ).center();
                    differences.emplace_back( nbCenter - elCenter, u[ neighborIdx ] - u[ elIndex ] );
                }

                if( onb.size() < dimension )
                {
                    GlobalCoordinate dx = differences.back().first;
                    for( const GlobalCoordinate &v : onb )
                        dx.axpy( -(dx*v), v );

                    const auto dxNorm = dx.two_norm();
                    if( dxNorm >= tolerance_ )
                        onb.push_back( dx /= dxNorm );
                }
            }
        }

        // reserve memory for constraints
        const std::size_t numConstraints = differences.size();
        constraints.resize( 2u*numConstraints );
        Optim::ActiveIndexMapper< SmallObjectAllocator< unsigned int > > active( GlobalCoordinate::dimension, constraints.size() );

        GlobalCoordinate negGradient( 0 );
        for( std::size_t i = 0u; i < numConstraints; ++i )
        {
            const double sign = (differences[ i ].second >= 0 ? 1 : -1);

            negGradient.axpy( sign, differences[ i ].first );

            constraints[ 2*i ].normal() = differences[ i ].first;
            constraints[ 2*i ].normal() *= sign;
            constraints[ 2*i ].rhs() = sign*differences[ i ].second;

            constraints[ 2*i+1 ].normal() = constraints[ 2*i ].normal();
            constraints[ 2*i+1 ].normal() *= -1;
            constraints[ 2*i+1 ].rhs() = 0;
        }

        // activate GlobalCoordinate::dimension constraints active in the origin
        active.clear();
        for( int i = 0; i < dimension; ++i )
            active.activate( 2*select[ i ]+1 );

        // solve
        du = 0;
        lp_( negGradient, constraints, du, active );

        //auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        //gridView().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
    }

    const GridView &gridView () const { return gridView_; }

private:
    GridView gridView_;
    BoundaryValue boundaryValue_;
    Real tolerance_;
    LP lp_;
    std::vector< std::vector< unsigned int > > faceAxes_;

    mutable std::vector< std::pair< GlobalCoordinate, StateVector > > differences_;
    mutable Constraints constraints_;
};



// lpReconstruction
// ----------------

template< class SV, class GV, class BV >
inline static LPReconstructionLocal< GV, SV, BV > lpReconstructionLocal ( const GV &gridView, BV boundaryValue, typename FieldTraits< SV >::real_type tolerance )
{
    return LPReconstructionLocal< GV, SV, BV >( gridView, std::move( boundaryValue ), std::move( tolerance ) );
}


} // namespace FV

} // namespace Dune

#endif // #ifndef DUNE_FV_LPRECONSTRUCTION_HH
