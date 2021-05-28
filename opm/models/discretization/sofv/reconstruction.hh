#ifndef DUNE_FEM_LIMITER_HH
#define DUNE_FEM_LIMITER_HH

//#include <dune/fem/space/finitevolume.hh>
//#include <dune/fem/io/parameter.hh>

#include <opm/material/densead/Evaluation.hpp>

#include <dune/fv/leastsquaresreconstruction.hh>
#include <dune/fv/lpreconstruction.hh>
//#include <dune/fv/qpreconstruction.hh>

#include "limiterutility.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Limited reconstruction.
   *
   * \ingroup PassBased
   */
  template <class TypeTag>
  class LimitedReconstruction
  {
    typedef LimitedReconstruction< TypeTag >     ThisType;
  public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid)             GridType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView)         GridViewType;
    typedef typename GET_PROP_TYPE(TypeTag, GridPart)         GridPartType;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper)        DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation)       Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager)   ThreadManager;


    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange = PrimaryVariables::dimension}; //  };
    enum { numDeri = Evaluation :: numVars };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar)           DomainFieldType; //should be Scalar

    typedef Dune::FieldVector< Evaluation, dimRange >         EvalDimVector;
    typedef Opm::MathToolbox<Evaluation> Toolbox;


    typedef LimiterUtility< DomainFieldType, Evaluation, dimDomain, dimRange >  LimiterUtilityType;

    typedef typename LimiterUtilityType :: RangeFieldType   RangeFieldType;
    typedef typename LimiterUtilityType :: RangeType        RangeType;

    typedef typename LimiterUtilityType :: DomainType       DomainType;

    typedef std::vector< Dune::FieldVector< Evaluation, dimRange > >      SolutionVector;


    typedef typename LimiterUtilityType :: GradientEvalType             GradientEvalType;
    typedef typename LimiterUtilityType :: SingleGradientEvalType       SingleGradientEvalType;
    typedef typename LimiterUtilityType :: CheckType         CheckType;
    typedef typename LimiterUtilityType :: ComboSetType      ComboSetType;
    typedef typename LimiterUtilityType :: KeyType           KeyType;

    typedef std::map< int, ComboSetType > ComboSetMapType ;

    typedef MinModLimiter< Evaluation > LimiterFunctionType;
    //typedef SuperBeeLimiter< Field > LimiterFunctionType;

    typedef typename GridPartType :: template Codim< 0 > :: EntityType   EntityType;
    typedef typename EntityType :: Geometry                              Geometry;

    typedef typename LimiterUtilityType::MatrixStorage MatrixCacheEntry;
    typedef std::map< KeyType, MatrixCacheEntry > MatrixCacheType;

    static const bool StructuredGrid     = false; //GridPartCapabilities::isCartesian< GridPartType >::v;

    //! type of cartesian grid checker
    typedef CheckCartesian< GridPartType >  CheckCartesianType;

    enum ReconstructionScheme
    {
      firstOrder = 0, // first order scheme
      secondRecon = 1, // 2nd order selected reconstruction
      secondLeast = 2, // 2nd order least squares
      secondLP    = 3,  // 2nd order linear programming
      secondQP    = 4   // 2nd order quadratic programming
    };

    static ReconstructionScheme getReconstructionSchemeId()
    {
      static const std::string reconstructions[] = { "1st", "2nd-R", "2nd-LS", "2nd-LP", "2nd-QP" };
      return (ReconstructionScheme) Dune::Fem::Parameter::getEnum( "finitevolume.reconstruction", reconstructions );
    }

    static const std::string& schemeName (const int scheme)
    {
      static const std::string reconstructions[] = { "1st", "2nd-R", "2nd-LS", "2nd-LP", "2nd-QP" };
      return reconstructions[ scheme ];
    }

    struct BoundaryValue
    {
      const ThisType& op_;
      BoundaryValue( const ThisType& op ) : op_( op ) {}

      RangeType operator () ( const typename GridViewType::Intersection &i,
                              const DomainType &x,
                              const DomainType &n,
                              const RangeType &uIn ) const
      {
        return uIn;
        //return op_.model().problem().boundaryValue( x, op_.time() );
      }
    };

    typedef Dune::FV::LimitedLeastSquaresReconstruction< GridViewType, RangeType, BoundaryValue > LSRecon;
    typedef Dune::FV::LimitedLeastSquaresReconstructionLocal< GridViewType, Evaluation, BoundaryValue > LSReconLocal;

    typedef Dune::FV::LPReconstruction< GridViewType, RangeType, BoundaryValue > LinearProgramming;
    typedef Dune::FV::LPReconstructionLocal< GridViewType, Evaluation, BoundaryValue > LinearProgrammingLocal;

    //typedef Dune::FV::QPReconstruction< GridViewType, RangeType, BoundaryValue > QuadraticProgramming;

  public:
    LimitedReconstruction( const GridPartType &gridPart, const DofMapper& dofMapper)
      : gridPart_( gridPart )
      , dofMapper_( dofMapper )
      , lsRecon_( static_cast< GridViewType > (gridPart_), BoundaryValue( *this ) )
      , lsReconLocal_( static_cast< GridViewType > (gridPart_), BoundaryValue( *this ) )
      , linProg_( static_cast< GridViewType > (gridPart_), BoundaryValue( *this )
                  , Dune::Fem::Parameter::getValue<double>("finitevolume.linearprogramming.tol", 1e-8 ) )
      //, quadProg_( static_cast< GridViewType > (gridPart_), BoundaryValue( *this )
      //            , Dune::Fem::Parameter::getValue<double>("finitevolume.linearprogramming.tol", 1e-8 ) )
      , limiterFunction_()
      , storedComboSets_()
      , matrixCacheVec_( gridPart_.grid().maxLevel() + 1 )
      , cartesianGrid_( CheckCartesianType::check( gridPart_ ) )
    {
      time_ = 0;

      for( unsigned int i=0; i<ThreadManager::maxThreads(); ++i )
      {
        linProgLocal_.push_back(
            LinearProgrammingLocal( static_cast< GridViewType > (gridPart_), BoundaryValue( *this )
                                  , Dune::Fem::Parameter::getValue<double>("finitevolume.linearprogramming.tol", 1e-8 ) ) );
      }
    }

    double time () const { return time_; }

    //! calculate internal reconstruction based on mean values of U
    void update( //const double time,
                 const SolutionVector & U,
                 const Simulator& simulator,
                 const int schemeId = secondLP,
                 const bool recompute = true )
    {
      // set time
      time_ = 0;//time;

      ElementContext elemCtx(simulator);
      const size_t sizeGrid = gridPart_.indexSet().size( 0 ) ;
      elemCounter_.resize(sizeGrid);
      size_t size = 0; //gridPart_.indexSet().size( 0 ) ;
      const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
      for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
      {
        const auto& entity = *it ;
        elemCtx.updateStencil(entity);
        const auto& stencil = elemCtx.stencil(0);
        size_t numDof = stencil.numDof();
        const unsigned int entityIndex = dofMapper_.index( entity );
        elemCounter_[entityIndex] = size;
        size += numDof;
      }

      if( values_.size() != sizeGrid )
      {
        // resize gradient vector
        values_.resize( sizeGrid );
        centers_.resize( sizeGrid );
        numbers_.resize( sizeGrid );
        factor_.resize( sizeGrid );

        gradients_.resize( size );

        const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
        for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
        {
          const auto& entity = *it ;
          const unsigned int entityIndex = dofMapper_.index( entity );
          const auto center = entity.geometry().center();
          for( int i=0; i<dimDomain; ++i )
            centers_[ entityIndex ][ i ] = center[ i ] ;
        }
      }

      if( schemeId == secondLeast )
      {
//        // helper class for evaluation of average value of discrete function
        // helper class for evaluation of average value of discrete function
        EvalAverage average( U, dofMapper_ );

        const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
        for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
        {
          const auto& entity = *it ;
          const unsigned int entityIndex = dofMapper_.index( entity );
          average.evaluate( entity, values_[ entityIndex ] );
        }

        lsRecon_( gridPart_.indexSet(), elemCtx, values_, gradients_, factor_, recompute );
      }
      else if( schemeId == secondLP )
      {
        // helper class for evaluation of average value of discrete function
        EvalAverage average( U, dofMapper_ );

        const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
        for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
        {
          const auto& entity = *it ;
          const unsigned int entityIndex = dofMapper_.index( entity );
          average.evaluate( entity, values_[ entityIndex ] );
        }

        linProg_( gridPart_.indexSet(), elemCtx, values_, gradients_ ); //, factor_, reco);
       //   std::cout << " gradients " << size << std::endl;
      }
      /*
      else if( schemeId == secondQP )
      {
#if HAVE_DUNE_OPTIM
        // helper class for evaluation of average value of discrete function
        EvalAverage average( U, dofMapper_ );

        const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
        for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
        {
          const auto& entity = *it ;
          const unsigned int entityIndex = dofMapper_.index( entity );
          average.evaluate( entity, values_[ entityIndex ] );
        }

        quadProg_( gridPart_.indexSet(), values_, gradients_ ); //, factor_, reco);
#else
        std::cerr << "dune-optim needed for quadratic programming reconstruction" << std::endl;
        std::abort();
#endif
      }
*/
      else
      {
        // not implemented anymore
        assert( false );
        std::abort();

        /*
        assert( schemeId == secondRecon );
        if( recompute )
        {
          for( size_t i = 0; i<size; ++i ) numbers_[ i ].clear();
        }

        const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
        for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
        {
          //applyLocal( *it, U, recompute );
        }
        */
      }
    }

    //! calculate internal reconstruction based on mean values of U
    template <class ElementContext>
    void updateLocal(const ElementContext& elemCtx,
                     unsigned int upstreamDofIdx,
                     unsigned int timeIdx,
                     std::vector< Evaluation > & U,
                     const DomainType& point,
                     Evaluation& value,
                     const int schemeId = secondLP,
                     const bool recompute = true )
    {
      // set time
      time_ = 0; // maybe timeIdx??? ;//time;

      const auto& stencil = elemCtx.stencil(timeIdx);
      const size_t numDof = stencil.numDof();

      FieldVector<Evaluation, dimDomain> gradient;
      //Evaluation number;

      unsigned focusDofIdx = elemCtx.focusDofIndex();

      // set constant values for non-focused dofs
      for (unsigned int dofIdx = 0; dofIdx < numDof; ++dofIdx)
      {
          if (focusDofIdx != dofIdx)
              U[ dofIdx ] = Toolbox::createConstant(Toolbox::value(U[ dofIdx ]));
      }

      if( schemeId == secondLeast )
      {
        static bool firstCall = true ;
        if( firstCall )
        {
          std::cout << "Reconstruction: Using Least Squares approach!" << std::endl<<std::endl;
          firstCall = false;
        }

        Evaluation factor;
        lsReconLocal_( elemCtx, upstreamDofIdx, timeIdx, U, gradient, factor, recompute );
      }
      else if( schemeId == secondLP )
      {
        static bool firstCall = true ;
        if( firstCall )
        {
          std::cout << "Reconstruction: Using Linear Programming approach!" << std::endl<<std::endl;
          firstCall = false;
        }

        linProgLocal_[ ThreadManager::threadId() ]( elemCtx, upstreamDofIdx, timeIdx, U, gradient );
      }
      else
      {
        std::abort();
        assert( schemeId == secondRecon );
        //applyLocal( globalToLocal, elemCtx, elIndex, values, gradient, factor, number, recompute );
      }

      // Higher order reconstruction
      // f(x) = u + grad * ( x - center )
      value = U[ upstreamDofIdx ];
      DomainType x( point );
      x -= stencil.subControlVolume( upstreamDofIdx ).center();

      value += gradient * x ;
    }

    template <class ElementContext>
    void applyLocal(const std::map<size_t, size_t>& globalToLocal, const ElementContext& elemCtx, const unsigned int entityIndex,
                    std::vector<Evaluation> values, FieldVector<Evaluation, dimDomain> gradient, Evaluation factor, Evaluation number, const bool recompute )
    {

      // helper class for evaluation of average value of discrete function
      //EvalAverage average( U, dofMapper_ );
      // get geometry
      const auto& stencil = elemCtx.stencil(0);
      const auto& entity = stencil.element(entityIndex);
      const Geometry& geo = entity.geometry();

      // cache geometry type
      const GeometryType geomType = entity.type();

      // get bary center of element
      const DomainType& entityCenter = centers_[ entityIndex ];

      // evaluate mean value on current entity
      // average.evaluate( entity, values_[ entityIndex ] );

      // boundary is true if boundary segment was found
      // nonconforming is true if entity has at least one non-conforming intersections
      // cartesian is true if the grid is cartesian and no nonconforming refinement present
      typename LimiterUtilityType::Flags flags( cartesianGrid_, true );


      std::vector< DomainType > barysFull;
      std::vector< RangeType >  nbValsFull;

      // setup neighbors barycenter and mean value for all neighbors
      LimiterUtilityType::setupNeighborValues( globalToLocal, gridPart_, elemCtx, entityIndex, values, entityCenter, values[ entityIndex ],
                                               centers_, StructuredGrid, flags, barys_, nbVals_, barysFull, nbValsFull);

      // mark entity as finished, even if not limited everything necessary was done
      //assert( entityIndex  < visited_.size() );
      //visited_[ entityIndex ] = true ;

      ComboSetType& comboSet = storedComboSets_[ nbVals_.size() ];
      if( comboSet.empty() )
      {
        // create combination set
        LimiterUtilityType::buildComboSet( nbVals_.size(), comboSet );
      }
      ComboSetType& comboSetFull = storedComboSets_[ nbValsFull.size() ];
      //std::cout << "Neigh " << nbValsFull.size() << " combinations " << comboSetFull.size() << std::endl;

      if( comboSetFull.empty() )
      {
        // create combination set
        LimiterUtilityType::buildComboSet( nbValsFull.size(), comboSetFull );
      }

      // reset values
      localGradients_.clear();
      comboVec_.clear();

      // level is only needed for Cartesian grids to access the matrix caches
      const int matrixCacheLevel = ( flags.cartesian ) ? entity.level() : 0 ;
      assert( matrixCacheLevel < (int) matrixCacheVec_.size() );
      MatrixCacheType& matrixCache = matrixCacheVec_[ matrixCacheLevel ];

      // LimiterUtilityType linear functions, stored in localGradients_ and comboVec_
      LimiterUtilityType::
        calculateLinearFunctions( comboSet, geomType, flags,
                                  barys_, nbVals_,
                                  matrixCache,
                                  localGradients_,
                                  comboVec_);

      // functor for checking the physicallity of reconstructed functions
      CheckPhysical checkPhysical( geo, values_[ entityIndex ], entityCenter );

      if( recompute )
      {
        // Limiting
        LimiterUtilityType::limitFunctions( limiterFunction_, checkPhysical,
                                            comboVec_, barysFull, nbValsFull, localGradients_, factor );
        // std::cout << "nbvals: " << nbVals_.size() << std::endl;
      }

      assert( gradients_.size() > 0 );
      // take maximum of limited functions and store in gradients
      LimiterUtilityType::getMaxFunction(localGradients_, gradient,
                                         factor, number, factor );

      /*
      if( ! recompute )
      {
        for( int r=0; r<RangeType::dimension; ++r )
          gradients_[ entityIndex ][ r ] *= factor_[ entityIndex ][ r ];
      }
      */
    }

  private:
    LimitedReconstruction( const LimitedReconstruction& );

  protected:
    struct CheckPhysical
    {
      std::vector< DomainType > points_;
      const RangeType& value_;

      CheckPhysical( const Geometry& geometry,
                     const RangeType& entityValue,
                     const DomainType& entityCenter )
        : points_( geometry.corners() ),
          value_( entityValue )
      {
        const int nCorners = points_.size();
        for( int i=0; i<nCorners; ++i )
        {
          const auto corner = geometry.corner( i );
          for( int d=0; d<dimDomain; ++d )
          {
            points_[ i ][ d ]  = corner[ d ];
          }
          points_[ i ] -= entityCenter;
        }
      }

      bool operator() ( const int r, const DomainType& gradient ) const
      {
        std::cout << " Solution is NOT checked to be physical in operator()!! "  << std::endl;
        return true;
        /*
        const int nCorners = points_.size();
        for( int i=0; i<nCorners; ++i )
        {
          const double functionValue = value_[ r ] + ( points_[ i ] * gradient );
          if( ! model_.physical( functionValue ) )
            return false ;
        }
        return true;
        */
      }
    };

    struct EvalAverage {
      //static const int dimRange = RangeType::dimension;
      const SolutionVector& U_;
      const DofMapper& dofMapper_;
      EvalAverage( const SolutionVector& U, const DofMapper& dofMapper )
        : U_( U )
        , dofMapper_( dofMapper )
      {}

      bool evaluate( const EntityType& entity, RangeType& value ) const
      {
        const RangeType& uEn = U_[ dofMapper_.index( entity ) ];
        for( int i=0; i<dimRange; ++i )
        {
          value[ i ] = uEn[ i ];
        }
        return true;
      }

      template <class IntersectionType, class IntersectionGeometry>
      bool boundaryValue( const EntityType& entity,
                          const IntersectionType& intersection,
                          const IntersectionGeometry& interGeo,
                          const DomainType& globalPoint,
                          const RangeType& entityValue,
                          RangeType& neighborValue ) const
      {
        return false;
      }
    };

    class LocalFunction
    {
      //static const int dimRange = RangeType::dimension;
      const Geometry       geometry_;
      const DomainType&    center_;
      const RangeType     value_;
      const GradientEvalType&  gradient_;
    public:
      LocalFunction( const EntityType& entity, const DomainType& center,
                     const RangeType& value, const GradientEvalType& gradient )
        : geometry_( entity.geometry() ),
          center_( center ),
          value_( value ),
          gradient_( gradient )
      {}

      const Geometry& geometry () const { return geometry_; }

      template <class LocalPoint>
      void evaluate( const LocalPoint& local, RangeType& value ) const
      {
        // compute point of evaluation
        DomainType x = geometry_.global( Fem::coordinate( local ) );
        evaluateGlobal( x, value );
      }

      void evaluateGlobal( const DomainType& point, RangeType& value ) const
      {
        // compute point of evaluation
        DomainType x( point );
        x -= center_;

        // set u
        value = value_;

        //std::cout<<"value "<<value<<std::endl;

        for( int i=0; i<dimRange; ++i )
        {
          // f(x) = u + grad * ( x - center )
          value[i] += gradient_[ i ] * x ;

          //std::cout<<"Gradient[i] " << i <<gradient_[ i ]<<std::endl;
          //std::cout<<"Gradient[i] "<<gradient_[ i ].derivative[0]<<std::endl;
          //std::cout<<"Gradient[i] "<<gradient_[ i ].derivative[1]<<std::endl;

          //std::cout<<"value[i] "<<value[ i ]<<std::endl;

        }
      }

//      void evaluateGlobal( const DomainType& point, EvalDimVector& value ) const
//      {
//        // compute point of evaluation
//        DomainType x( point );
//        x -= center_;

//        for( int i=0; i<dimRange; ++i )
//        {
//          //int idx = i * (numDeri+1);
//          // set u
//          //value = value_;

//          // f(x) = u + grad * ( x - center )
//          Field val = value_[ i ] + ( gradient_[ i ] * x );
//          //++ idx ;
//          value[ i ].setValue( val );
//          //for( int d=0; d<numDeri; ++d, ++idx )
//          //{
//          //  Field val = value_[ idx ]; // + ( gradient_[ idx ] * x );
//          //  value[ i ].setDerivative( d, val );
//          //}
//          //std::cout<<"Gradient[i] "<<gradient_[ i ]<<std::endl;
//        }
//      }
    };

    class SimpleLocalFunction
    {
      //static const int dimRange = RangeType::dimension;
      const DomainType&    center_;
      const Evaluation&    value_;
      const SingleGradientEvalType& gradient_;
    public:
      SimpleLocalFunction( const DomainType& center,
                           const Evaluation& value,
                           const SingleGradientEvalType& gradient)
        : center_( center ),
          value_( value ),
          gradient_( gradient )
      {}

      template <class Result>
      void evaluateGlobal( const DomainType& point, Result& value ) const
      {
        // compute point of evaluation
        DomainType x( point );
        x -= center_;

        // set u
        value = value_;
        // f(x) = u + grad * ( x - center )
        value += gradient_ * x ;
      }

    };


    void applyLocal( const EntityType& entity, const SolutionVector& U, const bool recompute )
    {
      // helper class for evaluation of average value of discrete function
      EvalAverage average( U, dofMapper_ );

      const unsigned int entityIndex = dofMapper_.index( entity );

      // get geometry
      const Geometry& geo = entity.geometry();

      // cache geometry type
      const GeometryType geomType = entity.type();

      // get bary center of element
      const DomainType& entityCenter = centers_[ entityIndex ];

      // evaluate mean value on current entity
      average.evaluate( entity, values_[ entityIndex ] );

      // boundary is true if boundary segment was found
      // nonconforming is true if entity has at least one non-conforming intersections
      // cartesian is true if the grid is cartesian and no nonconforming refinement present
      typename LimiterUtilityType::Flags flags( cartesianGrid_, true );


      std::vector< DomainType > barysFull;
      std::vector< RangeType >  nbValsFull;

      // setup neighbors barycenter and mean value for all neighbors
      LimiterUtilityType::setupNeighborValues( gridPart_, entity, average, entityCenter, values_[ entityIndex ],
                                               centers_, StructuredGrid, flags, barys_, nbVals_, barysFull, nbValsFull);

      // mark entity as finished, even if not limited everything necessary was done
      //assert( entityIndex  < visited_.size() );
      //visited_[ entityIndex ] = true ;

      ComboSetType& comboSet = storedComboSets_[ nbVals_.size() ];
      if( comboSet.empty() )
      {
        // create combination set
        LimiterUtilityType::buildComboSet( nbVals_.size(), comboSet );
      }
      ComboSetType& comboSetFull = storedComboSets_[ nbValsFull.size() ];
      //std::cout << "Neigh " << nbValsFull.size() << " combinations " << comboSetFull.size() << std::endl;

      if( comboSetFull.empty() )
      {
        // create combination set
        LimiterUtilityType::buildComboSet( nbValsFull.size(), comboSetFull );
      }

      // reset values
      localGradients_.clear();
      comboVec_.clear();

      // level is only needed for Cartesian grids to access the matrix caches
      const int matrixCacheLevel = ( flags.cartesian ) ? entity.level() : 0 ;
      assert( matrixCacheLevel < (int) matrixCacheVec_.size() );
      MatrixCacheType& matrixCache = matrixCacheVec_[ matrixCacheLevel ];

      // calculate linear functions, stored in localGradients_ and comboVec_
      LimiterUtilityType::
        calculateLinearFunctions( comboSet, geomType, flags,
                                  barys_, nbVals_,
                                  matrixCache,
                                  localGradients_,
                                  comboVec_);

      // functor for checking the physicallity of reconstructed functions
      CheckPhysical checkPhysical( geo, values_[ entityIndex ], entityCenter );

      std::vector< RangeType > factors;
      if( recompute )
      {
        // Limiting
        LimiterUtilityType::limitFunctions( limiterFunction_, checkPhysical,
                                            comboVec_, barysFull, nbValsFull, localGradients_, factors );
        // std::cout << "nbvals: " << nbVals_.size() << std::endl;
      }

      assert( gradients_.size() > 0 );
      // take maximum of limited functions and store in gradients
      LimiterUtilityType::getMaxFunction(localGradients_, gradients_[ entityIndex ],
                                         factor_[ entityIndex ], numbers_[ entityIndex ], factors );

      /*
      if( ! recompute )
      {
        for( int r=0; r<RangeType::dimension; ++r )
          gradients_[ entityIndex ][ r ] *= factor_[ entityIndex ][ r ];
      }
      */
    }



  public:
    typedef LocalFunction  LocalFunctionType;

    //! return local reconstruction
    LocalFunctionType localFunction( const EntityType& entity )
    {
      const unsigned int entityIndex = gridPart_.indexSet().index( entity );
      return LocalFunctionType( entity, centers_[ entityIndex ],
                                values_[ entityIndex ], gradients_[ entityIndex ] );
    }

    //! return local reconstruction
    const LocalFunctionType localFunction( const EntityType& entity ) const
    {
      const unsigned int entityIndex = gridPart_.indexSet().index( entity );
      return LocalFunctionType( entity, centers_[ entityIndex ],
                                values_[ entityIndex ], gradients_[ entityIndex ] );
    }

    //! return local reconstruction
    const LocalFunctionType localFunction( const EntityType& entity, const size_t dofIdx ) const
    {
      return LocalFunctionType( entity, centers_[ dofIdx ],
                                values_[ dofIdx ], gradients_[ dofIdx ] );
    }

    //! return local reconstruction
    const LocalFunctionType localFunctionLocal( const EntityType& entity,
                                                const size_t upstreamDofIdx,
                                                const size_t focusIdx,
                                                const size_t newFocusIdx ) const
    {
      const unsigned int entityIndex = gridPart_.indexSet().index( entity );

      size_t gradientIndex = elemCounter_[ entityIndex ] + newFocusIdx;
      // when upstream is not focus the derivatives need to be removed
      RangeType value = values_[ entityIndex ];
      if( upstreamDofIdx != focusIdx )
      {
        for( int i=0; i<RangeType::dimension; ++i )
        {
          value[ i ] = Toolbox::createConstant( Toolbox::value(value[ i ]));
        }
      }
      if (newFocusIdx < 0) {
          GradientEvalType gradient = gradients_[ elemCounter_[ entityIndex ] ];
          for( int i=0; i<RangeType::dimension; ++i )
              for( int j=0; j<dimDomain; ++j )
          {
              gradient[ i ][ j ] = Toolbox::createConstant( Toolbox::value(gradient[ i ][ j ]));
          }

          return LocalFunctionType( entity, centers_[ entityIndex ],
                                    value, gradient );
      }
      return LocalFunctionType( entity, centers_[ entityIndex ],
                                value, gradients_[ gradientIndex ] );
    }

    //! return local reconstruction
    template <class Result>
    void evaluate( const size_t entityIndex,
                   const int upstreamDofIdx,
                   const int focusIdx,
                   const int newFocusIdx,
                   const int component,
                   const DomainType& point,
                   Result& result ) const
    {
      // when upstream is not focus the derivatives need to be removed
      Evaluation value = values_[ entityIndex ][ component ];
      if( upstreamDofIdx != focusIdx )
      {
        value = Toolbox::createConstant( Toolbox::value( value ) );
      }

      if (newFocusIdx < 0)
      {
          SingleGradientEvalType gradient = gradients_[ elemCounter_[ entityIndex ] ][ component ];
          for( int j=0; j<dimDomain; ++j )
          {
            gradient[ j ] = Toolbox::createConstant( Toolbox::value(gradient[ j ]));
          }

          SimpleLocalFunction( centers_[ entityIndex ], value, gradient ).evaluateGlobal( point, result );
      }
      else
      {
        size_t gradientIndex = elemCounter_[ entityIndex ] + newFocusIdx;
        SimpleLocalFunction( centers_[ entityIndex ], value, gradients_[ gradientIndex ][ component] ).evaluateGlobal( point, result );
      }
    }

  protected:
    bool checkPhysical( const LocalFunctionType& localRecon ) const
    {
      /*const auto& geometry = localRecon.geometry();
      const int corners = geometry.corners();
      RangeType value;
      for( int i = 0; i < corners; ++ i )
      {
        localRecon.evaluateGlobal( geometry.corner( i ), value );
        if( ! model_.problem().physical( value ) )
          return false ;
      }*/
      std::cout << " Solution is NOT checked to be physical in checkPhysical func!! "  << std::endl;
      return true;

    }

    const GridPartType&     gridPart_;
    const DofMapper&        dofMapper_;

    const LSRecon lsRecon_;
    const LSReconLocal lsReconLocal_;
    const LinearProgramming linProg_;
    std::vector< LinearProgrammingLocal > linProgLocal_;

    //const QuadraticProgramming quadProg_;

    LimiterFunctionType limiterFunction_;

    mutable ComboSetMapType storedComboSets_;

    mutable std::vector< GradientEvalType > localGradients_;
    mutable std::vector< CheckType >    comboVec_;

    mutable std::vector< DomainType > barys_;
    mutable std::vector< RangeType >  nbVals_;
    mutable std::vector< MatrixCacheType > matrixCacheVec_;

    std::vector< RangeType >    values_;
    std::vector< DomainType >   centers_;
    std::vector< GradientEvalType > gradients_;
    std::vector< std::vector< int > > numbers_;
    std::vector< RangeType >  factor_;
    std::vector< size_t> elemCounter_;

    double time_;

    const bool cartesianGrid_;
  };

} // end namespace Fem

} // end namespace Dune

#endif // DUNE_FEM_LIMITER_HH
