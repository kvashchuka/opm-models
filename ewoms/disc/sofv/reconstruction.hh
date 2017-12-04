#ifndef EWOMS_RECONSTRUCTION_HH
#define EWOMS_RECONSTRUCTION_HH

#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/io/parameter.hh>

#include "limiterutility.hh"
#include "limitermodel.hh"

namespace Ewoms
{
//namespace Fem
//{

  /**
   * \brief Limited reconstruction.
   *
   * \ingroup PassBased
   */
  template <class TypeTag> //Model, class DiscreteFunction>
  class LimitedReconstruction
  {
  public:
//    typedef DiscreteFunction  DiscreteFunctionType;
//    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
//      DiscreteFunctionSpaceType;
//
//    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
//    typedef typename DiscreteFunctionSpaceType :: RangeType  RangeType;
      typedef typename GET_PROP_TYPE(TypeTag, Grid) GridType;
      typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPartType;
      typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
      typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
      //typedef typename Dune::Fem::BasicGridFunctionAdapter::GridPartType GridPartType;
      typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
      enum { dimDomain = GridType::dimensionworld };
      enum { dimRange  = PrimaryVariables::dimension };

      typedef typename GET_PROP_TYPE(TypeTag, Scalar) DomainFieldType; //should be Scalar
      typedef double Field; //can be Scalar
      typedef Field RangeFieldType;
      typedef Field FieldType;

      typedef Dune::FieldVector<typename GridType::ctype, dimDomain> DomainType;
      typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

    typedef LimiterUtility< TypeTag >      LimiterUtilityType;
    typedef typename LimiterUtilityType :: GradientType      GradientType;
    typedef typename LimiterUtilityType :: CheckType         CheckType;
    typedef typename LimiterUtilityType :: ComboSetType      ComboSetType;
    typedef typename LimiterUtilityType :: KeyType           KeyType;

    typedef LimiterModel< TypeTag > Model;
    typedef typename Model :: Traits :: LimiterFunctionType  LimiterFunctionType;

   // typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
   // typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
   // typedef typename GridPartType :: GridType                        GridType;

   // typedef typename DiscreteFunctionSpaceType :: EntityType         EntityType;
   typedef typename GridView::template Codim<0>::Entity EntityType;
    typedef typename EntityType :: Geometry   Geometry;
      //typedef typename EntityType :: GeometryType GeometryType;

    typedef typename LimiterUtilityType::MatrixStorage MatrixCacheEntry;
    typedef std::map< KeyType, MatrixCacheEntry > MatrixCacheType;

    static const bool StructuredGrid     = Dune::Fem::GridPartCapabilities::isCartesian< GridPartType >::v;

    //! type of cartesian grid checker
    typedef Dune::Fem::CheckCartesian< GridPartType >  CheckCartesianType;

  public:
    LimitedReconstruction( const GridPartType &gridPart, const Model& limiterModel, const DofMapper& dofMapper)
      : //space_( space )
        gridPart_( gridPart )
      , model_( limiterModel )
      , limiterFunction_()
      , conformingComboSet_()
      , comboSet_()
      , matrixCacheVec_( gridPart_.grid().maxLevel() + 1 )
      , cartesianGrid_( CheckCartesianType::check( gridPart_ ) )
    {}

    //! calculate internal reconstruction based on mean values of U
    void update( const SolutionVector & U, const DofMapper& dofMapper )
    {
      const size_t size = gridPart_.indexSet().size( 0 ) ;
      // resize gradient vector
      values_.resize( size );
      gradients_.resize( size );

      const auto end = gridPart_.template end< 0, Dune::Interior_Partition >();
      for( auto it = gridPart_.template begin< 0, Dune::Interior_Partition >(); it != end; ++it )
      {
        applyLocal( *it, U, dofMapper );
      }
    }

  private:
    LimitedReconstruction( const LimitedReconstruction& );

  protected:
    struct CheckPhysical
    {
      std::vector< DomainType > points_;
      const Model& model_;
      const RangeType& value_;

      CheckPhysical( const Model& model,
                     const Geometry& geometry,
                     const RangeType& entityValue,
                     const DomainType& entityCenter )
        : points_( geometry.corners() ),
          model_( model ),
          value_( entityValue )
      {
        const int nCorners = points_.size();
        for( int i=0; i<nCorners; ++i )
        {
          points_[ i ]  = geometry.corner( i );
          points_[ i ] -= entityCenter;
        }
      }

      bool operator() ( const int r, const DomainType& gradient ) const
      {
        const int nCorners = points_.size();
        for( int i=0; i<nCorners; ++i )
        {
          const double functionValue = value_[ r ] + ( points_[ i ] * gradient );
          if( ! model_.physical( functionValue ) )
            return false ;
        }
        return true;
      }
    };

    struct EvalAverage {
      enum { dimRange  = PrimaryVariables::dimension };
      const SolutionVector& U_;
      const DofMapper dofMapper_;
     // typedef typename SolutionVector :: LocalFunctionType LocalFunctionType;
      EvalAverage( const SolutionVector & U, const DofMapper& dofMapper )
        : U_( U )
          , dofMapper_ (dofMapper)
      {}

      bool evaluate( const EntityType& entity, RangeType& value ) const
      {
        //const LocalFunctionType& uEn = U_.localFunction( entity );
        unsigned dofIdx = dofMapper_.index(entity);
        for( int i=0; i<dimRange; ++i )
        {
          value[ i ] = U_[ dofIdx][i];
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
        enum { dimRange  = PrimaryVariables::dimension };
      const Geometry       geometry_;
      const DomainType     center_;
      const RangeType&     value_;
      const GradientType&  gradient_;
    public:
      LocalFunction( const EntityType& entity, const RangeType& value, const GradientType& gradient )
        : geometry_( entity.geometry() ),
          center_( geometry_.center() ),
          value_( value ),
          gradient_( gradient )
      {}

      const Geometry& geometry () const { return geometry_; }

      template <class LocalPoint>
      void evaluate( const LocalPoint& local, RangeType& value ) const
      {
        // compute point of evaluation
        DomainType x = geometry_.global( Dune::Fem::coordinate( local ) );
        evaluateGlobal( x, value );
      }

      void evaluateGlobal( const DomainType& point, RangeType& value ) const
      {
        // compute point of evaluation
        DomainType x( point );
        x -= center_;

        // set u
        value = value_;

        for( int i=0; i<dimRange; ++i )
        {
          // f(x) = u + grad * ( x - center )
          value[ i ] += gradient_[ i ] * x;
          //std::cout<<"Gradient[i] "<<gradient_[ i ]<<std::endl;
        }
      }
    };


    void applyLocal( const EntityType& entity, const SolutionVector & U, const DofMapper & dofMapper)
    {
      // helper class for evaluation of average value of discrete function
      EvalAverage average( U, dofMapper );

     // const unsigned int entityIndex = U.space().gridPart().indexSet().index( entity );
      unsigned entityIndex = dofMapper.index(entity); //dofIdx usually

      // get geometry
      const Geometry& geo = entity.geometry();

      // cache geometry type
     // const GeometryType geomType = geo.type();

      // get bary center of element
      const DomainType entityCenter = geo.center();

      // evaluate mean value on current entity
      average.evaluate( entity, values_[ entityIndex ] );

      // boundary is true if boundary segment was found
      // nonconforming is true if entity has at least one non-conforming intersections
      // cartesian is true if the grid is cartesian and no nonconforming refinement present
      typename LimiterUtilityType::Flags flags( cartesianGrid_, true );

      // setup neighbors barycenter and mean value for all neighbors
      LimiterUtilityType::setupNeighborValues( gridPart_, entity, average, entityCenter, values_[ entityIndex ],
                                               StructuredGrid, flags, barys_, nbVals_ );

      // mark entity as finished, even if not limited everything necessary was done
      //assert( entityIndex  < visited_.size() );
      //visited_[ entityIndex ] = true ;

      // create combination set
      //TODO multipleGeometryTypes uses space_, which we don't have. Make sure to check geometry types in some other way
      const ComboSetType& comboSet = LimiterUtilityType::
        setupComboSet( nbVals_.size(), flags.nonConforming, /*space_.multipleGeometryTypes()*/ true,
                       conformingComboSet_, comboSet_ );

      // reset values
      localGradients_.clear();
      comboVec_.clear();

      // level is only needed for Cartesian grids to access the matrix caches
      const int matrixCacheLevel = ( flags.cartesian ) ? entity.level() : 0 ;
      assert( matrixCacheLevel < (int) matrixCacheVec_.size() );
      MatrixCacheType& matrixCache = matrixCacheVec_[ matrixCacheLevel ];

      // calculate linear functions, stored in localGradients_ and comboVec_
      LimiterUtilityType::calculateLinearFunctions( comboSet, /*geomType,*/ flags,
                                barys_, nbVals_,
                                matrixCache,
                                localGradients_,
                                comboVec_ );

      // functor for checking the physicallity of reconstructed functions
      CheckPhysical checkPhysical( model_, geo, values_[ entityIndex ], entityCenter );

      // Limiting
      LimiterUtilityType::limitFunctions( limiterFunction_, checkPhysical,
                                          comboVec_, barys_, nbVals_, localGradients_);

      // take maximum of limited functions and store in gradients
      LimiterUtilityType::getMaxFunction(localGradients_, gradients_[ entityIndex ] );

      /*
      // if reconstructed solution is not physical, drop it
      if( ! checkPhysical( localFunction( entity ) ) )
          gradients_[ entityIndex ] = GradientType( DomainType( 0 ) );
          */
    }

  public:
    typedef LocalFunction  LocalFunctionType;

    //! return local reconstruction
    LocalFunctionType localFunction( const EntityType& entity )
    {
      const unsigned int entityIndex = gridPart_.indexSet().index( entity );
      return LocalFunctionType( entity, values_[ entityIndex ], gradients_[ entityIndex ] );
    }

    //! return local reconstruction
    const LocalFunctionType localFunction( const EntityType& entity ) const
    {
      const unsigned int entityIndex = gridPart_.indexSet().index( entity );
      return LocalFunctionType( entity, values_[ entityIndex ], gradients_[ entityIndex ] );
    }

      const GridPartType GridPart() const
      {
          return gridPart_;
      };

      const GridPartType& gridPart_;

  protected:
    bool checkPhysical( const LocalFunctionType& localRecon ) const
    {
      const auto& geometry = localRecon.geometry();
      const int corners = geometry.corners();
      RangeType value;
      for( int i = 0; i < corners; ++ i )
      {
        localRecon.evaluateGlobal( geometry.corner( i ), value );
        if( value[ 0 ] > 1.0 || value[ 0 ] < 0.0 )
          return false ;
      }
      return true;
    }

   // const DiscreteFunctionSpaceType& space_;
  // const GridPartType& gridPart_;
    const Model& model_;

    LimiterFunctionType limiterFunction_;

    mutable ComboSetType conformingComboSet_;
    mutable ComboSetType comboSet_;

    mutable std::vector< GradientType > localGradients_;
    mutable std::vector< CheckType >    comboVec_;

    mutable std::vector< DomainType > barys_;
    mutable std::vector< RangeType >  nbVals_;
    mutable std::vector< MatrixCacheType > matrixCacheVec_;

    std::vector< RangeType >    values_;
    std::vector< GradientType > gradients_;

    const bool cartesianGrid_;
  };

//} // end namespace Fem

} // end namespace Dune

#endif // DUNE_FEM_LIMITER_HH
