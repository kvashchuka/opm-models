#ifndef LIMITER_MODEL_HH
#define LIMITER_MODEL_HH

#include <cmath>
#include <dune/common/fvector.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/fem/io/parameter.hh>


namespace Ewoms
{

//namespace Fem
//{

template <class TypeTag> //, class ProblemType>
class LimiterModel;

template <class TypeTag>
class LimiterModelTraits
{
public:
  //typedef Problem  ProblemType;
  typedef double Field; //can be Scalar
  typedef Field RangeFieldType;
  typedef Field FieldType;
  //typedef GridPart GridPartType;
  //typedef typename GridPart::GridType GridType;
  typedef typename GET_PROP_TYPE(TypeTag, Grid) GridType;
  typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPart;
  typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange  = PrimaryVariables::dimension };
  enum { dimGradRange = dimRange * dimDomain };
  typedef FieldVector<typename GridType::ctype, dimDomain> DomainType;
  typedef FieldVector<typename GridType::ctype, dimDomain-1> FaceDomainType;
  typedef FieldVector<RangeFieldType, dimRange> RangeType;
  typedef FieldVector<RangeFieldType, dimGradRange> GradientType;
  typedef typename GridPart::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType :: Intersection IntersectionType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef MinModLimiter< FieldType > LimiterFunctionType ;

//   typedef SuperBeeLimiter< FieldType > LimiterFunctionType ;
//   typedef VanLeerLimiter< FieldType > LimiterFunctionType ;

  typedef LimiterModel< TypeTag > ModelType;
};

/** \brief LimiterModel for second order extension by RK */
template <class TypeTag>
class LimiterModel
{
public:
  //typedef Problem ProblemType;
  typedef LimiterModelTraits< TypeTag > Traits;
  enum { dimDomain = Traits::dimensionworld };
  enum { dimGrid = Traits::dimension };
  enum { dimRange = Traits :: dimRange };
  typedef typename Traits::FieldType FieldType;
  typedef typename Traits::GridType GridType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::FaceDomainType FaceDomainType;
    typedef typename Traits::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType :: Intersection IntersectionType;
  typedef typename Traits::EntityType EntityType;

protected:
  typedef LimiterModel<TypeTag> ThisType;

public:
  LimiterModel( /*const ProblemType& problem*/ )
    //:problem_(problem)
  {}


  // has flux
  bool hasFlux() const
  {
    return true;
  }

  // has source term
  bool hasSource() const
  {
    return false;
  }

  //! needed for limitation to determine inflow and outflow intersections
    // needs to be rewritten, probably in FilterVelocity from DarcyFlux
  inline void velocity( const EntityType& entity,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        DomainType& velocity) const
  {
    velocity = 0.0;
    //velocity = problem_.velocity();
    //problem_.velocity( x, time, u, velocity );
  }


  /** \brief adjust average values, e.g. transform to primitive or something similar */
  void adjustAverageValue( const EntityType& entity,
                           const DomainType& xLocal,
                           RangeType& u ) const
  {
  }

  inline void reAdjust( RangeType &u ) const
  {
  }

  // we have physical check for this model
  bool hasPhysical() const
  {
    return true;
  }



  // calculate jump between left and right value
  inline bool physical( const double value ) const
  {
    return true; //(value >= 0.0) && (value <= 1.0);
  }

  // calculate jump between left and right value
  inline bool physical(  const EntityType& entity,
                          const DomainType& xLocal,
                          const RangeType& u ) const
  {
    bool result = true;
    //for( int i=0; i<dimRange; ++ i )
    //  result = result && physical( u[ i ] );
    return result;
  }

  inline void jump(
        const IntersectionType& it,
        const double time,
        const FaceDomainType& x,
        const RangeType& uLeft,
        const RangeType& uRight,
        RangeType &shockIndicator ) const
  {
    shockIndicator[0] = (uLeft[0] - uRight[0])/(0.5 * (uLeft[0] + uRight[0]));
  }

  // calculate jump between left and right value
  inline void adaptationIndicator(
                   const IntersectionType& it,
                   const double time,
                   const FaceDomainType& x,
                   const RangeType& uLeft,
                   const RangeType& uRight,
                   RangeType& indicator) const
  {
    // use jump of values
    jump( it, time, x, uLeft, uRight, indicator );
  }

  inline bool hasBoundaryValue(const IntersectionType& it,
                               const double time,
                               const FaceDomainType& x) const
  {
    return false;
  }

  template <class LocalEvaluation>
  inline void boundaryValue(const LocalEvaluation& local,
                            const RangeType& uLeft,
                            RangeType& uRight) const
  {
    // copy value (for limitation)
    uRight = uLeft;
  }

protected:
    //probably darcyflux problme here
    //const ProblemType& problem_;
};

//} // end namespace Fem
} // end namespace Dune

#endif
