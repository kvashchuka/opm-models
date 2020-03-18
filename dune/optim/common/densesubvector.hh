#ifndef DUNE_OPTIM_COMMON_DENSESUBVECTOR_HH
#define DUNE_OPTIM_COMMON_DENSESUBVECTOR_HH

#include <functional>
#include <utility>

#include <dune/common/densevector.hh>

#include <dune/optim/common/integerrange.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class V, class Indices >
  class DenseSubVector;



  // DenseMatVecTraits for DenseSubVector
  // ------------------------------------

  template< class V, class Indices >
  struct DenseMatVecTraits< DenseSubVector< V, Indices > >
  {
    typedef DenseSubVector< V, Indices > derived_type;
    typedef V container_type;

    typedef typename V::value_type value_type;
    typedef typename std::decay< decltype( std::ref( std::declval< Indices >() ).get() ) >::type::size_type size_type;
  };



  // FieldTraits for DenseSubVector
  // ------------------------------

  template< class V, class Indices >
  struct FieldTraits< DenseSubVector< V, Indices > >
  {
    typedef typename FieldTraits< typename V::value_type >::field_type field_type;
    typedef typename FieldTraits< typename V::value_type >::real_type real_type;
  };


  // DenseSubVector
  // --------------

  template< class V, class Indices >
  class DenseSubVector
    : public DenseVector< DenseSubVector< V, Indices > >
  {
    typedef DenseSubVector< V, Indices > This;
    typedef DenseVector< DenseSubVector< V, Indices > > Base;

  public:
    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;

    DenseSubVector ( V &v, Indices indices ) : v_( v ), indices_( std::move( indices ) ) {}

    DenseSubVector ( const This &other ) : v_( other.v_ ), indices_( other.inices_ ) {}
    DenseSubVector ( This &&other ) : v_( other.v_ ), indices_( std::move( other.indices_ ) ) {}

    using Base::operator=;

    size_type size () const { return std::ref( indices_ ).get().size(); }

    const value_type &operator[] ( size_type i ) const { return v_[ std::ref( indices_ ).get()[ i ] ]; }
    value_type &operator[] ( size_type i ) { return v_[ std::ref( indices_ ).get()[ i ] ]; }

  private:
    typename DenseMatVecTraits< This >::container_type &v_;
    Indices indices_;
  };



  // denseSubVector
  // --------------

  template< class V, class Indices >
  inline static DenseSubVector< V, Indices > denseSubVector ( V &v, Indices indices )
  {
    return DenseSubVector< V, Indices >( v, std::move( indices ) );
  }

  template< class V >
  inline static auto denseSubVector ( V &v, std::size_t size, std::size_t offset = 0 ) noexcept
    -> decltype( denseSubVector( v, integerRange( size, offset ) ) )
  {
    return denseSubVector( v, integerRange( size, offset ) );
  }

  template< class V >
  inline static auto denseSubVector ( V &v, std::size_t size, std::size_t offset, std::size_t stride ) noexcept
    -> decltype( denseSubVector( v, integerRange( size, offset, stride ) ) )
  {
    return denseSubVector( v, integerRange( size, offset, stride ) );
  }

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_COMMON_DENSESUBVECTOR_HH
