#ifndef DUNE_COMMON_GETREFERENCE_HH
#define DUNE_COMMON_GETREFERENCE_HH

#include <functional>
#include <utility>

namespace Dune
{

  // getReference
  // ------------

  template< class T >
  inline static T &getReference ( T &value )
  {
    return value;
  }

  template< class T >
  inline static T &getReference ( const std::reference_wrapper< T > &value )
  {
    return value;
  }

  template< class T >
  inline static T &getReference ( std::reference_wrapper< T > &value )
  {
    return value;
  }



  // getReferredType
  // ---------------

  template< class T >
  using getReferredType = typename std::remove_reference< decltype( getReference( std::declval< T >() ) ) >::type;

} // namespace Dune

#endif // #ifndef DUNE_COMMON_GETREFERENCE_HH
