#ifndef DUNE_OPTIM_STD_ITERATOR_HH
#define DUNE_OPTIM_STD_ITERATOR_HH

#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

namespace Dune
{

  namespace Std
  {

    // iterator_t
    // ----------

    namespace __iterator_t
    {

      using std::begin;

      template< class I >
      decltype( begin( std::declval< I & >() ) ) begin ( I & );

      template< class I >
      decltype( begin( std::declval< I & >() ) ) begin ( std::reference_wrapper< I > );


      template< class I >
      decltype( begin( std::declval< const I & >() ) ) cbegin ( I & );

      template< class I >
      decltype( begin( std::declval< const I & >() ) ) cbegin ( std::reference_wrapper< I > );


      template< class I >
      decltype( std::declval< I & >().rbegin() ) rbegin ( I & );

      template< class I >
      decltype( std::declval< I & >().rbegin() ) rbegin ( std::reference_wrapper< I > );


      template< class I >
      decltype( std::declval< const I & >().rbegin() ) crbegin ( I & );

      template< class I >
      decltype( std::declval< const I & >().rbegin() ) crbegin ( std::reference_wrapper< I > );

    } // namespace __iterator_t

    template< class I >
    using iterator_t = decltype( __iterator_t::begin( std::declval< I >() ) );

    template< class I >
    using const_iterator_t = decltype( __iterator_t::cbegin( std::declval< I >() ) );

    template< class I >
    using reverse_iterator_t = decltype( __iterator_t::rbegin( std::declval< I >() ) );

    template< class I >
    using const_reverse_iterator_t = decltype( __iterator_t::crbegin( std::declval< I >() ) );

  } // namespace Std

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_STD_ITERATOR_HH
