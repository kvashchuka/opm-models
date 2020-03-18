#ifndef DUNE_OPTIM_STD_SUBARRAY_HH
#define DUNE_OPTIM_STD_SUBARRAY_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/optim/std/iterator.hh>

namespace Dune
{

  namespace Std
  {

    // subarray_iterator
    // -----------------

    template< class A, class I, class R, class P >
    struct subarray_iterator
    {
      typedef typename std::remove_const< A >::type array_type;

      friend class subarray_iterator< array_type, I, typename array_type::reference, typename array_type::pointer >;
      friend class subarray_iterator< const array_type, I, typename array_type::const_reference, typename array_type::const_pointer >;

      subarray_iterator ( A &array, I indices_iterator )
        : array_( &array ), indices_iterator_( indices_iterator )
      {}

      subarray_iterator ( const subarray_iterator< array_type, I, typename array_type::reference, typename array_type::pointer > &other )
        : array_( other.array_ ), indices_iterator_( other.indices_iterator_ )
      {}

      subarray_iterator ( const subarray_iterator< const array_type, I, typename array_type::const_reference, typename array_type::const_pointer > &other )
        : array_( other.array_ ), indices_iterator_( other.indices_iterator_ )
      {}

      subarray_iterator &operator= ( const subarray_iterator< array_type, I, typename array_type::reference, typename array_type::pointer > &other )
      {
        array_ = other.array_;
        indices_iterator_ = other.indices_iterator_;
        return *this;
      }

      subarray_iterator &operator= ( const subarray_iterator< const array_type, I, typename array_type::const_reference, typename array_type::const_pointer > &other )
      {
        array_ = other.array_;
        indices_iterator_ = other.indices_iterator_;
        return *this;
      }

      bool operator== ( const subarray_iterator< array_type, I, typename array_type::reference, typename array_type::pointer > &other ) const
      {
        return (array_ == other.array_) && (indices_iterator_ == other.indices_iterator_);
      }

      bool operator== ( const subarray_iterator< const array_type, I, typename array_type::const_reference, typename array_type::const_pointer > &other ) const
      {
        return (array_ == other.array_) && (indices_iterator_ == other.indices_iterator_);
      }

      bool operator!= ( const subarray_iterator< array_type, I, typename array_type::reference, typename array_type::pointer > &other ) const
      {
        return (array_ != other.array_) || (indices_iterator_ != other.indices_iterator_);
      }

      bool operator!= ( const subarray_iterator< const array_type, I, typename array_type::const_reference, typename array_type::const_pointer > &other ) const
      {
        return (array_ != other.array_) || (indices_iterator_ != other.indices_iterator_);
      }

      subarray_iterator &operator++ () { ++indices_iterator_; return *this; }
      subarray_iterator operator++ ( int ) { subarray_iterator other( *this ); ++(*this); return other; }

      P operator-> () const { return &(*array_)[ *indices_iterator_ ]; }
      R operator* () const { return (*array_)[ *indices_iterator_ ]; }

    private:
      A *array_;
      I indices_iterator_;
    };



    // subarray
    // --------

    template< class A, class I >
    struct subarray
    {
      typedef A array_type;
      typedef std::decay_t< decltype( std::ref( std::declval< const I & >() ).get() ) > indices_type;

      typedef typename array_type::value_type value_type;

      typedef typename array_type::reference reference;
      typedef typename array_type::const_reference const_reference;

      typedef typename array_type::pointer pointer;
      typedef typename array_type::const_pointer const_pointer;

      typedef typename indices_type::size_type size_type;

      typedef subarray_iterator< array_type, const_iterator_t< I >, reference, pointer > iterator;
      typedef subarray_iterator< const array_type, const_iterator_t< I >, const_reference, const_pointer > const_iterator;
      typedef subarray_iterator< array_type, const_reverse_iterator_t< I >, reference, pointer > reverse_iterator;
      typedef subarray_iterator< const array_type, const_reverse_iterator_t< I >, const_reference, const_pointer > const_reverse_iterator;

      subarray ( array_type &array, I indices )
        : array_( array ), indices_( std::move( indices ) )
      {}

      // iterator

      const_iterator begin () const { return { array_, indices().begin() }; }
      iterator begin () { return { array_, indices().begin() }; }
      const_iterator end () const { return { array_, indices().end() }; }
      iterator end () { return { array_, indices().end() }; }

      const_reverse_iterator rbegin () const { return { array_, indices().begin() }; }
      reverse_iterator rbegin () { return { array_, indices().begin() }; }
      const_reverse_iterator rend () const { return { array_, indices().end() }; }
      reverse_iterator rend () { return { array_, indices().end() }; }

      const_iterator cbegin () const { return { array_, indices().begin() }; }
      const_iterator cend () const { return { array_, indices().end() }; }

      const_reverse_iterator crbegin () const { return { array_, indices().begin() }; }
      const_reverse_iterator crend () const { return { array_, indices().end() }; }

      // capacity

      bool empty () const noexcept { return indices().emtpy(); }
      size_type size () const noexcept { return indices().size(); }

      // element access

      const_reference operator[] ( size_type n ) const { return array_[ indices()[ n ] ]; }
      reference operator[] ( size_type n ) { return array_[ indices()[ n ] ]; }

      const_reference at ( size_type n ) const { return array_[ indices().at( n ) ]; }
      reference at ( size_type n ) { return array_[ indices().at( n ) ]; }

      const_reference front () const { return array_[ indices().front() ]; }
      reference front () { return array_[ indices().back() ]; }
      const_reference back () const { return array_[ indices().front() ]; }
      reference back () { return array_[ indices().back() ]; }

      const indices_type &indices () const { return std::ref( indices_ ).get(); }
      indices_type &indices () { return std::ref( indices_ ).get(); }

    private:
      array_type &array_;
      I indices_;
    };



    // make_subarray
    // -------------

    template< class A, class I >
    inline static subarray< A, I > make_subarray ( A &array, I indices ) noexcept
    {
      return subarray< A, I >( array, std::move( indices ) );
    }

  } // namespace Std

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_STD_SUBARRAY_HH
