#ifndef DUNE_OPTIM_COMMON_INTEGERRANGE_HH
#define DUNE_OPTIM_COMMON_INTEGERRANGE_HH

#include <type_traits>

namespace Dune
{

  // IntegerRange
  // ------------

  template< class T, class S = T >
  struct IntegerRange;

  template< class T >
  struct IntegerRange< T, T >
  {
    typedef T value_type;
    typedef T size_type;

    explicit IntegerRange ( size_type size, size_type offset = 0, size_type stride = 1 )
      : size_( size ), offset_( offset ), stride_( stride )
    {}

    value_type operator[] ( size_type i ) const { return offset_ + i*stride_; }

    bool empty () const { return size_ == size_type( 0 ); }
    size_type size () const { return size_; }

  private:
    size_type size_, offset_, stride_;
  };

  template< class T, T stride >
  struct IntegerRange< T, std::integral_constant< T, stride > >
  {
    typedef T value_type;
    typedef T size_type;

    explicit IntegerRange ( size_type size, size_type offset = 0 )
      : size_( size ), offset_( offset )
    {}

    operator IntegerRange< T > () const { return IntegerRange< T >( size_, offset_, stride ); }

    value_type operator[] ( size_type i ) const { return offset_ + i*stride; }

    bool empty () const { return size_ == size_type( 0 ); }
    size_type size () const { return size_; }

  private:
    size_type size_, offset_;
  };



  // integerRange
  // ------------

  template< class T >
  inline static IntegerRange< T, std::integral_constant< T, 1 > > integerRange ( T size, T offset = 0 ) noexcept
  {
    return IntegerRange< T, std::integral_constant< T, 1 > >( size, offset );
  }

  template< class T >
  inline static IntegerRange< T > integerRange ( T size, T offset, T stride ) noexcept
  {
    return IntegerRange< T >( size, offset, stride );
  }

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_COMMON_INTEGERRANGE_HH
