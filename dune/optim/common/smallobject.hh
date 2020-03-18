#ifndef DUNE_OPTIM_COMMON_SMALLOBJECT_HH
#define DUNE_OPTIM_COMMON_SMALLOBJECT_HH

#include <cassert>
#include <memory>
#include <new>
#include <utility>

namespace Dune
{

  // SmallObjectPool
  // ---------------

  class SmallObjectPool
  {
    union Block
    {
      Block *next;
      unsigned int blocks;
    };

  public:
    static const std::size_t blockSize = sizeof( Block );
    static const std::size_t maxBlocks = (std::size_t( 1 ) << 10) - std::size_t( 1 );
    static const std::size_t maxSize = maxBlocks * blockSize;

    static void *allocate ( unsigned int size )
    {
      const unsigned int blocks = (size + (blockSize-1)) / blockSize;
      if( blocks >= maxBlocks )
        return nullptr;
      Block *&next = list( blocks );
      Block *current = next;
      if( current )
        next = current->next;
      else
        current = new Block[ blocks+1 ];
      current->blocks = blocks;
      return current + 1;
    }

    static void free ( void *p )
    {
      if( p )
      {
        Block *current = reinterpret_cast< Block * >( p ) - 1;
        const unsigned int blocks = current->blocks;
        Block *&next = list( blocks );
        current->next = next;
        next = current;
      }
    }

  private:
    SmallObjectPool ()
    {
      for( std::size_t i = 0; i < maxBlocks; ++i )
        list_[ i ] = nullptr;
    }

    ~SmallObjectPool ()
    {
      for( std::size_t i = 0; i < maxBlocks; ++i )
      {
        for( Block *next = list_[ i ]; next; )
        {
          Block *current = next;
          next = current->next;
          delete[] current;
        }
      }
    }

    static SmallObjectPool &instance ()
    {
      static SmallObjectPool inst;
      return inst;
    }

    static Block *&list ( unsigned int blocks )
    {
      assert( blocks < maxBlocks );
      return instance().list_[ blocks ];
    }

    Block *list_[ maxBlocks ];
  };



  // SmallObject
  // -----------

  struct SmallObject
  {
    void *operator new ( size_t size )
    {
      return SmallObjectPool::allocate( size );
    }

    void operator delete ( void *p )
    {
      SmallObjectPool::free( p );
    }
  };



  // SmallObjectAllocator
  // --------------------

  template< class T >
  struct SmallObjectAllocator;

  template<>
  struct SmallObjectAllocator< void >
  {
    typedef void value_type;

    typedef void *pointer;
    typedef const void *const_pointer;

    template< class U > struct rebind { typedef SmallObjectAllocator< U > other; };
  };

  template< class T >
  struct SmallObjectAllocator
  {
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T *pointer;
    typedef const T *const_pointer;

    typedef T &reference;
    typedef const T &const_reference;

    template< class U > struct rebind { typedef SmallObjectAllocator< U > other; };

    SmallObjectAllocator () noexcept {}
    template< class U > SmallObjectAllocator ( const SmallObjectAllocator< U > & ) noexcept {}
    ~SmallObjectAllocator () noexcept {}

    pointer address ( reference r ) const noexcept { return &r; }
    const_pointer address ( const_reference r ) const noexcept { return &r; }

    pointer allocate ( size_type n, SmallObjectAllocator< void >::const_pointer hint = 0 );
    void deallocate ( pointer p, size_type n ) { SmallObjectPool::free( p ); }

    template< class ...Args >
    void construct ( pointer p, Args &&...args ) { ::new( p ) T( std::forward< Args >( args )...); }
    void destroy ( pointer p ) { p->~T(); }

    size_type max_size () const noexcept { return (SmallObjectPool::maxSize / sizeof( T )); }
  };



  // Implementation of SmallObjectAllocator
  // --------------------------------------

  template< class T >
  inline typename SmallObjectAllocator< T >::pointer
  SmallObjectAllocator< T >::allocate ( size_type n, SmallObjectAllocator< void >::const_pointer hint )
  {
    return static_cast< pointer >( SmallObjectPool::allocate( n * sizeof( T ) ) );
  }


  template< class T >
  bool operator== ( const SmallObjectAllocator< T > &, const SmallObjectAllocator< T > & ) throw()
  {
    return true;
  }


  template< class T >
  bool operator!= ( const SmallObjectAllocator< T > &, const SmallObjectAllocator< T > & ) throw()
  {
    return false;
  }

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_COMMON_SMALLOBJECT_HH
