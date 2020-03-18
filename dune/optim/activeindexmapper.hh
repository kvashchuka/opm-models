#ifndef DUNE_OPTIM_ACTIVEINDEXMAPPER_HH
#define DUNE_OPTIM_ACTIVEINDEXMAPPER_HH

#include <algorithm>
#include <limits>
#include <memory>

namespace Dune
{

  namespace Optim
  {

    // ActiveIndexMapper
    // -----------------

    // Implements a mapping to select some active subindices
    template< class A = std::allocator< unsigned int > >
    struct ActiveIndexMapper
    {
      typedef A Allocator;

      typedef unsigned int size_type;

      typedef const size_type *const_iterator;

      struct InactiveIterator;

      static std::mutex m;

      ActiveIndexMapper ( size_type maxSize, size_type baseSize, const Allocator &allocator = Allocator() )
        : maxSize_( maxSize ),
          baseSize_( baseSize ),
          size_( 0 ),
          allocator_( allocator )
      {
        m.lock();
        map_ = allocator_.allocate( maxSize_ + baseSize_ );
        std::fill( invmap(), invmap() + baseSize_, -1 );
        m.unlock();
      }

      ActiveIndexMapper ( const ActiveIndexMapper & ) = delete;
      ActiveIndexMapper ( ActiveIndexMapper && ) = delete;

      ~ActiveIndexMapper ()
      {
        m.lock();
        allocator_.deallocate( map_, maxSize_ + baseSize_ );
        m.unlock();
      }

      ActiveIndexMapper &operator= ( const ActiveIndexMapper & ) = delete;
      ActiveIndexMapper &operator= ( ActiveIndexMapper && ) = delete;

      size_type operator[] ( size_type index ) const
      {
        assert( index < size_ );
        return map_[ index ];
      }

      const_iterator begin () const { return map_; }
      const_iterator end () const { return map_ + size_; }

      std::reverse_iterator< const_iterator > rbegin () const { return end(); }
      std::reverse_iterator< const_iterator > rend () const { return begin(); }

      size_type activate ( size_type baseIndex )
      {
        assert( baseIndex < baseSize_ );
        if( invmap()[ baseIndex ] == std::numeric_limits< size_type >::max() )
        {
          assert( size_ < maxSize_ );
          invmap()[ baseIndex ] = size_;
          map_[ size_ ] = baseIndex;
          return size_++;
        }
        else
          return invmap()[ baseIndex ];
      }

      void clear ()
      {
        while( size_ > 0 )
          invmap()[ map()[ --size_ ] ] = std::numeric_limits< size_type >::max();
      }

      InactiveIterator beginInactive () const;
      InactiveIterator endInactive () const;

      bool isActive ( size_type baseIndex ) const
      {
        assert( baseIndex < baseSize_ );
        return (invmap()[ baseIndex ] < std::numeric_limits< size_type >::max());
      }

      void update ( size_type index, size_type baseIndex )
      {
        assert( index < size() );
        assert( !isActive( baseIndex ) );
        invmap()[ map()[ index ] ] = std::numeric_limits< size_type >::max();
        invmap()[ baseIndex ] = index;
        map()[ index ] = baseIndex;
      }

      size_type range () const { return baseSize_; }

      void release ( size_type index )
      {
        assert( index < size_ );
        invmap()[ map()[ --size_ ] ] = index;
        invmap()[ map()[ index ] ] = std::numeric_limits< size_type >::max();
        map()[ index ] = map()[ size_ ];
      }

      size_type size () const { return size_; }

    private:
      size_type *map () const { return map_; }
      size_type *invmap () const { return map_ + maxSize_; }

      const size_type maxSize_;
      const size_type baseSize_;
      size_type size_;
      Allocator allocator_;
      size_type *map_;
    };



    // ActiveIndexMapper::InactiveIterator
    // -----------------------------------

    template< class A >
    struct ActiveIndexMapper< A >::InactiveIterator
    {
      InactiveIterator ( size_type index, size_type baseSize, const size_type *invmap )
        : baseSize_( baseSize ),
          invmap_( invmap )
      {
        index_ = nextStop( index );
      }

      const InactiveIterator &operator++ ()
      {
        index_ = nextStop( index_+1 );
        return *this;
      }

      size_type operator* () const { return index_; }

      bool operator== ( const InactiveIterator &other ) const
      {
        return (index_ == other.index_) && (invmap_ == other.invmap_);
      }

      bool operator!= ( const InactiveIterator &other ) const
      {
        return (index_ != other.index_) || (invmap_ != other.invmap_);
      }

    private:
      size_type nextStop ( size_type index ) const
      {
        while( (index < baseSize_) && (invmap_[ index ] < std::numeric_limits< size_type >::max()) )
          ++index;
        return index;
      }

      size_type index_, baseSize_;
      const size_type *invmap_;
    };



    // Implementation of ActiveIndexMapper
    // -----------------------------------

    template< class A >
    typename ActiveIndexMapper< A >::InactiveIterator
    ActiveIndexMapper< A >::beginInactive () const
    {
      return InactiveIterator( 0, baseSize_, invmap() );
    }

    template< class A >
    typename ActiveIndexMapper< A >::InactiveIterator
    ActiveIndexMapper< A >::endInactive () const
    {
      return InactiveIterator( baseSize_, baseSize_, invmap() );
    }   

    template< class A >
    std::mutex ActiveIndexMapper< A >::m;
  } // namespace Optim


} // namespace Dune

#endif // #ifndef DUNE_OPTIM_ACTIVEINDEXMAPPER_HH
