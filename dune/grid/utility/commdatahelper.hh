#ifndef DUNE_GRID_UTILITY_COMMDATAHELPER_HH
#define DUNE_GRID_UTILITY_COMMDATAHELPER_HH

#include <type_traits>
#include <utility>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class K, int SIZE >
  class FieldVector;

  template< class K, int ROWS, int COLS >
  class FieldMatrix;



  // CommDataHelper
  // --------------

  template< class T >
  struct CommDataHelper
  {
    static_assert( std::is_trivial< T >::value, "Default CommDataHelper is only defined for trivial data types." );

    typedef T Data;

    template< class MessageBuffer >
    static T read ( MessageBuffer &buffer )
    {
      T value;
      buffer.read( value );
      return std::move( value );
    }

    template< class MessageBuffer >
    static void write ( MessageBuffer &buffer, const T &value )
    {
      buffer.write( value );
    }

    static std::size_t size () noexcept { return 1u; }
  };



  // CommDataHelper for FieldVector
  // ------------------------------

  template< class K, int SIZE >
  struct CommDataHelper< FieldVector< K, SIZE > >
  {
    typedef typename CommDataHelper< K >::Data Data;

    template< class MessageBuffer >
    static FieldVector< K, SIZE > read ( MessageBuffer &buffer )
    {
      FieldVector< K, SIZE > value;
      for( int i = 0; i < SIZE; ++i )
        value[ i ] = CommDataHelper< K >::read( buffer );
      return std::move( value );
    }

    template< class MessageBuffer >
    static void write ( MessageBuffer &buffer, const FieldVector< K, SIZE > &value )
    {
      for( int i = 0; i < SIZE; ++i )
        CommDataHelper< K >::write( buffer, value[ i ] );
    }

    static std::size_t size () noexcept { return SIZE; }
  };



  // CommDataHelper for FieldMatrix
  // ------------------------------

  template< class K, int ROWS, int COLS >
  struct CommDataHelper< FieldMatrix< K, ROWS, COLS > >
  {
    typedef typename CommDataHelper< K >::Data Data;

    template< class MessageBuffer >
    static FieldMatrix< K, ROWS, COLS > read ( MessageBuffer &buffer )
    {
      FieldMatrix< K, ROWS, COLS > value;
      for( int i = 0; i < ROWS; ++i )
        for( int j = 0; j < COLS; ++j )
          value[ i ][ j ] = CommDataHelper< K >::read( buffer );
      return std::move( value );
    }

    template< class MessageBuffer >
    static void write ( MessageBuffer &buffer, const FieldMatrix< K, ROWS, COLS > &value )
    {
      for( int i = 0; i < ROWS; ++i )
        for( int j = 0; j < COLS; ++j )
          CommDataHelper< K >::write( buffer, value[ i ][ j ] );
    }

    static std::size_t size () noexcept { return ROWS * COLS; }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_UTILITY_COMMDATAHELPER_HH
