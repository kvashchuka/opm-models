#ifndef DUNE_GRID_UTILITY_VECTORCOMMDATAHANDLE_HH
#define DUNE_GRID_UTILITY_VECTORCOMMDATAHANDLE_HH

#include <cassert>
#include <cstddef>

#include <type_traits>
#include <utility>

#include <dune/geometry/type.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/utility/commdatahelper.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class, template< int > class >
  class MultipleCodimMultipleGeomTypeMapper;

  template< class, int >
  class SingleCodimSingleGeomTypeMapper;



  // Internal Forward Declarations
  // -----------------------------

  template< class Mapper, class Vector, class Function >
  class VectorCommDataHandle;



  namespace __VectorCommDataHandle
  {

    template< class Vector >
    using CommDataHelper = Dune::CommDataHelper< typename std::decay< typename Vector::value_type >::type >;

    template< class Mapper, class Vector, class Function >
    using CommDataHandleIF = Dune::CommDataHandleIF< VectorCommDataHandle< Mapper, Vector, Function >, typename CommDataHelper< Vector >::Data >;

  } // namespace __VectorCommDataHandle



  // VectorCommDataHandle
  // --------------------

  template< class GV, int cd, class Vector, class Function >
  class VectorCommDataHandle< SingleCodimSingleGeomTypeMapper< GV, cd >, Vector, Function >
    : public __VectorCommDataHandle::CommDataHandleIF< SingleCodimSingleGeomTypeMapper< GV, cd >, Vector, Function >
  {
    typedef VectorCommDataHandle< SingleCodimSingleGeomTypeMapper< GV, cd >, Vector, Function > This;

    typedef SingleCodimSingleGeomTypeMapper< GV, cd > Mapper;

    typedef __VectorCommDataHandle::CommDataHelper< Vector > CommDataHelper;

  public:
    VectorCommDataHandle ( const Mapper &mapper, Vector &vector, Function function = Function() )
      : mapper_( mapper ), vector_( vector ), function_( function )
    {}

    bool contains ( int dim, int codim ) const { return (codim == cd); }
    bool fixedsize ( int dim, int codim ) const { return true; }

    template< class Entity >
    std::size_t size ( const Entity &entity ) const
    {
      typename Mapper::Index index;
      return (mapper_.contains( entity, index ) ? CommDataHelper::size() : 0);
    }

    template< class Buffer, class Entity >
    void gather ( Buffer &buffer, const Entity &entity ) const
    {
      typename Mapper::Index index;
      if( mapper_.contains( entity, index ) )
        CommDataHelper::write( buffer, vector_[ index ] );
    }

    template< class Buffer, class Entity >
    void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
    {
      assert( n == size( entity ) );
      typename Mapper::Index index;
      if( mapper_.contains( entity, index ) )
        vector_[ index ] = function_( vector_[ index ], CommDataHelper::read( buffer ) );
    }

  private:
    const Mapper &mapper_;
    Vector &vector_;
    Function function_;
  };

  template< class GV, template< int > class Layout, class Vector, class Function >
  class VectorCommDataHandle< MultipleCodimMultipleGeomTypeMapper< GV, Layout >, Vector, Function >
    : public __VectorCommDataHandle::CommDataHandleIF< MultipleCodimMultipleGeomTypeMapper< GV, Layout >, Vector, Function >
  {
    typedef VectorCommDataHandle< MultipleCodimMultipleGeomTypeMapper< GV, Layout >, Vector, Function > This;

    typedef MultipleCodimMultipleGeomTypeMapper< GV, Layout > Mapper;

    typedef __VectorCommDataHandle::CommDataHelper< Vector > CommDataHelper;

  public:
    VectorCommDataHandle ( const Mapper &mapper, Vector &vector, Function function = Function() )
      : mapper_( mapper ), vector_( vector ), function_( function )
    {}

    bool contains ( int dim, int codim ) const
    {
      return codim == 0;

      // mapper_.layout is private, so assume it is default-constructed
      Layout< GV::dimension > layout;

      // a codim is contained, if the layout contains any geometry type of
      // dimension dim - codim
      const int mydim = dim - codim;
      for( unsigned int id = 1; id < 2u*(1u << ((mydim-1)-1)); id += 2 )
      {
        if( layout.contains( GeometryType( id, mydim ) ) )
          return true;
      }
      return layout.contains( GeometryType( GeometryType::none, mydim ) );
    }

    bool fixedsize ( int dim, int codim ) const
    {
      // if a codim is not contained, it has fixed size (namely zero)
      if( !contains( dim, codim ) )
        return true;

      // mapper_.layout is private, so assume it is default-constructed
      Layout< GV::dimension > layout;

      // a codim has fixed size, if the layout contains all geometry types of
      // dimension dim - codim
      const int mydim = dim - codim;
      for( unsigned int id = 1; id < 2u*(1u << ((mydim-1)-1)); id += 2 )
      {
        if( !layout.contains( GeometryType( id, mydim ) ) )
          return false;
      }
      return layout.contains( GeometryType( GeometryType::none, mydim ) );
    }

    template< class Entity >
    std::size_t size ( const Entity &entity ) const
    {
      typename Mapper::Index index;
      return (mapper_.contains( entity, index ) ? CommDataHelper::size() : 0);
    }

    template< class Buffer, class Entity >
    void gather ( Buffer &buffer, const Entity &entity ) const
    {
      typename Mapper::Index index;
      if( mapper_.contains( entity, index ) )
        CommDataHelper::write( buffer, vector_[ index ] );
    }

    template< class Buffer, class Entity >
    void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
    {
      assert( n == size( entity ) );
      typename Mapper::Index index;
      if( mapper_.contains( entity, index ) )
        vector_[ index ] = function_( vector_[ index ], CommDataHelper::read( buffer ) );
    }

  private:
    const Mapper &mapper_;
    Vector &vector_;
    Function function_;
  };



  // vectorCommDataHandle
  // --------------------

  template< class Mapper, class Vector, class Function >
  inline static VectorCommDataHandle< Mapper, Vector, Function >
  vectorCommDataHandle( const Mapper &mapper, Vector &vector, Function function )
  {
    return VectorCommDataHandle< Mapper, Vector, Function >( mapper, vector, std::move( function ) );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_UTILITY_VECTORCOMMDATAHANDLE_HH
