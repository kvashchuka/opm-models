#ifndef DUNE_OPTIM_SOLVER_ROWREFVECTOR_HH
#define DUNE_OPTIM_SOLVER_ROWREFVECTOR_HH

#include <dune/common/densevector.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class K >
  class RowReferenceVector;



  // External Forward Declarations
  // -----------------------------

  template< class K, class A >
  class FullyDynamicMatrix;



  // DenseMatVecTraits for RowReferenceVector
  // ----------------------------------------

  template< class K >
  struct DenseMatVecTraits< RowReferenceVector< K > >
  {
    typedef RowReferenceVector< K > derived_type;
    typedef K value_type;
    typedef std::size_t size_type;
  };



  // FieldTraits for RowReferenceVector
  // ----------------------------------

  template< class K >
  struct FieldTraits< RowReferenceVector< K > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };



  template< class K >
  struct const_reference< RowReferenceVector< K > >
  {
    typedef RowReferenceVector< const K > type;
  };

  template< class K >
  struct const_reference< RowReferenceVector< const K > >
  {
    typedef RowReferenceVector< const K > type;
  };


  template< class K >
  struct mutable_reference< RowReferenceVector< K > >
  {
    typedef RowReferenceVector< K > type;
  };

  template< class K >
  struct mutable_reference< RowReferenceVector< const K > >
  {
    typedef RowReferenceVector< K > type;
  };



  // RowReferenceVector
  // ------------------

  template< class K >
  class RowReferenceVector
    : public DenseVector< RowReferenceVector< K > >
  {
    typedef DenseVector< RowReferenceVector< K > > Base;

  public:
    template< class, class >
    friend class FullyDynamicMatrix;

    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;

  protected:
    // only FullDynamicMatrix should be abled to call ctor
    RowReferenceVector ( K *data, size_type size )
      : _data( data ),
        _size( size )
    {}

  public:
    // allow copying
    RowReferenceVector ( const RowReferenceVector &other )
      : _data( other._data ),
        _size( other._size )
    {}

    using Base::operator=;

    size_type size () const { return _size; }
    const K &operator[] ( size_type i ) const { return _data[ i ]; }
    K &operator[] ( size_type i ) { return _data[ i ]; }

  private:
    K *_data;
    size_type _size;
  };

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_SOLVER_ROWREFVECTOR_HH
