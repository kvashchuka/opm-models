#ifndef DUNE_OPTIM_SOLVER_FULLYDYNMATRIX_HH
#define DUNE_OPTIM_SOLVER_FULLYDYNMATRIX_HH

#include <dune/common/densematrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/optim/solver/rowrefvector.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class K, class A >
  class FullyDynamicMatrix;



  // DenseMatVecTraits for FullyDynamicMatrix
  // ----------------------------------------

  template< class K, class A >
  struct DenseMatVecTraits< FullyDynamicMatrix< K, A > >
  {
    typedef FullyDynamicMatrix< K, A > derived_type;

    typedef DynamicVector< K > row_type;

    typedef RowReferenceVector< K > row_reference;
    typedef RowReferenceVector< const K > const_row_reference;

    typedef K value_type;
    typedef typename A::size_type size_type;
  };



  // FieldTraits for FullyDynamicMatrix
  // ----------------------------------

  template< class K, class A >
  struct FieldTraits< FullyDynamicMatrix< K, A > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };



  // FullyDynamicMatrix
  // ------------------

  template< class K, class A = std::allocator< K > >
  class FullyDynamicMatrix
    : public DenseMatrix< FullyDynamicMatrix< K, A > >
  {
    typedef FullyDynamicMatrix< K, A > This;
    typedef DenseMatrix< FullyDynamicMatrix< K, A > > Base;

  public:
    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;
    typedef typename Base::row_type row_type;
    typedef typename Base::row_reference row_reference;
    typedef typename Base::const_row_reference const_row_reference;

    typedef A allocator_type;
    typedef std::pair< size_type, size_type > dimensions;

  protected:
    typedef typename A::template rebind< K >::other real_allocator_type;
    typedef typename real_allocator_type::pointer pointer;

  public:
    /**
     * \name Construction and Destruction
     * \{
     */

    //! \brief Default constructor
    explicit FullyDynamicMatrix ( const allocator_type &allocator = allocator_type() )
      : size_( std::make_pair( 0, 0 ) ),
        capacity_( std::make_pair( 0, 0 ) ),
        allocator_( allocator )
    {}

    /** \brief constructor default initializing all matrix values */
    FullyDynamicMatrix ( size_type r, size_type c, const allocator_type &allocator = allocator_type() )
      : size_( std::make_pair( 0, 0 ) ),
        capacity_( std::make_pair( 0, 0 ) ),
        allocator_( allocator )
    {
      resize( r, c );
    }

    //! \brief Constructor initializing the whole matrix with a scalar
    FullyDynamicMatrix ( size_type r, size_type c, const value_type &v = value_type(),
                         const allocator_type &allocator = allocator_type() )
      : size_( std::make_pair( 0, 0 ) ),
        capacity_( std::make_pair( 0, 0 ) ),
        allocator_( allocator )
    {
      resize( r, c, v );
    }

    //! \brief Constructor initializing the whole matrix with a scalar
    FullyDynamicMatrix ( dimensions size, const value_type &v = value_type(),
                         const allocator_type &allocator = allocator_type() )
      : size_( std::make_pair( 0, 0 ) ),
        capacity_( std::make_pair( 0, 0 ) ),
        allocator_( allocator )
    {
      resize( size, v );
    }

    /** \brief copy constructor */
    FullyDynamicMatrix ( const This &other )
      : size_( other.size_ ),
        capacity_( size_ ),
        allocator_( other.allocator_ )
    {
      const size_type numEntries = capacity_.first * capacity_.second;
      values_ = (numEntries > 0 ?allocator_.allocate( numEntries ) : nullptr);

      // copy values
      for( size_type row = 0; row < size_.first; ++row )
        for( size_type col = 0; col < size_.second; ++col )
          allocator_.construct( values_ + entry( row, col, capacity_ ), *( other.values_ + entry( row, col, other.capacity_ ) ) );
    }

    /** \brief move constructor */
    FullyDynamicMatrix ( This &&other )
      : values_( other.values_ ),
        size_( other.size_ ),
        capacity_( other.capacity_ ),
        allocator_( std::move( other.allocator_ ) )
    {
      other.values_ = nullptr;
      other.size_ = std::make_pair( 0, 0 );
      other.capacity_ = std::make_pair( 0, 0 );
    }

    /** \brief destructor */
    ~FullyDynamicMatrix ()
    {
      clear();
      shrink_to_fit();
    }

    /** \} */

    //===== assignment
    using Base::operator=;

    This &operator= ( const This &other )
    {
      // just paranoia to do nothing in case of self assignment
      if( this == &other )
        return *this;

      clear();
      if( capacity_ != other.size_ )
        doReserve( other.size_ );

      // copy values
      size_ = other.size_;
      for( size_type row = 0; row < size_.first; ++row )
        for( size_type col = 0; col < size_.second; ++col )
          allocator_.construct( values_ + entry( row, col, capacity_ ), *( other.values_ + entry( row, col, other.capacity_ ) ) );

      return *this;
    }

    This &operator= ( This &&other )
    {
      swap( other );
      return *this;
    }

    /**
     * \brief Size and Capacity
     * \{
     */

    void resize ( size_type r, size_type c, const value_type &v = value_type() )
    {
      resize( std::make_pair( r, c ), v );
    }

    void resize ( dimensions size, const value_type &v = value_type() )
    {
      if( size.first == size_.first && size.second == size_.second )
        return;

      if( size.first > capacity_.first || size.second > capacity_.second )
        reserve( std::make_pair( std::max( 2 * capacity_.first, size.first ), std::max( 2 * capacity_.second, size.second ) ) );

      // delete extra rows
      for( size_type row = size.first; row < size_.first; ++row )
        for( size_type col = 0; col < size_.second; ++col )
          allocator_.destroy( values_ + entry( row, col, capacity_ ) );

      size_.first = std::min( size_.first, size.first );

      for( size_type col = size.second; col < size_.second; ++col )
        for( size_type row = 0; row < size_.first; ++row )
          allocator_.destroy( values_ + entry( row, col, capacity_ ) );

      size_.second = std::min( size_.second, size.second );


      for( size_type row = size_.first; row < size.first; ++row )
        for( size_type col = 0; col < size_.second; ++col )
          allocator_.construct( values_ + entry( row, col, capacity_ ), v );

      size_.first = std::max( size_.first, size.first );

      for( size_type col = size_.second; col < size.second; ++col )
        for( size_type row = 0; row < size_.first; ++row )
          allocator_.construct( values_ + entry( row, col, capacity_ ), v );

      size_.second = std::max( size_.second, size.second );
    }

    void reserve ( size_type r, size_type c )
    {
      reserve( std::make_pair( r, c ) );
    }

    void reserve ( dimensions capacity  )
    {
      capacity = std::make_pair( std::max( capacity.first, capacity_.first ), std::max( capacity.second, capacity_.second ) );
      if( capacity.first > capacity_.first || capacity.second > capacity_.second )
        doReserve( capacity );
    }

    void clear () { resize( std::make_pair( 0, 0 ) ); }

    bool empty () { return (size_.first * size_.second == 0);  }

    dimensions size () const { return size_; }
    dimensions capacity () const { return capacity_; }

    void shrink_to_fit ()
    {
      if( capacity_ != size_ )
        doReserve( size_ );
    }

    /** \} */

    // make this thing a matrix
    size_type mat_rows () const { return size_.first; }
    size_type mat_cols () const { return size_.second; }

    row_reference mat_access ( size_type i ) { return row_reference( values_ + i * capacity_.second, size_.second ); }
    const_row_reference mat_access ( size_type i ) const { return const_row_reference( values_ + i * capacity_.second, size_.second ); }

    /**
     * \name Modifiers
     * \{
     */

    template< class X >
    void push_row ( const DenseVector< X > &v )
    {
      resize( size_.first+1, std::max( size_.second, v.size() ) );

      const size_type row = size_.first-1;
      for( size_type col = 0; col < v.size(); ++col )
        *(values_ + entry( row, col, capacity_ )) = v[ col ];
    }

    template< class X >
    void push_col ( const DenseVector< X > &v )
    {
      resize( std::max( size_.first, v.size() ), size_.second+1 );

      const size_type col = size_.second-1;
      for( size_type row = 0; row < v.size(); ++row )
        *(values_ + entry( row, col, capacity_ )) = v[ row ];
    }

    void swap ( This &other )
    {
      using std::swap;
      swap( values_, other.values_ );
      swap( size_, other.size_ );
      swap( capacity_, other.capacity_ );
      // swap allocators
    }

    /** \} */

    allocator_type getallocator_ () const { return allocator_; }

  private:
    static size_type entry ( size_type r, size_type c, dimensions capacity )
    {
      return r * capacity.second + c;
    }

    void doReserve ( dimensions capacity )
    {
      const size_type numEntries = capacity.first * capacity.second;
      pointer values = (numEntries > 0 ?allocator_.allocate( numEntries ) : nullptr);

      // move matrix
      for( size_type row = 0; row < size_.first; ++row )
        for( size_type col = 0; col < size_.second; ++col )
        {
          allocator_.construct( values + entry( row, col, capacity ), *( values_ + entry( row, col, capacity_ ) ) );
          allocator_.destroy( values_ + entry( row, col, capacity_ ) );
        }

      if( values_ )
        allocator_.deallocate( values_, capacity_.first * capacity_.second );

      values_ = values;
      capacity_ = capacity;
    }

    pointer values_ = nullptr;
    dimensions size_, capacity_;
    real_allocator_type allocator_;
  };

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_SOLVER_FULLYDYNMATRIX_HH
