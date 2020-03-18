#ifndef DUNE_OPTIM_COMMON_MATRIX_HH
#define DUNE_OPTIM_COMMON_MATRIX_HH

#include <dune/optim/common/vector.hh>

namespace Dune
{

  //! This interface provides uniform access to all kinds of matrices.
  template< class ct, class Imp >
  class MatrixInterface
  {
    typedef MatrixInterface< ct, Imp > This;

  public:
    //! Export the field type
    typedef ct Field;

    enum
    {
      //! Declare any derived class to implement MatrixInterface
      _MatrixInterface_ = true
    };

    //! Returns a const reference to the matrix element indexed by row and column
    const Field &operator() ( int row, int column ) const
    {
      return asImp().operator()( row, column );
    }

    //! Returns a reference to thhe matrix element indexed by row and column
    Field &operator() ( int row, int column )
    {
      return asImp().operator()( row, column );
    }

    //! Initialize the matrix with a scalar
    Imp &operator= ( const Field s )
    {
      return asImp().operator=( s );
    }

    //! Returns the number of columns
    const unsigned int columns () const
    {
      return asImp().columns();
    }

    //! Multiply this matrix by any vector implementing VectorInterface
    template< class DomainType, class RangeType >
    RangeType &multiply ( const DomainType &x, RangeType &y ) const
    {
      return asImp().multiply( x, y );
    }

    //! Multiply the transposed of this matrix by any vector implementing VectorInterface
    template< class DomainType, class RangeType >
    RangeType &multiplyTransposed ( const DomainType &x, RangeType &y ) const
    {
      return asImp().multiply( x, y );
    }

    //! Returns the number of rows
    const unsigned int rows () const
    {
      return asImp().rows();
    }

  protected:
    //! Barton-Nackman trick
    const Imp &asImp () const
    {
      return static_cast< const Imp & >( *this );
    }

    //! Barton-Nackman trick
    Imp &asImp ()
    {
      return static_cast< Imp & >( *this );
    }
  };



  //! Default implementation of MatrixInterface
  template< class ct, class Imp >
  class MatrixDefault
    : public MatrixInterface< ct, Imp >
  {
    typedef MatrixDefault< ct, Imp > This;
    typedef MatrixInterface< ct, Imp > Base;

  protected:
    using Base::asImp;

  public:
    //! Export the field type
    typedef typename Base::Field Field;

    //! Initialize the matrix with a scalar
    Imp &operator= ( const Field s )
    {
      const unsigned int rows = Base::rows();
      const unsigned int cols = Base::cols();

      for( unsigned int i = 0; i < rows; ++i )
        for( unsigned int j = 0; j < cols; ++j )
          (*this)( i, j ) = s;
      return asImp();
    }

    //! Multiply this matrix by any vector implementing VectorInterface
    template< class DomainVector, class RangeVector >
    RangeVector &multiply ( const DomainVector &x, RangeVector &y ) const
    {
      const unsigned int rows = Base::rows();
      const unsigned int cols = Base::columns();
      assert( cols == x.size() );
      assert( rows == y.size() );

      for( unsigned int i = 0; i < rows; ++i )
      {
        Field &yi = y[ i ];
        yi = 0;
        for( unsigned int j = 0; j < cols; ++j )
          yi += (*this)( i, j ) * x[ j ];
      }
      return y;
    }

    //! Multiply the transposed of this matrix by any vector implementing VectorInterface
    template< class RangeVector, class DomainVector >
    DomainVector &multiplyTransposed ( const RangeVector &x, DomainVector &y ) const
    {
      const unsigned int rows = Base::rows();
      const unsigned int cols = Base::columns();
      assert( rows == x.size() );
      assert( cols == y.size() );

      for( unsigned int i = 0; i < cols; ++i )
      {
        Field &yi = y[ i ];
        yi = 0;
        for( unsigned int j = 0; j < rows; ++j )
          yi += (*this)( j, i ) * x[ j ];
      }
      return y;
    }

    template< class RangeVector, class DomainVector >
    Field multiplyBoth ( const RangeVector &x, const DomainVector &y ) const
    {
      const unsigned int rows = Base::rows();
      const unsigned int cols = Base::columns();
      assert( rows == x.size() );
      assert( cols == y.size() );

      Field z( 0 );
      for( unsigned int i = 0; i < rows; ++i )
      {
        Field zi( 0 );
        for( unsigned int j = 0; j < cols; ++j )
          zi += (*this)( i, j ) * x[ j ];
        z += y[ i ] * zi;
      }
      return z;
    }
  };



  //! Readonly default implementation of MatrixInterface
  template< class ct, class Imp >
  class MatrixDefaultReadOnly
    : public MatrixDefault< ct, Imp >
  {
    typedef MatrixDefaultReadOnly< ct, Imp > This;
    typedef MatrixDefault< ct, Imp > Base;

  public:
    //! Export the field type
    typedef typename Base::Field Field;

    //! Initialize the matrix with a scalar
    Imp &operator= ( const Field s )
    {
      DUNE_THROW( MathError, "Trying to assign to a readonly matrix." );
    }

    Field &operator() ( int row, int column )
    {
      DUNE_THROW( MathError, "Trying to assign to a readonly matrix." );
    }
  };



  //! Group 4 matrices together to a block matrix (without doing any copying)
  template< class Matrix11Type, class Matrix12Type, class Matrix21Type, class Matrix22Type >
  class BlockMatrix
    : public MatrixDefault< typename Matrix11Type::Field, BlockMatrix< Matrix11Type, Matrix12Type, Matrix21Type, Matrix22Type > >
  {
  public:
    //! Export the field type
    typedef typename Matrix11Type::Field Field;

  private:
    typedef BlockMatrix< Matrix11Type, Matrix12Type, Matrix21Type, Matrix22Type > ThisType;
    typedef MatrixDefault< Field, ThisType > BaseType;

    const Matrix11Type &matrix11_; //!< Reference to the upper left submatrix
    const Matrix12Type &matrix12_; //!< Reference to the upper right submatrix
    const Matrix21Type &matrix21_; //!< Reference to the lower left submatrix
    const Matrix22Type &matrix22_; //!< Reference to the lower right submatrix

    const unsigned int columns11_; //!< Number of columns in the two left submatrices
    const unsigned int rows11_;    //!< Number of rows in the upper two submatrices
    const unsigned int columns_;   //!< Number of columns in the block matrix
    const unsigned int rows_;      //!< Number of rows in the block matrix

  public:
    //! Constructor taking 4 matrices, grouping them to a block matrix
    BlockMatrix ( const Matrix11Type &matrix11, const Matrix12Type &matrix12, const Matrix21Type &matrix21, const Matrix22Type &matrix22 )
      : matrix11_( matrix11 ), matrix12_( matrix12 ), matrix21_( matrix21 ), matrix22_( matrix22 ),
        columns11_( matrix11_.columns() ), rows11_( matrix11_.rows() ),
        columns_( columns11_ + matrix12_.columns() ), rows_( rows11_ + matrix21_.rows() )
    {
      assert( matrix21_.columns() == columns11_ );
      assert( matrix22_.columns() == columns_ - columns11_ );
      assert( matrix12_.rows() == rows11_ );
      assert( matrix22_.rows() == rows_ - rows11_ );
    }

    //! Returns a const reference to thhe matrix element indexed by row and column
    const Field &operator() ( unsigned int row, unsigned int column ) const
    {
      const int row2 = row - rows11_;
      const int column2 = column - columns11_;

      if( row2 < 0 )
      {
        if( column2 < 0 )
          return matrix11_( row, column );
        else
          return matrix12_( row, column2 );
      }
      else {
        if( column2 < 0 )
          return matrix21_( row2, column );
        else
          return matrix22_( row2, column2 );
      }
    }

    //! Returns a const reference to thhe matrix element indexed by row and column
    Field &operator() ( unsigned int row, unsigned int column )
    {
      const int row2 = row - rows11_;
      const int column2 = column - columns11_;

      if( row2 < 0 )
      {
        if( column2 < 0 )
          return matrix11_( row, column );
        else
          return matrix12_( row, column2 );
      }
      else {
        if( column2 < 0 )
          return matrix21_( row2, column );
        else
          return matrix22_( row2, column2 );
      }
    }

    //! Initialize the matrix with a scalar
    ThisType &operator= ( const Field s )
    {
      matrix11_ = s;
      matrix12_ = s;
      matrix21_ = s;
      matrix22_ = s;
      return *this;
    }

    //! Returns the number of columns
    const unsigned int columns () const
    {
      return columns_;
    }

    //! Returns the number of rows
    const unsigned int rows () const
    {
      return rows_;
    }
  };



  // This class is capable of wrapping a DUNE FieldMatrix
  template< class M >
  class FieldMatrixWrapper
    : public MatrixDefault< typename M::field_type, FieldMatrixWrapper< M > >
  {
    typedef FieldMatrixWrapper< M > This;
    typedef MatrixDefault< typename M::field_type, This > Base;

  public:
    typedef typename Base::Field Field;

    typedef M MatrixType;

    FieldMatrixWrapper ( const MatrixType &matrix )
      : matrix_( matrix )
    {}

    const Field &operator() ( unsigned int row, unsigned int column ) const
    {
      return matrix_[ row ][ column ];
    }

    //! Initialize the matrix with a scalar
    This &operator= ( const Field s )
    {
      return ((Base *)this)->operator=( s );
    }

    Field &operator() ( int row, int column )
    {
      return matrix_[ row ][ column ];
    }

    const unsigned int columns () const
    {
      return MatrixType::cols;
    }

    const unsigned int rows () const
    {
      return MatrixType::rows;
    }

  private:
    const MatrixType &matrix_;
  };



  //! A representation of the null matrix with just one element
  template< class ct >
  class NullMatrix
    : public MatrixDefaultReadOnly< ct, NullMatrix< ct > >
  {
    typedef NullMatrix< ct > This;
    typedef MatrixDefaultReadOnly< ct, This > Base;

  public:
    //! Export the field type
    typedef typename Base::Field Field;

    //! Constructor requiring the number of rows and columns
    NullMatrix ( unsigned int columns, unsigned int rows )
      : columns_( columns ),
        rows_( rows ),
        null_( 0 )
    {}

    //! Returns a const reference to thhe matrix element indexed by row and column
    const Field &operator() ( unsigned int row, unsigned int column ) const
    {
      assert( row < rows_ );
      assert( column < columns_ );
      return null_;
    }

    //! Initialize the matrix with a scalar
    This &operator= ( const Field s )
    {
      return ((Base *)this)->operator=( s );
    }

    //! Returns the number of columns
    unsigned int columns () const { return columns_; }

    //! Multiply this matrix by any vector implementing VectorInterface
    template< class DomainType, class RangeType >
    RangeType &multiply ( const DomainType &x, RangeType &y ) const
    {
      static_assert( Fem::SupportsVectorInterface< DomainType >::value, "DomainType must be a vector" );
      static_assert( Fem::SupportsVectorInterface< RangeType >::value, "RangeType must be a vector" );

      y.assign( 0 );
      return y;
    }

    //! Multiply the transposed of this matrix by any vector implementing VectorInterface
    template< class DomainType, class RangeType >
    RangeType &multiplyTransposed ( const DomainType &x, RangeType &y ) const
    {
      static_assert( Fem::SupportsVectorInterface< DomainType >::value, "DomainType must be a vector" );
      static_assert( Fem::SupportsVectorInterface< RangeType >::value, "RangeType must be a vector" );

      y.assign( 0 );
      return y;
    }

    //! Returns the number of rows
    unsigned int rows () const { return rows_; }

  private:
    const unsigned int columns_;
    const unsigned int rows_;
    const Field null_;
  };



  //! This class transposes a matrix without actually copying it. It simply wraps the matrix.
  template< class M >
  class TransposedMatrix
    : MatrixDefault< typename M::Field, TransposedMatrix< M > >
  {
    typedef TransposedMatrix< M > This;
    typedef MatrixDefault< typename M::Field, This > Base;

  public:
    //! Export the field type
    typedef typename Base::Field Field;

    typedef M MatrixType;

    //! Constructor wrapping a matrix
    TransposedMatrix ( const MatrixType &matrix )
      : matrix_( matrix )
    {}

    //! Returns a const reference to thhe matrix element indexed by row and column
    const Field &operator() ( int row, int column ) const
    {
      return matrix_( column, row );
    }

    //! Returns a reference to thhe matrix element indexed by row and column
    Field &operator() ( int row, int column )
    {
      return matrix_( column, row );
    }

    //! Initialize the matrix with a scalar
    This &operator= ( const Field s )
    {
      matrix_ = s;
      return *this;
    }

    //! Returns the number of columns
    const unsigned int columns () const
    {
      return matrix_.rows();
    }

    //! Multiply this matrix by any vector implementing VectorInterface
    template< class DomainType, class RangeType >
    RangeType &multiply ( const DomainType &x, RangeType &y ) const
    {
      return matrix_.multiplyTransposed( x, y );
    }

    //! Multiply the transposed of this matrix by any vector implementing VectorInterface
    template< class DomainType, class RangeType >
    RangeType &multiplyTransposed ( const DomainType &x, RangeType &y ) const
    {
      return matrix_.multiply( x, y );
    }

    //! returns the number of rows
    const unsigned int rows () const
    {
      return matrix_.columns();
    }

  private:
    const MatrixType &matrix_; //!< Reference to the wrapped matrix
  };



  //! Print any matrix implementing MatrixInterface
  template< class M >
  inline std::ostream &
  operator<< ( std::ostream &out, const MatrixInterface< typename M::Field, M > &m )
  {
    const unsigned int rows = m.rows();
    const unsigned int cols = m.columns();

    out << rows << "x" << cols << " matrix"<< std::endl;
    for( unsigned int i = 0; i < rows; ++i )
    {
      out << "|";
      for( unsigned int j = 0; j < cols; ++j )
      {
        char field[ 32 ];
        sprintf( field, " %8.3f ", m( i, j ) );
        out << field;
      }
      out << "|" << std::endl;
    }
    return (out << std::endl);
  }

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_COMMON_MATRIX_HH
