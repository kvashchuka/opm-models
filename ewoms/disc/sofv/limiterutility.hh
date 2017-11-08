#ifndef EWOMS_LIMITERUTILITY_HH
#define EWOMS_LIMITERUTILITY_HH

#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/common/typeindexedtuple.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/pass/localdg.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>

#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/misc/compatibility.hh>

//*************************************************************
namespace Ewoms
{
//namespace Fem
//{

  template< class TypeTag >
  struct LimiterUtility
  {
//    typedef DiscreteFunctionSpace DiscreteFunctionSpaceType ;
//    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;
//    typedef typename FunctionSpaceType :: DomainType      DomainType;
//    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
//    typedef typename FunctionSpaceType :: RangeType       RangeType;
//    typedef typename FunctionSpaceType :: RangeFieldType  RangeFieldType;
//
//    static const int dimRange  = FunctionSpaceType :: dimRange;
//    static const int dimDomain = FunctionSpaceType :: dimDomain;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) DomainFieldType; //should be Scalar
      typedef double Field; //can be Scalar
      typedef Field RangeFieldType;
      typedef Field FieldType;

      typedef typename GET_PROP_TYPE(TypeTag, Grid) GridType;
      typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPart;
      typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
      enum { dimDomain = GridType::dimensionworld };
      enum { dimRange  = PrimaryVariables::dimension };
      //? dimWorld = dimDomain? enum { dimWorld = GridView::dimensionworld };

      typedef Dune::FieldVector<typename GridType::ctype, dimDomain> DomainType;
      typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

      typedef Dune::FieldVector<DomainFieldType, dimDomain> GradientType;
      //  typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector; this is a range gradient type, here we need domain gradient type
      typedef Dune::FieldMatrix<DomainFieldType, dimDomain, dimDomain> MatrixType;

      typedef typename GridView::template Codim<0>::Entity EntityType;
      typedef typename EntityType :: GeometryType GeometryType;


    //typedef FieldVector< DomainType , dimRange > GradientType;
    //typedef FieldMatrix< DomainFieldType, dimDomain , dimDomain > MatrixType;

   // static const int dimGrid = DiscreteFunctionSpaceType::GridType::dimension;
      static const int dimGrid = GridType::dimension;

    typedef Dune::DGFEntityKey<int> KeyType;
    typedef std::vector<int> CheckType;
    typedef std::pair< KeyType, CheckType > VectorCompType;
    typedef std::set< VectorCompType > ComboSetType;

    struct Flags
    {
      bool boundary;
      bool nonConforming;
      bool cartesian;
      bool limiter;
      Flags( bool cart, bool lim )
        : boundary( false ), nonConforming( false ), cartesian( cart ), limiter( lim ) {}
    };

    static bool functionIsPhysical( const DomainType& xGlobal,
                                    const RangeType& value,
                                    const GradientType& gradient )
    {
      RangeType res = value ;
      for ( int i=0; i<dimRange; ++ i )
      {
        res[ i ] += xGlobal * gradient[ i ];
        if( res[ i ] < 0.0 || res[ i ] > 1.0 )
          return false ;
      }
      return true;
    }

    //! limit all functions
    template <class LimiterFunction, class CheckPhysical, class CheckSet >
    static void limitFunctions(const LimiterFunction& limiterFunction,
                               const CheckPhysical& checkPhysical,
                               const std::vector< CheckSet >& comboVec,
                               const std::vector< DomainType >& barys,
                               const std::vector< RangeType  >& nbVals,
                               std::vector< GradientType >& gradients)
    {
      // get accuracy threshold
      const double limitEps = limiterFunction.epsilon();

      const size_t numFunctions = gradients.size();
      // for all functions check with all values
      for(size_t j=0; j<numFunctions; ++j)
      {
        const std::vector<int> & v = comboVec[j];

        // loop over dimRange
        for(int r=0; r<RangeType::dimension; ++r)
        {
          RangeFieldType minimalFactor = 1;
          DomainType& D = gradients[ j ][ r ];

          const auto endit = v.end();
          for(auto it = v.begin(); it != endit ; ++it )
          {
            // get current number of entry
            const size_t k = *it;

            // evaluate values for limiter function
            const DomainFieldType d = nbVals[ k ][ r ];
            const DomainFieldType g = D * barys[ k ];

            const DomainFieldType length2 =  barys[ k ].two_norm2();

            // if length is to small then the grid is corrupted
            assert( length2 > 1e-14 );

            // if the gradient in direction of the line
            // connecting the barycenters is very small
            // then neglect this direction since it does not give
            // valuable contribution to the linear function
            // call limiter function
            // g = grad L ( w_E,i - w_E ) ,  d = u_E,i - u_E
            DomainFieldType localFactor = limiterFunction( g, d );

            const DomainFieldType factor = (g*g) / length2 ;
            if( localFactor < 1.0 &&  factor  < limitEps )
            {
              localFactor = 1.0 - std::tanh( factor / limitEps );
            }

            // take minimum
            minimalFactor = std::min( localFactor , minimalFactor );
          }

          // scale linear function
          D *= minimalFactor;

          // if non-physical then set gradient to zero
          if( ! checkPhysical( r, D ) )
          {
            D = 0;
          }
        }
      }
    }

    // chose function with maximal gradient
    static void getMaxFunction(const std::vector< GradientType >& gradients,
                               GradientType& maxGradient)
    {
        enum { dimRange  = PrimaryVariables::dimension };

      const size_t numFunctions = gradients.size();

      const int startFunc = 0;
      std::vector< size_t > number(dimRange, startFunc);

      RangeType max (0);
      for(int r=0; r<dimRange; ++r)
      {
        max[r] = gradients[ startFunc ][ r ].two_norm2();
      }

      for(size_t l=1; l<numFunctions; ++l)
      {
        for(int r=0; r<dimRange; ++r)
        {
          RangeFieldType D_abs = gradients[ l ][ r ].two_norm2();
          if( D_abs > max[r] )
          {
            number[r] = l;
            max[r] = D_abs;
          }
        }
      }

      for(int r=0; r<dimRange; ++r)
      {
        maxGradient[r] = gradients[ number[r] ][r];
      }
    }


    // fill combination vector recursive
    template< class SetType, int i >
    struct FillVector
    {
      static void fill(const int neighbors,
                       const int start,
                       SetType& comboSet,
                       std::vector<int>& v)
      {
        assert( (int) v.size() > dimGrid-i );
        for(int n = start; n<neighbors; ++n)
        {
          v[ dimGrid - i ] = n;
          FillVector< SetType, i-1 >::fill( neighbors, n + 1, comboSet, v);
        }
      }
    };

    // termination of fill combination vector
    template< class SetType >
    struct FillVector< SetType, 1 >
    {
      static void fill(const int neighbors,
                       const int start,
                       SetType& comboSet,
                       std::vector<int>& v)
      {
        typedef typename SetType :: value_type VectorCompType;
        assert( (int) v.size() == dimGrid );
        for(int n = start; n<neighbors; ++n)
        {
          v[ dimGrid-1 ] = n;
          comboSet.insert( VectorCompType( v, std::vector<int> () ) );
        }
      }
    };

    // build combo set
    static void buildComboSet(const int neighbors,
                              ComboSetType& comboSet)
    {
      // clear set
      comboSet.clear();

      // maximal number of neighbors
      std::vector<int> v(dimGrid,0);
      FillVector< ComboSetType, dimGrid >::fill( neighbors, 0, comboSet, v );

      // create set containing all numbers
      std::set<int> constNumbers;
      typedef typename std::set<int> :: iterator el_iterator;
      for(int i = 0; i<neighbors; ++i)
      {
        constNumbers.insert(i);
      }

      const int checkSize = neighbors - dimGrid ;
      typedef typename ComboSetType :: iterator iterator;
      const iterator endit = comboSet.end();
      for(iterator it = comboSet.begin(); it != endit; ++it)
      {
        const KeyType& v = (*it).first;
        CheckType& check = const_cast<CheckType&> ((*it).second);

        // reserve memory
        check.reserve ( checkSize );

        // get set containing all numbers
        std::set<int> numbers (constNumbers);

        // remove the key values
        for(int i=0; i<dimGrid; ++i)
        {
          el_iterator el = numbers.find( v[i] );
          numbers.erase( el );
        }

        // generate check vector
        el_iterator endel = numbers.end();
        for(el_iterator elit = numbers.begin();
            elit != endel; ++ elit)
        {
          check.push_back( *elit );
        }
      }
    }

    // setup set storing combinations of linear functions
    static const ComboSetType& setupComboSet(const int neighbors,
                                             const bool nonConforming,
                                             const bool multipleGeomTypes,
                                             ComboSetType& conformingComboSet,
                                             ComboSetType& comboSet )
    {
      // check for new build (this set is constant)
      if( conformingComboSet.size() == 0 )
      {
        buildComboSet(neighbors, conformingComboSet);
      }

      // in case of non-conforming grid or grid with more than one
      // element type this has to be re-build on every element
      if( nonConforming || multipleGeomTypes )
      {
        buildComboSet(neighbors, comboSet);
        return comboSet;
      }
      else
      {
        return conformingComboSet;
      }
    }


    template <class GridPart, class Entity, class EvalAverage>
    static void
    setupNeighborValues( const GridPart& gridPart,
                         const Entity& entity,
                         const EvalAverage& average,
                         const DomainType& entityCenter,
                         const RangeType&  entityValue,
                         const bool StructuredGrid,
                         Flags& flags,
                         std::vector< DomainType >& baryCenters,
                         std::vector< RangeType  >& neighborValues )
    {
      typedef Entity EntityType;

      //! type of cartesian grid checker
      typedef Dune::Fem::CheckCartesian< GridPart >  CheckCartesianType;

      // get local references
      std::vector< DomainType >& barys  = baryCenters;
      std::vector< RangeType >&  nbVals = neighborValues;

      barys.reserve( dimGrid * dimGrid );
      nbVals.reserve( dimGrid * dimGrid );

      // clear old values
      barys.clear();
      nbVals.clear();

      typedef typename GridPart :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType :: Intersection  IntersectionType;

      const auto endit = gridPart.iend( entity );

      const bool cartesianGrid = flags.cartesian;

      // loop over all neighbors
      for (auto it = gridPart.ibegin( entity ); it != endit; ++it )
      {
        const IntersectionType& intersection = *it;

        const bool hasBoundary = intersection.boundary();
        const bool hasNeighbor = intersection.neighbor();

        // check cartesian
        if( !StructuredGrid && flags.cartesian )
        {
          // check whether this element is really cartesian
          flags.cartesian &= ( ! CheckCartesianType::checkIntersection(intersection) );
        }

        DomainType lambda( 1 );
        RangeType neighborValue( 0 );

        /////////////////////////////////////
        //  if we have a neighbor
        /////////////////////////////////////
        if ( hasNeighbor )
        {
          // check all neighbors
          const EntityType& neighbor = intersection.outside();

          // nonConforming case
          flags.nonConforming |= (! intersection.conforming() );

          // this is true in the periodic case
          if( ! hasBoundary )
          {
            // get barycenter of neighbor
            lambda = neighbor.geometry().center();
            // calculate difference
            lambda -= entityCenter;
          }

          // evaluate average value on neighbor
          flags.limiter |= average.evaluate( neighbor, neighborValue );

          // calculate difference
          neighborValue -= entityValue;

        } // end neighbor

        ////////////////////////////
        // --boundary
        ////////////////////////////
        // use ghost cell approach for limiting,
        if( intersection.boundary() )
        {
          // we have entity with boundary intersections
          flags.boundary = true ;

          typedef typename IntersectionType :: Geometry LocalGeometryType;
          const LocalGeometryType& interGeo = intersection.geometry();

          /////////////////////////////////////////
          // construct bary center of ghost cell
          /////////////////////////////////////////

          // get unit normal
          lambda = intersection.centerUnitOuterNormal();

          // get one point of intersection
          DomainType point ( interGeo.corner( 0 ) );
          point -= entityCenter;

          const double length = (point * lambda);
          const double factorLength = (cartesianGrid) ? 2.0 * length : length ;
          lambda *= factorLength;

          assert( lambda.two_norm () > 0 );

          /////////////////////////////////////////////////
          /////////////////////////////////////////////////
          // only when we don't have a neighbor
          // this can be true in periodic case
          if( ! hasNeighbor )
          {
            // check for boundary Value
            const DomainType pointOnBoundary = lambda + entityCenter;

            // evaluate data on boundary
            if( average.boundaryValue( entity, intersection, interGeo, pointOnBoundary, entityValue, neighborValue ) )
            {
              neighborValue -= entityValue;
            }
          }

        } //end boundary

        // store difference of mean values
        nbVals.push_back(neighborValue);

        // store difference between bary centers
        barys.push_back(lambda);

      } // end intersection iterator
    }

    // matrix assemblers for the reconstruction matrices
    template<int dimension,int dimensionworld = dimension, class DomainFieldType = double>
    struct MatrixAssemblerForLinearReconstruction
    {
      typedef Dune::DGFEntityKey<int> KeyType;
      typedef Dune::FieldMatrix< DomainFieldType, dimensionworld , dimensionworld > MatrixType;
      typedef Dune::FieldVector< DomainFieldType , dimensionworld > DomainType;

      static inline
      void assembleMatrix(const KeyType& v,
                          const std::vector< DomainType >& barys,
                          MatrixType& matrix)
      {
        assert(dimension == dimensionworld);
        const int size = v.size();
        for(int i=0; i<size; ++i)
        {
          matrix[ i ] = barys[ v[ i ] ];
        }
      }
    };


    template<class DomainFieldType>
    struct MatrixAssemblerForLinearReconstruction<1, 2, DomainFieldType>
    {
      typedef Dune::DGFEntityKey<int> KeyType;
      typedef Dune::FieldMatrix< DomainFieldType, 2 , 2 > MatrixType;
      typedef Dune::FieldVector< DomainFieldType , 2 > DomainType;

      static inline
      void assembleMatrix(const KeyType& v,
                          const std::vector< DomainType >& barys,
                          MatrixType& matrix)
      {
        const int size = v.size();
        assert( size==1 );

        for(int i=0; i<size; ++i)
        {
          matrix[ i ] = barys[ v[ i ] ];
        }
        //take a vector that is orthogonal to the first one
        matrix[ 1 ] [ 0 ]  = (-1.) * barys[ v[ 0 ] ] [ 1 ];
        matrix[ 1 ] [ 1 ]  = barys[ v[ 0 ] ] [ 0 ];
      }
    };

    template<class DomainFieldType>
    struct MatrixAssemblerForLinearReconstruction<2, 3, DomainFieldType>
    {
      typedef Dune::DGFEntityKey<int> KeyType;
      typedef Dune::FieldMatrix< DomainFieldType, 3 , 3 > MatrixType;
      typedef Dune::FieldVector< DomainFieldType , 3 > DomainType;

      static inline
      void assembleMatrix(const KeyType& v,
                          const std::vector< DomainType >& barys,
                          MatrixType& matrix)
      {
        const int size = v.size();
        assert( size==2 );

        for(int i=0; i<size; ++i)
        {
          matrix[ i ] = barys[ v[ i ] ];
        }
        //take the cross product of first and second vector as the third vector
        matrix[ 2 ] [ 0 ]  = barys[ v[ 0 ] ] [ 1 ] *  barys[ v[ 1 ] ] [ 2 ] - barys[ v[ 0 ] ] [ 2 ] *  barys[ v[ 1 ] ] [ 1 ];
        matrix[ 2 ] [ 1 ]  = barys[ v[ 0 ] ] [ 2 ] *  barys[ v[ 1 ] ] [ 0 ] - barys[ v[ 0 ] ] [ 0 ] *  barys[ v[ 1 ] ] [ 2 ];
        matrix[ 2 ] [ 2 ]  = barys[ v[ 0 ] ] [ 0 ] *  barys[ v[ 1 ] ] [ 1 ] - barys[ v[ 0 ] ] [ 1 ] *  barys[ v[ 1 ] ] [ 0 ];
      }
    };

    struct MatrixIF
    {
      virtual bool apply( const KeyType& v,
                          const std::vector< DomainType >& barys,
                          const std::vector< RangeType >& nbVals,
                          GradientType& dM ) = 0;
      virtual MatrixIF* clone() const = 0;

      virtual ~MatrixIF () {}

      static constexpr double detEps_ = 1e-12 ;
    };

    struct RegularMatrix : public MatrixIF
    {
      MatrixType inverse_ ;
      bool inverseCalculated_ ;

      RegularMatrix() : inverseCalculated_( false ) {}

      MatrixIF* clone() const { return new RegularMatrix( *this ); }

      bool inverse(const KeyType& v, const std::vector< DomainType >& barys )
      {
        if( ! inverseCalculated_ )
        {
          MatrixType matrix;
          // setup matrix
          MatrixAssemblerForLinearReconstruction<dimGrid,dimDomain, DomainFieldType>
            :: assembleMatrix(v, barys, matrix);
          // invert matrix
          RangeFieldType det = Dune :: FMatrixHelp :: invertMatrix( matrix, inverse_ );
          if( std::abs( det ) > MatrixIF :: detEps_ )
          {
            inverseCalculated_ = true ;
          }
        }
        return inverseCalculated_ ;
      }

      bool apply( const KeyType& v,
                  const std::vector< DomainType >& barys,
                  const std::vector< RangeType >& nbVals,
                  GradientType& dM )
      {
        // if matrix is regular
        if( inverse( v, barys ) )
        {
          DomainType rhs ;
          const int vSize = v.size();
          // calculate D
          for(int r=0; r<dimRange; ++r)
          {
            for(int i=0; i<vSize; ++i)
            {
              rhs[ i ] = nbVals[ v[ i ] ][ r ];
            }

            //if we have surface grid, then we need additional row to make the matrix quadratic. Claim that
            //derivative in normal (to entity) direction is zero
            // whats with 1,2 ???
            if( (dimGrid == 2) && (dimDomain == 3) )
            {
              rhs[vSize] = 0;
            }

            // get solution
            inverse_.mv( rhs, dM[r] );
          }
          return true;
        }
        return false;
      }
    };

    struct LeastSquaresMatrix : public MatrixIF
    {
      // dimension
      enum { dim = dimGrid };
      // new dimension is dim + 1
      enum { newDim = dim + 1 };

      // apply least square by adding another point
      // this should make the linear system solvable

      // need new matrix type containing one row more
      typedef Dune::FieldMatrix<DomainFieldType, newDim , dimDomain> NewMatrixType;
      typedef Dune::FieldVector<DomainFieldType, newDim > NewVectorType;

      MatrixType inverse_ ;
      // new matrix
      NewMatrixType A_ ;

      bool inverseCalculated_ ;

      LeastSquaresMatrix() : inverseCalculated_( false ) {}

      MatrixIF* clone() const { return new LeastSquaresMatrix( *this ); }

      bool inverse(const KeyType& nV, const std::vector< DomainType >& barys )
      {
        if( ! inverseCalculated_ )
        {
          assert( (int) nV.size() == newDim );

          // create matrix
          for(int k=0; k<newDim; ++k)
          {
            A_[k] = barys[ nV[k] ];
          }

          MatrixType matrix ;

          // matrix = A^T * A
          multiply_AT_A(A_, matrix);

            // invert matrix
          RangeFieldType det = Dune :: FMatrixHelp :: invertMatrix(matrix, inverse_ );

          if( std::abs( det ) > MatrixIF :: detEps_ )
          {
            inverseCalculated_ = true ;
          }
        }
        return inverseCalculated_ ;
      }

      bool apply( const KeyType& nV,
                  const std::vector< DomainType >& barys,
                  const std::vector< RangeType >& nbVals,
                  GradientType& dM )
      {
        if( inverse( nV, barys ) )
        {
          // need new right hand side
          NewVectorType newRhs;
          DomainType rhs ;

          // calculate D
          for(int r=0; r<dimRange; ++r)
          {
            // get right hand side
            for(int i=0; i<newDim; ++i)
            {
              newRhs[i] = nbVals[ nV[i] ][r];
            }

            // convert newRhs to matrix
            A_.mtv(newRhs, rhs);

            // get solution
            inverse_.mv( rhs, dM[r] );
          }
          return true ;
        }
        return false ;
      }

      // matrix = A^T * A
      template <class NewMatrixType, class MatrixType>
      void multiply_AT_A(const NewMatrixType& A, MatrixType& matrix) const
      {
        assert( (int) MatrixType :: rows == (int) NewMatrixType :: cols );

        for(int row=0; row< MatrixType :: rows; ++row)
        {
          for(int col=0; col< MatrixType :: cols; ++col)
          {
            matrix[row][col] = 0;
            for(int k=0; k<NewMatrixType :: rows;  ++k)
            {
              matrix[row][col] += A[k][row] * A[k][col];
            }
          }
        }
      }

    };

    class MatrixStorage
    {
    protected:
      MatrixIF* matrix_;
    public:
      MatrixStorage() : matrix_( 0 ) {}
      explicit MatrixStorage( const MatrixIF& matrix )
       :  matrix_( matrix.clone() )
      {}

      MatrixStorage( const MatrixStorage& other ) : matrix_( 0 )
      {
        assign( other );
      }

      MatrixStorage& operator = ( const MatrixStorage& other )
      {
        removeObj();
        assign( other );
        return *this;
      }

      ~MatrixStorage() { removeObj(); }

      MatrixIF* matrix() { assert( matrix_ ); return matrix_ ; }

    protected:
      void removeObj()
      {
        delete matrix_; matrix_ = 0;
      }

      void assign( const MatrixStorage& other )
      {
        matrix_ = ( other.matrix_ ) ? (other.matrix_->clone() ) : 0 ;
      }
    };


    // get linear function from reconstruction of the average values
    template <class MatrixCacheType>
    static void calculateLinearFunctions(const ComboSetType& comboSet,
                                         const GeometryType& geomType,
                                         const Flags& flags,
                                         const std::vector< DomainType >& baryCenters,
                                         const std::vector< RangeType  >& neighborValues,
                                         MatrixCacheType& matrixCache,
                                         std::vector< GradientType >& deoMods,
                                         std::vector< CheckType >&  comboVec )
    {
      // initialize combo vecs
      const size_t comboSize = comboSet.size();
      deoMods.reserve( comboSize );
      comboVec.reserve( comboSize );

      deoMods.clear();
      comboVec.clear();

      // assert( !StructuredGrid || cartesian );
      // use matrix cache in case of structured grid
      const bool useCache = flags.cartesian
                            && ! flags.nonConforming
                            && ! flags.boundary;

      static const int dim = dimGrid;

      typedef typename ComboSetType :: iterator iterator;

      // calculate linear functions
      // D(x) = U_i + D_i * (x - w_i)
      const iterator endit = comboSet.end();
      for(iterator it = comboSet.begin(); it != endit; ++it)
      {
        // get tuple of numbers
        const KeyType& v = (*it).first;

        RegularMatrix regInverse ;

        MatrixIF* inverse = & regInverse ;

        if( useCache )
        {
          typedef typename MatrixCacheType :: iterator iterator ;
          iterator matrixEntry = matrixCache.find( v );
          if( matrixEntry != matrixCache.end() )
          {
            inverse = matrixEntry->second.matrix();
          }
          else
          {
            if( regInverse.inverse( v, baryCenters ) )
            {
              matrixCache[ v ] = MatrixStorage( regInverse );
            }
          }
        }

        // create new instance of limiter coefficients
        GradientType dM;

        // if applied is not true the inverse is singular
        const bool applied = inverse->apply( v, baryCenters, neighborValues, dM );

        if( applied )
        {
          // store linear function
          deoMods.push_back( dM );
          comboVec.push_back( (*it).second );
        }
        else
        {
          // apply least square by adding another point
          // this should make the linear system solvable

          // creare vector with size = dim+1
          // the first dim components are equal to v
          CheckType nV( dim+1 );
          for(int i=0; i<dim; ++i) nV[i] = v[i];

          // take first point of list of points to check
          CheckType check ( (*it).second );
          assert( check.size() > 0 );

          // get check iterator
          typedef typename CheckType :: iterator CheckIteratorType;
          const CheckIteratorType checkEnd = check.end();

          // take first entry as start value
          CheckIteratorType checkIt = check.begin();

          // matrix should be regular now
          if( checkIt == checkEnd )
          {
            //TODO unccomment this check!
            // should work for 1 or 2 otherwise error
            //DUNE_THROW(InvalidStateException,"Check vector has no entries!");
          }

          CheckType checkNums ( check );
          checkNums.reserve( checkNums.size() + dim );
          for(int i=0; i<dim; ++i)
          {
            checkNums.push_back( v[i] );
          }

          for( ; checkIt != checkEnd; ++ checkIt )
          {
            ////////////////////////////////////////////
            // apply least squares approach here
            ////////////////////////////////////////////
            LeastSquaresMatrix lsInverse ;

            // assign last element
            nV[dim] = *checkIt ;

            KeyType newV ( nV );

            bool matrixNotSingular = lsInverse.apply( newV, baryCenters, neighborValues, dM ) ;

            // if matrix was valid add to functions
            if( matrixNotSingular )
            {
              // store linear function
              deoMods.push_back( dM );

              // store vector with points to check
              comboVec.push_back( checkNums );
            }
          }
        }
      } // end solving
    }

  };


  inline std::ostream& operator << (std::ostream& out, const Dune::DGFEntityKey<int>& key )
  {
    out << "(";
    for(int i=0; i<key.size(); ++i )
      out << key[i] << " ";
    out << ")";
    return out;
  }

  // base class for Limiters, mainly reading limitEps
  struct LimiterFunctionBase
  {
    const double limitEps_;
    LimiterFunctionBase()
      : limitEps_( getEpsilon() )
    {}

    //! get tolerance for shock detector
    static double getEpsilon()
    {
      // default value is 1e-8
     // return Parameter::getValue("femdg.limiter.limiteps", double(1e-8) );
      return double(1e-8);
    }

    //! return epsilon for limiting
    double epsilon () const { return limitEps_; }

  protected:
    void printInfo(const std::string& name ) const
    {
      //TODO what with parameter file in ewoms? what would be the substitute?
      if( /*Parameter::verbose()*/ true )
      {
        std::cout << "LimiterFunction: " << name << " with limitEps = " << limitEps_ << std::endl;
      }
    }
  };

  template <class Field>
  struct MinModLimiter : public LimiterFunctionBase
  {
    using LimiterFunctionBase :: limitEps_ ;
    typedef Field FieldType;

    MinModLimiter() : LimiterFunctionBase()
    {
      printInfo( "minmod" );
    }

    FieldType operator ()(const FieldType& g, const FieldType& d ) const
    {
      const FieldType absG = std::abs( g );

      // avoid rounding errors
      // gradient of linear functions is nerly zero
      if ( absG < 1e-10 ) return 1;

      // if product smaller than zero set localFactor to zero
      // otherwise d/g until 1
      const FieldType localFactor =
            ( (d * g) < 0.0 ) ? 0.0 : std::min( 1.0 , ( d / g ) );

      return localFactor;
    }
  };

  template <class Field>
  struct SuperBeeLimiter : public LimiterFunctionBase
  {
    using LimiterFunctionBase :: limitEps_ ;
    typedef Field FieldType;

    SuperBeeLimiter() : LimiterFunctionBase()
    {
      printInfo( "superbee" );
    }

    FieldType operator ()(const FieldType& g, const FieldType& d ) const
    {
      // avoid rounding errors
      // gradient of linear functions is nerly zero
      if ( std::abs( g ) < limitEps_ ) return 1;

      const FieldType r = d / g ;

      // if product smaller than zero set localFactor to zero
      // otherwise d/g until 1
      const FieldType localFactor =
            ( (g * d) < 0.0 ) ? 0.0 :
            ( std::max( std::min( 2.0 * r , 1.0 ), std::min( r, 2.0 ) ) );

      return localFactor;
    }
  };

  template <class Field>
  struct VanLeerLimiter : public LimiterFunctionBase
  {
    using LimiterFunctionBase :: limitEps_ ;
    typedef Field FieldType;

    VanLeerLimiter() : LimiterFunctionBase()
    {
      printInfo( "vanLeer" );
    }

    FieldType operator ()(const FieldType& g, const FieldType& d ) const
    {
      const FieldType absG = std::abs( g );
      const FieldType absD = std::abs( d );

      // avoid rounding errors
      if ( (absG < limitEps_) && (absD < limitEps_) ) return 1;
      if ( absG < limitEps_ ) return 1;

      const FieldType r = d / g ;
      const FieldType absR = std::abs( r );

      return ( absR + r ) / ( 1.0 + absR );
    }
  };

//} // end namespace Fem
} // end namespace Dune

#endif // DUNE_LIMITERUTILITY_HH
