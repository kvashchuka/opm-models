#ifndef DUNE_OPTIM_CONSTRAINT_LINEAR_HH
#define DUNE_OPTIM_CONSTRAINT_LINEAR_HH

namespace Dune
{

  namespace Optim
  {

    // LinearConstraint
    // ----------------

    template< class V >
    class LinearConstraint
    {
      typedef LinearConstraint< V > This;

    public:
      typedef V Vector;

      typedef typename Vector::value_type Field;

      typedef Field ctype;
      typedef Vector NormalType;

      LinearConstraint () = default;

      explicit LinearConstraint ( const Vector &normal ) : normal_( normal ) {}

      template< class X >
      ctype evaluate ( const X &x ) const
      {
        return x*normal() - rhs();
      }

      const Vector &normal () const { return normal_; }
      Vector &normal () { return normal_; }

      const ctype &rhs () const { return rhs_; }
      ctype &rhs () { return rhs_; }

      unsigned int size () const { return normal_.size(); }

    private:
      Vector normal_;
      Field rhs_;
    };

    // LinearConstraint
    // ----------------

    template< class V, class Field >
    class LinearConstraintEval
    {
      typedef LinearConstraintEval< V, Field > This;

    public:
      typedef V Vector;

      //typedef typename Vector::value_type Field;

      typedef Field ctype;
      typedef Vector NormalType;

      LinearConstraintEval () = default;

      explicit LinearConstraintEval ( const Vector &normal ) : normal_( normal ) {}

      template< class X >
      ctype evaluate ( const X &x ) const
      {
        return x*normal() - rhs();
      }

      const Vector &normal () const { return normal_; }
      Vector &normal () { return normal_; }

      const ctype &rhs () const { return rhs_; }
      ctype &rhs () { return rhs_; }

      unsigned int size () const { return normal_.size(); }

    private:
      Vector normal_;
      Field rhs_;
    };

  } // namespace Optim

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_CONSTRAINT_LINEAR_HH
