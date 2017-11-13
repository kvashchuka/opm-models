// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::SofvDiscretization
 */
#ifndef EWOMS_SOFV_DISCRETIZATION_HH
#define EWOMS_SOFV_DISCRETIZATION_HH

#include <opm/material/densead/Math.hpp>

#include "sofvproperties.hh"
#include "sofvstencil.hh"
#include "sofvgridcommhandlefactory.hh"
#include "sofvbaseoutputmodule.hh"

#include <ewoms/linear/elementborderlistfromgrid.hh>
#include <ewoms/disc/common/fvbasediscretization.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/finitevolume.hh>
#include "reconstruction.hh"
#include "limitermodel.hh"
#include "limiterutility.hh"
#endif

namespace Ewoms {
template <class TypeTag>
class SofvDiscretization;
}

namespace Ewoms {
namespace Properties {
//! Set the stencil
SET_PROP(SofvDiscretization, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Ewoms::SofvStencil<Scalar, GridView> type;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(SofvDiscretization, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(SofvDiscretization, Discretization, Ewoms::SofvDiscretization<TypeTag>);

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
SET_TYPE_PROP(SofvDiscretization, DiscBaseOutputModule,
              Ewoms::SofvBaseOutputModule<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(SofvDiscretization, GridCommHandleFactory,
              Ewoms::SofvGridCommHandleFactory<TypeTag>);

#if HAVE_DUNE_FEM
//! Set the DiscreteFunctionSpace
SET_PROP(SofvDiscretization, DiscreteFunctionSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPart;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::Fem::FunctionSpace<typename GridPart::GridType::ctype,
                                     Scalar,
                                     GridPart::GridType::dimensionworld,
                                     numEq> FunctionSpace;
public:
    typedef Dune::Fem::FiniteVolumeSpace< FunctionSpace, GridPart, 0 > type;
};
        //SET_BOOL_PROP(SofvDiscretization, higherOrder_, true);
        //TODO Maybe put it into property system
        const bool higherOrder_ = true;
#endif



//! Set the border list creator for to the one of an element based
//! method
SET_PROP(SofvDiscretization, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Ewoms::Linear::ElementBorderListFromGrid<GridView, ElementMapper> type;
};

//! For the element centered finite volume method, ghost and overlap elements must be
//! assembled to calculate the fluxes over the process boundary faces of the local
//! process' grid partition
SET_BOOL_PROP(SofvDiscretization, LinearizeNonLocalElements, true);

//! locking is not required for the element centered finite volume method because race
//! conditions cannot occur since each matrix/vector entry is written exactly once
SET_BOOL_PROP(SofvDiscretization, UseLinearizationLock, false);

} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
/*!
 * \ingroup SofvDiscretization
 *
 * \brief The base class for the element-centered finite-volume discretization scheme.
 */
template<class TypeTag>
class SofvDiscretization : public FvBaseDiscretization<TypeTag>
{
    typedef FvBaseDiscretization<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    typedef LimiterModel<TypeTag> LimiterModelType;
    typedef LimitedReconstruction< TypeTag > ReconstructionType;


public:
    SofvDiscretization(Simulator& simulator)
        : ParentType(simulator),
          limiterModel_(),
          reconstruction_(simulator.gridManager().gridPart(), limiterModel_, dofMapper())
    { }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return "sofv"; }

    /*!
 * \brief Called by the update() method before it tries to
 *        apply the newton method. This is primary a hook
 *        which the actual model can overload.
 */
    void updateBegin()
    {
        if( higherOrder_ )
        {
            // compute linear reconstructions
            reconstruction_.update( asImp_().solution(/*timeIdx=*/0), dofMapper() );
        }
    }

    /*!
     * \brief Returns the number of global degrees of freedom (DOFs) due to the grid
     */
    size_t numGridDof() const
    { return static_cast<size_t>(this->gridView_.size(/*codim=*/0)); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper& dofMapper() const
    { return this->elementMapper(); }

    /*!
     * \brief Syncronize the values of the primary variables on the
     *        degrees of freedom that overlap with the neighboring
     *        processes.
     *
     * For the Element Centered Finite Volume discretization, this
     * method retrieves the primary variables corresponding to
     * overlap/ghost elements from their respective master process.
     */
    void syncOverlap()
    {
        // syncronize the solution on the ghost and overlap elements
        typedef GridCommHandleGhostSync<PrimaryVariables,
                                        SolutionVector,
                                        DofMapper,
                                        /*commCodim=*/0> GhostSyncHandle;

        auto ghostSync = GhostSyncHandle(this->solution(/*timeIdx=*/0),
                                         asImp_().dofMapper());
        this->gridView().communicate(ghostSync,
                                     Dune::InteriorBorder_All_Interface,
                                     Dune::ForwardCommunication);
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter& res)
    { res.template serializeEntities</*codim=*/0>(asImp_(), this->gridView_); }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        res.template deserializeEntities</*codim=*/0>(asImp_(), this->gridView_);
        this->solution(/*timeIdx=*/1) = this->solution(/*timeIdx=*/0);
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
    mutable ReconstructionType reconstruction_;
    LimiterModelType limiterModel_;
    const bool higherOrder_ = true;
};
} // namespace Ewoms

#endif
