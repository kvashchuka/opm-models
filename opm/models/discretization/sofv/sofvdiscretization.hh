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
 * \copydoc Opm::SofvDiscretization
 */
#ifndef EWOMS_SOFV_DISCRETIZATION_HH
#define EWOMS_SOFV_DISCRETIZATION_HH

#include <opm/material/densead/Math.hpp>

#include "sofvproperties.hh"
#include "sofvstencil.hh"
#include "sofvgridcommhandlefactory.hh"
#include "sofvbaseoutputmodule.hh"

#include <opm/simulators/linalg/elementborderlistfromgrid.hh>
#include <opm/models/discretization/common/fvbasediscretization.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>

#include <opm/models/io/vtkblackoilsolventmodule.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/finitevolume.hh>
#include "reconstruction.hh"
#include "limiterutility.hh"

namespace Opm {
template <class TypeTag>
class SofvDiscretization;
}

BEGIN_PROPERTIES

//! Set the stencil
SET_PROP(SofvDiscretization, Stencil)
{
private:

public:
    typedef Opm::SofvStencil<TypeTag> type;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(SofvDiscretization, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(SofvDiscretization, Discretization, Opm::SofvDiscretization<TypeTag>);

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
SET_TYPE_PROP(SofvDiscretization, DiscBaseOutputModule,
              Opm::SofvBaseOutputModule<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(SofvDiscretization, GridCommHandleFactory,
              Opm::SofvGridCommHandleFactory<TypeTag>);

//! Set the DiscreteFunctionSpace
/*
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
*/

//SET_BOOL_PROP(SofvDiscretization, EnableHigherOrder, false);

//! Set the border list creator for to the one of an element based
//! method
SET_PROP(SofvDiscretization, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Opm::Linear::ElementBorderListFromGrid<GridView, ElementMapper> type;
};

//! For the element centered finite volume method, ghost and overlap elements must be
//! assembled to calculate the fluxes over the process boundary faces of the local
//! process' grid partition
SET_BOOL_PROP(SofvDiscretization, LinearizeNonLocalElements, true);

//! locking is not required for the element centered finite volume method because race
//! conditions cannot occur since each matrix/vector entry is written exactly once
SET_BOOL_PROP(SofvDiscretization, UseLinearizationLock, false);


SET_BOOL_PROP(SofvDiscretization, EnableSolvent, false);
SET_BOOL_PROP(SofvDiscretization, EnablePolymer, false);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup SofvDiscretization
 *
 * \brief The base class for the element-centered finite-volume discretization scheme.
 */
template<class TypeTag>
class SofvDiscretization : public FvBaseDiscretization<TypeTag>
{
    typedef FvBaseDiscretization<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) GridType;
    //typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPartType;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, DiscreteFunctionSpace) DiscreteFunctionSpaceType;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;

    //enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)         Scalar;


    typedef typename Opm::MathToolbox<Evaluation> Toolbox;
    const bool localRecons = EWOMS_GET_PARAM(TypeTag, bool, EnableLocalReconstruction);
    const int  schemeId_   = EWOMS_GET_PARAM(TypeTag, int , ReconstructionSchemeId);

    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange = PrimaryVariables::dimension}; //  };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    typedef Dune::FieldVector< Evaluation, dimRange >      EvalVector;
    typedef Dune::FieldVector< Scalar, dimRange >      Vector;

    // intersection iterator type
    typedef typename GridView::IntersectionIterator IntersectionIteratorType;
    // intersection type
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    // geometry of intersection
    typedef typename IntersectionType::Geometry IntersectionGeometryType;

    typedef Dune::Fem::LimitedReconstruction< TypeTag > ReconstructionType;

    typedef typename ReconstructionType :: DomainType         DomainType;
    typedef typename ReconstructionType :: DomainFieldType    DomainFieldType;
    typedef typename ReconstructionType :: RangeType          RangeType;
    typedef typename ReconstructionType :: RangeFieldType     RangeFieldType;

    typedef typename ReconstructionType::LocalFunctionType  ReconstructedLocalFunctionType;

    typedef typename GridView::template Codim<0>::Entity EntityType;
    // geometry type
    typedef typename EntityType::Geometry GeometryType;
    // global coordinates
    typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;
    // local coordinates
    typedef typename GeometryType::LocalCoordinate LocalCoordinateType;


public:
    using ParentType::simulator_;
    using ParentType::enableHigherOrder;
    using ParentType::reconstructOnlySolventOrPolymer;

    SofvDiscretization(Simulator& simulator)
        : ParentType(simulator),
          reconstruction_(simulator.vanguard().gridPart(), dofMapper())
    {
        localTotalMobility_.resize( ThreadManager::maxThreads() );
    }

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
        updateReconstruction_();
    }

    void computeError()
    {
        ElementContext elemCtx(simulator_);
        //size_t numDof = asImp_().numGridDof();
        unsigned timeIdx = 0;
        Scalar time = simulator_.time() + simulator_.timeStepSize();
        //PrimaryVariables exactSol;// = asImp_().solution(timeIdx);
        SolutionVector& solution = asImp_().solution(timeIdx);

        double errorL1 = 0.0;
        double errorL1_Pressure = 0.0;
        double errorL2 = 0.0;
        auto count = 0;
        int total_count = 0;
        auto total_volume = 0.0;
        auto checkSatisPhysical = true;

        if (checkSatisPhysical) {
            const auto & gridView = simulator_.gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;

                // deal with the current element
                elemCtx.updatePrimaryStencil(elem);

                for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(timeIdx); dofIdx++) {
                    unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);

                    if (solution[globalIdx][1] > 1.0 || solution[globalIdx][1] < 0.0){
                        std::cout << "Unphysical solution detected! Value = " << solution[globalIdx][1] << ", global Index " << globalIdx;
                        std::cout << ", time " << time << std::endl;
                        total_count ++;
                    }

                }
            }
        }

        if (std::abs(time - simulator_.endTime()) < 1.0) {
            // iterate through the grid and evaluate error
            const auto & gridView = simulator_.gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;

                // deal with the current element
                elemCtx.updatePrimaryStencil(elem);

                Scalar exactWetSaturation = 0.0;
                Scalar exactWetPressure = 0.0;

                for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(timeIdx); dofIdx++) {
                    unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
                    exactWetSaturation += simulator_.problem().exactWetSat(elemCtx, dofIdx, time);
                    exactWetPressure += simulator_.problem().exactWetPressure(elemCtx, dofIdx, time);
                    Scalar volume = elemCtx.stencil(timeIdx).subControlVolume(dofIdx).volume(); //elemCtx.dofVolume(dofIdx, timeIdx);
                    total_volume += volume;
                    //take only saturation, ignore pressure

                    errorL1 += std::abs(exactWetSaturation - solution[globalIdx][1]) * volume;
                    errorL1_Pressure += std::abs(exactWetPressure - solution[globalIdx][0]) * volume;
                    errorL2 += std::abs(exactWetSaturation - solution[globalIdx][1]) * std::abs(exactWetSaturation - solution[globalIdx][1]) * volume;

                    if (solution[globalIdx][1] > 1.0 || solution[globalIdx][1] < 0.0)
                        count++;

                }
            }

            std::cout << " Error-L1 = " << errorL1 << std::endl;
            std::cout << " Error-L2 = " << std::sqrt(errorL2) << std::endl;
            std::cout << " Total number of unphysical saturation solutions obtained during the whole simulation " << total_count << std::endl;
            std::cout << " Number of unphysical saturation solutions on the final time step " << count << std::endl;
            std::cout << " Total volume " << total_volume << std::endl;
            std::cout << " pressure_error = " << (errorL1_Pressure*100)/1e5 << std::endl;
        }
    }

    void updateReconstruction_ ( )
    {
        if(localRecons)
            return;

        if( enableHigherOrder() )
        {
            //std::cout << "start update reconstruction" << std::endl;
            ElementContext elemCtx(simulator_);


            size_t numDof = asImp_().numGridDof();

            elementDofInfo_.resize( numDof );

            typedef Dune::FieldVector< Evaluation, dimRange >      Vector;
            std::vector<Vector> totalMobility(numDof);

            // iterate through the grid and evaluate the initial condition
            const auto & gridView = simulator_.gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;

                // deal with the current element
                // also compute stencil for setting up neighbor topo
                elemCtx.updateStencil(elem);
                elemCtx.updateIntensiveQuantities( 0 );

                // map the local degree of freedom index to the global one
                unsigned int elemIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);

                // setup neighbor topology
                elementDofInfo_[ elemIdx ].resize( elemCtx.numDof(/*timeIdx=*/0) );
                for (unsigned int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx)
                {
                    elementDofInfo_[ elemIdx ][ dofIdx ] = elemCtx.globalSpaceIndex( dofIdx, /*timeIdx=*/0);
                }

                // Only update primary dofs
                //totalMobility[globalIdx].resize(elemCtx.numDof(/*timeIdx=*/0));
                // loop over all element vertices, i.e. sub control volumes
                for (unsigned int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx)
                {
                    // map the local degree of freedom index to the global one
                    unsigned int globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                    //std::cout << "dofIdx = " << dofIdx << std::endl;

                    totalMobility[globalIdx] = 0;

                    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    {
                        if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                            continue;
                        }
                        //totalMobility[globalIdx][phaseIdx] = elemCtx.intensiveQuantities(dofIdx,  /*timeIdx=*/0).mobility(phaseIdx).value();
                        //std::cout << elemCtx.intensiveQuantities(dofIdx,0).mobility(phaseIdx) << " mob base"  << std::endl;
                        //totalMobility[globalIdx][phaseIdx] = elemCtx.intensiveQuantities(dofIdx,  /*timeIdx=*/0).mobility(phaseIdx).value();
                        // for the case when AD is switched OFF
                        //totalMobility[globalIdx][phaseIdx] = elemCtx.intensiveQuantities(dofIdx,/*timeIdx=*/0).mobility(phaseIdx);
                        const auto& mob = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).mobility(phaseIdx);
                        //int idx = phaseIdx*(numDeri+1) ;
                        totalMobility[globalIdx][ phaseIdx ] = mob ;
                        //++idx;
                        //for(int d=0; d<numDeri; ++d, ++idx )
                        //    totalMobility[globalIdx][ idx ] = mob.derivative( d );
                        //totalMobility[globalIdx][phaseIdx] = Toolbox::value(elemCtx.intensiveQuantities(dofIdx,  /*timeIdx=*/0).mobility(phaseIdx));
                    }
                    if (GET_PROP_VALUE(TypeTag, EnableSolvent)) {
                        const auto& mob = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).solventMobility();
                        //int idx = numPhases; //*(numDeri+1) ;
                        totalMobility[globalIdx][ numPhases ] = mob;
                        //++idx;
                        //for(int d=0; d<numDeri; ++d, ++idx )
                        //    totalMobility[globalIdx][ idx ] = mob.derivative( d );
                    }


                }
            }

            // compute linear reconstructions
            reconstruction_.update( totalMobility, simulator_, schemeId_ );
            //std::cout << "Function from sofvdiscretization updateBegin() got called " << std::endl;
            //std::cout << "end update reconstruction" << std::endl;
        }

    }

    void updateReconstructionLocal_ (const ElementContext& elemCtx, unsigned upstreamDofIdx, unsigned phaseIdx,
                                     const GlobalCoordinateType& point,
                                     Evaluation& reconstructedMobility ) const
    {
        if(!localRecons)
            return;

#warning hack
        const bool isSolvent = (phaseIdx == FluidSystem::numPhases) && GET_PROP_VALUE(TypeTag, EnableSolvent);
        const bool isPolymer = (phaseIdx == FluidSystem::numPhases) && GET_PROP_VALUE(TypeTag, EnablePolymer);

        // if we have more then one phase and
        const bool onlySolvent = reconstructOnlySolventOrPolymer();
        if( enableHigherOrder() )
        {
            size_t timeIdx = 0;
            const auto& stencil = elemCtx.stencil(timeIdx);
            const size_t numDof = stencil.numDof();
            auto& totalMobility = localTotalMobility_[ ThreadManager::threadId() ];

            totalMobility.resize( numDof );
            // loop over all element vertices, i.e. sub control volumes
            for (unsigned int dofIdx = 0; dofIdx < numDof; ++dofIdx)
            {
                if (isSolvent) {
                    const auto& mob = elemCtx.intensiveQuantities(dofIdx, timeIdx).solventMobility();
                    totalMobility[dofIdx] = mob;
                }
                else if (isPolymer)
                {
                   const auto& c = elemCtx.intensiveQuantities(dofIdx, timeIdx).polymerConcentration();
                   totalMobility[dofIdx] = c;
                }

                else {
                    const auto& mob = elemCtx.intensiveQuantities(dofIdx, timeIdx).mobility(phaseIdx);
                    totalMobility[dofIdx] = mob;
                }
            }

            if (onlySolvent && !(isSolvent || isPolymer) ) {
                if (elemCtx.focusDofIndex() == upstreamDofIdx )
                    reconstructedMobility = totalMobility[ upstreamDofIdx ];
                else
                    reconstructedMobility = Opm::scalarValue(totalMobility[upstreamDofIdx]);
            }
            else
            {
                // compute linear reconstructions
                reconstruction_.updateLocal(elemCtx, upstreamDofIdx, timeIdx, totalMobility, point, reconstructedMobility, schemeId_ );
            }
        }
    }

    void evaluateReconstruction(const ElementContext& elemCtx,
                                const int upstreamDofIdx,
                                const int phaseIdx,
                                const GlobalCoordinateType& point,
                                Evaluation& upResult) const
    {
        if( enableHigherOrder() )
        {
            if( localRecons )
            {
                updateReconstructionLocal_(elemCtx, upstreamDofIdx, phaseIdx, point, upResult);
            }
            else
            {
                unsigned int focusIdx = elemCtx.focusDofIndex();
                unsigned int globalFocusIdx = elemCtx.globalSpaceIndex( focusIdx, /*timeIdx=*/0);

                // TODO: This information can be pre-computed
                //const auto& upstreamElement = elemCtx.stencil( 0 ).element( upstreamDofIdx );
                unsigned int upstreamGlobalIdx = elemCtx.globalSpaceIndex( upstreamDofIdx, /*timeIdx=*/0 );

                const auto& upstreamInfo = elementDofInfo_[ upstreamGlobalIdx ];

                // default for newFocusIdx is the old focusIdx
                int newFocusIdx = -1; // TODO: Tor Harald please check! focusIdx;
                const unsigned int size = upstreamInfo.size();
                for( unsigned int dof = 0; dof < size; ++dof )
                {
                    if( globalFocusIdx == upstreamInfo[ dof ] )
                    {
                        newFocusIdx = dof;
                        break ;
                    }
                }

              //#warning is this sufficient?
                //int newFocusIdx = focusIdx;
                //if (upstreamDofIdx > 0 )
                //    newFocusIdx = 0;

                reconstruction_.evaluate( upstreamGlobalIdx, upstreamDofIdx, focusIdx, newFocusIdx, phaseIdx, point, upResult );

//                if( upstreamDofIdx > 0 )
//                {
//                  for( int i=0; i<RangeType::dimension; ++i )
//                  {
//                   interiorResult[ i ] = Toolbox::createConstant( Toolbox::value(interiorResult[ i ]));
//                  }
//                }
            }

//            if (upstreamDofIdx == 0) {
                //reconstruction_.localFunction( elemCtx.element(), upstreamDofIdx).evaluateGlobal( point, upResult );
                //            } else {
//                const auto& elem = elemCtx.stencil( 0 ).element( upstreamDofIdx );
//               ElementContext elemCtxUp(simulator_);
//                elemCtxUp.updateStencil(elem);
//                unsigned focusDofIdx = elemCtx.focusDofIndex();
//                elemCtxUp.setFocusDofIndex(focusDofIdx);
//                elemCtxUp.updateAllIntensiveQuantities();
//                updateReconstructionLocal_(elemCtxUp, 0, phaseIdx);
//                reconstruction_.localFunction( elemCtxUp.element(), 0).evaluateGlobal( point, upResult );
//            }
        }
        else
        {
            std::abort();
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

    double eps_ = 1e-6;

    mutable std::vector< std::vector< Evaluation > > localTotalMobility_;

    std::vector< std::vector< unsigned int > > elementDofInfo_;

};

} // namespace Opm

#else // #if HAVE_DUNE_FEM

// fallback to EcfvDiscretization in case dune-fem is not available
#define SofvDiscretization EcfvDiscretization

//namespace Opm {
//    template <class TypeTag>
//    using SofvDiscretization = EcfvDiscretization< TypeTag >;
//} // namespace Opm

#endif // #if HAVE_DUNE_FEM


#endif
