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
 * \copydoc Opm::ImmiscibleIntensiveQuantities
 */
#ifndef EWOMS_IMMISCIBLE_INTENSIVE_QUANTITIES_HH
#define EWOMS_IMMISCIBLE_INTENSIVE_QUANTITIES_HH

#include "immiscibleproperties.hh"

#include <opm/models/common/energymodule.hh>

#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {
/*!
 * \ingroup ImmiscibleModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume for the immiscible multi-phase model.
 */
template <class TypeTag>
class ImmiscibleIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public EnergyIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { dimWorld = GridView::dimensionworld };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldVector<Evaluation, numPhases> EvalPhaseVector;

    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;
    typedef Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy> EnergyIntensiveQuantities;
    typedef Opm::ImmiscibleFluidState<Evaluation, FluidSystem,
                                      /*storeEnthalpy=*/enableEnergy> FluidState;

public:
    ImmiscibleIntensiveQuantities()
    { }

    ImmiscibleIntensiveQuantities(const ImmiscibleIntensiveQuantities& other) = default;

    ImmiscibleIntensiveQuantities& operator=(const ImmiscibleIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        // material law parameters
        const auto& problem = elemCtx.problem();
        const typename MaterialLaw::Params& materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        Opm::Valgrind::CheckDefined(priVars);

        Evaluation sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            const Evaluation& Salpha = priVars.makeEvaluation(saturation0Idx + phaseIdx, timeIdx);
            fluidState_.setSaturation(phaseIdx, Salpha);
            sumSat += Salpha;
        }
        fluidState_.setSaturation(numPhases - 1, 1 - sumSat);

        EvalPhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(pC);

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(relativePermeability_);

        const Evaluation& p0 = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, p0 + (pC[phaseIdx] - pC[0]));

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updateAll(fluidState_);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            // compute and set the density
            const Evaluation& rho = FluidSystem::density(fluidState_, paramCache, phaseIdx);
            fluidState_.setDensity(phaseIdx, rho);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx]/mu;
        }

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // exact wetting saturation
        //exactWetSat_ = problem.exactWetSat(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the intrinsic permeability tensor a degree of freedom.
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \brief Returns the relative permeability of a given phase
     *        within the control volume.
     *
     * \copydetails Doxygen::phaseIdxParam
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \copydetails Doxygen::phaseIdxParam
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    const Evaluation& porosity() const
    { return porosity_; }

    //const Evaluation& exactWetSat() const
    //{ return exactWetSat_; }

    const Evaluation& solventMobility() const
    { throw std::runtime_error("solventMobility() called but solvents are disabled"); }

    const Evaluation& polymerConcentration() const
    { throw std::runtime_error("polymerConcentration() called but polymer is disabled"); }




protected:
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation exactWetSat_;
    DimMatrix intrinsicPerm_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Opm

#endif
