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
 * \copydoc Opm::TransportProblem
 */
#ifndef EWOMS_TRANSPORT_PROBLEM_HH
#define EWOMS_TRANSPORT_PROBLEM_HH

#include <opm/models/io/dgfvanguard.hh>
#include <opm/models/immiscible/immiscibleproperties.hh>
#include <opm/models/discretization/common/fvbaseadlocallinearizer.hh>

#include <opm/material/fluidmatrixinteractions/PureTransport.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/components/Artificial_liquid.hpp>
#include <opm/material/components/Dnapl.hpp>
//#include <opm/material/common/Unused.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <sstream>
#include <string>
#include <iostream>


#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

#if HAVE_DUNE_POLYGONGRID
#include <dune/polygongrid/grid.hh>
#include <dune/polygongrid/dgf.hh>
#endif

#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>

#include <opm/grid/CpGrid.hpp>


#define USE_DUNE_POLYGONGRID 0
#define USE_DUNE_CPGRID 0
#define USE_DUNE_ALUGRID 0
#define USE_DUNE_YASPGRID_3D 0

namespace Opm {
template <class TypeTag>
class TransportProblem;
}

BEGIN_PROPERTIES

#if USE_DUNE_CPGRID
NEW_TYPE_TAG(TransportBaseProblem, INHERITS_FROM(EclCpGridVanguard));
#else
NEW_TYPE_TAG(TransportBaseProblem);
#endif

// declare the properties specific for the lens problem
NEW_PROP_TAG(TransportLowerLeftX);
NEW_PROP_TAG(TransportLowerLeftY);
NEW_PROP_TAG(TransportLowerLeftZ);
NEW_PROP_TAG(TransportUpperRightX);
NEW_PROP_TAG(TransportUpperRightY);
NEW_PROP_TAG(TransportUpperRightZ);

SET_TYPE_PROP(TransportBaseProblem, Problem, Opm::TransportProblem<TypeTag>);
//SET_BOOL_PROP(TransportBaseProblem, EnableEclOutput, false);

#if USE_DUNE_ALUGRID
SET_TYPE_PROP(TransportBaseProblem, Grid, Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > );
        SET_STRING_PROP(TransportBaseProblem, GridFile, "../tests/data/groundwater_2d.dgf");
#elif USE_DUNE_POLYGONGRID
SET_TYPE_PROP(TransportBaseProblem, Grid, Dune::PolygonGrid< double > );
#elif USE_DUNE_CPGRID
 // already added
#elif USE_DUNE_YASPGRID_3D
        SET_TYPE_PROP(TransportBaseProblem, Grid, Dune::YaspGrid<3>);
        SET_STRING_PROP(TransportBaseProblem, GridFile, "../../tests/data/groundwater_3d.dgf");
#else
// Use Dune-grid's YaspGrid
SET_TYPE_PROP(TransportBaseProblem, Grid, Dune::YaspGrid<1>);
        SET_STRING_PROP(TransportBaseProblem, GridFile, "../../tests/data/groundwater_1d.dgf");
#endif


// Set the problem property

//SET_STRING_PROP(TransportBaseProblem, Grid, EclDeckFileName, "data/CPGRID_INCL.DATA");

//SET_TYPE_PROP(TransportBaseProblem, Grid, Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > );
//SET_TYPE_PROP(TransportBaseProblem, Grid, Dune::PolyhedralGrid< 3, 3 >);

// Set the GridManager property
//SET_TYPE_PROP(TransportBaseProblem, GridManager, Opm::DgfGridManager<TypeTag>);

// The default DGF file to load
//SET_STRING_PROP(TransportBaseProblem, GridFile, "../tests/data/outflow.dgf");
//SET_STRING_PROP(TransportBaseProblem, GridFile, "../tests/data/simplex2.dgf");
//SET_STRING_PROP(TransportBaseProblem, GridFile, "../tests/data/groundwater_3d.dgf");
//SET_STRING_PROP(TransportBaseProblem, GridFile, "../tests/data/dbls_20.msh.dgf");
// Set the wetting phase
SET_PROP(TransportBaseProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::Artificial_liquid<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TransportBaseProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::Artificial_liquid<Scalar> > type;
};

// Set the material Law
SET_PROP(TransportBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx> Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::PureTransport<Traits> EffectiveLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffectiveLaw> type;
};

// Write the solutions of individual newton iterations?
SET_BOOL_PROP(TransportBaseProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(TransportBaseProblem, NumericDifferenceMethod, +1);

// Enable gravity
SET_BOOL_PROP(TransportBaseProblem, EnableGravity, false);

// define the properties specific for the lens problem
SET_SCALAR_PROP(TransportBaseProblem, TransportLowerLeftX, 0.0);
SET_SCALAR_PROP(TransportBaseProblem, TransportLowerLeftY, 0.0);
SET_SCALAR_PROP(TransportBaseProblem, TransportLowerLeftZ, 0.0);
SET_SCALAR_PROP(TransportBaseProblem, TransportUpperRightX, 1.0);
SET_SCALAR_PROP(TransportBaseProblem, TransportUpperRightY, 1.0);
SET_SCALAR_PROP(TransportBaseProblem, TransportUpperRightZ, 1.0);

//SET_SCALAR_PROP(TransportBaseProblem, DomainSizeX, 1.0);
//SET_SCALAR_PROP(TransportBaseProblem, DomainSizeY, 1.0);
//SET_SCALAR_PROP(TransportBaseProblem, DomainSizeZ, 1.0);

//SET_INT_PROP(TransportBaseProblem, CellsX, 50);
//SET_INT_PROP(TransportBaseProblem, CellsY, 50);
//SET_INT_PROP(TransportBaseProblem, CellsZ, 50);

// The default for the end time of the simulation
SET_SCALAR_PROP(TransportBaseProblem, EndTime, 8.772e3);
int a = 1000;
// The default for the initial time step size of the simulation
SET_SCALAR_PROP(TransportBaseProblem, InitialTimeStepSize, a);
//SET_SCALAR_PROP(TransportBaseProblem, MaxTimeStepSize, a);
//SET_SCALAR_PROP(TransportBaseProblem, MinTimeStepSize, a);

// By default, include the intrinsic permeability tensor to the VTK output files
SET_BOOL_PROP(TransportBaseProblem, VtkWriteIntrinsicPermeabilities, true);

// enable the storage cache by default for this problem
SET_BOOL_PROP(TransportBaseProblem, EnableStorageCache, false);

// enable the cache for intensive quantities by default for this problem
SET_BOOL_PROP(TransportBaseProblem, EnableIntensiveQuantityCache, false);

END_PROPERTIES


namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the model, the depth of the domain is implicitly
 * assumed to be 1 m everywhere.
 *
 * On the top and the bottom of the domain no-flow boundary conditions
 * are used, while free-flow conditions apply on the left and right
 * boundaries; DNAPL is injected at the top boundary from 3m to 4m at
 * a rate of 0.04 kg/(s m^2).
 *
 * At the boundary on the left, a free-flow condition using the
 * hydrostatic pressure scaled by a factor of 1.125 is imposed, while
 * on the right, it is just the hydrostatic pressure. The DNAPL
 * saturation on both sides is zero.
 */
template <class TypeTag>
class TransportProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nonWettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    TransportProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-9;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C
 /*       lensLowerLeft_[0] = EWOMS_GET_PARAM(TypeTag, Scalar, TransportLowerLeftX);
        lensLowerLeft_[1] = EWOMS_GET_PARAM(TypeTag, Scalar, TransportLowerLeftY);
        lensUpperRight_[0] = EWOMS_GET_PARAM(TypeTag, Scalar, TransportUpperRightX);
        lensUpperRight_[1] = EWOMS_GET_PARAM(TypeTag, Scalar, TransportUpperRightY);

        if (dimWorld == 3) {
            lensLowerLeft_[2] = EWOMS_GET_PARAM(TypeTag, Scalar, TransportLowerLeftZ);
            lensUpperRight_[2] = EWOMS_GET_PARAM(TypeTag, Scalar, TransportUpperRightZ);
        }
*/
        // residual saturations
        //lensMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.18);
        //lensMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);
        outerMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.05);
        outerMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);

        // parameters for the Van Genuchten law: alpha and n
        //lensMaterialParams_.setVgAlpha(0.00045);
        //lensMaterialParams_.setVgN(7.3);
        //outerMaterialParams_.setVgAlpha(0.0037);
        //outerMaterialParams_.setVgN(4.7);

        //lensMaterialParams_.finalize();
        outerMaterialParams_.finalize();

    //    lensK_ = this->toDimMatrix_(9.05e-12);
        outerK_ = this->toDimMatrix_(4.6e-10);

    /*    if (dimWorld == 3) {
            this->gravity_ = 0;
            this->gravity_[1] = -9.81;
        }*/
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, TransportLowerLeftX,
                             "The x-coordinate of the lens' lower-left corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, TransportLowerLeftY,
                             "The y-coordinate of the lens' lower-left corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, TransportUpperRightX,
                             "The x-coordinate of the lens' upper-right corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, TransportUpperRightY,
                             "The y-coordinate of the lens' upper-right corner "
                             "[m].");

        if (dimWorld == 3) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, TransportLowerLeftZ,
                                 "The z-coordinate of the lens' lower-left "
                                 "corner [m].");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, TransportUpperRightZ,
                                 "The z-coordinate of the lens' upper-right "
                                 "corner [m].");
        }
    }

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        //const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);

        return outerK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return 0.4; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        //const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);

        return outerMaterialParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    //! \}

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        typedef typename GET_PROP_TYPE(TypeTag, LocalLinearizerSplice) LLS;

        bool useAutoDiff = std::is_same<LLS, TTAG(AutoDiffLocalLinearizer)>::value;

        std::ostringstream oss;
        oss << "Transport_" << Model::name()
            << "_" << Model::discretizationName()
            << "_" << (useAutoDiff?"ad":"fd");
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::beginTimeStep
     */
    void beginTimeStep()
    { }

    /*!
     * \copydoc FvBaseProblem::beginIteration
     */
    void beginIteration()
    { }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        //this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl << std::flush;
        }
#endif // NDEBUG
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
           // Sw=1, wetting phase enters from the left boundary

            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;
            Sw = 1.0;

            Scalar p0, p1;
            p0 = 1e5;
            p1 = 5e4;
            pw = p0;

            // specify a full fluid state using pw and Sw
           // const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);

            // free flow boundary. we assume incompressible fluids
            Scalar densityW = WettingPhase::density(temperature_, Scalar(pw));
            Scalar densityN = NonwettingPhase::density(temperature_, Scalar(pw));

            Opm::ImmiscibleFluidState<Scalar, FluidSystem,
                    /*storeEnthalpy=*/false> fs;
            fs.setSaturation(wettingPhaseIdx, Sw);
            fs.setSaturation(nonWettingPhaseIdx, 1.0 - Sw);
            fs.setTemperature(T);

            //Scalar pC[numPhases];
            //MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wettingPhaseIdx, pw);
            fs.setPressure(nonWettingPhaseIdx, pw); // + pC[nonWettingPhaseIdx] - pC[wettingPhaseIdx]);

            fs.setDensity(wettingPhaseIdx, densityW);
            fs.setDensity(nonWettingPhaseIdx, densityN);

            fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(temperature_, fs.pressure(wettingPhaseIdx)));
            fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(temperature_, fs.pressure(nonWettingPhaseIdx)));

            // impose an freeflow boundary condition
            values.setInFlow(context, spaceIdx, timeIdx, fs);

        } else
        if (onRightBoundary_(pos)) {
            //lower pressure p1, Sw=0.1 same as in initial state
            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;
            Sw = 0.1;

            Scalar p0, p1;
            p0 = 1e5;
            p1 = 5e4;
            pw = p1;

            // specify a full fluid state using pw and Sw
            //const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);

            // free flow boundary. we assume incompressible fluids
            Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/Scalar(pw));
            Scalar densityN = NonwettingPhase::density(temperature_, /*pressure=*/Scalar(pw));

            Opm::ImmiscibleFluidState<Scalar, FluidSystem,
            /*storeEnthalpy=*/false> fs;
            fs.setSaturation(wettingPhaseIdx, Sw);
            fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);
            fs.setTemperature(T);

            //Scalar pC[numPhases];
            //MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wettingPhaseIdx, pw);
            fs.setPressure(nonWettingPhaseIdx, pw); // + pC[nonWettingPhaseIdx] - pC[wettingPhaseIdx]);

            fs.setDensity(wettingPhaseIdx, densityW);
            fs.setDensity(nonWettingPhaseIdx, densityN);

            fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(temperature_, fs.pressure(wettingPhaseIdx)));
            fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(temperature_, fs.pressure(nonWettingPhaseIdx)));

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else {
            // no flow boundary
            values.setNoFlow();
        }

    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        //Scalar depth = this->boundingBoxMax()[1] - pos[1];

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

{
                Scalar p0, p1, pw;
                p0 = 1e5;
                p1 = 5e4;
                pw = p0 - (p0-p1)*pos[0];
                fs.setPressure(wettingPhaseIdx, /*pressure=*/pw);
                fs.setPressure(nonWettingPhaseIdx, /*pressure=*/pw);

                fs.setTemperature(temperature_);

                Scalar Sw;
                Sw = 0.1;
                fs.setSaturation(wettingPhaseIdx, Sw);
                fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updatePhase(fs, wettingPhaseIdx);
            Scalar densityW = FluidSystem::density(fs, paramCache, wettingPhaseIdx);
            Scalar densityN = FluidSystem::density(fs, paramCache, nonWettingPhaseIdx);

    fs.setDensity(wettingPhaseIdx, densityW);
    fs.setDensity(nonWettingPhaseIdx, densityN);

    fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(temperature_, fs.pressure(wettingPhaseIdx)));
    fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(temperature_, fs.pressure(nonWettingPhaseIdx)));

//            // hydrostatic pressure (assuming incompressibility)
//            Scalar pw = 1e5;// - densityW * this->gravity()[1] * depth;
//
//            // calculate the capillary pressure
//            const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
//            Scalar pC[numPhases];
//            MaterialLaw::capillaryPressures(pC, matParams, fs);
//
//            // make a full fluid state
//            fs.setPressure(wettingPhaseIdx, pw);
//            fs.setPressure(nonWettingPhaseIdx, pw + (pC[wettingPhaseIdx] - pC[nonWettingPhaseIdx]));

            // assign the primary variables
            values.assignNaive(fs);
        }


    }

    template <class Context>
    Scalar exactWetPressure (const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        Scalar p0, p1, pw;
        p0 = 1e5;
        p1 = 5e4;
        pw = p0 - (p0-p1)*pos[0];

        return pw;

    }

    template <class Context>
    Scalar exactWetSat (const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        //Scalar end_time = 8.772e3;

        Scalar time = this->simulator().time();
        Scalar oldDt = this->simulator().timeStepSize();
        time += oldDt;

        Scalar frontPosition = velocity_ * time; //end_time;

        Scalar exactSw;


        if (pos[0] < frontPosition + 1e-6 )
            exactSw = 1.0;
        else
            exactSw = 0.1;

        return exactSw;
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = Scalar(0.0); }

    //! \}

private:

    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onUpperRightCorner_ (const GlobalPosition& pos) const
    {
        return ((pos[0] > this->boundingBoxMax()[0] - well_radius - eps_) && (pos[1] > this->boundingBoxMax()[1] - well_radius - eps_) && (pos[2] > this->boundingBoxMax()[2] - well_radius - eps_));
    }

    bool onLowerLeftCorner_ (const GlobalPosition pos) const
    {
        return ((pos[1] < this->boundingBoxMin()[1] + well_radius + eps_) && (pos[0] < this->boundingBoxMin()[0] + well_radius + eps_) && (pos[2] < this->boundingBoxMin()[2] + well_radius + eps_));
    }

    bool onInlet_(const GlobalPosition& pos) const
    {
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar lambda = (this->boundingBoxMax()[0] - pos[0]) / width;
        return onUpperBoundary_(pos) && 0.5 < lambda && lambda < 2.0 / 3.0;
    }

bool inSquare_(const GlobalPosition& pos) const
{
    return ((pos[0] > 0.3 && pos[0] < 0.5) && (pos[1] > 0.3 && pos[1] < 0.5));
}


    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    double well_radius = 0.02;

    DimMatrix lensK_;
    DimMatrix outerK_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;

    Scalar temperature_;
    Scalar eps_;
    Scalar velocity_ = 5.75e-5;
};

} // namespace Opm

#endif
