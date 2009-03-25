/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH
#define DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/2p2c/2p2ctraits.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

//// forward declaration of the 2p2c box model
//template<class ProblemT, class TwoPTwoCTraitsT>
//class TwoPTwoCBoxModel;


///////////////////////////////////////////////////////////////////////////
// TwoPTwoCBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief 2P-2C specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxJacobian for the 2P-2C twophase flow.
 */
template<class ProblemT,
         class BoxTraitsT,
         class TwoPTwoCTraitsT,
         class Implementation>
class TwoPTwoCBoxJacobianBase : public BoxJacobian<ProblemT,
                                                   BoxTraitsT,
                                                   Implementation >
{
protected:
//    friend class TwoPTwoCBoxModel<ProblemT, TwoPTwoCTraitsT>;

    typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                    BoxTraitsT,
                                    TwoPTwoCTraitsT,
                                    Implementation>              ThisType;
    typedef BoxJacobian<ProblemT, BoxTraitsT, Implementation>    ParentType;
    typedef ProblemT                                Problem;
    typedef typename Problem::DomainTraits          DomTraits;
    typedef BoxTraitsT                              BoxTraits;
    typedef TwoPTwoCTraitsT                         TwoPTwoCTraits;
    typedef Dune::CollectiveCommunication<Problem>  CollectiveCommunication;


    enum {
        dim              = DomTraits::dim,
        dimWorld         = DomTraits::dimWorld,

        numEq            = BoxTraits::numEq,
        numPhases        = TwoPTwoCTraits::numPhases,
        numComponents    = TwoPTwoCTraits::numComponents,

        pressureIdx      = TwoPTwoCTraits::pressureIdx,
        switchIdx        = TwoPTwoCTraits::switchIdx,

        wPhase           = TwoPTwoCTraits::wPhase,
        nPhase           = TwoPTwoCTraits::nPhase,

        wComp            = TwoPTwoCTraits::wComp,
        nComp            = TwoPTwoCTraits::nComp,

        wPhaseOnly       = TwoPTwoCTraits::wPhaseOnly,
        nPhaseOnly       = TwoPTwoCTraits::nPhaseOnly,
        bothPhases       = TwoPTwoCTraits::bothPhases
    };
    static const int formulation  = TwoPTwoCTraits::formulation;
    enum {
        pWsN             = TwoPTwoCTraits::pWsN,
        pNsW             = TwoPTwoCTraits::pNsW,
    };


    typedef typename DomTraits::Scalar                Scalar;
    typedef typename DomTraits::CoordScalar           CoordScalar;
    typedef typename DomTraits::Grid                  Grid;
    typedef typename DomTraits::Element               Element;
    typedef typename DomTraits::ElementIterator       ElementIterator;
    typedef typename Element::EntityPointer           ElementPointer;
    typedef typename DomTraits::LocalPosition         LocalPosition;
    typedef typename DomTraits::GlobalPosition        GlobalPosition;
    typedef typename DomTraits::VertexIterator        VertexIterator;

    typedef typename BoxTraits::SolutionVector      SolutionVector;
    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;

    typedef typename TwoPTwoCTraits::PhasesVector        PhasesVector;
    typedef typename TwoPTwoCTraits::VariableVertexData  VariableVertexData;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    /*!
     * \brief Cached data for the each vertex of the element.
     */
    struct ElementData
    {
        VariableVertexData vertex[BoxTraits::ShapeFunctionSetContainer::maxsize];
    };

    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVertexData {
        int phaseState;
        int oldPhaseState;
    };

    /*!
     * \brief Function to update variable data of the vertices of the
     *        the current element (essentially secondary variables)
     */
    void updateVarVertexData_(VariableVertexData &vertDat,
                              const SolutionVector &vertSol,
                              int phaseState,
                              const Element &element,
                              int localIdx,
                              Problem &problem,
                              Scalar temperature) const
    {
        const GlobalPosition &global = element.geometry().corner(localIdx);
        const LocalPosition &local =
            DomTraits::referenceElement(element.type()).position(localIdx,
                                                                 dim);

        if (formulation == pWsN)
            {
                vertDat.pW = vertSol[pressureIdx];
                if (phaseState == bothPhases) vertDat.satN = vertSol[switchIdx];
                else if (phaseState == wPhaseOnly) vertDat.satN = 0.0;
                else if (phaseState == nPhaseOnly) vertDat.satN = 1.0;
                else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                vertDat.satW = 1.0 - vertDat.satN;
                vertDat.pC = problem.materialLaw().pC(vertDat.satW,
                                                      global,
                                                      element,
                                                      local);

                vertDat.pN = vertDat.pW + vertDat.pC;
            }
        else if (formulation == pNsW)
            {
                vertDat.pN = vertSol[pressureIdx];
                if (phaseState == bothPhases) vertDat.satW = vertSol[switchIdx];
                else if (phaseState == wPhaseOnly) vertDat.satW = 1.0;
                else if (phaseState == nPhaseOnly) vertDat.satW = 0.0;
                else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                vertDat.satN = 1.0 - vertDat.satW;
                vertDat.pC = problem.materialLaw().pC(vertDat.satW,
                                                      global,
                                                      element,
                                                      local);

                vertDat.pW = vertDat.pN - vertDat.pC;
            }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        // Solubilities of components in phases
        if (phaseState == bothPhases) {
            vertDat.massfrac[wComp][nPhase] = problem.multicomp().xWN(vertDat.pN, temperature);
            vertDat.massfrac[nComp][wPhase] = problem.multicomp().xAW(vertDat.pN, temperature);
        }
        else if (phaseState == wPhaseOnly) {
            vertDat.massfrac[wComp][nPhase] = 0.0;
            vertDat.massfrac[nComp][wPhase] = vertSol[switchIdx];
        }
        else if (phaseState == nPhaseOnly){
            vertDat.massfrac[wComp][nPhase] = vertSol[switchIdx];
            vertDat.massfrac[nComp][wPhase] = 0.0;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

        vertDat.massfrac[wComp][wPhase] = 1.0 - vertDat.massfrac[nComp][wPhase];
        vertDat.massfrac[nComp][nPhase] = 1.0 - vertDat.massfrac[wComp][nPhase];
        //                    vertDat.phaseState = phaseState;

        vertDat.density[wPhase] = problem.wettingPhase().density(temperature,
                                                                 vertDat.pW,
                                                                 vertDat.massfrac[nComp][wPhase]);
        vertDat.density[nPhase] = problem.nonwettingPhase().density(temperature,
                                                                    vertDat.pN,
                                                                    vertDat.massfrac[wComp][nPhase]);

        // Mobilities
        vertDat.mobility[wPhase] = problem.materialLaw().mobW(vertDat.satW,
                                                              global,
                                                              element,
                                                              local,
                                                              temperature,
                                                              vertDat.pW);
        vertDat.mobility[nPhase] = problem.materialLaw().mobN(vertDat.satN,
                                                              global,
                                                              element,
                                                              local,
                                                              temperature,
                                                              vertDat.pN);

        // diffusion coefficients
        vertDat.diffCoeff[wPhase] = problem.wettingPhase().diffCoeff(temperature, vertDat.pW);
        vertDat.diffCoeff[nPhase] = problem.nonwettingPhase().diffCoeff(temperature, vertDat.pN);
    }

public:
    TwoPTwoCBoxJacobianBase(ProblemT &problem)
        : ParentType(problem),
          staticVertexDat_(problem.numVertices())
    {
        switchFlag_ = false;
    };

    /*!
     * \brief Set the current grid element.
     */
    void setCurrentElement(const Element &element)
    {
        ParentType::setCurrentElement_(element);
    };

    /*!
     * \brief Set the parameters for the calls to the remaining
     *        members.
     */
    void setParams(const Element &element, LocalFunction &curSol, LocalFunction &prevSol)
    {
        setCurrentElement(element);

        // TODO: scheme which allows not to copy curSol and
        // prevSol all the time
        curSol_ = curSol;
        updateElementData_(curElemDat_, curSol_, false);
        curSolDeflected_ = false;

        prevSol_ = prevSol;
        updateElementData_(prevElemDat_, prevSol_, true);
    };

    /*!
     * \brief Vary a single component of a single vertex of the
     *        local solution for the current element.
     *
     * This method is a optimization, since if varying a single
     * component at a degree of freedom not the whole element cache
     * needs to be recalculated. (Updating the element cache is very
     * expensive since material laws need to be evaluated.)
     */
    void deflectCurSolution(int vert, int component, Scalar value)
    {
        // make sure that the original state can be restored
        if (!curSolDeflected_) {
            curSolDeflected_ = true;

            curSolOrigValue_ = curSol_[vert][component];
            curSolOrigVarData_ = curElemDat_.vertex[vert];
        }

        int globalIdx = ParentType::problem_.vertexIdx(ParentType::curElement_(),
                                                       vert);

        curSol_[vert][component] = value;
        asImp_()->updateVarVertexData_(curElemDat_.vertex[vert],
                                       curSol_[vert],
                                       staticVertexDat_[globalIdx].phaseState,
                                       this->curElement_(),
                                       vert,
                                       this->problem_,
                                       Implementation::temperature_(curSol_[vert]));
    }

    /*!
     * \brief Restore the local jacobian to the state before
     *        deflectCurSolution() was called.
     *
     * This only works if deflectSolution was only called with
     * (vertex, component) as arguments.
     */
    void restoreCurSolution(int vert, int component)
    {
        curSolDeflected_ = false;
        curSol_[vert][component] = curSolOrigValue_;
        curElemDat_.vertex[vert] = curSolOrigVarData_;
    };

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element in the 2P-2C
     *        model.
     *
     * This function should not include the source and sink terms.
     */
    void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
    {
        result = Scalar(0);

        // if flag usePrevSol is set, the solution from the previous time step is used,
        // otherwise the current solution is used. The secondary variables are used accordingly.
        // This computes the derivative of the storage term.
        const LocalFunction &sol   = usePrevSol ? this->prevSol_ : this->curSol_;
        const ElementData &elementCache = usePrevSol ? prevElemDat_  : curElemDat_;

        Scalar satN = elementCache.vertex[scvId].satN;
        Scalar satW = elementCache.vertex[scvId].satW;

        // assume porosity defined at vertices
        Scalar porosity =
            this->problem_.porosity(this->curElement_(), scvId);

        // storage of component water
        result[wComp] =
            porosity*(elementCache.vertex[scvId].density[wPhase]*
                      satW*
                      elementCache.vertex[scvId].massfrac[wComp][wPhase]
                      + elementCache.vertex[scvId].density[nPhase]*
                      satN*
                      elementCache.vertex[scvId].massfrac[wComp][nPhase]);

        // storage of component air
        result[nComp] =
            porosity*(elementCache.vertex[scvId].density[nPhase]*
                      satN*
                      elementCache.vertex[scvId].massfrac[nComp][nPhase]
                      + elementCache.vertex[scvId].density[wPhase]*
                      satW*
                      elementCache.vertex[scvId].massfrac[nComp][wPhase]);

        // storage of energy (if nonisothermal model is used)
        asImp_()->heatStorage(result, scvId, sol, elementCache);
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(SolutionVector &flux, int faceId) const
    {
        // set flux vector to zero
        int i = this->curElementGeom_.subContVolFace[faceId].i;
        int j = this->curElementGeom_.subContVolFace[faceId].j;

        // normal vector, value of the area of the scvf
        const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

        // get global coordinates of verts i,j
        const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
        const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

        // get local coordinates of verts i,j
        const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
        const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

        const ElementData &elemDat = this->curElemDat_;
        const VariableVertexData &vDat_i = elemDat.vertex[i];
        const VariableVertexData &vDat_j = elemDat.vertex[j];

        GlobalPosition pGrad[numPhases];
        GlobalPosition xGrad[numPhases];
        GlobalPosition tempGrad(0.0);
        for (int phase = 0; phase < numPhases; ++phase) {
            pGrad[phase] = Scalar(0);
            xGrad[phase] = Scalar(0);
        }

        GlobalPosition tmp(0.0);
        PhasesVector pressure(0.0), massfrac(0.0);
        PhasesVector densityIJ(0.);

        // calculate FE gradient (grad p for each phase)
        for (int idx = 0; idx < this->curElementGeom_.numVertices; idx++) // loop over adjacent vertices
            {
                // FEGradient at vertex idx
                const LocalPosition &feGrad = this->curElementGeom_.subContVolFace[faceId].grad[idx];

                pressure[wPhase] = elemDat.vertex[idx].pW;
                pressure[nPhase] = elemDat.vertex[idx].pN;

                // compute sum of pressure gradients for each phase
                for (int phase = 0; phase < numPhases; phase++)
                    {
                        // the pressure gradient
                        tmp = feGrad;
                        tmp *= pressure[phase];
                        pGrad[phase] += tmp;
                        densityIJ[phase] += elemDat.vertex[idx].density[phase] *
                            this->curElementGeom_.subContVolFace[faceId].shapeValue[idx];
                    }

                // the concentration gradient of the non-wetting
                // component in the wetting phase
                tmp = feGrad;
                tmp *= elemDat.vertex[idx].massfrac[nComp][wPhase];
                xGrad[wPhase] += tmp;

                // the concentration gradient of the wetting component
                // in the non-wetting phase
                tmp = feGrad;
                tmp *= elemDat.vertex[idx].massfrac[wComp][nPhase];
                xGrad[nPhase] += tmp;

                // temperature gradient
                asImp_()->updateTempGrad(tempGrad, feGrad, this->curSol_, idx);
            }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phase=0; phase < numPhases; phase++)
            {
                tmp = this->problem_.gravity();
                tmp *= densityIJ[phase];
                pGrad[phase] -= tmp;
            }

        // calculate the permeability tensor
        Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        harmonicMeanK_(K, Kj);

        // magnitute of darcy velocity of each phase projected
        // on the normal of the sub-control volume's face
        PhasesVector vDarcyOut;
        // temporary vector for the Darcy velocity
        GlobalPosition vDarcy;
        for (int phase=0; phase < numPhases; phase++)
            {
                K.mv(pGrad[phase], vDarcy);  // vDarcy = K * grad p
                vDarcyOut[phase] = vDarcy*normal;
            }

        // find upsteam and downstream verts
        const VariableVertexData *upW = &vDat_i;
        const VariableVertexData *dnW = &vDat_j;
        const VariableVertexData *upN = &vDat_i;
        const VariableVertexData *dnN = &vDat_j;

        if (vDarcyOut[wPhase] > 0) {
            std::swap(upW, dnW);
        };
        if (vDarcyOut[nPhase] > 0)  {
            std::swap(upN, dnN);
        };

        // Upwind parameter
        Scalar alpha = 1.0; // -> use only the upstream vertex

        ////////
        // advective flux of the wetting component
        ////////

        // flux in the wetting phase
        flux[wComp] =  vDarcyOut[wPhase] * (
                                            alpha* // upstream verts
                                            (  upW->density[wPhase] *
                                               upW->mobility[wPhase] *
                                               upW->massfrac[wComp][wPhase])
                                            +
                                            (1-alpha)* // downstream verts
                                            (  dnW->density[wPhase] *
                                               dnW->mobility[wPhase] *
                                               dnW->massfrac[wComp][wPhase]));

        // flux in the non-wetting phase
        flux[wComp] += vDarcyOut[nPhase] * (
                                            alpha* // upstream vert
                                            (  upN->density[nPhase] *
                                               upN->mobility[nPhase] *
                                               upN->massfrac[wComp][nPhase])
                                            +
                                            (1-alpha)* // downstream vert
                                            (  dnN->density[nPhase] *
                                               dnN->mobility[nPhase] *
                                               dnN->massfrac[wComp][nPhase]) );

        ////////
        // advective flux of the non-wetting component
        ////////

        // flux in the wetting phase
        flux[nComp]   = vDarcyOut[nPhase] * (
                                             alpha * // upstream verts
                                             (  upN->density[nPhase] *
                                                upN->mobility[nPhase] *
                                                upN->massfrac[nComp][nPhase])
                                             +
                                             (1-alpha) * // downstream vert
                                             (  dnN->density[nPhase] *
                                                dnN->mobility[nPhase] *
                                                dnN->massfrac[nComp][nPhase]) );

        // flux in the non-wetting phase
        flux[nComp]  += vDarcyOut[wPhase] * (
                                             alpha * // upstream vert
                                             (  upW->density[wPhase] *
                                                upW->mobility[wPhase] *
                                                upW->massfrac[nComp][wPhase])
                                             +
                                             (1-alpha) * // downstream vert
                                             (  dnW->density[wPhase] *
                                                dnW->mobility[wPhase] *
                                                dnW->massfrac[nComp][wPhase]) );

        ////////
        // advective flux of energy
        ////////
        asImp_()->advectiveHeatFlux(flux, vDarcyOut, alpha, upW, dnW, upN, dnN);

        /////////////////////////////
        // DIFFUSION
        /////////////////////////////

        SolutionVector normDiffGrad;

        Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of liquid
        Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of gas
        SolutionVector avgDensity;

        // Diffusion coefficient

        // calculate tortuosity at the nodes i and j needed
        // for porous media diffusion coefficient
        Scalar tauW_i, tauW_j;
        Scalar tauN_i, tauN_j;

        Scalar porosity_i = this->problem_.porosity(this->curElement_(), i);
        Scalar porosity_j = this->problem_.porosity(this->curElement_(), j);

        tauW_i = pow(porosity_i * vDat_i.satW,
                     7.0/3) / porosity_i;
        tauW_j = pow(porosity_j * vDat_j.satW,
                     7.0/3) / porosity_j;
        tauN_i = pow(porosity_i * vDat_i.satN,
                     7.0/3) / porosity_i;
        tauN_j = pow(porosity_j * vDat_j.satN,
                     7.0/3) / porosity_j;

        // arithmetic mean of porous media diffusion coefficient
        Scalar Dwn, Daw;

        // approximate the effective cross sections for
        // diffusion by the harmonic mean of the volume
        // occupied by the phases
        /*
          Dwn = harmonicMean_(vDat_i.satN * tauN_i * vDat_i.diffCoeff[nPhase],
          vDat_j.satN * tauN_j * vDat_j.diffCoeff[nPhase]);
          Daw = harmonicMean_(vDat_i.satW * tauW_i * vDat_i.diffCoeff[wPhase],
          vDat_j.satW * tauW_j * vDat_j.diffCoeff[wPhase]);

          //                std::cout << "Daw: " << Daw << ", Dwn: " << Dwn << "\n";
          */

        // arithmetic mean (2e-7)
        Dwn = 1./2*(porosity_i * vDat_i.satN * tauN_i * vDat_i.diffCoeff[nPhase] +
                    porosity_j * vDat_j.satN * tauN_j * vDat_j.diffCoeff[nPhase]);
        Daw = 1./2*(porosity_i * vDat_i.satW * tauW_i * vDat_i.diffCoeff[wPhase] +
                    porosity_j * vDat_j.satW * tauW_j * vDat_j.diffCoeff[wPhase]);
        if (vDat_i.satN == 0 || vDat_j.satN == 0)
            Dwn = 0;
        if (vDat_i.satW == 0 || vDat_j.satW == 0)
            Daw = 0;

        // calculate the arithmetic mean of densities
        avgDensity[wPhase] = 0.5*(vDat_i.density[wPhase] + vDat_j.density[wPhase]);
        avgDensity[nPhase] = 0.5*(vDat_i.density[nPhase] + vDat_j.density[nPhase]);

        diffusionAW = Daw * avgDensity[wPhase] * (xGrad[wPhase] * normal);
        diffusionWW = - diffusionAW;
        diffusionWN = Dwn * avgDensity[nPhase] * (xGrad[nPhase] * normal);
        diffusionAN = - diffusionWN;

        // add diffusion of water to water flux
        flux[wComp] += diffusionWW + diffusionWN;

        // add diffusion of air to air flux
        flux[nComp] += diffusionAN + diffusionAW;

        ////////
        // diffusive flux of energy (only for non-isothermal
        // models)
        ////////
        asImp_()->diffusiveHeatFlux(flux, faceId, tempGrad);
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(SolutionVector &q, int localVertexIdx)
    {
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }


    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * Called by TwoPTwoCBoxModel::initial()
     */
    void initStaticData()
    {
        setSwitched(false);

        VertexIterator it = this->problem_.vertexBegin();
        VertexIterator endit = this->problem_.vertexEnd();
        for (; it != endit; ++it)
            {
                int globalIdx = this->problem_.vertexIdx(*it);
                const GlobalPosition &globalPos = it->geometry().corner(0);

                // initialize phase state
                staticVertexDat_[globalIdx].phaseState =
                    this->problem_.initialPhaseState(*it, globalIdx, globalPos);
                staticVertexDat_[globalIdx].oldPhaseState =
                    staticVertexDat_[globalIdx].phaseState;
            }
    }

    /*!
     * \brief Update the static data of a single vert and do a
     *        variable switch if necessary.
     */
    void updateStaticData(SpatialFunction &curGlobalSol, SpatialFunction &oldGlobalSol)
    {
        bool wasSwitched = false;

        VertexIterator it = this->problem_.vertexBegin();
        for (; it != this->problem_.vertexEnd(); ++it)
            {
                int globalIdx = this->problem_.vertexIdx(*it);
                const GlobalPosition &global = it->geometry().corner(0);

                wasSwitched = primaryVarSwitch_(curGlobalSol,
                                                globalIdx,
                                                global)
                    || wasSwitched;
            }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        wasSwitched = this->problem_.grid().comm().max(wasSwitched);

        setSwitched(wasSwitched);
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhaseState()
    {
        int numVertices = this->problem_.numVertices();
        for (int i = 0; i < numVertices; ++i)
            staticVertexDat_[i].oldPhaseState = staticVertexDat_[i].phaseState;
    }

    /*!
     * \brief Reset the current phase state of all verts to the old one after an update failed
     */
    void resetPhaseState()
    {
        int numVertices = this->problem_.numVertices();
        for (int i = 0; i < numVertices; ++i)
            staticVertexDat_[i].phaseState = staticVertexDat_[i].oldPhaseState;
    }

    /*!
     * \brief Return true if the primary variables were switched
     *        after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Set whether there was a primary variable switch after in the last
     *        timestep.
     */
    void setSwitched(bool yesno)
    {
        switchFlag_ = yesno;
    }

    /*!
     * \brief Compute Darcy velocity for output
     *
     */
    //    template<class PressureFunction, class Vector, class Grid, class Problem>
    void calculateDarcyVelocity(GlobalPosition &velocity,
                                const int faceId)
    {
        velocity = 0.0;
        return ; // TODO;

        int i = this->curElementGeom_.subContVolFace[faceId].i;
        int j = this->curElementGeom_.subContVolFace[faceId].j;

        // normal vector of value of the area of the scvf
        const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

        // get global coordinates of verts i,j
        const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
        const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

        // get local coordinates of verts i,j
        const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
        const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

        const ElementData &elemDat = this->curElemDat_;
        const VariableVertexData &vDat_i = elemDat.vertex[i];
        const VariableVertexData &vDat_j = elemDat.vertex[j];

        GlobalPosition pGrad[numPhases];
        GlobalPosition tempGrad(0.0);
        for (int phase = 0; phase < numPhases; ++phase) {
            pGrad[phase] = Scalar(0);
        }

        GlobalPosition tmp(0.0);
        PhasesVector pressure(0.0);
        PhasesVector densityIJ(0.);

        // calculate FE gradient (grad p for each phase)
        for (int idx = 0; idx < this->curElementGeom_.numVertices; idx++) // loop over vertices
            {
                // FEGradient at vertex idx
                const LocalPosition &feGrad = this->curElementGeom_.subContVolFace[faceId].grad[idx];
                pressure[wPhase] = elemDat.vertex[idx].pW;
                pressure[nPhase] = elemDat.vertex[idx].pN;

                // compute sum of pressure gradients for each phase
                for (int phase = 0; phase < numPhases; phase++)
                    {
                        // the pressure gradient
                        tmp = feGrad;
                        tmp *= pressure[phase];

                        pGrad[phase] += tmp;
                        densityIJ[phase] += elemDat.vertex[idx].density[phase] *
                            this->curElementGeom_.subContVolFace[faceId].shapeValue[idx];
                    }
            }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phase=0; phase < numPhases; phase++)
            {
                tmp = this->problem_.gravity();
                tmp *= densityIJ[phase];
                pGrad[phase] -= tmp;
            }

        // calculate the permeability tensor
        Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        harmonicMeanK_(K, Kj);

        // magnitude of Darcy velocity of each phase projected
        // on the normal of the sub-control volume's face
        PhasesVector vDarcyOut;
        // temporary vector for the Darcy velocity
        GlobalPosition vDarcy;
        for (int phase=0; phase < numPhases; phase++)
            {
                K.mv(pGrad[phase], vDarcy);  // vDarcy = K * grad p
                vDarcyOut[phase] = vDarcy*normal;
            }

        // find upsteam and downstream verts
        const VariableVertexData *upW = &vDat_i;
        const VariableVertexData *dnW = &vDat_j;
        const VariableVertexData *upN = &vDat_i;
        const VariableVertexData *dnN = &vDat_j;

        if (vDarcyOut[wPhase] > 0) {
            std::swap(upW, dnW);
        };
        if (vDarcyOut[nPhase] > 0)  {
            std::swap(upN, dnN);
        };

        // Upwind parameter
        Scalar alpha = 1.0; // -> use only the upstream vertex

        ////////
        // advective flux of the wetting component
        ////////

        // flux in the wetting phase
        velocity[wPhase] =  vDarcyOut[wPhase] * (
                                                 alpha* // upstream verts
                                                 (  upW->density[wPhase] *
                                                    upW->mobility[wPhase])
                                                 +
                                                 (1-alpha)* // downstream vert
                                                 (  dnW->density[wPhase] *
                                                    dnW->mobility[wPhase]));

        // flux in the non-wetting phase
        velocity[nPhase] =  vDarcyOut[nPhase] * (
                                                 alpha* // upstream verts
                                                 (  upN->density[nPhase] *
                                                    upN->mobility[nPhase])
                                                 +
                                                 (1-alpha)* // downstream vert
                                                 (  dnN->density[nPhase] *
                                                    dnN->mobility[nPhase]));
    }

    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SpatialFunction &globalSol, Dune::FieldVector<Scalar, 4> &mass)
    {
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        unsigned numVertices = this->problem_.numVertices();
        LocalFunction curSol(numVertices);
        ElementData   elemDat;
        VariableVertexData tmp;
        int state;
        Scalar vol, poro, rhoN, rhoW, satN, satW, xAW, xWW, xWN, xAN, pW, Te;
        Scalar massNComp(0.), massNCompNPhase(0.), massWComp(0.), massWCompWPhase(0.);

        enum
        {   gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase state

        mass = 0;
        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
        Scalar minTe = 1e100;
        Scalar maxTe = -1e100;
        Scalar minX = 1e100;
        Scalar maxX = -1e100;

        // Loop over elements
        for (; elementIt != endit; ++elementIt)
            {

                setCurrentElement(*elementIt);
                this->restrictToElement(curSol, globalSol);
                updateElementData_(elemDat, curSol, false);
                // get geometry type

                int numLocalVerts = elementIt->template count<dim>();

                // Loop over element vertices
                for (int i = 0; i < numLocalVerts; ++i)
                    {
                        int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                        vol = this->curElementGeom_.subContVol[i].volume;

                        state =  staticVertexDat_[globalIdx].phaseState;
                        poro = this->problem_.porosity(this->curElement_(), i);
                        rhoN = elemDat.vertex[i].density[nPhase];
                        rhoW = elemDat.vertex[i].density[wPhase];
                        satN = elemDat.vertex[i].satN;
                        satW = elemDat.vertex[i].satW;
                        xAW = elemDat.vertex[i].massfrac[nComp][wPhase];
                        xWW = elemDat.vertex[i].massfrac[wComp][wPhase];
                        xWN = elemDat.vertex[i].massfrac[wComp][nPhase];
                        xAN = elemDat.vertex[i].massfrac[nComp][nPhase];
                        pW = elemDat.vertex[i].pW;
                        Te = Implementation::temperature_((*globalSol)[globalIdx]);
                        massNComp = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);
                        massNCompNPhase = vol * poro * satN * rhoN * xAN;
                        massWComp = vol * poro * (satW * rhoW * xWW + satN * rhoN * xWN);
                        massWCompWPhase = vol * poro * satW * rhoW * xWW;

                        // get minimum and maximum values of primary variables
                        minSat = std::min(minSat, satN);
                        maxSat = std::max(maxSat, satN);
                        minP = std::min(minP, pW);
                        maxP = std::max(maxP, pW);
                        minX = std::min(minX, xAW);
                        maxX = std::max(maxX, xAW);
                        minTe = std::min(minTe, Te);
                        maxTe = std::max(maxTe, Te);

                        // IF PARALLEL: check all processors to get minimum and maximum
                        //values of primary variables
                        // also works for sequential calculation
                        minSat = collectiveCom_.min(minSat);
                        maxSat = collectiveCom_.max(maxSat);
                        minP = collectiveCom_.min(minP);
                        maxP = collectiveCom_.max(maxP);
                        minX = collectiveCom_.min(minX);
                        maxX = collectiveCom_.max(maxX);
                        minTe = collectiveCom_.min(minTe);
                        maxTe = collectiveCom_.max(maxTe);

                        // calculate total mass
                        mass[0] += massNComp;       // total mass of nonwetting component
                        mass[1] += massNCompNPhase; // mass of nonwetting component in nonwetting phase
                        mass[2] += massWComp;       // total mass of wetting component
                        mass[3] += massWCompWPhase; // mass of wetting component in wetting phase

                        // IF PARALLEL: calculate total mass including all processors
                        // also works for sequential calculation
                        mass = collectiveCom_.sum(mass);
                    }

            }
        if(collectiveCom_.rank() == 0) // IF PARALLEL: only print by processor with rank() == 0
            {
                // print minimum and maximum values
                std::cout << "nonwetting phase saturation: min = "<< minSat
                          << ", max = "<< maxSat << std::endl;
                std::cout << "wetting phase pressure: min = "<< minP
                          << ", max = "<< maxP << std::endl;
                std::cout << "mass fraction nComp: min = "<< minX
                          << ", max = "<< maxX << std::endl;
                std::cout << "temperature: min = "<< minTe
                          << ", max = "<< maxTe << std::endl;
            }
    }
    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        unsigned numElements = this->problem_.numElements();
        ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pN =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoN =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobN =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracAinW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracAinN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracWinW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracWinN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *temperature  = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *phaseState   = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityX    = writer.template createField<Scalar, 1>(numElements);
        ScalarField *velocityY    = writer.template createField<Scalar, 1>(numElements);
        ScalarField *velocityZ    = writer.template createField<Scalar, 1>(numElements);

        LocalFunction curSol(numVertices);
        ElementData   elemDat;
        VariableVertexData tmp;

        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();

        for (; elementIt != endit; ++elementIt)
            {
                int numLocalVerts = elementIt->template count<dim>();

                setCurrentElement(*elementIt);
                this->restrictToElement(curSol, globalSol);
                updateElementData_(elemDat, curSol, false);

                for (int i = 0; i < numLocalVerts; ++i)
                    {
                        int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                        (*pW)[globalIdx] = elemDat.vertex[i].pW;
                        (*pN)[globalIdx] = elemDat.vertex[i].pN;
                        (*pC)[globalIdx] = elemDat.vertex[i].pC;
                        (*Sw)[globalIdx] = elemDat.vertex[i].satW;
                        (*Sn)[globalIdx] = elemDat.vertex[i].satN;
                        (*rhoW)[globalIdx] = elemDat.vertex[i].density[wPhase];
                        (*rhoN)[globalIdx] = elemDat.vertex[i].density[nPhase];
                        (*mobW)[globalIdx] = elemDat.vertex[i].mobility[wPhase];
                        (*mobN)[globalIdx] = elemDat.vertex[i].mobility[nPhase];
                        (*massfracAinW)[globalIdx] = elemDat.vertex[i].massfrac[nComp][wPhase];
                        (*massfracAinN)[globalIdx] = elemDat.vertex[i].massfrac[nComp][nPhase];
                        (*massfracWinW)[globalIdx] = elemDat.vertex[i].massfrac[wComp][wPhase];
                        (*massfracWinN)[globalIdx] = elemDat.vertex[i].massfrac[wComp][nPhase];
                        (*temperature)[globalIdx] = Implementation::temperature_((*globalSol)[globalIdx]);
                        (*phaseState)[globalIdx] = staticVertexDat_[globalIdx].phaseState;
                    };

                // Vector containing the velocity at the element
                GlobalPosition velocity[numPhases];
                GlobalPosition elementVelocity[numPhases];

                // loop over the phases
                for (int phase=0; phase < numPhases; phase++)
                    {
                        elementVelocity[phase] = 0;

                        int elementIdx = this->problem_.elementIdx(*elementIt);
                        for (int faceId = 0; faceId< this->curElementGeom_.numEdges; faceId++)
                            {
                                velocity[phase] = 0;
                                asImp_()->calculateDarcyVelocity(velocity[phase],
                                                                 faceId);
                                elementVelocity[phase] += velocity[phase];
                            }
                        elementVelocity[phase] *= 1.0/this->curElementGeom_.numEdges;
                        (*velocityX)[elementIdx] = elementVelocity[0][phase];
                        if (dim >= 2)
                            (*velocityY)[elementIdx] = elementVelocity[1][phase];
                        if (dim == 3)
                            (*velocityZ)[elementIdx] = elementVelocity[2][phase];
                    }
            }

        writer.addVertexData(pW, "pW");
        writer.addVertexData(pN, "pN");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(Sw, "SW");
        writer.addVertexData(Sn, "SN");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(rhoN, "rhoN");
        writer.addVertexData(mobW, "mobW");
        writer.addVertexData(mobN, "mobN");
        writer.addVertexData(massfracAinW, "XaW");
        writer.addVertexData(massfracAinN, "XaN");
        writer.addVertexData(massfracWinW, "XwW");
        writer.addVertexData(massfracWinN, "XwN");
        writer.addVertexData(temperature, "T");
        writer.addVertexData(phaseState, "phase state");
        writer.addCellData(velocityX, "Vx");
        if (dim >= 2)
            writer.addCellData(velocityY, "Vy");
        if (dim == 3)
            writer.addCellData(velocityZ, "Vz");
    }


protected:
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }


    void updateElementData_(ElementData &dest, const LocalFunction &sol, bool isOldSol)
    {
        int phaseState;
        int numVertices = this->curElement_().template count<dim>();
        for (int i = 0; i < numVertices; i++) {
            int iGlobal = ParentType::problem_.vertexIdx(ParentType::curElement_(), i);
            if (isOldSol)
                phaseState = staticVertexDat_[iGlobal].oldPhaseState;
            else
                phaseState = staticVertexDat_[iGlobal].phaseState;
            asImp_()->updateVarVertexData_(dest.vertex[i],
                                           sol[i],
                                           phaseState,
                                           this->curElement_(),
                                           i,
                                           this->problem_,
                                           Implementation::temperature_(sol[i]));
        }
    }


    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SpatialFunction &globalSol,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        int phaseState    = staticVertexDat_[globalIdx].phaseState;
        int newPhaseState = phaseState;
        Scalar pW = 0;
        Scalar pN = 0;
        Scalar satW = 0;
        Scalar satN = 0;

        if (formulation == pWsN)
            {
                pW = (*globalSol)[globalIdx][pressureIdx];
                if      (phaseState == bothPhases) satN = (*globalSol)[globalIdx][switchIdx];
                else if (phaseState == wPhaseOnly) satN = 0.0;
                else if (phaseState == nPhaseOnly) satN = 1.0;

                Scalar satW = 1 - satN;

                Scalar pC = this->problem_.pC(satW, globalIdx, globalPos);
                pN = pW + pC;
            }

        // Evaluate saturation and pressures
        else if (formulation == pNsW)
            {
                pN = (*globalSol)[globalIdx][pressureIdx];
                satW = 0.0;
                if      (phaseState == bothPhases) satW = (*globalSol)[globalIdx][switchIdx];
                else if (phaseState == wPhaseOnly) satW = 1.0;
                else if (phaseState == nPhaseOnly) satW = 0.0;

                satN = 1 - satW;

                Scalar pC = this->problem_.pC(satW, globalIdx, globalPos);
                pW = pN - pC;
            }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        Scalar temperature = Implementation::temperature_((*globalSol)[globalIdx]);

        // phase state is checked and a switch is performed
        if (phaseState == nPhaseOnly)
            {
                Scalar xWN = (*globalSol)[globalIdx][switchIdx];
                Scalar xWNmax = this->problem_.multicomp().xWN(pN, temperature);

                if (xWN > xWNmax)
                    {
                        // wetting phase appears
                        std::cout << "wetting phase appears at vertex " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        newPhaseState = bothPhases;
                        if (formulation == pNsW) (*globalSol)[globalIdx][switchIdx] = 1e-3;
                        else if (formulation == pWsN) (*globalSol)[globalIdx][switchIdx] = 1 - 1e-3;
                    };
            }
        else if (phaseState == wPhaseOnly)
            {
                Scalar xAW = (*globalSol)[globalIdx][switchIdx];
                Scalar xAWmax = this->problem_.multicomp().xAW(pN, temperature);

                if (xAW > xAWmax)
                    {
                        // non-wetting phase appears
                        std::cout << "Non-wetting phase appears at vertex " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        if (formulation == pNsW)
                            (*globalSol)[globalIdx][switchIdx] = 1 - 1e-3;
                        else if (formulation == pWsN)
                            (*globalSol)[globalIdx][switchIdx] = 1e-3;
                        newPhaseState = bothPhases;
                    }
            }
        else if (phaseState == bothPhases) {
            if (satN < 0) {
                // non-wetting phase disappears
                std::cout << "Non-wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                (*globalSol)[globalIdx][switchIdx]
                    = this->problem_.multicomp().xAW(pN, temperature)*(1 - 1e-2);
                newPhaseState = wPhaseOnly;
            }
            else if (satW < 0) {
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                (*globalSol)[globalIdx][switchIdx]
                    = this->problem_.multicomp().xWN(pN, temperature)*(1 - 1e-2);
                newPhaseState = nPhaseOnly;
            }
        }

        staticVertexDat_[globalIdx].phaseState = newPhaseState;

        return phaseState != newPhaseState;
    }

    // harmonic mean of the permeability computed directly.  the
    // first parameter is used to store the result.
    static void harmonicMeanK_(Tensor &Ki, const Tensor &Kj)
    {
        for (int kx=0; kx < Tensor::rows; kx++){
            for (int ky=0; ky< Tensor::cols; ky++){
                if (Ki[kx][ky] != Kj[kx][ky]) {
                    Ki[kx][ky] = harmonicMean_(Ki[kx][ky], Kj[kx][ky]);
                }
            }
        }
    }

    // returns the harmonic mean of two scalars
    static Scalar harmonicMean_(Scalar x, Scalar y)
    {
        if (x == 0 || y == 0)
            return 0;
        return (2*x*y)/(x + y);
    };

    // parameters given in constructor
    std::vector<StaticVertexData> staticVertexDat_;
    bool                          switchFlag_;
    int                           formulation_;

    // current solution
    LocalFunction      curSol_;
    ElementData        curElemDat_;

    // needed for restoreCurSolution()
    bool               curSolDeflected_;
    Scalar             curSolOrigValue_;
    VariableVertexData curSolOrigVarData_;

    // previous solution
    LocalFunction      prevSol_;
    ElementData        prevElemDat_;
    CollectiveCommunication collectiveCom_;
};


} // end namepace

#endif
