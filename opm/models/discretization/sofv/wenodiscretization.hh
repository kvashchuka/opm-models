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
 * \copydoc Opm::WenoDiscretization
 */
#ifndef EWOMS_WENO_DISCRETIZATION_HH
#define EWOMS_WENO_DISCRETIZATION_HH

#include <opm/material/densead/Math.hpp>

#include "sofvproperties.hh"
#include "sofvstencil.hh"
#include "sofvgridcommhandlefactory.hh"
#include "sofvbaseoutputmodule.hh"

#include <opm/simulators/linalg/elementborderlistfromgrid.hh>
#include <opm/models/discretization/common/fvbasediscretization.hh>

#include <unordered_set>

namespace Opm {
template <class TypeTag>
class WenoDiscretization;
}

BEGIN_PROPERTIES

//! Set the stencil
SET_PROP(WenoDiscretization, Stencil)
{
public:
    typedef Opm::SofvStencil<TypeTag> type;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(WenoDiscretization, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(WenoDiscretization, Discretization, Opm::WenoDiscretization<TypeTag>);

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
SET_TYPE_PROP(WenoDiscretization, DiscBaseOutputModule,
              Opm::SofvBaseOutputModule<TypeTag>);

//! Set the DiscreteFunctionSpace
/*
SET_PROP(WenoDiscretization, DiscreteFunctionSpace)
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

//! The class to create grid communication handles
SET_TYPE_PROP(WenoDiscretization, GridCommHandleFactory,
              Opm::SofvGridCommHandleFactory<TypeTag>);

//! Set the border list creator for to the one of an element based
//! method
SET_PROP(WenoDiscretization, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Opm::Linear::ElementBorderListFromGrid<GridView, ElementMapper> type;
};

//! For the element centered finite volume method, ghost and overlap elements must be
//! assembled to calculate the fluxes over the process boundary faces of the local
//! process' grid partition
SET_BOOL_PROP(WenoDiscretization, LinearizeNonLocalElements, true);

//! locking is not required for the element centered finite volume method because race
//! conditions cannot occur since each matrix/vector entry is written exactly once
SET_BOOL_PROP(WenoDiscretization, UseLinearizationLock, false);

SET_BOOL_PROP(WenoDiscretization, EnableHigherOrder, true);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup WenoDiscretization
 *
 * \brief The base class for the element-centered finite-volume discretization scheme.
 */
template<class TypeTag>
class WenoDiscretization : public FvBaseDiscretization<TypeTag>
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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename Opm::MathToolbox<Evaluation> Toolbox;
    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange  = PrimaryVariables::dimension };
    typedef Dune::FieldVector< Scalar, dimDomain > RangeType;
    typedef Dune::FieldVector<Evaluation, dimDomain> EvalDimVector;
    typedef typename GridView::ctype CoordScalar;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename GridView::template Codim<0>::Entity EntityType;
    // geometry type
    typedef typename EntityType::Geometry GeometryType;
    // global coordinates
    typedef typename GeometryType::GlobalCoordinate GlobalCoordinateType;
    // local coordinates
    typedef typename GeometryType::LocalCoordinate LocalCoordinateType;

    // Vector and Matrix used to weno reconstruction
    typedef Dune::FieldMatrix<CoordScalar, 3, 3> Matrix;
    typedef Dune::FieldVector<CoordScalar, 3> Vector;


public:
    using ParentType::simulator_;
    using ParentType::enableHigherOrder;

    WenoDiscretization(Simulator& simulator)
        : ParentType(simulator)
    {
        //std::vector<RangeType> totalMobility(numDof);
        // iterate through the grid and evaluate the initial condition
        const auto & gridView = simulator_.gridManager().gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        unsigned numElements = gridView.size(/*codim=*/0);
        triplets_.resize(numElements);
        invMatrix_.resize(numElements);
        //TODO get this from element!!
        unsigned globalIdx = 0;
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            updateVertexNeighbourSets_(elem, globalIdx);
            globalIdx++;
        }

    }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return "weno"; }

    /*!
 * \brief Called by the update() method before it tries to
 *        apply the newton method. This is primary a hook
 *        which the actual model can overload.
 */
    void updateBegin()
    {
        updateReconstruction_();
    }

    // The reconstruction is pre-computed
    void updateReconstruction_ ( )
    {
    }

    void evaluateReconstruction(const ElementContext& elemCtx,
                                const int exteriorDofIdx,
                                const GlobalCoordinateType& point,
                                EvalDimVector& interiorResult,
                                EvalDimVector& exteriorResult ) const
    {

        unsigned focusDofIdx = elemCtx.focusDofIndex();
        // Add the first order part
        const auto& interior = elemCtx.intensiveQuantities(0, /*timeIdx=*/0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (0 != static_cast<int>(focusDofIdx))
                interiorResult[phaseIdx] = Toolbox::value(interior.mobility(phaseIdx));
            else
                interiorResult[phaseIdx] = interior.mobility(phaseIdx);
        }
        const auto& exterior = elemCtx.intensiveQuantities(exteriorDofIdx, /*timeIdx=*/0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (exteriorDofIdx != static_cast<int>(focusDofIdx))
                exteriorResult[phaseIdx] = Toolbox::value(exterior.mobility(phaseIdx));
            else
                exteriorResult[phaseIdx] = exterior.mobility(phaseIdx);
        }

        if( enableHigherOrder() )
        {
            unsigned globalIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                typedef Dune::FieldMatrix<Scalar, 3, 3> DimMatrix3;
                typedef Dune::FieldVector<Evaluation, 3> EvalDimVector3;
                int numPoly = triplets_[globalIdx].size();
                std::vector<EvalDimVector> sigmas;
                // compute the sigmas (gradient)
                for (int polyIdx = 0; polyIdx < numPoly; ++polyIdx) {
                    const std::tuple<int,int, int>& elementIntPairs = triplets_[globalIdx][polyIdx];

                    if (std::get<0>(elementIntPairs) != exteriorDofIdx && std::get<1>(elementIntPairs) != exteriorDofIdx && std::get<2>(elementIntPairs) != exteriorDofIdx ) {
                        // only add contribution from triplets including the exterior cell
                        continue;
                    }
                    const DimMatrix3& C = invMatrix_[globalIdx][polyIdx];

                    // Use initial values? (x, 1)
                    if (dimDomain == 2) {
                    const auto& int0 = elemCtx.intensiveQuantities(0, 0);
                    const auto& int1 = elemCtx.intensiveQuantities(std::get<0>(elementIntPairs), 0);
                    const auto& int2 = elemCtx.intensiveQuantities(std::get<1>(elementIntPairs), 0);
                    Evaluation mob0 = int0.mobility(phaseIdx);
                    Evaluation mob1 = int1.mobility(phaseIdx);
                    Evaluation mob2 = int2.mobility(phaseIdx);
                    if (0 != static_cast<int>(focusDofIdx))
                        mob0 = Toolbox::value(mob0);
                    if (std::get<0>(elementIntPairs) != static_cast<int>(focusDofIdx))
                        mob1 = Toolbox::value(mob1);
                    if (std::get<1>(elementIntPairs) != static_cast<int>(focusDofIdx))
                        mob2 = Toolbox::value(mob2);

                    EvalDimVector3 p = {mob0, mob1, mob2};
                    EvalDimVector sigma(0.0);

                    C.mv(p, sigma);

                    sigmas.emplace_back (sigma);
                    } else if (dimDomain ==3) {
                    const auto& int0 = elemCtx.intensiveQuantities(0, 0);
                    const auto& int1 = elemCtx.intensiveQuantities(std::get<0>(elementIntPairs), 0);
                    const auto& int2 = elemCtx.intensiveQuantities(std::get<1>(elementIntPairs), 0);
                    const auto& int3 = elemCtx.intensiveQuantities(std::get<2>(elementIntPairs), 0);
                    Evaluation mob0 = int0.mobility(phaseIdx);
                    Evaluation mob1 = int1.mobility(phaseIdx);
                    Evaluation mob2 = int2.mobility(phaseIdx);
                    Evaluation mob3 = int3.mobility(phaseIdx);
                    if (0 != static_cast<int>(focusDofIdx))
                        mob0 = Toolbox::value(mob0);
                    if (std::get<0>(elementIntPairs) != static_cast<int>(focusDofIdx))
                        mob1 = Toolbox::value(mob1);
                    if (std::get<1>(elementIntPairs) != static_cast<int>(focusDofIdx))
                        mob2 = Toolbox::value(mob2);
                    if (std::get<2>(elementIntPairs) != static_cast<int>(focusDofIdx))
                        mob3 = Toolbox::value(mob3);

                    EvalDimVector3 p = {mob1-mob0, mob2-mob0, mob3-mob0};
                    EvalDimVector sigma(0.0);

                    C.mv(p, sigma);

                    sigmas.emplace_back (sigma);

                    }
                }

                // use unit volume for now.
                // we may want to weight contribution according to cell volume
                Scalar vol = 1; //stencil.subControlVolume(upstreamIdx).volume();

                // Compute the betas used to compute the weights
                std::vector<Evaluation> betas;

                const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);

                //std::cout << sigmas.size() <<  " " << stencil.numDof() <<std::endl;
                if (int(stencil.numDof()) < (dimDomain*2 + 0) )
                    continue;

                //if (sigmas.size() < (dimDomain-1)*2 + 1)
                //    continue;

                Evaluation sumbeta(0.0);
                //std::vector<Scalar> gammas = {0.8, 0.1,0.1}
                for (size_t polyIdx = 0; polyIdx < sigmas.size(); ++polyIdx) {
                    Scalar gamma = 1.0 / sigmas.size();
                    //Scalar gamma = gammas[polyIdx];
                    Evaluation beta = gamma / pow(1e-6 + sigmas[polyIdx].two_norm2()*vol,2);
                    sumbeta +=beta;
                    betas.emplace_back (beta);
                }

                EvalDimVector weno(0.0);
                for (size_t polyIdx = 0; polyIdx < betas.size(); ++polyIdx) {
                    const Evaluation weight = betas[polyIdx] / sumbeta;
                    for (int i = 0; i < dimDomain; ++i) {
                        weno[i] += sigmas[polyIdx][i] * weight;
                    }
                }

                // Add the second order term
                const auto& posIn = elemCtx.pos(0, /*timeIdx=*/0);
                RangeType distVecIn(point);
                distVecIn -= posIn;
                for (int i = 0; i < dimDomain; ++i) {
                    interiorResult[phaseIdx] += weno[i] * distVecIn[i];
                }
                const auto& posEx = elemCtx.pos(exteriorDofIdx, /*timeIdx=*/0);
                RangeType distVecEx(point);
                distVecEx -= posEx;
                for (int i = 0; i < dimDomain; ++i) {
                    exteriorResult[phaseIdx] += weno[i] * distVecEx[i];
                }
            }
        }
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
        double errorL2 = 0.0;

        if (std::abs(time - simulator_.endTime()) < 1.0) {
            // iterate through the grid and evaluate error
            const auto & gridView = simulator_.gridManager().gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;

                // deal with the current element
                elemCtx.updatePrimaryStencil(elem);

                Scalar exactWetSaturation = 0.0;

                for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(timeIdx); dofIdx++) {
                    unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
                    exactWetSaturation += simulator_.problem().exactWetSat(elemCtx, dofIdx, time);
                    Scalar volume = elemCtx.stencil(timeIdx).subControlVolume(dofIdx).volume(); //elemCtx.dofVolume(dofIdx, timeIdx);
                    //take only saturation, ignore pressure
                   // if (std::abs(solution[globalIdx][1]) < 1.0 + eps_){
                        errorL1 += std::abs(exactWetSaturation - solution[globalIdx][1]) * volume;
                        errorL2 += std::abs(exactWetSaturation - solution[globalIdx][1]) * std::abs(exactWetSaturation - solution[globalIdx][1]) * volume;
                    //}
                    //std::cout << " exact solution on each element  " << exactWetSaturation  << std::endl;
                }
            }

            std::cout << " Error-L1 = " << errorL1 << std::endl;
            std::cout << " Error-L2 = " << std::sqrt(errorL2) << std::endl;
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

    std::vector<std::vector<std::tuple<int, int, int> > > triplets_;
    std::vector<std::vector<Matrix>> invMatrix_;

    void updateVertexNeighbourSets_(const Element& element, unsigned globalIdx) {

        ElementContext elemCtx(simulator_);
        // deal with the current element
        elemCtx.updateStencil(element);
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);

        // Create the set of verticies accociated to this element.
        std::set<int> vertexIndices;
        const auto & gridView = simulator_.gridManager().gridView();
        const auto& endIsIt = gridView.iend(element);

        // add the corners of the element
        const int numVerticesInElement = element.subEntities(dimDomain);
        for (int i = 0; i < numVerticesInElement; ++i) {
            const int vertexIdx = gridView.indexSet().subIndex ( element, i, dimDomain);
            vertexIndices.insert(vertexIdx);
        }

        // add the hanging nodes
        auto isIt = gridView.ibegin(element);
        for (; isIt != endIsIt; ++isIt) {
            const auto& intersection = *isIt;
            if (!intersection.neighbor())
                continue;

            const Element& outside = intersection.outside();
             if (!intersection.conforming ()) {
                const auto& intersectionGeometry = intersection.geometry();
                const auto& outsideGeometry = outside.geometry();

                const int numCornersOutsideElment = outsideGeometry.corners();
                const int numCornersIntersection = intersectionGeometry.corners();
                for (int i = 0; i < numCornersOutsideElment; ++i) {
                    for (int j = 0; j < numCornersIntersection; ++j) {
                        if (intersectionGeometry.corner(j) == outsideGeometry.corner(i) ) {
                            const int vertexIdx = gridView.indexSet().subIndex ( outside, i, dimDomain);
                            if (vertexIndices.find(vertexIdx) == vertexIndices.end())
                            {
                                //std::cout << "hangingNodesIndices "<< vertexIdx << std::endl;
                            }
                            vertexIndices.insert(vertexIdx);
                            continue;
                        }
                    }
                }
             }
        }

        // Find the elements sharing a vertex
        const int numVertices = vertexIndices.size();
        std::vector<std::vector<int>> vertexNeighbourSets(numVertices);

        std::set<int>::iterator it;
        int count = 0;
        for (it = vertexIndices.begin(); it != vertexIndices.end(); ++it)
        {
            const int vertexIdx = *it; // Note the "*" here

            int neighborIdx = 1;
            auto isIt = gridView.ibegin(element);
            for (; isIt != endIsIt; ++isIt) {
                const auto& intersection = *isIt;

                if (!intersection.neighbor())
                    continue;

                const Element& outside = intersection.outside();
                const int numVerticesOutside = outside.subEntities(dimDomain);
                for (int j = 0; j < numVerticesOutside; ++j) {
                    const int vertexIdxOutside = gridView.indexSet().subIndex ( outside, j, dimDomain);
                    // add neigbouring elements if they share a node
                    if(vertexIdx == vertexIdxOutside)
                        vertexNeighbourSets[count].emplace_back(neighborIdx);
                }

                ++ neighborIdx;
            }
            count ++;
        }

        // Store the triplets.
        // The current element always has index 0 and is not added to the list
        // i.e. we only store two elements int the tuple.
        std::unordered_set< std::tuple<int, int, int>, tuple_hash > unordered;
        for (int setIdx = 0; setIdx < numVertices; ++setIdx) {
            for (size_t i = 0; i < vertexNeighbourSets[setIdx].size(); i++) {
                for (size_t j = i+1; j < vertexNeighbourSets[setIdx].size(); j++) {
                    if (dimDomain == 2) {
                        std::tuple<int,int, int> tmp = std::make_tuple(vertexNeighbourSets[setIdx][i], vertexNeighbourSets[setIdx][j], -1);
                        if (std::get<0>(tmp) != std::get<1>(tmp) && std::get<0>(tmp) != std::get<2>(tmp) && std::get<1>(tmp) != std::get<2>(tmp))
                            unordered.insert(tmp);
                    } else if (dimDomain == 3) {
                        for (size_t k = j+1; k < vertexNeighbourSets[setIdx].size(); k++) {
                            std::tuple<int,int, int> tmp = std::make_tuple(vertexNeighbourSets[setIdx][i], vertexNeighbourSets[setIdx][j], vertexNeighbourSets[setIdx][k]);
                            if (std::get<0>(tmp) != std::get<1>(tmp) && std::get<0>(tmp) != std::get<2>(tmp) && std::get<1>(tmp) != std::get<2>(tmp))
                                unordered.insert(tmp);
                        }
                    }
                }
            }
        }
        triplets_[globalIdx] = std::vector<std::tuple<int,int,int>>(unordered.begin(), unordered.end());

        int numberOfTriplets = triplets_[globalIdx].size();
        invMatrix_[globalIdx].resize(numberOfTriplets);
        for (int i = 0; i < numberOfTriplets; ++i) {
            std::tuple<int,int,int> triplet = triplets_[globalIdx][i];
            Vector vec0;
            Vector vec1;
            Vector vec2;
            if (dimDomain == 2) {

            const auto& el0 = stencil.subControlVolume(0);
            const auto& el1 = stencil.subControlVolume(std::get<0>(triplet));
            const auto& el2 = stencil.subControlVolume(std::get<1>(triplet));


            for (int j = 0; j < dimDomain; ++j) {
                vec0[j] = el0.center()[j];
                vec1[j] = el1.center()[j];
                vec2[j] = el2.center()[j];
            }
            // add unit 3rd coordinate for 2d problems.
            if (dimDomain == 2) {
                vec0[2] = 1;
                vec1[2] = 1;
                vec2[2] = 1;
            }

            } else if (dimDomain == 3) {
                const auto& el0 = stencil.subControlVolume(0);
                const auto& el1 = stencil.subControlVolume(std::get<0>(triplet));
                const auto& el2 = stencil.subControlVolume(std::get<1>(triplet));
                const auto& el3 = stencil.subControlVolume(std::get<2>(triplet));
                for (int j = 0; j < dimDomain; ++j) {
                    vec0[j] = el1.center()[j]-el0.center()[j];
                    vec1[j] = el2.center()[j]-el0.center()[j];
                    vec2[j] = el3.center()[j]-el0.center()[j];
                }
            }

            Matrix C = Matrix( { vec0, vec1, vec2 } );

            try {
                C.invert(); }
            catch( ... ) {
                // If the centers are on the same line, the matrix is
                // not invertible, pertube to make it invertible.
                vec0[0] *=1.00001;
                vec1[1] *=0.99999;
                Matrix C = Matrix( { vec0, vec1, vec2 } );
                C.invert();
            }
            invMatrix_[globalIdx][i] = C;
        }
    }

    struct tuple_hash {
        inline std::size_t operator()(const std::tuple<int,int,int> & v) const {
            return std::get<0>(v)*31 + std::get<1>(v)*5+std::get<2>(v);
        }
    };

};
} // namespace Opm

#endif
