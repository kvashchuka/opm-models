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
 * \copydoc Opm::SofvStencil
 */
#ifndef EWOMS_SOFV_STENCIL_HH
#define EWOMS_SOFV_STENCIL_HH

#include <opm/models/utils/quadraturegeometries.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/geometry/type.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Opm {
/*!
 * \ingroup SofvDiscretization
 *
 * \brief Represents the stencil (finite volume geometry) of a single
 *        element in the SOFV discretization.
 *
 * The SOFV discretization is a element centered finite volume
 * approach. This means that each element corresponds to a control
 * volume.
 */
template <class TypeTag>
class SofvStencil
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dimWorld = GridView::dimensionworld };

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::Entity Element;
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
#endif

    //typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
    //                                                  Dune::MCMGElementLayout> ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef DofMapper  ElementMapper;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldVector<Scalar, dimWorld> WorldVector;

    const bool localRecons = EWOMS_GET_PARAM(TypeTag, bool, EnableLocalReconstruction);


public:
    typedef Element        Entity;
    typedef ElementMapper  Mapper;

    typedef typename Element::Geometry LocalGeometry;

    /*!
     * \brief Represents a sub-control volume.
     *
     * For element centered finite volumes, this is equivalent to the
     * element, in the vertex centered finite volume approach, this
     * corresponds to the intersection of a finite volume and the
     * grid element.
     */
    class SubControlVolume
    {
    public:
        // default construct an uninitialized object.
        // this is only here because std::vector needs it...
        SubControlVolume()
        {}

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        SubControlVolume(const Element& element)
            : element_(element)
#else
        SubControlVolume(const ElementPointer& elementPtr)
            : elementPtr_(elementPtr)
#endif
        { update(); }

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        void update(const Element& element)
        { element_ = element; }
#else
        void update(const ElementPointer& elementPtr)
        { elementPtr_ = elementPtr; }
#endif

        void update()
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            const auto& geometry = element_.geometry();
#else
            const auto& geometry = elementPtr_->geometry();
#endif
            centerPos_ = geometry.center();
            volume_ = geometry.volume();
        }

        /*!
         * \brief The global position associated with the sub-control volume
         */
        const GlobalPosition& globalPos() const
        { return centerPos_; }

        /*!
         * \brief The center of the sub-control volume
         */
        const GlobalPosition& center() const
        { return centerPos_; }

        /*!
         * \brief The volume [m^3] occupied by the sub-control volume
         */
        Scalar volume() const
        { return volume_; }

        /*!
         * \brief The geometry of the sub-control volume.
         */
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        const LocalGeometry geometry() const
        { return element_.geometry(); }
#else
        const LocalGeometry geometry() const
        { return elementPtr_->geometry(); }
#endif

    private:
        GlobalPosition centerPos_;
        Scalar volume_;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        Element element_;
#else
        ElementPointer elementPtr_;
#endif
    };

    /*!
     * \brief Represents a face of a sub-control volume.
     */
    class SubControlVolumeFace
    {
    public:
        SubControlVolumeFace()
        {}

        SubControlVolumeFace(const Intersection& intersection, unsigned localNeighborIdx)
        {
            exteriorIdx_ = static_cast<unsigned short>(localNeighborIdx);

            normal_ = intersection.centerUnitOuterNormal();

            const auto& geometry = intersection.geometry();
            integrationPos_ = geometry.center();
            area_ = geometry.volume();
        }

        /*!
         * \brief Returns the local index of the degree of freedom to
         *        the face's interior.
         */
        unsigned short interiorIndex() const
        {
            // The local index of the control volume in the interior
            // of a face of the stencil in the element centered finite
            // volume discretization is always the "central"
            // element. In this implementation this element always has
            // index 0....
            return 0;
        }

        /*!
         * \brief Returns the local index of the degree of freedom to
         *        the face's outside.
         */
        unsigned short exteriorIndex() const
        { return exteriorIdx_; }

        /*!
         * \brief Returns the global position of the face's
         *        integration point.
         */
        const GlobalPosition& integrationPos() const
        { return integrationPos_; }

        /*!
         * \brief Returns the outer unit normal at the face's
         *        integration point.
         */
        const WorldVector& normal() const
        { return normal_; }

        /*!
         * \brief Returns the area [m^2] of the face
         */
        Scalar area() const
        { return area_; }

    private:
        GlobalPosition integrationPos_;
        WorldVector normal_;

        Scalar area_;

        unsigned short exteriorIdx_;
    };

    SofvStencil(const GridView& gridView, const Mapper& mapper)
        : gridView_(gridView)
        , elementMapper_(mapper)
    { }

    void updateTopology(const Element& element)
    {
        auto isIt = gridView_.ibegin(element);
        const auto& endIsIt = gridView_.iend(element);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        // add the "center" element of the stencil
        subControlVolumes_.clear();
        subControlVolumes_.emplace_back(/*SubControlVolume(*/element/*)*/);
        elements_.clear();
        elements_.emplace_back(element);
#else
        // add the "center" element of the stencil
        ElementPointer ePtr(element);
        subControlVolumes_.clear();
        subControlVolumes_.emplace_back(/*SubControlVolume(*/ePtr/*)*/);
        elements_.clear();
        elements_.emplace_back(ePtr);
#endif
        interiorFaces_.clear();
        boundaryFaces_.clear();

        for (; isIt != endIsIt; ++isIt) {
            const auto& intersection = *isIt;
            // if the current intersection has a neighbor, add a
            // degree of freedom and an internal face, else add a
            // boundary face
            if (intersection.neighbor()) {
                elements_.emplace_back( intersection.outside() );
                subControlVolumes_.emplace_back(/*SubControlVolume(*/elements_.back()/*)*/);
                interiorFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, subControlVolumes_.size() - 1/*)*/);
            }
            else {
                boundaryFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, - 10000/*)*/);
            }
        }

        numPrimaryDof_ = elements_.size();
        if (localRecons) {
            // add layer 2
            //numPrimaryDof_ = elements_.size();
            for (size_t dofIdx = 1; dofIdx < numPrimaryDof_; ++dofIdx) {
                const Element& element = elements_[dofIdx];
                auto isIt = gridView_.ibegin(element);
                const auto& endIsIt = gridView_.iend(element);
                for (; isIt != endIsIt; ++isIt) {
                    const auto& intersection = *isIt;
                    // if the current intersection has a neighbor, add a
                    // degree of freedom and an internal face, else add a
                    // boundary face
                    if (intersection.neighbor()) {
                        if (newElement(intersection.outside())) {
                            elements_.emplace_back( intersection.outside() );
                            subControlVolumes_.emplace_back(/*SubControlVolume(*/elements_.back()/*)*/);
                            //interiorFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, subControlVolumes_.size() - 1/*)*/);
                        }
                    }
                    else {
                        //boundaryFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, - 10000/*)*/);
                    }
                }
            }
        } else {
            //numPrimaryDof_ = 1;
        }

        globalToLocal_.clear();
        const int nDofs = numDof();
        for( int dof = 0; dof<nDofs; ++dof )
        {
            globalToLocal_.insert( std::make_pair( globalSpaceIndex( dof ), dof ) );
        }

        /*
        differences_.clear();
        const auto& center = subControlVolumes_[ 0 ].center();

        const int nBnd = boundaryFaces_.size();
        differences_.reserve( nDofs + nBnd );
        for( int dof = 1; dof<nDofs; ++dof )
        {
            // store neighbor center - element center
            differences_.push_back( subControlVolumes_[ dof ].center() );
            differences_.back() -= center;
        }

        for( int bnd = 0; bnd < nBnd; ++bnd )
        {
            differences_.push_back( boundaryFaces_[ bnd ].integrationPos() );
            differences_.back() -= center;
        }
        */
    }

    void updatePrimaryTopology(const Element& element)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        // add the "center" element of the stencil
        subControlVolumes_.clear();
        subControlVolumes_.emplace_back(/*SubControlVolume(*/element/*)*/);
        elements_.clear();
        elements_.emplace_back(element);
#else
        // add the "center" element of the stencil
        ElementPointer ePtr(element);
        subControlVolumes_.clear();
        subControlVolumes_.emplace_back(/*SubControlVolume(*/ePtr/*)*/);
        elements_.clear();
        elements_.emplace_back(ePtr);
#endif

//        interiorFaces_.clear();
//        boundaryFaces_.clear();

//        auto isIt = gridView_.ibegin(element);
//        const auto& endIsIt = gridView_.iend(element);
//        for (; isIt != endIsIt; ++isIt) {
//            const auto& intersection = *isIt;
//            // if the current intersection has a neighbor, add a
//            // degree of freedom and an internal face, else add a
//            // boundary face
//            if (intersection.neighbor()) {
//                elements_.emplace_back( intersection.outside() );
//                subControlVolumes_.emplace_back(/*SubControlVolume(*/elements_.back()/*)*/);
//                interiorFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, subControlVolumes_.size() - 1/*)*/);
//            }
//            else {
//                boundaryFaces_.emplace_back(/*SubControlVolumeFace(*/intersection, - 10000/*)*/);
//            }
//        }
    }

    void update(const Element& element)
    {
        updateTopology(element);
    }

    void updateCenterGradients()
    {
        assert(false); // not yet implemented
    }

    /*!
     * \brief Return the element to which the stencil refers.
     */
    const Element& element() const
    { return element( 0 ); }

    /*!
     * \brief Return the grid view of the element to which the stencil
     *        refers.
     */
    const GridView& gridView() const
    { return *gridView_; }

    /*!
     * \brief Returns the number of degrees of freedom which the
     *        current element interacts with.
     */
    size_t numDof() const
    { return subControlVolumes_.size(); }

    /*!
     * \brief Returns the number of degrees of freedom which are contained
     *        by within the current element.
     *
     * Primary DOFs are always expected to have a lower index than
     * "secondary" DOFs.
     *
     * For element centered finite elements, this is only the central DOF.
     */
    size_t numPrimaryDof() const
    { return 1; }

    size_t numPrimaryDof2() const
    { return numPrimaryDof_; }

    /*!
     * \brief Return the global space index given the index of a degree of
     *        freedom.
     */
    unsigned globalSpaceIndex(unsigned dofIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof());

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        return static_cast<unsigned>(elementMapper_.index(element(dofIdx)));
#else
        return static_cast<unsigned>(elementMapper_.map(element(dofIdx)));
#endif
    }

    unsigned globalSpaceIndex(const Element& element) const
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        return static_cast<unsigned>(elementMapper_.index(element));
#else
        return static_cast<unsigned>(elementMapper_.map(element));
#endif
    }

    /*!
     * \brief Return partition type of a given degree of freedom
     */
    Dune::PartitionType partitionType(unsigned dofIdx) const
    { return elements_[dofIdx]->partitionType(); }

    /*!
     * \brief Return the element given the index of a degree of
     *        freedom.
     *
     * If no degree of freedom index is passed, the element which was
     * passed to the update() method is returned...
     */
    const Element& element(unsigned dofIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof());

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        return elements_[dofIdx];
#else
        return *elements_[dofIdx];
#endif
    }

    /*!
     * \brief Return the entity given the index of a degree of
     *        freedom.
     */
    const Entity& entity(unsigned dofIdx) const
    {
        return element( dofIdx );
    }

    /*!
     * \brief Returns the sub-control volume object belonging to a
     *        given degree of freedom.
     */
    const SubControlVolume& subControlVolume(unsigned dofIdx) const
    { return subControlVolumes_[dofIdx]; }

    /*!
     * \brief Returns the number of interior faces of the stencil.
     */
    size_t numInteriorFaces() const
    { return interiorFaces_.size(); }

    /*!
     * \brief Returns the face object belonging to a given face index
     *        in the interior of the domain.
     */
    const SubControlVolumeFace& interiorFace(unsigned bfIdx) const
    { return interiorFaces_[bfIdx]; }

    /*!
     * \brief Returns the number of boundary faces of the stencil.
     */
    size_t numBoundaryFaces() const
    { return boundaryFaces_.size(); }

    /*!
     * \brief Returns the boundary face object belonging to a given
     *        boundary face index.
     */
    const SubControlVolumeFace& boundaryFace(unsigned bfIdx) const
    { return boundaryFaces_[bfIdx]; }

    size_t globalToLocal( const size_t globalSpaceIdx ) const
    {
        //auto item = globalToLocal_.find( globalSpaceIdx );
        //assert( item != globalToLocal_.end() );
        //return item->second;
        return globalToLocal_.at( globalSpaceIdx );
    }

    size_t globalToLocal( const Element& element ) const
    {
        return globalToLocal( globalSpaceIndex( element ) );
    }

    /*
    const std::vector<GlobalPosition>& differences() const {
        return differences_;
    }
    */

private:

    bool newElement(const Element& element) const
    {
        for (auto& elem : elements_)
            if(globalSpaceIndex(element) == globalSpaceIndex(elem))
                return false;

        return true;
    }

protected:
    const GridView&       gridView_;
    const ElementMapper&  elementMapper_;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
    std::vector<Element> elements_;
#else
    std::vector<ElementPointer> elements_;
#endif
    std::map< size_t, size_t > globalToLocal_;

    std::vector<SubControlVolume>      subControlVolumes_;
    std::vector<SubControlVolumeFace>  interiorFaces_;
    std::vector<SubControlVolumeFace>  boundaryFaces_;

    //std::vector<GlobalPosition> differences_;

    size_t  numPrimaryDof_;
};

} // namespace Opm


#endif

