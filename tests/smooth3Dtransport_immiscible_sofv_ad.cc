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
 * \brief Two-phase test for the immiscible model which uses the element-centered finite
 *        volume discretization in conjunction with automatic differentiation
 !*/

#include "config.h"

#include <opm/models/utils/start.hh>

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/discretization/sofv/sofvdiscretization.hh>
#include <opm/models/discretization/sofv/wenodiscretization.hh>

#include "problems/smooth3Dtransportproblem.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(SmoothTransportProblemSofvAd, INHERITS_FROM(ImmiscibleTwoPhaseModel, SmoothTransportBaseProblem));

// use the element centered finite volume spatial discretization
//SET_TAG_PROP(TransportProblemSofvAd, SpatialDiscretizationSplice, WenoDiscretization);
SET_TAG_PROP(SmoothTransportProblemSofvAd, SpatialDiscretizationSplice, SofvDiscretization);

// use automatic differentiation for this simulator
SET_TAG_PROP(SmoothTransportProblemSofvAd, LocalLinearizerSplice, AutoDiffLocalLinearizer);

END_PROPERTIES

int main(int argc, char **argv)
{
    typedef TTAG(SmoothTransportProblemSofvAd) ProblemTypeTag;
    return Opm::start<ProblemTypeTag>(argc, argv);
}

