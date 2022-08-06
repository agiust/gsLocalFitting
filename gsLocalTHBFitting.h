/** @file gsLocalTHBFitting.h

    @brief Adaptive local fitting using (truncated) hierarchical splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Giust
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsLocalFitting/gsLocalFitting.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsLocalFitting/gsLocalFittingUtils.h>

namespace gismo {

template<int d, typename T>
class gsLocalTHBFitting : public gsLocalFitting<d,T> {

/**
    \brief
    This class applies local hierarchical fitting to a given function.

    \tparam T coefficient type

    \ingroup HSplines
*/

public:

    // Default constructor
    gsLocalTHBFitting();

    //Constructor
    gsLocalTHBFitting(gsHTensorBasis<d,T>& basis, const gsFunctionExpr<T>& f);

    //Compute bounding boxes used by the quasi-interpolant
    void computeBBoxes(gsMatrix<int>& lcMatrix, gsMatrix<int>& ucMatrix);

    //Mark the element with error above a fixed threshold
    std::vector<bool> buildMarkVector(T threshold);

    //Refine the marked with a suitable extension
    void doRefinement(std::vector<bool> markedElements, int extension);

    //Perform the iterative fitting
    void iterativeFitting(int numIter, T threshold);


protected:

    //Initialize the local tensor product knot vector needed for each local basis
    void initLocKnotVector(
        int i,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix,
        std::vector<gsKnotVector<T>>& knotVectorContainer
        );

    //Initialize local basis
    typename gsBasis<T>::uPtr initLocBas(
        int i,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix,
        std::vector<gsKnotVector<T>>& knotVectorContainer
        );

    //Build local basis
    void buildLocBas(
        int i,
        gsTHBSplineBasis<d,T>& localBasis,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix
        );

    //Set coefficient for the solution
    void setCoeff(
        int i,
        gsMatrix<T>& x,
        gsBasis<T>& localBasis,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix
        );


}; // class gsLocalTHBFitting


}// namespace gismo


