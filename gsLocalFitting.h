/** @file gsLocalFitting.h

    @brief Local fitting using tensor-product splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Giust
*/


#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo {

template<int d, typename T>
class gsLocalFitting {

public:

/**
    \brief
    This class applies local fitting to a given function.

    \tparam T coefficient type
*/

    //Default constructor
    gsLocalFitting(): bas(nullptr), geoResult(nullptr){
    }

    //Constructor
    gsLocalFitting(gsBasis<T>& basis, const gsFunctionExpr<T>& f);

    //Destructor
    virtual ~gsLocalFitting();

    //Compute bounding boxes used in the quasi-interpolant step
    virtual void computeBBoxes(gsMatrix<int>& lcMatrix, gsMatrix<int>& ucMatrix);

    //Perform local fitting
    void doLocalFitting();

    /// Gives back the computed approximation
    gsGeometry<T>* result() const {
        return geoResult; }

    //Compute L_inf error
    T computeError() const;

    //Iterative fitting
    virtual void iterativeFitting(int numIter, T threshold);

protected:

    //Initialize the local tensor product knot vector needed for each local basis
    virtual void initLocKnotVector(
        int i,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix,
        std::vector<gsKnotVector<T>>& knotVectorContainer
        );

    //Initialize local basis
    virtual typename gsBasis<T>::uPtr initLocBas(
        int i,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix,
        std::vector<gsKnotVector<T>>& knotVectorContainer
        );

    //Compute points in the physical space from the parametric one
    void computePoints(gsBasis<T>& basis, gsMatrix<T>& parP, gsMatrix<T>& phyP) const;

    //Assemble local linear system
    void assembleSystem(
        gsMatrix<T>& A_mat,
        gsMatrix<T>& B,
        gsBasis<T>& localBasis,
        const gsMatrix<T>& parP,
        const gsMatrix<T>& phyP
        );

    //Set computed coefficient
    virtual void setCoeff(
        int i,
        gsMatrix<T>& x,
        gsBasis<T>& localBasis,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix
        );

    //Perform all the numerical computations
    void compute(
        int i,
        gsBasis<T>& localBasis,
        const gsMatrix<T>& parP,
        const gsMatrix<T>& phyP,
        const gsMatrix<int>& lcMatrix,
        const gsMatrix<int>& ucMatrix
        );

    //Basis
    gsBasis<T>* bas;

    //The solution vector
    gsMatrix<T> solV;

    //Function
    gsFunctionExpr<T> func;

    //Resulting geometry
    gsGeometry<T>* geoResult;

};// class gsLocalFitting


}// namespace gismo


