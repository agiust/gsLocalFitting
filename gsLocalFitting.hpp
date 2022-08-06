/** @file gsLocalFitting.hpp

    @brief Provides implementation of the local fitting algorithm with tensor-product splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Giust

*/


#include <gsCore/gsBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsMatrix/gsMatrix.h>


namespace gismo{

template<int d,typename T>
gsLocalFitting<d,T>::gsLocalFitting(gsBasis<T>& basis, const gsFunctionExpr<T>& f):
    bas(&basis), func(f), geoResult(nullptr){solV.setZero(basis.size(),d+1);};


template<int d,typename T>
gsLocalFitting<d,T>::~gsLocalFitting()
{
    if ( geoResult )
        delete geoResult;
};


template<int d,typename T>
void gsLocalFitting<d,T>::computeBBoxes(gsMatrix<int>& lcMatrix, gsMatrix<int>& ucMatrix) {

    gsTensorBSplineBasis<d, T>* basis = dynamic_cast< gsTensorBSplineBasis<d,T>* > (this->bas);

    for (int i = 0; i < lcMatrix.cols(); ++i) {
        for (int j = 0; j <d ; ++j) {
            lcMatrix(j,i)=basis->elementSupport(i)(j,0);
            ucMatrix(j,i)=basis->elementSupport(i)(j,1);
        }
    }

}


template<int d,typename T>
void gsLocalFitting<d,T>::doLocalFitting() {

    gsMatrix<int> lcMatrix(d, bas->size());
    gsMatrix<int> ucMatrix(d, bas->size());
    computeBBoxes(lcMatrix,ucMatrix);
    solV.setZero(bas->size(),d+1);

    for (int i = 0; i <bas->size() ; ++i) {

        //compute local knot vector
        std::vector<gsKnotVector<T>>knotVectorContainer;

        //init local basis
        typename gsBasis<T>::uPtr localBasis=initLocBas(i, lcMatrix, ucMatrix, knotVectorContainer);

        //compute local points
        gsMatrix<T>parP(d, localBasis->numElements()*pow(localBasis->degree(0)+1,d));
        gsMatrix<T>phyP(parP.cols(), parP.rows()+1);
        computePoints(*localBasis,parP,phyP);

        //solve and set coeff
        compute(i, *localBasis, parP, phyP, lcMatrix, ucMatrix);

    }

    //set result
    gsMatrix<T> xFinal;
    xFinal=solV;

    geoResult = bas->makeGeometry( give(xFinal) ).release();

};


template<int d,typename T>
void gsLocalFitting<d,T>::initLocKnotVector(int i, const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix,
    std::vector<gsKnotVector<T>>& knotVectorContainer) {

    knotVectorContainer.clear();

    for (int j = 0; j <d ; ++j) {

        gsTensorBSplineBasis<d, T>* basis = dynamic_cast< gsTensorBSplineBasis<d,T>* > (this->bas);

        gsKnotVector<T> knotV(basis->knot(j, basis->degree(0) + lcMatrix(j, i) ),
                              basis->knot(j, basis->degree(0) + ucMatrix(j, i) ),
                              ucMatrix(j, i)  - lcMatrix(j, i) - 1,
                              basis->degree(0) + 1);

        knotVectorContainer.push_back(knotV);

    }

};


template<int d,typename T>
typename gsBasis<T>::uPtr gsLocalFitting<d,T>::initLocBas(int i, const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix,
    std::vector<gsKnotVector<T>>& knotVectorContainer){

    //initialize local knot vector
    initLocKnotVector(i, lcMatrix, ucMatrix,knotVectorContainer);

    //initialize local tp basis
    gsTensorBSplineBasis<d,T>newBasis(knotVectorContainer);

    //get basis
    typename gsBasis<T>::uPtr newLocalBasis(new gsTensorBSplineBasis<d,T>(newBasis));

    return newLocalBasis;

}


template<int d,typename T>
void gsLocalFitting<d,T>::computePoints(gsBasis<T>& basis, gsMatrix<T>& parP, gsMatrix<T>& phyP) const{

    typename gsBasis<T>::domainIter domIter = basis.makeDomainIterator();

    gsGaussRule<T> gr(basis, 1, 1);
    gsVector<T> v;
    gsMatrix<T> pts,vectorPoints(d, gr.numNodes() * basis.numElements());

    int counter = 0;
    for (; domIter->good(); domIter->next())
    {
        gr.mapTo(domIter->lowerCorner(), domIter->upperCorner(), pts, v);
        for (index_t i = 0; i != pts.cols(); ++i)
            vectorPoints.col(counter++) = pts.col(i);
    }

    parP=vectorPoints;
    for (int l = 0; l < parP.rows(); ++l) {
        for (int i = 0; i <parP.cols() ; ++i) {
            phyP(i,l)=parP(l,i);
        }
    }

    gsMatrix<T>point(d,1);
    for (int j = 0; j <phyP.rows() ; ++j) {
        for (int i = 0; i <d ; ++i) {
            point(i,0)=phyP(j,i);
        }

        phyP(j,d)=func.eval(point)(0,0);

    }

};


template<int d,typename T>
void gsLocalFitting<d,T>::assembleSystem(gsMatrix<T>& A_mat, gsMatrix<T>& B, gsBasis<T>& localBasis,
    const gsMatrix<T>& parP, const gsMatrix<T>& phyP) {

    //Assemble local least squares
    const int num_points = phyP.rows();

    gsMatrix<T> value, curr_point;
    gsMatrix<index_t> actives;

    for(index_t k = 0; k != num_points; ++k)
    {
        curr_point = parP.col(k);
        localBasis.eval_into(curr_point, value);
        localBasis.active_into(curr_point, actives);
        const index_t numActive = actives.rows();

        for (index_t i = 0; i != numActive; ++i)
        {
            const index_t ii = actives.at(i);
            B.row(ii) += value.at(i) * (phyP.row(k));
        }
    }

    gsMatrix<T>A(phyP.rows(),localBasis.size());
    gsMatrix<T>point(d,1);

    for (int j = 0; j <phyP.rows() ; ++j) {
        for (int k = 0; k <d ; ++k) {
            point(k,0)=phyP(j,k);
        }
        for (int i = 0; i <localBasis.size() ; ++i) {
            A(j,i)=localBasis.evalSingle(i,point)(0,0);
        }
    }

    A_mat=A.transpose()*A;

};


template<int d,typename T>
void gsLocalFitting<d,T>::compute(int i, gsBasis<T>& localBasis, const gsMatrix<T>& parP, const gsMatrix<T>& phyP, const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix)  {

    const int num_basis=localBasis.size();
    const int dimension=phyP.cols();

    //left side matrix
    gsMatrix<T> A_mat(num_basis , num_basis );
    A_mat.setZero();

    //right side vector
    gsMatrix<T> m_B(num_basis , dimension);
    m_B.setZero();

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b
    assembleSystem(A_mat, m_B,localBasis,parP,phyP);

    // solves for many right hand side columns
    gsMatrix<T> x;
    x = A_mat.partialPivLu() .solve(m_B);

    //set coeff
    setCoeff(i,x,localBasis,lcMatrix,ucMatrix);

}


template<int d,typename T>
void gsLocalFitting<d,T>::setCoeff(int i, gsMatrix<T>& x, gsBasis<T>& localBasis, const gsMatrix<int>& lcMatrix, const gsMatrix<int>& matrixOfTPUC) {

    gsTensorBSplineBasis<d, T>* basis = dynamic_cast< gsTensorBSplineBasis<d,T>* > (this->bas);
    gsTensorBSplineBasis<d, T>* newLocalBasis = dynamic_cast< gsTensorBSplineBasis<d,T>* > (&localBasis);

    //Determine the correct coefficient for the selected basis function
    gsMatrix<int> functionSupport(d,2);

    for (int k = 0; k <d ; ++k) {
        for (int j = 0; j <2 ; ++j) {
            functionSupport(k,j)=basis->elementSupport(i)(k,j);
        }
    }

    int nCoeff=0;
    gsMatrix<T>anchorsGlobal=bas->anchors();
    gsMatrix<T>anchorsLocal=localBasis.anchors();

    for(int t=0;t<localBasis.size();t++){

        gsMatrix<int> functionSupportLocal(d,2);

        for (int k = 0; k <d ; ++k) {
            for (int j = 0; j <2 ; ++j) {
                functionSupportLocal(k,j)=newLocalBasis->elementSupport(t)(k,j)+lcMatrix(k,i);
            }
        }
        if (functionSupport==functionSupportLocal && (anchorsGlobal.col(i)-anchorsLocal.col(t)).norm()<1e-14){

            nCoeff=t;

            break;
        }
    }

    solV.row(i)= x.row(nCoeff);

}


template<int d,typename T>
real_t gsLocalFitting<d,T>::computeError() const {

    gsMatrix<T>parP(d, bas->numElements()*pow(bas->degree(0)+1,d));
    gsMatrix<T>phyP(parP.cols(), parP.rows()+1);

    computePoints(*bas,parP,phyP);

    gsMatrix<T>exactValues;
    geoResult->eval_into(parP, exactValues);
    gsMatrix<T> errorVector;

    errorVector=(exactValues.row(d)-phyP.col(d).transpose()).cwiseAbs();
    real_t error = errorVector.maxCoeff();
    std::cout<<"The max error is "<<error<<std::endl;

    return error;

}


template<int d,typename T>
void gsLocalFitting<d,T>::iterativeFitting(int numIter, real_t threshold) {

    for (int i = 0; i <numIter ; ++i) {

        doLocalFitting();
        T error=computeError();

        if (error<threshold){
            std::cout<<"you have reached the threshold"<<std::endl;
            break;
        }

        bas->uniformRefine();

    }

}


}// namespace gismo
