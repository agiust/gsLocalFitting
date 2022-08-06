/** @file gsLocalTHBFitting.hpp

    @brief Provides implementation of the local fitting algorithm with thb splines.

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
#include <gsLocalFitting/gsLocalFitting.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsLocalFitting/gsLocalFittingUtils.h>


namespace gismo{

template<int d,typename T>
gsLocalTHBFitting<d,T>::gsLocalTHBFitting(gsHTensorBasis<d,T>& basis, const gsFunctionExpr<T>& f):
    gsLocalFitting<d,T>::gsLocalFitting(basis,f){};


template<int d,typename T>
void gsLocalTHBFitting<d,T>::computeBBoxes(gsMatrix<int>& lcMatrix, gsMatrix<int>& ucMatrix) {

    gsHTensorBasis<d, T>* basis = dynamic_cast< gsHTensorBasis<d,T>* > (this->bas);

    for (int i = 0; i < lcMatrix.cols(); ++i) {
        for (int j = 0; j <d ; ++j) {
            lcMatrix(j,i)=
                (basis->tensorLevel(basis->levelOf(i)).elementSupport(basis->flatTensorIndexOf(i)))(j,0)*(pow(2,basis->maxLevel()-basis->levelOf(i)));
            ucMatrix(j,i)=
                (basis->tensorLevel(basis->levelOf(i)).elementSupport(basis->flatTensorIndexOf(i)))(j,1)*(pow(2,basis->maxLevel()-basis->levelOf(i)));
        }
    }

}


template<int d,typename T>
void gsLocalTHBFitting<d,T>::initLocKnotVector(int i, const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix, std::vector<gsKnotVector<T>>& knotVectorContainer) {

    knotVectorContainer.clear();
    gsHTensorBasis<d, T>* basis = dynamic_cast< gsHTensorBasis<d,T>* > (this->bas);
    int maxL=basis->maxLevel();//max level of the global thb basis
    int minL=basis->levelOf(i);;//min level which overlaps my selected basis function

    for (int j = 0; j <d ; ++j) {

        gsKnotVector<T> knotV((basis->getBases()[basis->levelOf(i)])->knot(j,basis->maxDegree()+lcMatrix(j,i)/((1<<(maxL-minL)))),
                              (basis->getBases()[basis->levelOf(i)])->knot(j,basis->maxDegree()+ucMatrix(j,i)/((1<<(maxL-minL)))),
                              ucMatrix(j,i)/((1<<(maxL-minL)))-lcMatrix(j,i)/((1<<(maxL-minL)))-1,
                              basis->maxDegree() + 1);
        knotVectorContainer.push_back(knotV);

    }

};


template<int d,typename T>
void gsLocalTHBFitting<d,T>::buildLocBas(int i, gsTHBSplineBasis<d,T>& localBas,
    const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix) {

    gsHTensorBasis<d, T>* basis = dynamic_cast< gsHTensorBasis<d, T>* > (this->bas);

    int counter = 0;
    std::vector<std::vector<int>> knotVectorContainerMap;
    std::vector<int> indMap;

    //In indMap I express all the indices of the local basis
    // in terms of the indices of the finest level of the global THB basis
    for (int k = 0; k <d ; ++k) {
        for (int j = 0; j != ucMatrix(k, i) - lcMatrix(k, i) + 1; j++) {
            indMap.push_back(lcMatrix(k, i) + counter);
            counter = counter + 1;
        }

        knotVectorContainerMap.push_back(indMap);
        counter=0;
        indMap.clear();
    }

    //get the boxes
    gsMatrix<int> bLC;
    gsMatrix<int> bUC;
    gsVector<int> selLevel;
    basis->tree().getBoxes(bLC, bUC, selLevel);

    bool intersection;

    gsMatrix<int> box2Insert(2, d);
    gsMatrix<int> box2InsertMap(2, d);
    unsigned level2Insert;

    std::vector<unsigned> boxIn;
    int minL = basis->levelOf(i);
    int maxL=basis->maxLevel();
    std::vector<int> indMapLevel;
    std::vector< std::vector<int> > knotVectorContainerMapLevel;

    for (int ii = 0; ii != bLC.rows(); ii++) {

        //check if the box and the local basis intersect
        if (doOverlap(knotVectorContainerMap,bLC,bUC,ii))
        {
            intersection = true;
        } else { intersection = false; }

        //if the selected box and my local basis intersect, insert the box.
        //box2Insert is expressed in terms of the indices of the global basis
        if (intersection == true) {

            level2Insert = selLevel(ii) - minL;
            for (int k = 0; k <d ; ++k) {
                box2Insert(0, k)=
                    std::max(knotVectorContainerMap[k][0], bLC(ii, k)) / (1 << (maxL - selLevel(ii, 0)));

                box2Insert(1, k) =
                    std::min(knotVectorContainerMap[k][knotVectorContainerMap[k].size() - 1], bUC(ii, k)) /
                    (1 << (maxL - selLevel(ii, 0)));
            }

            counter=0;
            knotVectorContainerMapLevel.clear();
            indMapLevel.clear();

            for (int k = 0; k <d ; ++k) {
                for (int j = 0; j != ucMatrix(k, i) / (1 << (maxL - selLevel(ii, 0))) -
                    lcMatrix(k, i) / (1 << (maxL - selLevel(ii, 0))) + 1; j++) {
                    indMapLevel.push_back(lcMatrix(k, i) / (1 << (maxL - selLevel(ii, 0))) + counter);
                    counter = counter + 1;
                }
                knotVectorContainerMapLevel.push_back(indMapLevel);
                indMapLevel.clear();
                counter=0;
            }


            for (int k = 0; k <d ; ++k) {
                for (int r = 0; r != knotVectorContainerMapLevel[k].size(); r++) {
                    if (knotVectorContainerMapLevel[k][r] == box2Insert(0, k)) {
                        box2InsertMap(0, k) = r;
                    }
                    if (knotVectorContainerMapLevel[k][r] == box2Insert(1, k)) {
                        box2InsertMap(1, k) = r;
                    }

                }
            }

            std::vector<int> boxIn;
            boxIn.push_back(level2Insert);
            for (int k = 0; k <2 ; ++k) {
                for (int j = 0; j <d ; ++j) {
                    boxIn.push_back(box2InsertMap(k, j));

                }

            }

            localBas.refineElements(boxIn);
            boxIn.clear();

        }

    }

}


template<int d,typename T>
typename gsBasis<T>::uPtr gsLocalTHBFitting<d,T>::initLocBas(int i, const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix,
    std::vector<gsKnotVector<T>>& knotVectorContainer){

    //initialize local knot vector
    initLocKnotVector(i,lcMatrix,ucMatrix,knotVectorContainer);

    //initialize local tp basis
    gsTensorBSplineBasis<d,T>oldBas(knotVectorContainer);

    //initialize local thb basis
    gsTHBSplineBasis<d,T> newBas(oldBas);

    //build local thb basis
    buildLocBas(i,newBas,lcMatrix,ucMatrix);

    //get final basis
    typename gsBasis<T>::uPtr localBas(new gsTHBSplineBasis<d,T>(newBas));

    return localBas;

}


template<int d,typename T>
void gsLocalTHBFitting<d,T>::setCoeff(int i, gsMatrix<T>& x, gsBasis<T>& localBasis,
    const gsMatrix<int>& lcMatrix, const gsMatrix<int>& ucMatrix) {

    //determine the correct coefficient for the selected basis function
    gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->bas);
    gsTHBSplineBasis<d, T>* newLocalBasis = static_cast< gsTHBSplineBasis<d,T>* > (&localBasis);

    int maxL=basis->maxLevel();
    int functionLevel=basis->levelOf(i);
    gsMatrix<int> functionSupport(d,2);

    //put in functionSupport the support of the function expressed in indices of the global thb basis
    gsMatrix<int>flat;
    flat=basis->tensorLevel(functionLevel).elementSupport(basis->flatTensorIndexOf(i));
    functionSupport=flat;

    int counter=0;
    std::vector<std::vector<int>>knotVectorContainerMapLevel;
    knotVectorContainerMapLevel.clear();
    std::vector<int>indMapLevel;
    indMapLevel.clear();

    for (int k = 0; k <d ; ++k) {
        for (int j = 0; j != ucMatrix(k, i) / (1 << (maxL - functionLevel)) -
            lcMatrix(k, i) / (1 << (maxL - functionLevel)) + 1; j++) {
            indMapLevel.push_back(lcMatrix(k, i) / (1 << (maxL -functionLevel)) + counter);
            counter = counter + 1;
        }
        knotVectorContainerMapLevel.push_back(indMapLevel);
        indMapLevel.clear();
        counter=0;
    }

    int xCoeff=0;
    gsMatrix<real_t>anchorsGlobal=basis->anchors();
    gsMatrix<real_t>anchorsLocal=newLocalBasis->anchors();

    for(int t=0;t!=newLocalBasis->size();t++){

        gsMatrix<int> functionSupportLocal(d,2);
        gsMatrix<int>flatLocal(d,2);
        flatLocal=
            (newLocalBasis->tensorLevel(newLocalBasis->levelOf(t)).elementSupport(newLocalBasis->flatTensorIndexOf(t)));
        //put in functionSupportLocal the support of the function expressed in indices of the local thb basis
        functionSupportLocal=flatLocal;

        gsMatrix<int>comp(d,2);
        for (int k = 0; k <d ; ++k) {
            comp(k,0)=knotVectorContainerMapLevel[k][functionSupportLocal(k, 0)];
            comp(k,1)=knotVectorContainerMapLevel[k][functionSupportLocal(k, 1)];
        }

        if ( (anchorsGlobal.col(i) - anchorsLocal.col( t)).norm() < 1e-14 && functionSupport==comp){
            xCoeff = t;
            break;
        }

    }

    this->solV.row(i)= x.row(xCoeff);

}


template<int d,typename T>
std::vector<bool> gsLocalTHBFitting<d,T>::buildMarkVector(double threshold){

    gsHTensorBasis<d, T>* basis = dynamic_cast< gsHTensorBasis<d,T>* > (this->bas);

    gsHDomainIterator<real_t,d> domIter(*basis);//Create domain iterator
    int numberOfElements = this->bas->numElements();//Find number of elements
    gsGaussRule<> gr(*basis,1,1);
    std::vector<bool> markedElements;
    markedElements.reserve(numberOfElements);

    for (; domIter.good(); domIter.next())
    {
        gsVector<> v;
        gsMatrix<> pts, vectorPoints(d, gr.numNodes());
        gr.mapTo(domIter.lowerCorner(), domIter.upperCorner(), pts, v);
        gsMatrix<> valuePointsFitting(1, vectorPoints.cols());
        gsMatrix<> valuePoints(1, vectorPoints.cols());
        gsMatrix<> errorVector(1, vectorPoints.cols());
        for(index_t i = 0; i!=pts.cols();++i) {
            vectorPoints.col(i) = pts.col(i);
        }

        //Evaluate my approximated function in the nodes of the current element
        this->geoResult->eval_into(vectorPoints, valuePointsFitting);
        //Evaluate my known function in the nodes of the current element
        this->func.eval_into(vectorPoints, valuePoints);
        //Evaluate error
        errorVector=(valuePoints-valuePointsFitting.row(d)).cwiseAbs();

        bool toMark= ( (errorVector.array() > threshold ).any() );

        markedElements.push_back( toMark );//marking
    }

    return markedElements;

}


template<int d,typename T>
void gsLocalTHBFitting<d,T>::doRefinement(std::vector<bool> markedElements, int extension){
    //Do the refinement
    // numMarked: Number of marked cells on current patch (1 patch in this case), also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<double_t> refBoxes;

    const int numEl = this->bas->numElements();
    numMarked = std::count_if(markedElements.begin() + poffset,
                              markedElements.begin() + poffset + numEl,
                              std::bind2nd(std::equal_to<bool>(), true));

    refBoxes.resize(d, 2 * numMarked);

    numMarked = 0;

    typename gsBasis<double_t>::domainIter domainIter = this->bas->makeDomainIterator();
    for (; domainIter->good(); domainIter->next()) {
        if (markedElements[globalCount++])
        {
            // Construct degenerate box by setting both
            // corners equal to the center
            refBoxes.col(2 * numMarked) = domainIter->centerPoint();
            refBoxes.col(2 * numMarked + 1) = domainIter->centerPoint();

            // Advance marked cells counter
            numMarked++;

        }

    }

    // Refine all the found refBoxes
    this->bas->refine(refBoxes, extension);

}


template<int d,typename T>
void gsLocalTHBFitting<d,T>::iterativeFitting(int numIter, real_t threshold) {

    for (int i = 0; i < numIter; ++i) {
        gsLocalFitting<d, T>::doLocalFitting();
        double error=gsLocalFitting<d, T>::computeError();
        if (error<threshold){
            std::cout<<"you have reached the threshold"<<std::endl;
            break;
        }
        std::vector<bool> markedElements = this->buildMarkVector(threshold);
        this->doRefinement(markedElements, 1);
    }

}


}// namespace gismo
