/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthImplexContact.cpp,v $
// $Revision: 1.0 $
// $Date: 2020-XX-XX XX:XX:XX $
                                                                        
// Written: Onur Deniz Akan         (onur.akan@iusspavia.it)
//          Dr. Massimo Petracca    (m.petracca@asdea.net)
//          Prof. Guido Camata      (g.camata@unich.it)
//
// Created: May 2020
//
// Description: This file contains the implementation for the ZeroLengthImplexContact class.


#include "ZeroLengthImplexContact.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>

#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>
#include <elementAPI.h>


// static data
const int ZeroLengthImplexContact::numNDS = 2;


// initialize the class wide variables
Matrix ZeroLengthImplexContact::zlcSM4(4,4);
Matrix ZeroLengthImplexContact::zlcSM6(6,6);
Matrix ZeroLengthImplexContact::zlcSM12(12,12);
Matrix ZeroLengthImplexContact::zlcTM4(4,4);
Matrix ZeroLengthImplexContact::zlcTM6(6,6);
Matrix ZeroLengthImplexContact::zlcTM12(12,12);
Vector ZeroLengthImplexContact::zlcRV4(4);
Vector ZeroLengthImplexContact::zlcRV6(6);
Vector ZeroLengthImplexContact::zlcRV12(12);
Vector ZeroLengthImplexContact::zlcDV4(4);
Vector ZeroLengthImplexContact::zlcDV6(6);
Vector ZeroLengthImplexContact::zlcDV12(12);

// Remaining Steps To Do:
/* List:

9) write function sendSelf

10) write function recvSelf

11) fix impl-ex

12) ask and correct the signs of residual and stiffness

*/

// class utilities
namespace utils {

    // 2D Tools
    inline Vector make2DVector(double x, double y) {

        Vector v(2);
        v(0) = x;
        v(1) = y;
        return v;
    }

    inline double norm2D(Vector& A) {

        double norm = sqrt(A(0)*A(0) + A(1) * A(1));
        return norm;
    }

    inline void normalize2D(Vector& A) {
        
        double norm = utils::norm2D(A);
        A(1) = A(1) / norm;
        A(2) = A(2) / norm;
    }

    inline double multiply2D(Vector& A, Vector& B) {

        double result = 0;
        double length = A.Size();

        for (int i = 0; i < length; i++) {
            result += A(i) * B(i);
        }
        return result;
    }

    // 3D Tools
    inline Vector make3DVector(double x, double y, double z) {

        Vector v(3);
        v(0) = x;
        v(1) = y;
        v(2) = z;
        return v;
    }

    inline void cross(const Vector& A, const Vector& B, Vector& C) {

        C(0) = A(1) * B(2) - A(2) * B(1);
        C(1) = A(2) * B(0) - A(0) * B(2);
        C(2) = A(0) * B(1) - A(1) * B(0);
    }

    // Simple math tools
    inline int sign(double& num) {

        if (num < 0.0)
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }
}


void* OPS_ZeroLengthImplexContact(void) {

    double SmallNumber = 1.0e-6;
    Element* theElement = 0;

    // some kudos
    static int counter = 0;
    if (++counter == 1)
        opserr << "ZeroLengthImplexContact element - Written: O.D.Akan (IUSS), M.Petracca (ASDEA), G.Camata (UNICH)\n";
    
    // model dimension
    int ndm = OPS_GetNDM();

    // a quick check on number of args
    if (OPS_GetNumRemainingInputArgs() < 7) {
	opserr << "ZeroLengthImplexContact: WARNING: too few arguments \n" <<
	    "want - element zeroLengthImplexContact eleTag? iNode? jNode? Kn? Kt? mu? C? <-orient $x1 $x2 $x3> <-intType type?>\n";
	return 0;
    }
    
    // start with mandatory inputs
    // read eleTag, iNode, jNode
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "ZeroLengthImplexContact: WARNING: invalid int inputs\n";
	return 0;
    }

    // read Kn, Kt, mu, C
    double ddata[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "ZeroLengthImplexContact: WARNING: invalid double inputs\n";
	return 0;
    }

    // continue with optional inputs
    Vector x_e(3); x_e(0) = 1.0; x_e(1) = 0.0; x_e(2) = 0.0;                        // initialize orientation vector
    int integrationType = 0;                                                        // implicit by default
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* inputstring = OPS_GetString();
        if (strcmp(inputstring, "-orient") == 0) {                                  // #1 read global element orientation
            if (ndm == 2) {
                if (OPS_GetNumRemainingInputArgs() < 2) {
                    opserr << "ZeroLengthImplexContact: WARNING: insufficient orient values in 2D\n";
                    return 0;
                }
                numdata = 3;
                if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
                    opserr << "ZeroLengthImplexContact: WARNING: invalid double input after -orient\n";
                    return 0;
                }
            }
            else if (ndm == 3) {
                if (OPS_GetNumRemainingInputArgs() < 3) {
                    opserr << "ZeroLengthImplexContact: WARNING: insufficient orient values in 3D\n";
                    return 0;
                }
                numdata = 3;
                if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
                    opserr << "ZeroLengthImplexContact: WARNING: invalid double input after -orient\n";
                    return 0;
                }
            }
            else {
                opserr << "ZeroLengthImplexContact: WARNING: -orient: model dimension is invalid! \n";
                return 0;
            }
        }
        else if (strcmp(inputstring, "-intType") == 0) {                             // #2 read type of integration 
            numdata = 1;
            if (OPS_GetIntInput(&numdata, &integrationType) < 0) {
                opserr << "ZeroLengthImplexContact: WARNING: invalid integer after -intType\n";
                return 0;
            }
        }
    }
    // input reading stage is complete

    // check integration type and revert back to implicit if neither 1 or 0
    if (integrationType != 1 && integrationType != 0) {                             
        opserr << "ZeroLengthImplexContact: WARNING: type of integration is set to IMPLICIT due to invalid flag\n";
        integrationType = 0;
    }
    // check the normal vector and normalize
    if (x_e.Norm() < SmallNumber) {
        opserr << "ZeroLengthImplexContact: WARNING: normal vector is NOT valid!: -orient $x1 $x2 $x3 cannot be < 0, 0, 0 >\n";
        return 0;
    }
    x_e.Normalize();   // crucial step for the computaion of geom. transformation matix: method computeRotMatrix()

    // finally, create the element
    theElement = new ZeroLengthImplexContact(idata[0], idata[1], idata[2], ddata[0], ddata[1], 
        ddata[2], ddata[3], ndm, integrationType, x_e[0], x_e[1], x_e[2]);

    if (theElement == 0) {
        opserr << "WARNING: out of memory: element zeroLengthImplexContact " << idata[0] <<
            " iNode? jNode? Kn? Kt? mu? C? <-orient $x1 $x2 $x3> <-intType type?>\n";
    }
    
    return theElement;
}


// public:
    // constructor
ZeroLengthImplexContact::ZeroLengthImplexContact(int tag, int Nd1, int Nd2,
                    double Kn, double Kt, double fcoeff, double C,
                    int ndm, bool itype, double xN, double yN, double zN)
    :Element(tag,ELE_TAG_ZeroLengthImplexContact),
    connectedExternalNodes(numNDS),
    Knormal(Kn), Kfriction(Kt), mu(fcoeff), cohesion(C),
    numDIM(ndm), numDOF(0), doImplEx(itype), Xorient(3), gap(0),
    theNodes(), theDispVector(0), theResidVector(0), theStiffMatrix(0), theTransMatrix(0) {

  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)
    opserr << "FATAL ZeroLengthImplexContact::fullConstructor() - failed to create an ID of correct size\n";

  for (int j = 0; j < numNDS; j++) {
      theNodes[j] = 0;
  }
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  Xorient(0) = xN;
  Xorient(1) = yN;
  Xorient(2) = zN;
  }


    // null constructor
ZeroLengthImplexContact::ZeroLengthImplexContact(void)
  :Element(0,ELE_TAG_ZeroLengthImplexContact),     
    connectedExternalNodes(numNDS),
    Knormal(0), Kfriction(0), mu(0), cohesion(0),
    numDIM(0), numDOF(0), doImplEx(0), Xorient(3), gap(0),
    theNodes(), theDispVector(0), theResidVector(0), theStiffMatrix(0), theTransMatrix(0) {

  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2)
    opserr << "FATAL ZeroLengthImplexContact::nullConstructor() - failed to create an ID of correct size\n";

  for (int j=0; j<numNDS; j++ ) {
    theNodes[j] = 0;
  }
}


    // destructor
ZeroLengthImplexContact::~ZeroLengthImplexContact(void) {

    // delete pointers
    for (int i = 0; i < numNDS; i++) {
        if (theNodes[i])
            delete theNodes[i];
    }

    if (theDispVector != 0)
        delete theDispVector;

    if (theResidVector != 0)
        delete theResidVector;

    if (theStiffMatrix != 0)
        delete theStiffMatrix;

    if (theTransMatrix != 0)
        delete theTransMatrix;
}


    // public methods to obtain information about dof & connectivity
int ZeroLengthImplexContact::getNumExternalNodes(void) const {

return numNDS;
}


const ID & ZeroLengthImplexContact::getExternalNodes(void)  {

  return connectedExternalNodes;
}


Node ** ZeroLengthImplexContact::getNodePtrs(void) {

  return theNodes;
}


int ZeroLengthImplexContact::getNumDOF(void) {

  return numDOF;
}


void ZeroLengthImplexContact::setDomain(Domain *theDomain) {

   // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
        if (theNodes[0] == 0) {
            opserr << "WARNING ZeroLengthImplexContact::setDomain() - Nd1: " << Nd1 <<
                " does not exist in the model\n";
        }
        else {
            opserr << "WARNING ZeroLengthImplexContact::setDomain() - Nd2: " << Nd2 <<
                " does not exist in the model\n";
        }
        // fill this in so don't segment fault later
        numDOF = 4;
        theStiffMatrix = &zlcSM4;
        theResidVector = &zlcRV4;
        theTransMatrix = &zlcTM4;
        theDispVector = &zlcDV4;
        return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if different dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
        opserr << "WARNING ZeroLengthImplexContact::setDomain(): nodes " << Nd1 << " and " << Nd2 << 
            "have different dof at element ends" << this->getTag() << endln;
        // fill this in so don't segment fault later
        numDOF = 4;
        theStiffMatrix = &zlcSM4;
        theResidVector = &zlcRV4;
        theTransMatrix = &zlcTM4;
        theDispVector = &zlcDV4;
        return;
    }	
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointer
    if (numDIM == 2 && dofNd1 == 2) {
        numDOF = 4;
        theStiffMatrix = &zlcSM4;
        theResidVector = &zlcRV4;
        theTransMatrix = &zlcTM4;
        theDispVector = &zlcDV4;
    }
    else if (numDIM == 2 && dofNd1 == 3) {
        numDOF = 6;
        theStiffMatrix = &zlcSM6;
        theResidVector = &zlcRV6;
        theTransMatrix = &zlcTM6;
        theDispVector = &zlcDV6;
    }
    else if (numDIM == 3 && dofNd1 == 3) {
        numDOF = 6;
        theStiffMatrix = &zlcSM6;
        theResidVector = &zlcRV6;
        theTransMatrix = &zlcTM6;
        theDispVector = &zlcDV6;
    }
    else if (numDIM == 3 && dofNd1 == 6) {
        numDOF = 12;
        theStiffMatrix = &zlcSM12;
        theResidVector = &zlcRV12;
        theTransMatrix = &zlcTM12;
        theDispVector = &zlcDV12;
    }
    else {
        opserr << "WARNING ZeroLengthImplexContact::setDomain cannot handle " << numDIM << " dofs at nodes in " <<
            dofNd1 << " problem\n";

        numDOF = 4;
        theStiffMatrix = &zlcSM4;
        theResidVector = &zlcRV4;
        theTransMatrix = &zlcTM4;
        theDispVector = &zlcDV4;
        return;
    }

    // finally, compute the initial gap vector and check if length is zero within tolerance
    Vector& Ndloc = *theDispVector;                                 // dump location data to disp. vector temporarily
    const Vector& Nsl = theNodes[0]->getCrds();                     // slave node coords.
    const Vector& Nms = theNodes[1]->getCrds();                     // master node coords.
    const Vector& UgS = theNodes[0]->getTrialDisp();                // slave node init. disp.
    const Vector& UgM = theNodes[1]->getTrialDisp();                // master node init. disp.
    // get node coords. + displacements and put them into a vector
    int slave = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
        Ndloc(i)        = Nms(i) + UgM(i);
        Ndloc(i+slave)  = Nsl(i) + UgS(i);
    }
    // compute global to local coord. sys. transformation matrix
    computeRotMatrix();
    // rotate displacement and node coords. from global to local coord. sys.
    Matrix& R = *theTransMatrix;
    Vector temp(numDOF);
    temp = Ndloc;
    Ndloc.addMatrixVector(0.0, R, temp, 1.0);
    // compute the initial gap in local coord. sys.
    Vector N2Ndist(3); // node-to-node distance in local coord. sys.
    for (int i = 0; i < numDIM; i++) {
        N2Ndist(i) = Ndloc(i+slave)-Ndloc(i);  //gap = slave - master
    }
    // save the initial gap in normal dir.
    gap = N2Ndist(0);
    // print the length of element if not zero
    double Len = N2Ndist.Norm();
    if (Len > LENTOL) {
        opserr << "ATTENTION ZeroLengthImplexContact::setDomain(): Element " << this->getTag() << " has L = " << Len << "\n";
    }
    // done
}   	 


    // public methods to set the state of the element
int ZeroLengthImplexContact::commitState(void) {
    
    // do the implicit correction if impl-ex
    if (doImplEx) {
        // get trial strain tensor at integration point in local crd. sys.
        computeElementStuff();
        // update material internal variables
        computeMaterialStuff(false, false);  // explicit_phase?, do_tangent?
    }
    // store converged material state variables [step (n) --> (n+1)]
    sv.eps_commit = sv.eps;
    sv.sig_commit = sv.sig;
    sv.beta_commit = sv.beta;
    sv.lambda_commit = sv.lambda;
    sv.dtime_n_commit = sv.dtime_n;
    // store previously committed variables [step (n-1) --> (n)]
    sv.beta_commit_old = sv.beta_commit;
    sv.lambda_commit_old = sv.lambda_commit;
    // done
    return 0;
}


int ZeroLengthImplexContact::revertToLastCommit(void) {

    // restore converged values [step (n+1) --> (n)]
    sv.eps = sv.eps_commit;
    sv.sig = sv.sig_commit;
    sv.beta = sv.beta_commit;
    sv.lambda = sv.lambda_commit;
    sv.dtime_n = sv.dtime_n_commit;
    //done
	return 0;
}


int ZeroLengthImplexContact::revertToStart(void) {   

    // reset state variables
    sv = StateVariables();
    // done
	return 0;
}


int ZeroLengthImplexContact::update(void) {

    // get trial strain tensor at integration point in local crd. sys.
    computeElementStuff();                                         
    // update material internal variables
    computeMaterialStuff(true,true);  // explicit_phase?, do_tangent?
    //done
    return 0;
}


    // public methods to obtain stiffness, mass, damping and residual information
const Matrix & ZeroLengthImplexContact::getTangentStiff(void) {

    Matrix& stiff = *theStiffMatrix;
    Matrix& R = *theTransMatrix;
    Matrix temp(numDOF, numDOF);
    
    // get trial strain tensor at integration point in local crd. sys.
    computeElementStuff();
    // integrate over g. pts. for residual and tangent (only one pt. in this case)
    computeMaterialStuff(true, true);       // explicit_phase?, do_tangent?
    // fill in tangent stiffness
    int idx = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            stiff(i, j) = sv.C(i, j);
            stiff(i, j + idx) = -1 * sv.C(i, j);
            stiff(i + idx, j) = -1 * sv.C(i, j);
            stiff(i + idx, j + idx) = sv.C(i, j);
        }
    }
    // rotate local stiff matrix and resid vector back to global cord. sys.
    temp = stiff;
    stiff.addMatrixTripleProduct(0.0, R, temp, 1.0);
    // done
    return stiff;
}


const Matrix & ZeroLengthImplexContact::getInitialStiff(void) {
   
    Matrix& stiff = *theStiffMatrix;
    Matrix& R = *theTransMatrix;
    Matrix temp(numDOF, numDOF);

    // get trial strain tensor at integration point in local crd. sys.
    computeElementStuff();
    // integrate over g. pts. for residual and tangent (only one pt. in this case)
    computeMaterialStuff(true, false);      // explicit_phase?, do_tangent? 
    // fill in initial stiffness
    int idx = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            stiff(i, j) = sv.C(i, j);
            stiff(i, j + idx) = -1 * sv.C(i, j);
            stiff(i + idx, j) = -1 * sv.C(i, j);
            stiff(i + idx, j + idx) = sv.C(i, j);
        }
    }
    // rotate local stiff matrix and resid vector back to global cord. sys.
    temp = stiff;
    stiff.addMatrixTripleProduct(0.0, R, temp, 1.0);
    // done
    return stiff;
}
    

const Matrix & ZeroLengthImplexContact::getDamp(void) {

    Matrix& damp = *theTransMatrix;
 	damp.Zero();                                              // no mass
	return damp;
}


const Matrix & ZeroLengthImplexContact::getMass(void) {
    
    Matrix& mass = *theTransMatrix;
 	mass.Zero();                                              // no mass 
	return mass;
}


void  ZeroLengthImplexContact::zeroLoad(void) {

                                                                                        // does nothing
}


int ZeroLengthImplexContact::addLoad(ElementalLoad *theLoad, double loadFactor) {

  return 0;                                                                             // meaningless to addLoad to a contact!
}


int ZeroLengthImplexContact::addInertiaLoadToUnbalance(const Vector &accel) {
  
  return 0;                                                                             // does nothing as element has no mass!
}


const Vector & ZeroLengthImplexContact::getResistingForce(void) {

    Vector& resid = *theResidVector;
    Matrix& R = *theTransMatrix;
    Vector temp(numDOF);

    // get trial strain tensor at integration point in local crd. sys.
    computeElementStuff();
    // integrate over g. pts. for residual and tangent (only one pt. in this case)
    computeMaterialStuff(true, false);  // explicit_phase?, do_tangent?         
    // fill the residual vector
    int idx = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
        resid(i) = sv.sig(i);
        resid(i + idx) = -1 * sv.sig(i);
    }
    // rotate local stiff matrix and resid vector back to global cord. sys.
    temp = resid;
    resid.addMatrixTransposeVector(0.0, R, temp, 1.0);
    return resid;
}


const Vector & ZeroLengthImplexContact::getResistingForceIncInertia(void) { 

    Vector& resid = *theResidVector;
    Matrix& R = *theTransMatrix;
    Vector temp(numDOF);

    // get trial strain tensor at integration point in local crd. sys.
    computeElementStuff();
    // integrate over g. pts. for residual and tangent (only one pt. in this case)
    computeMaterialStuff(true, false);  // explicit_phase?, do_tangent?      
    // fill the residual vector
    int idx = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
        resid(i) = sv.sig(i);
        resid(i + idx) = -1 * sv.sig(i);
    }
    // rotate local stiff matrix and resid vector back to global cord. sys.
    temp = resid;
    resid.addMatrixTransposeVector(0.0, R, temp, 1.0);
    return  resid;
}


    // public methods for element output
int ZeroLengthImplexContact::sendSelf(int commitTag, Channel &theChannel) {

    int res = 0;
    int dataTag = this->getDbTag();
    static Vector data(11);
    data(0)  = this->getTag();
    data(1)  = Knormal;
    data(2)  = Kfriction;
    data(3)  = mu;
    data(4)  = cohesion;
    data(5)  = doImplEx;
    data(6)  = Xorient(0);
    data(7)  = Xorient(1);
    data(8)  = Xorient(2);
    data(9)  = numDIM;
    data(10) = numDOF;

    res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING ZeroLengthImplexContact::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return -1;
    }
    res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
    if (res < 0) {
        opserr << "WARNING ZeroLengthImplexContact::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return -1;
    }
	return 0;
}


int ZeroLengthImplexContact::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {

    int res;
    int dataTag = this->getDbTag();
    static Vector data(11);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING ZeroLengthImplexContact::recvSelf() - failed to receive Vector\n";
        return -1;
    }
    this->setTag((int)data(0));
    Knormal         = data(1);
    Kfriction        = data(2);
    mu              = data(3);
    cohesion        = data(4);
    doImplEx        = (bool)data(5);
    Xorient(0)      = data(6);
    Xorient(1)      = data(7);
    Xorient(2)      = data(8);
    numDIM          = (int)data(9);
    numDOF          = (int)data(10);

    res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
    if (res < 0) {
        opserr << "WARNING ZeroLengthImplexContact::recvSelf() - failed to receive ID\n";
        return -1;
    }
    return 0;
}


int ZeroLengthImplexContact::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode) {
    
    // ensure setDomain() worked
    if (theNodes[0] == 0 || theNodes[1] == 0)
        return 0;

    static Vector v1(3);
    static Vector v2(3);
    float d1 = 1.0;

    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);
    return theViewer.drawPoint(v1, d1, 10);
}


void ZeroLengthImplexContact::Print(OPS_Stream &strm, int flag) {

    if (flag == 0) {                                                                                    // print everything
        strm << "Element: " << this->getTag();
        strm << " type: ZeroLengthImplexContact  iNode: " << connectedExternalNodes(0);
        strm << " jNode: " << connectedExternalNodes(1) << endln;
    } 
    else if (flag == 1) {
        strm << this->getTag() << endln;
    } 
}


Response* ZeroLengthImplexContact::setResponse(const char **argv, int argc, Information &eleInformation) {

    Vector& resid = *theResidVector;
    Matrix& stiff = *theStiffMatrix;

    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
        return new ElementResponse(this, 1, resid);
    }                                               
    else if (strcmp(argv[0], "stiff") == 0 || strcmp(argv[0], "stiffness") == 0) { 
        return new ElementResponse(this, 2, stiff);                                                     // tangent stiffness matrix
    }
    else {
        return 0;
    }
}


int ZeroLengthImplexContact::getResponse(int responseID, Information &eleInfo)
{
    if (responseID == 1) {
        return eleInfo.setVector(this->getResistingForce());
    }
    else if (responseID == 2) {
        return eleInfo.setMatrix(this->getTangentStiff());
    }
    else {
        return -1;
    }
}


int ZeroLengthImplexContact::updateParameter(int parameterID, double value)
{
    if (parameterID == 1) {
        // set user defined current time increment
        // this is useful for rate dependency in implicit mode and for the implex model
        // when using arc length or displacement control methods, where the pseudo time step
        // is actually the load factor.
        // if when this variable is first set or when it is set before the first commit
        // we set the committed variable to the same value
        sv.dtime_n = value;
        if (!sv.dtime_first_set) {
            sv.dtime_n_commit = sv.dtime_n;
            sv.dtime_first_set = true;
        }
        sv.dtime_is_user_defined = true;
    }
    // done
    return 0;
}


// protected:


// private:
void ZeroLengthImplexContact::computeRotMatrix(void) {

    Matrix& rotatn = *theTransMatrix;

    // initialize orthogonal vectors to Xorient
    Vector rY(3);
    Vector rZ(3);
    // create global Y and Z axes
    static Vector gY = utils::make3DVector(0.0, 1.0, 0.0);   //global Y axis
    static Vector gZ = utils::make3DVector(0.0, 0.0, 1.0);   //global Z axis
    
    // compute two orthonal vectors to Xorient
    if (fabs(Xorient ^ gY) < 0.99) {                 // assume X is not parallel to global Y, and assume Xorient is a unit vector
        utils::cross(Xorient, gY, rZ);
        rZ.Normalize();
        utils::cross(rZ, Xorient, rY);
        rY.Normalize();
    }
    else {                                      // if not, assume X is not parallel to global Z
        utils::cross(Xorient, gZ, rY);
        rY.Normalize();
        utils::cross(rY, Xorient, rZ);
        rZ.Normalize();
    }
    // assemble the rotation matrix
    int lim_i = numNDS;         // limit "i" to number of nodes 
    int lim_j = numDOF / 2;     // limit "j" to number of d.o.f. at one node
    if (numDIM == 2 && numDOF == 4) {
        for (int i = 0; i < lim_i; i++) {
            for (int j = 0; j < lim_j; j++) {
                rotatn( (lim_j*i)   , (lim_j*i+j) ) = Xorient(j);
                rotatn( (lim_j*i+1) , (lim_j*i+j) ) = rY(j);
            }
        }

    }
    else if (numDIM == 2 && numDOF == 6) {
        for (int i = 0; i < lim_i; i++) {
            for (int j = 0; j < lim_j; j++) {
                rotatn( (lim_j*i)   , (lim_j*i+j) )    = Xorient(j);
                rotatn( (lim_j*i+1) , (lim_j*i+j) )    = rY(j);
                rotatn( (lim_j*i+2) , (lim_j*i+j) )    = rZ(j);
            }
        }

    }
    else if (numDIM == 3 && numDOF == 6) {
        for (int i = 0; i < lim_i; i++) {
            for (int j = 0; j < lim_j; j++) {
                rotatn( (lim_j*i)   , (lim_j*i+j) )    = Xorient(j);
                rotatn( (lim_j*i+1) , (lim_j*i+j) )    = rY(j);
                rotatn( (lim_j*i+2) , (lim_j*i+j) )    = rZ(j);
            }
        }
    }
    else if (numDIM == 3 && numDOF == 12) {
        for (int i = 0; i < lim_i; i++) {
            for (int j = 0; j < lim_j; j++) {
                rotatn( (lim_j*i)   , (lim_j*i+j) )   = Xorient(j);
                rotatn( (lim_j*i+1) , (lim_j*i+j) )   = rY(j);
                rotatn( (lim_j*i+2) , (lim_j*i+j) )   = rZ(j);
                rotatn( (lim_j*i+3) , (lim_j*i+j) )   = Xorient(j);
                rotatn( (lim_j*i+4) , (lim_j*i+j) )   = rY(j);
                rotatn( (lim_j*i+5) , (lim_j*i+j) )   = rZ(j);
            }
        }
    }
    else {
        opserr << "WARNING ZeroLengthImplexContact::computeRotMatrix cannot compute rotation matrix\n";
        rotatn.Zero();
    }
}


void ZeroLengthImplexContact::computeMaterialStuff(bool e_phase, bool t_flag) {

    /* ************************************************************************************ **
    **                                                                                      **
    **   Description: IMPL-EX/Implicit coulomb friction contact material                    **
    **                                                                                      **
    ** ************************************************************************************ */
    
    Vector ST(2);                   // trial tangent stress
    Vector dir(2);                  // slip direction
    double dlambda = 0.0;           // slip multiplier
    double slipfunction = 0.0;      // slip function

    // get committed values
    sv.sig = sv.sig_commit;
    sv.beta = sv.beta_commit;
    sv.lambda = sv.lambda_commit;

    // time factor for explicit extrapolation
    double time_factor = 1.0;
    if (e_phase && doImplEx && (sv.dtime_n_commit > 0.0))
        time_factor = sv.dtime_n / sv.dtime_n_commit;
    // note: the implex method just wants the ratio of the new to the old time step
    // not the real time step, so it is just fine to assume it to 1.
    // otherwise we have to deal with the problem of the opensees pseudo-time step
    // being the load multiplier in continuation methods...
    time_factor = 1.0;

    // compute strain increment
    sv.deps = sv.eps - sv.eps_commit;

    // Contact Model
        // update beta
    if (e_phase && doImplEx) {
            // explicit extrapolation
        //sv.beta = sv.beta_commit + time_factor * (sv.beta_commit - sv.beta_commit_old);
        sv.beta = sv.beta_commit;
    }
    else {
            // implicit
        sv.beta = sv.beta_commit + sv.deps(0);
    }
        // compute contact stiffness
            // flush the old tangent
    sv.C.Zero();
    if (sv.beta <= 0.0) {
            // fill-in initial stiffness
        sv.C(0, 0) = Knormal;
        sv.C(1, 1) = Kfriction;
        sv.C(2, 2) = Kfriction;
    }
        // update normal stress
    sv.sig(0) = sv.C(0, 0) * sv.eps(0);

    // Friction Model
        // compute trial friction stress
    ST(0) = sv.sig_commit(1) + Kfriction * sv.deps(1);
    ST(1) = sv.sig_commit(2) + Kfriction * sv.deps(2);
        // compute stresses
    if (e_phase && doImplEx) {
            // explicit extrapolation
        sv.lambda = sv.lambda_commit + time_factor * (sv.lambda_commit - sv.lambda_commit_old);
            // compute slip multiplier
        dlambda = std::max(0.0, (sv.lambda - sv.lambda_commit));
            // update stresses
        sv.sig(1) = ST(0) / (1 + Kfriction * dlambda);
        sv.sig(2) = ST(1) / (1 + Kfriction * dlambda);
    }
    else {
            // compute shear strength
        double tau_y = cohesion - sv.sig(0) * mu;               // m-c yield surface (sigma_n is negative)
            // on tau12 - tau13 plane...
        double normST = utils::norm2D(ST);                      // magnitude of friction stress
        slipfunction = normST - tau_y;                          // compute slip function
        if (slipfunction < 0.0) {                               // stick
            // compute slip multipler
            dlambda = 0.0;
            // update stresses
            sv.sig(1) = ST(0);
            sv.sig(2) = ST(1);
        }
        else {                                                  //slip
            // compute slip multipler
            sv.lambda = slipfunction / (Kfriction* tau_y);
            dlambda = sv.lambda - sv.lambda_commit;
            // update stresses
            dir = ST; dir.Normalize();               // slip direction                     
            sv.sig(1) = tau_y * dir(0);
            sv.sig(2) = tau_y * dir(1);
        }
        // update slip multiplier
        sv.lambda = sv.lambda_commit + dlambda;
    }

    // compute tangent modulus tensor
    if (t_flag) {
        if (doImplEx) {
            // compute tangent modulus tensor
            sv.C(1, 1) = Kfriction / (1 + Kfriction * dlambda);
            sv.C(2, 2) = Kfriction / (1 + Kfriction * dlambda);
        }
        else {
            if (sv.beta <= 0.0) {
                if (slipfunction < 0.0) {                                           //stick
                    sv.C(1, 1) = Kfriction;
                    sv.C(2, 2) = Kfriction;
                }
                else {                                                              // slip
                    sv.C(1, 0) = mu * Knormal * dir(0) * utils::sign(sv.sig(0));
                    sv.C(2, 0) = mu * Knormal * dir(1) * utils::sign(sv.sig(0));
                }
            }
        }
    }
}


void ZeroLengthImplexContact::computeElementStuff(void) {

    /* ************************************************************************************ **
    **                                                                                      **
    **  Description: Implementation of a two-node zero-length element with one gauss point  **
    **                                                                                      **
    ** ************************************************************************************ */
    
    Vector& uloc = *theDispVector;
    Matrix& R = *theTransMatrix;
    Matrix theB(numDIM,numDOF);
    Vector u_trial(numDIM);
    Vector temp(numDOF);

    // step 1: compute global to local rotation matrix
    computeRotMatrix();

    // step 2: get global displacements, put them into a vector
    const Vector& UgM = theNodes[0]->getTrialDisp();                                // master node disp. 
    const Vector& UgS = theNodes[1]->getTrialDisp();                                // slave node disp.
    int slave = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
        uloc(i)         = UgM(i);
        uloc(i+slave)   = UgS(i);
    }

    // step 3: rotate displacement vector from global to local coord. sys.
    temp = uloc;
    uloc.addMatrixVector(0.0, R, temp, 1.0);  // vector_local = rotation_matrix * vector_global

    // step 4: interp. element disp. for each integration point using the B matrix (only one pt. in this case)
        // create the B-matrix (N1 and N2 are linear)
    for (int i = 0; i < numDIM; i++) {
        theB(i,i)          = 1;                    // derivative of N1
        theB(i,i+slave)    = -1;                     // derivative of N2
    }
    u_trial.addMatrixVector(0.0, theB, uloc, 1.0);  // displacement at the gauss point < 0
        // update current displacement jump
    for (int i = 0; i < numDIM; i++) {
        sv.eps(i) = u_trial(i);                     // keep displacement 3D regardless
    }
    sv.eps(0) += gap;                               // add initial gap
    // element stuff is done
}