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
#include <map>

/*
TODO:
- add initial disp
*/

namespace 
{
     class GlobalStorage
     {
     public:
         int size = 0;
         Matrix K; // stiffness
         Matrix K0; // initial stiffness
         Matrix M; // mass 
         Matrix D; // damping
         Vector U; // displacement
         Vector R; // residual

     public:
         GlobalStorage() = default;
         GlobalStorage& resize(int N) {
             if (N != size) {
                 K.resize(N, N);
                 K0.resize(N, N);
                 M.resize(N, N);
                 D.resize(N, N);
                 U.resize(N);
                 R.resize(N);
             }
             return *this;
         }
     };

     static GlobalStorage& getGlobalStorage(int N)
     {
         static std::map<int, GlobalStorage> gsmap;
         return gsmap[N].resize(N);
     }
}

// Remaining Steps To Do:
/* List:

9) write function sendSelf

10) write function recvSelf

11) fix impl-ex

12) ask and correct the signs of residual and stiffness

*/

// class utilities
namespace utils {

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

}

void* OPS_ZeroLengthImplexContact(void) {

    double SmallNumber = 1.0e-6;
    Element* theElement = nullptr;

    // some kudos
    static int counter = 0;
    if (++counter == 1)
        opserr << "ZeroLengthImplexContact element - Written: O.D.Akan (IUSS), M.Petracca (ASDEA), G.Camata (UNICH)\n";

    // model dimension
    int ndm = OPS_GetNDM();
    if (ndm < 2 || ndm > 3) {
        opserr << "ZeroLengthImplexContact: Unsupported NDM (" << ndm << "). It should be 2 or 3\n";
        return theElement;
    }

    // a quick check on number of args
    if (OPS_GetNumRemainingInputArgs() < 7) {
        opserr << "ZeroLengthImplexContact: WARNING: too few arguments \n" <<
            "want - element zeroLengthImplexContact eleTag? iNode? jNode? Kn? Kt? mu? C? <-orient $x1 $x2 $x3> <-intType type?>\n";
        return theElement;
    }

    // start with mandatory inputs
    // read eleTag, iNode, jNode
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "ZeroLengthImplexContact: WARNING: invalid int inputs\n";
        return theElement;
    }

    // read Kn, Kt, mu, C
    double ddata[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
        opserr << "ZeroLengthImplexContact: WARNING: invalid double inputs\n";
        return theElement;
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
                    return theElement;
                }
                numdata = 3;
                if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
                    opserr << "ZeroLengthImplexContact: WARNING: invalid double input after -orient\n";
                    return theElement;
                }
            }
            else if (ndm == 3) {
                if (OPS_GetNumRemainingInputArgs() < 3) {
                    opserr << "ZeroLengthImplexContact: WARNING: insufficient orient values in 3D\n";
                    return theElement;
                }
                numdata = 3;
                if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
                    opserr << "ZeroLengthImplexContact: WARNING: invalid double input after -orient\n";
                    return theElement;
                }
            }
            else {
                opserr << "ZeroLengthImplexContact: WARNING: -orient: model dimension is invalid! \n";
                return theElement;
            }
        }
        else if (strcmp(inputstring, "-intType") == 0) {                             // #2 read type of integration 
            numdata = 1;
            if (OPS_GetIntInput(&numdata, &integrationType) < 0) {
                opserr << "ZeroLengthImplexContact: WARNING: invalid integer after -intType\n";
                return theElement;
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
        return theElement;
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

ZeroLengthImplexContact::ZeroLengthImplexContact(int tag, int Nd1, int Nd2,
    double Kn, double Kt, double fcoeff, double C,
    int ndm, bool itype, double xN, double yN, double zN)
    : Element(tag, ELE_TAG_ZeroLengthImplexContact)
    , connectedExternalNodes(2)
    , Knormal(Kn)
    , Kfriction(Kt)
    , mu(fcoeff)
    , cohesion(C)
    , numDIM(ndm)
    , doImplEx(itype)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    Xorient[0] = xN;
    Xorient[1] = yN;
    Xorient[2] = zN;
}

ZeroLengthImplexContact::ZeroLengthImplexContact(void)
    : Element(0, ELE_TAG_ZeroLengthImplexContact)
{
}

ZeroLengthImplexContact::~ZeroLengthImplexContact(void)
{
}

int ZeroLengthImplexContact::getNumExternalNodes(void) const
{
    return 2;
}

const ID& ZeroLengthImplexContact::getExternalNodes(void)
{
    return connectedExternalNodes;
}

Node** ZeroLengthImplexContact::getNodePtrs(void)
{
    return theNodes.data();
}

int ZeroLengthImplexContact::getNumDOF(void)
{
    return numDOF[0] + numDOF[1];
}

void ZeroLengthImplexContact::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = nullptr;
        theNodes[1] = nullptr;
        return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    // if can't find both - send a warning message
    if ((theNodes[0] == nullptr) || (theNodes[1] == nullptr)) {
        opserr << 
            "FATAL ERROR ZeroLengthImplexContact::setDomain() - Nd1: " << Nd1 <<
            " and/or Nd2: " << Nd2 << " do not exist in the model.\n";
        exit(-1);
    }

    // now determine the number of dof and the dimension
    numDOF[0] = theNodes[0]->getNumberDOF();
    numDOF[1] = theNodes[1]->getNumberDOF();

    // check dofs
    if (numDIM == 2) {
        for (int i = 0; i < 2; ++i) {
            int idof = numDOF[i];
            if (idof < 2 || idof > 3) {
                opserr <<
                    "FATAL ERROR ZeroLengthImplexContact::setDomain() - #DOFs (" 
                    << idof << ") at node " << i + 1 
                    << " is not supported! it can be either 2 or 3\n";
                exit(-1);
            }
        }
    }
    else {
        for (int i = 0; i < 2; ++i) {
            int idof = numDOF[i];
            if ((idof != 3) && (idof != 4) && (idof != 6)) {
                opserr <<
                    "FATAL ERROR ZeroLengthImplexContact::setDomain() - #DOFs ("
                    << idof << ") at node " << i + 1
                    << " is not supported! it can be either 3, 4 or 6\n";
                exit(-1);
            }
        }
    }

    // call the base class method
    DomainComponent::setDomain(theDomain);
}

int ZeroLengthImplexContact::commitState(void) 
{
    // do the implicit correction if impl-ex
    if (doImplEx) {
        // get trial strain tensor at integration point in local crd. sys.
        computeStrain();
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

int ZeroLengthImplexContact::revertToLastCommit(void) 
{
    // restore converged values [step (n+1) --> (n)]
    sv.eps = sv.eps_commit;
    sv.sig = sv.sig_commit;
    sv.beta = sv.beta_commit;
    sv.lambda = sv.lambda_commit;
    sv.dtime_n = sv.dtime_n_commit;
    //done
	return 0;
}

int ZeroLengthImplexContact::revertToStart(void) 
{   
    // reset state variables
    sv = StateVariables();
    // done
	return 0;
}

int ZeroLengthImplexContact::update(void)
{
    computeStrain();
    updateInternal(true, true);
    return 0;
}

const Matrix & ZeroLengthImplexContact::getTangentStiff(void)
{
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    auto& K = gs.K;
    const auto& C = sv.C;
    computeGenericStiffness(C, K);
    return K;
}

const Matrix & ZeroLengthImplexContact::getInitialStiff(void) 
{
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    auto& K = gs.K0;

    static Matrix C0(3, 3);
    C0.Zero();
    const Vector& LGap = computeLocalInitialGap();
    double Un = LGap(0);
    if (Un < 0.0) {
        C0(0, 0) = Knormal;
        C0(1, 1) = C0(2, 2) = Kfriction;
    }
    
    computeGenericStiffness(C0, K);
    return K;
}

const Matrix & ZeroLengthImplexContact::getDamp(void)
{
    // get global storage for gloabl DOFset
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    gs.D.Zero();
    return gs.D;
}

const Matrix & ZeroLengthImplexContact::getMass(void)
{
    // get global storage for global DOFset
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    gs.M.Zero();
    return gs.M;
}

const Vector & ZeroLengthImplexContact::getResistingForce(void) {

    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);

    // get trial strain tensor at integration point in local crd. sys.
    computeStrain();
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

const Vector & ZeroLengthImplexContact::getResistingForceIncInertia(void)
{
    return getResistingForce();
}

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
    data(10) = numDOF[0];

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
    numDOF[0]       = (int)data(10);

    res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
    if (res < 0) {
        opserr << "WARNING ZeroLengthImplexContact::recvSelf() - failed to receive ID\n";
        return -1;
    }
    return 0;
}

int ZeroLengthImplexContact::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
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

Response* ZeroLengthImplexContact::setResponse(const char **argv, int argc, Information &eleInformation)
{
    // TODO: record contact forces (material stress)
    // : displacement jump, 
}

int ZeroLengthImplexContact::getResponse(int responseID, Information &eleInfo)
{
    // TODO: record contact forces (material stress)
    // : displacement jump,
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

const Matrix& ZeroLengthImplexContact::computeRotMatrix1()
{
    static Matrix T(3, 3);

    // initialize orthogonal vectors to Xorient
    static Vector rY(3);
    static Vector rZ(3);
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

    for (int j = 0; j < 3; j++) { // for each DOF (local DOFset)
        T(0, j) = Xorient(j);
        T(1, j) = rY(j);
        T(2, j) = rZ(j);
    }

    return T;
}

const Matrix& ZeroLengthImplexContact::computeRotMatrix2(void)
{
    static Matrix T2(6, 6);
    T2.Zero();

    const Matrix& T1 = computeRotMatrix1();

    for (int q = 0; q < 2; ++q) {
        int index = q * 3;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                T2(index + i, index + j) = T1(i, j);
    }

    return T2;
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
        double normST = ST.Norm();                     // magnitude of friction stress
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

/*
void ZeroLengthImplexContact::computeElementStuff(void) {

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
*/

const Vector& ZeroLengthImplexContact::computeLocalInitialGap()
{
    const Vector& PM = theNodes[0]->getCrds();
    const Vector& PS = theNodes[1]->getCrds();
    static Vector GGap(3);
    GGap.Zero();
    for (int i = 0; i < numDIM; ++i)
        GGap(i) = PS(i) - PM(i);
    // make it local
    static Vector LGap(3);
    const Matrix& T1 = computeRotMatrix1();
    LGap.addMatrixVector(0.0, T1, GGap, 1.0);
    return LGap;
}

void ZeroLengthImplexContact::computeStrain()
{
    // get global displacement for global DOFset
    const Vector& UgM = theNodes[0]->getTrialDisp();
    const Vector& UgS = theNodes[1]->getTrialDisp();

    // set global displacement for local DOFset 
    static Vector UGL(6);
    for (int i = 0; i < numDIM; i++) {
        UGL(i) = UgM(i);
        UGL(i + 3) = UgS(i);
    }

    // compute local displacement for local DOFset
    static Vector UL(6);
    const Matrix& T2 = computeRotMatrix2();
    UL.addMatrixVector(0.0, T2, UGL, 1.0);

    // compute diplacement jump US(2) - UM(1) (will be called strain) in local coordinates
    static Matrix B(3, 6);
    for (int i = 0; i < 3; i++) {
        B(i, i + 3) = 1;
        B(i, i) = -1;
    }
    sv.eps.addMatrixVector(0.0, B, UL, 1.0);  // displacement jump

    // add the initial gap (if any)
    const Vector& LGap = computeLocalInitialGap();
    sv.eps.addVector(1.0, LGap, 1.0);
}

void ZeroLengthImplexContact::updateInternal(bool do_implex, bool do_tangent)
{

}

void ZeroLengthImplexContact::computeGenericStiffness(const Matrix& C, Matrix& K)
{

}
