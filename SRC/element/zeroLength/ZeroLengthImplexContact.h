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

// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthImplexContact.h,v $
// $Revision: 1.0 $
// $Date: 2020-XX-XX XX:XX:XX $

#ifndef ZeroLengthImplexContact_h
#define ZeroLengthImplexContact_h

// Written: Onur Deniz Akan         (onur.akan@iusspavia.it)
//          Dr. Massimo Petracca    (m.petracca@asdea.net)
//          Prof. Guido Camata      (g.camata@unich.it)
//
// Created: May 2020
//

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                         ZeroLengthImplexContact element                        |
 +                                                                                +
 |--------------------------------------------------------------------------------|
 |                                                                                |
 +         Authors: Onur Deniz Akan (IUSS), Dr. Massimo Petracca (ASDEA)          +
 |                  and  Prof. Guido Camata (UNICH)                               |
 |                                                                                |
 +      Istituto Universitario di Studi Superiori di Pavia          (IUSS)        +
 |      Advanced Structural Design & Analysis Software Technology   (ASDEA)       |
 |		Universita degli Studi 'G. d'Annunzio' Chieti - Pescara	    (UNICH)       |
 +			                                                                      +
 |                                                                                |
 |                    Email: onur.akan@iusspavia.it (O.D.A.)                      |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- May 2020                                                     |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

 /* command

  element zeroLengthImplexContact $eleID $sNdID $mNdID $Kn $Kt $mu $c <-orient $x1 $x2 $x3> <-intType $type>

  where:
   $eleID   : element ID of this contact element
   $sNdID   : slave node ID
   $mNdID   : master node ID
   $Kn      : penalty in normal directions
   $Kt      : penalty in tangential directions
   $mu      : friction coefficient [mu = tan(Phi) where Phi is the angle of internal friction]
   $c       : cohesion
   -orient  : vector components in global coordinates defining local x-axis in 2D or 3D (optional)
                default = <1,0> in 2D and <1,0,0> in 3D
   -intType : type of integration scheme (optional)
                rFlag = 0 implicit, backward-euler integration scheme (default)
                rFlag = 1 impl-ex, implicit/explicit integration scheme


  Description: This file contains the class definition for ZeroLengthImplexContact.
  (1) A ZeroLengthImplexContact element is defined by two nodes in R^2 (x,y) or R^3 (x,y,z).
  (2) Penalty (Kn,Kt) method is used to constrain displacements in the normal direction,
      and in the tangential direction along with a Mohr-Coulomb frictional yield surface.
  (3) Optionally, the normal vector of the contact plane, orienting the contact axis can be
      specified by the user. Otherwise, the global X vector is chosen by default.
  (4) Optionally, the user can prefer using the IMPL-EX scheme over the standard Backward Euler
      return mapping scheme for the implicit-explicit integration of the contact material residual
      and tangent modulus.
  (5) Currently, only small deformations are considered, however an update on large deformations
      may be implemented in the near future.
  (6) In IMPL-EX scheme, the material tangent and residual at the current step (n+1) are computed
      by lineary extrapolating the material parameters commited at steps (n) and (n-1). Then, after
      the convergence criteria is achieved, extrapolated material parameters are corrected with a 
      step of traditional implicit computation from the last commited step (n). Finally, corrected 
      parameters are saved as the current step (n+1) commited variables.
  (7) Along with the Newton-Raphson method, IMPL-EX results in a step-linear solution of the material 
      non-linearity, and a symmetric, semi-positive definite tangent tensor regardless of the analytical
      tangent, hence improving the robustness of the solution. The order of accuracy is the same with
      the implicit (backward) euler integration scheme, however the commited error at a time step is larger.
      An automatic time-stepping algorithm is employed to a priori control this error. (Oliver et al, 2008)

   References:
   Armero, F. and Petocz, E. "A New Dissipative Time-Stepping Algorithm for Frictional
        Contact Problems: Formulation and Analysis", Comp. Methods Appl. Mech. Eng., 179,
        151-178 (1999)

   Laursen,T.A. and Simo, J.C. "A Continuum-based Finite Element Formulation for the
        Implicit Solution of Multibody, Large Deformation Frictional Contact Problems".
        Int. J. Numer. Methods Eng. , 36, 3451-3485 (1993)

   Oliver, J., Huespe, A.E. and Cante, J.C. "An Implicit/Explicit Integration Scheme to
        Increase Computability of Non-linear Material and Contact/Friction Problems", Comp.
        Methods Appl. Mech. Eng., 197, 1865-1889 (2008)
 */

#include <Element.h>
#include <Matrix.h>
#include <array>

#define	LENTOL 1.0e-6                                   // tolerance for zero-length of the element

class Node;
class Channel;
class Response;

class ZeroLengthImplexContact : public Element {

    public:
        class StateVariables {
            public:
                // material strain and stress [(0 = normal) (1 = tangent_1) (2 = tangent_2)]
                Vector eps              = Vector(3);        // material strain 
                Vector eps_commit       = Vector(3);        // commited material strain
                Vector deps             = Vector(3);        // strain increment
                Vector sig              = Vector(3);        // contact stress
                Vector sig_commit       = Vector(3);        // commited contact stress
                double beta             = 0.0;              // plastic multiplier (contact)
                double beta_commit      = 0.0;              // commited plastic multiplier
                double lambda           = 0.0;              // slip multiplier
                double lambda_commit    = 0.0;              // commited slip multiplier
                // state variables for implex
                double beta_commit_old  = 0.0;              // prev. commited plastic multiplier
                double lambda_commit_old = 0.0;             // prev. commited slip multiplier
                double dtime_n = 0.0;                       // time factor
                double dtime_n_commit = 0.0;                // commited time factor
                bool dtime_is_user_defined = false;
                bool dtime_first_set = false;
                // moduli
                Matrix C = Matrix(3, 3);                  // tangent modulus matrix
                // constructor
                StateVariables(void) = default;
                StateVariables(const StateVariables&) = default;
                StateVariables& operator = (const StateVariables&) = default;
        };

    public:
        // constructor
        ZeroLengthImplexContact(int tag, int Nd1, int Nd2,
            double Kn, double Kt, double fs, double c,
            int ndm, bool itype, double xN, double yN, double zN);

        // null constructor & destructor
        ZeroLengthImplexContact(void);
        ~ZeroLengthImplexContact(void);

        // public methods to obtain information about dof & connectivity
        int getNumExternalNodes(void) const;
        const ID& getExternalNodes(void);
        Node** getNodePtrs(void);
        int getNumDOF(void);
        void setDomain(Domain* theDomain);

        // public methods to set the state of the element
        int commitState(void);
        int revertToLastCommit(void);
        int revertToStart(void);
        int update(void);

        // public methods to obtain stiffness, mass, damping and residual information
        const Matrix& getTangentStiff(void);
        const Matrix& getInitialStiff(void);
        const Matrix& getDamp(void);
        const Matrix& getMass(void);

        const Vector& getResistingForce(void);
        const Vector& getResistingForceIncInertia(void);

        // public methods for element output
        int sendSelf(int commitTag, Channel& theChannel);
        int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
        int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
        void Print(OPS_Stream& s, int flag = 0);

        Response* setResponse(const char** argv, int argc, Information& eleInformation);
        int getResponse(int responseID, Information& eleInformation);
        int updateParameter(int parameterID, double value);

    private:
        // compute rotation matrix
        const Matrix& computeRotMatrix();                                    

        // residual and tangent computation
        void computeMaterialStuff(bool e_phase, bool t_flag);           // compute contact material response
                                                                                 
        // compute strain
        void computeStrain();
        // compute material response
        void updateInternal(bool do_implex, bool do_tangent);

    private:
        // element info
        ID connectedExternalNodes = ID(2); // contains the tags of the end nodes
        double Knormal = 0.0; // normal penalty
        double Kfriction = 0.0; // tangential penalty
        double mu = 0.0; // friction coefficient
        double cohesion = 0.0; // cohension
        int numDIM = 0; // model dimension
        std::array<int, 2> numDOF = { { 0, 0 } }; // no. of DOF at 1st and 2nd node
        bool doImplEx = false; // integration type flag
        Vector Xorient = Vector(3); // contact axis orientation in global coords.
        std::array<std::array<double, 3>, 2> U0 = { {{0.0,0.0,0.0}, {0.0,0.0,0.0}} }; // initial displacement at node 1 and 2

        // element pointers
        std::array<Node*, 2> theNodes = { { nullptr, nullptr } }; // node pointers

        // state variables
        StateVariables sv; // element & material state variables

};
#endif