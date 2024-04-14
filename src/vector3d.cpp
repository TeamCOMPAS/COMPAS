#include <iostream>
#include <math.h>
#include "vector3d.h"
#include "constants.h"


// Default constructor
Vector3d::Vector3d() {

    m_ObjectId = globalObjectId++; 

    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}

// Regular constructor - with parameters for x, y, and z
Vector3d::Vector3d(const double p_x, const double p_y, const double p_z) {

    m_ObjectId = globalObjectId++; 

    m_x = p_x;
    m_y = p_y;
    m_z = p_z;
}

// Regular constructor - initialise from std:vector
Vector3d::Vector3d(const DBL_VECTOR p_Vec) {

    m_ObjectId = globalObjectId++; 
    
    int numValuesSupplied = p_Vec.size();

    SHOW_WARN_IF(numValuesSupplied != 3, ERROR::EXPECTED_3D_VECTOR);

    m_x = numValuesSupplied >= 1 ? p_Vec[0] : 0.0;
    m_y = numValuesSupplied >= 2 ? p_Vec[1] : 0.0;
    m_z = numValuesSupplied >= 3 ? p_Vec[2] : 0.0;
}


/*
 * Redefine a vector from one coordinate basis to another using Euler Angles
 * 
 * For a vector defined in a new coordinate basis, (X',Y',Z'), we want to 
 * find it's values in a previous coordinate basis (X,Y,Z) so that we can 
 * add the new vector to vectors defined in the previous basis.
 *
 * For a change of basis from (X,Y,Z)->(X',Y',Z')
 * 
 *     ThetaE [0, pi] is the angle between Z and Z',
 *     Vector N := Z x Z' (cross product)
 *     PhiE [0, 2pi) is the angle between X and N
 *     PsiE [0, 2pi) is the angle between X' and N
 * 
 * These angles uniquely determine the change of basis, which is applied 
 * in the form of a rotation matrix R as a function of these angles.
 * 
 * V_{X,Y,Z} = R * V_{X',Y',Z'} gives the vector in the original coordinates.
 *
 * For details, see:
 * https://en.wikipedia.org/wiki/Euler_angles
 * https://en.wikipedia.org/wiki/Change_of_basis
 *
 * 
 * Vector3d::RotateVector(const double p_ThetaE, const double p_PhiE, const double p_PsiE)
 * 
 * @param   [IN]   p_ThetaE                    Euler angle Theta (rad) 
 * @param   [IN]   p_PhiE                      Euler angle Phi   (rad) 
 * @param   [IN]   p_PsiE                      Euler angle Psi   (rad) 
 * @return                                     Vector in previous basis
 */
Vector3d Vector3d::RotateVector(const double p_ThetaE, const double p_PhiE, const double p_PsiE) {
// Replace for convenience, undefined below
#define cTheta cos(p_ThetaE)
#define cPhi   cos(p_PhiE)
#define cPsi   cos(p_PsiE)
#define sTheta sin(p_ThetaE)
#define sPhi   sin(p_PhiE)
#define sPsi   sin(p_PsiE)

    Vector3d result = *this;        // default return is this vector

    // Define the Rotation Matrix     
    std::vector<DBL_VECTOR> rotationMatrix = {
        { cPhi * cPsi - sPhi * cTheta * sPsi ,  -cPhi * sPsi - sPhi * cTheta * cPsi ,  sTheta * sPhi },
        { sPhi * cPsi + cPhi * cTheta * sPsi ,  -sPhi * sPsi + cPhi * cTheta * cPsi , -sTheta * cPhi },
        { sTheta * sPsi                      ,  sTheta * cPsi                       ,  cTheta        }
    };

    // Apply rotation
    for (size_t row = 0; row < 3; row++) {
        for (size_t col = 0; col < 3; col++) {
            result[row] += result[col] * rotationMatrix[row][col];
        }
    }

    return result;

#undef cTheta
#undef cPhi
#undef cPsi
#undef sTheta
#undef sPhi
#undef sPsi
}
