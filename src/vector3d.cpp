#include <iostream>
#include <math.h>

#include "vector3d.h"
#include "constants.h"


Vector3d::Vector3d() {

    m_ObjectId = globalObjectId++; 

    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}

// Regular constructor - with parameters for x, y, and z
Vector3d::Vector3d(const double p_x, const double p_y, const double p_z) {

    // Initialize member variables
    
    m_ObjectId = globalObjectId++; 

    m_x = p_x;
    m_y = p_y;
    m_z = p_z;
}

// Regular constructor - initialise from std:vector
Vector3d::Vector3d(const DBL_VECTOR p_Vec) {

    m_ObjectId = globalObjectId++; 
    
    size_t numValuesSupplied = p_Vec.size();

    THROW_ERROR_IF(numValuesSupplied != 3, ERROR::EXPECTED_3D_VECTOR);  // this is a coding error

    m_x = p_Vec[0];
    m_y = p_Vec[1];
    m_z = p_Vec[2];
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
 * Vector3d::ChangeBasis(const double p_ThetaE, const double p_PhiE, const double p_PsiE)
 * 
 * @param   [IN]   p_ThetaE                    Euler angle Theta (rad) 
 * @param   [IN]   p_PhiE                      Euler angle Phi   (rad) 
 * @param   [IN]   p_PsiE                      Euler angle Psi   (rad) 
 * @return                                     Vector in previous basis
 */
Vector3d Vector3d::ChangeBasis(const double p_ThetaE, const double p_PhiE, const double p_PsiE) {
// For convenience, undefined below
#define cTheta cos(p_ThetaE)
#define cPhi   cos(p_PhiE)
#define cPsi   cos(p_PsiE)
#define sTheta sin(p_ThetaE)
#define sPhi   sin(p_PhiE)
#define sPsi   sin(p_PsiE)

    // Define the Rotation Matrix     
    std::vector<DBL_VECTOR> rotationMatrix = {
        { cPhi * cPsi - sPhi * cTheta * sPsi ,  -cPhi * sPsi - sPhi * cTheta * cPsi ,  sTheta * sPhi },
        { sPhi * cPsi + cPhi * cTheta * sPsi ,  -sPhi * sPsi + cPhi * cTheta * cPsi , -sTheta * cPhi },
        { sTheta * sPsi                      ,  sTheta * cPsi                       ,  cTheta        }
    };

    Vector3d oldVector = *this;             // input
    Vector3d newVector = Vector3d(0,0,0);   // output

    // Apply rotation
    for (size_t row = 0; row < 3; row++) {
        for (size_t col = 0; col < 3; col++) {
            newVector[row] += oldVector[col] * rotationMatrix[row][col];
        }
    }

    return newVector;

#undef cTheta
#undef cPhi
#undef cPsi
#undef sTheta
#undef sPhi
#undef sPsi
}


/*
 * Right multiply a matrix by a vector
 *
 * Vector3d MatrixMult(const std::vector<DBL_VECTOR>& p_matrix, const Vector3d& p_vector) 
 *
 * @param   [IN]   p_Matrix                     matrix 
 * @param   [IN]   p_Vec                        vector
 * @return                                      angle between them, in radians
 */
Vector3d Vector3d::MatrixMult(const std::vector<DBL_VECTOR>& p_Matrix, const Vector3d& p_Vec) {

    Vector3d result = Vector3d(0.0, 0.0, 0.0);

    size_t numRowsSupplied = p_Matrix.size();
    THROW_ERROR_IF(numRowsSupplied != 3, ERROR::EXPECTED_3D_VECTOR);        // this is a code defect
    for (size_t row = 0; row < numRowsSupplied; row++) {
        size_t numColsSupplied = p_Matrix[row].size();
        THROW_ERROR_IF(numColsSupplied != 3, ERROR::EXPECTED_3D_VECTOR);    // this is a code defect
        for (size_t col = 0; col < numColsSupplied; col++) {
            result[row] += p_Matrix[row][col] * p_Vec[col];
        }
    }

    return result;
}


/*
 * Rotate a vector about the X axis.
 *
 * Vector3d RotateVectorAboutX(const double p_Theta)
 *
 * @param   [IN]   p_Theta                     Rotation angle (rad) 
 * @return                                     Vector after rotation
 */
Vector3d Vector3d::RotateVectorAboutX(const double p_Theta) {
// For convenience, undefined below
#define cTheta cos(p_Theta)
#define sTheta sin(p_Theta)

    // Define the Rotation Matrix
    std::vector<DBL_VECTOR> RotationMatrix = {     
        { 1.0,  0.0,     0.0    },
        { 0.0,  cTheta, -sTheta },
        { 0.0,  sTheta,  cTheta }};

    // Apply rotation and return result
    return MatrixMult(RotationMatrix, *this);

#undef cTheta
#undef sTheta
}


/*
 * Rotate a vector about the Y axis.
 *
 * Vector3d RotateVectorAboutY( const double p_Theta)
 *
 * @param   [IN]   p_Theta                     Rotation angle (rad) 
 * @return                                     Vector after rotation
 */
Vector3d Vector3d::RotateVectorAboutY( const double p_Theta) {
// For convenience, undefined below
#define cTheta cos(p_Theta)
#define sTheta sin(p_Theta)

    // Define the Rotation Matrix
    std::vector<DBL_VECTOR> RotationMatrix = {
        {  cTheta, 0.0, sTheta },
        {  0.0,    1.0, 0.0    },
        { -sTheta, 0.0, cTheta }};

    // Apply rotation and return result
    return MatrixMult(RotationMatrix, *this);

#undef cTheta
#undef sTheta
}


/*
 * Rotate a vector about the Z axis.
 *
 * Vector3d RotateVectorAboutZ( const double p_Theta)
 *
 * @param   [IN]   p_Theta                     Rotation angle (rad) 
 * @return                                     Vector after rotation
 */
Vector3d Vector3d::RotateVectorAboutZ( const double p_Theta) {
// For convenience, undefined below
#define cTheta cos(p_Theta)
#define sTheta sin(p_Theta)

    // Define the Rotation Matrix
    std::vector<DBL_VECTOR> RotationMatrix = {
        { cTheta, -sTheta,  0.0 },
        { sTheta,  cTheta,  0.0 },
        { 0.0,     0.0,     1.0 }};

    // Apply rotation and return result
    return MatrixMult(RotationMatrix, *this);

#undef cTheta
#undef sTheta
}
