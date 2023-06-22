#include <iostream>
#include <math.h>
#include "vector3d.h"
#include "constants.h"

///////////////////////////////////////////////////////////
// 
// Initialize Vector3d object
// 
///////////////////////////////////////////////////////////

// Default constructor
Vector3d::Vector3d() {

    // Initialize member variables

    m_ObjectId = globalObjectId++; 

    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;

}

// Initialize from values 

Vector3d::Vector3d(double p_x, double p_y, double p_z) {

    m_ObjectId = globalObjectId++; 

    m_x = p_x;
    m_y = p_y;
    m_z = p_z;
}

Vector3d::Vector3d(DBL_VECTOR p_v) {

    m_ObjectId = globalObjectId++; 
    
    m_x = p_v[0];
    m_y = p_v[1];
    m_z = p_v[2];
}

// Update after initialization
void Vector3d::UpdateVector(const double p_x, const double p_y, const double p_z) {
    m_x = p_x;
    m_y = p_y;
    m_z = p_z;
}

///////////////////////////////
// Overload operators

// Indexing Operator
double& Vector3d::operator[] (size_t p_i) {
    switch (p_i) {
        case 0: return m_x;
        case 1: return m_y;
        case 2: return m_z;
        default: throw "Not an index!\n";       // JR: we should fix this one day - the error should be handled rather than just stopping the program with an exception
    }
}

// Overload setting operator
void Vector3d::operator =(Vector3d p_SourceVec) {
    UpdateVector(p_SourceVec[0], p_SourceVec[1], p_SourceVec[2]);
}

// Overload vector addition operator
Vector3d Vector3d::operator +(Vector3d p_Vec) {
    Vector3d newVect;

    for (size_t i=0; i<3; i++) {
        newVect[i] = (*this)[i] + p_Vec[i];
    }
    return newVect;
}

// Overload vector subtraction operator
Vector3d Vector3d::operator -(Vector3d p_Vec) {
    Vector3d newVect;

    for (size_t i=0; i<3; i++) {
        newVect[i] = (*this)[i] - p_Vec[i];
    }
    return newVect;
}

// Overload scalar multiplication operator...
Vector3d Vector3d::operator *(const double p_Scalar) {
    Vector3d newVect;

    for (size_t i=0; i<3; i++) { 
        newVect[i] = (*this)[i] * p_Scalar; 
    }
    return newVect;
}
// ...and in reverse order (need the reverser to be a free function)
Vector3d operator *(double p_Scalar, Vector3d p_Vec) { return p_Vec * p_Scalar; }

// Overload scalar division operator
Vector3d Vector3d::operator /(const double p_Scalar) {
    Vector3d newVect;

    for (size_t i=0; i<3; i++) { 
        newVect[i] = (*this)[i] / p_Scalar; 
    }
    return newVect;
}


// Overload output << operator for Vector3d
std::ostream &operator <<(std::ostream &os, Vector3d const p_Vec) {
    return os << "{" << p_Vec[0] << ", " << p_Vec[1] << ", " << p_Vec[2] << "}";
}


//////////////////////////////////
// Add in common vector calculations

/*
 * Calculate the magnitude of a velocity vector, the speed.
 * 
 *
 * @return                                       The magnitude of the velocity vector (speed)
 */
double Vector3d::Magnitude() const {
    // Straightforward application of pythagorean theorem 
    return std::sqrt(linalg::dot((*this), (*this)));
}


/*
 * Convert the Vector3d to a DBL_VECTOR
 * 
 *
 * @return                                       The analogous DBL_VECTOR
 */
DBL_VECTOR Vector3d::asDBL_VECTOR() {

    Vector3d v = (*this);
    DBL_VECTOR vFinal = { v[0], v[1], v[2] };

    return vFinal;

}


/*
 * Redefines a vector from one coordinate basis to another using Euler Angles
 * 
 * For a vector defined in a new coordinate basis, (X',Y',Z'), we want to 
 * find it's values in a previous coordinate basis (X,Y,Z) so that we can 
 * add the new vector to vectors defined in the previous basis.
 *
 * For a change of basis from (X,Y,Z)->(X',Y',Z'),
 * ThetaE [0, pi] is the angle between Z and Z',
 * Vector N := Z x Z' (cross product)
 * PhiE [0, 2pi) is the angle between X and N
 * PsiE [0, 2pi) is the angle between X' and N
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
 * @param   [IN]   p_ThetaE                    Euler angle Theta (rad) 
 * @param   [IN]   p_PhiE                      Euler angle Phi   (rad) 
 * @param   [IN]   p_PsiE                      Euler angle Psi   (rad) 
 * @return                                     Vector in previous basis
 */
Vector3d Vector3d::RotateVector(const double p_ThetaE, const double p_PhiE, const double p_PsiE) {

    // Replace for convenience, undefine below
    #define cTheta cos(p_ThetaE)
    #define cPhi   cos(p_PhiE)
    #define cPsi   cos(p_PsiE)
    #define sTheta sin(p_ThetaE)
    #define sPhi   sin(p_PhiE)
    #define sPsi   sin(p_PsiE)

    // Define the Rotation Matrix     
    std::vector<DBL_VECTOR> RotationMatrix = {
        { cPhi*cPsi - sPhi*cTheta*sPsi ,  -cPhi*sPsi - sPhi*cTheta*cPsi ,  sTheta*sPhi },
        { sPhi*cPsi + cPhi*cTheta*sPsi ,  -sPhi*sPsi + cPhi*cTheta*cPsi , -sTheta*cPhi },
        { sTheta*sPsi                  ,  sTheta*cPsi                   ,  cTheta      }};

    #undef cTheta
    #undef cPhi
    #undef cPsi
    #undef sTheta
    #undef sPhi
    #undef sPsi

    // Multiply RotationMatrix * p_oldVector
    Vector3d oldVector = (*this);
    Vector3d newVector = Vector3d(0.0, 0.0, 0.0);

    for (size_t i=0; i< 3; i++) {
        for (size_t j=0; j<3; j++) {
            newVector[i] += RotationMatrix[i][j] * oldVector[j];
        }
    }
    return newVector;
}


/*
 * Returns the unit vector in the direction of the given vector
 * 
 *
 * @return                                     Unit vector of input
 */
Vector3d Vector3d::UnitVector() {

    Vector3d newVector = (*this)/(*this).Magnitude(); 

    return newVector;
}


///////////////////////////////////////////////////////////
// 
// Linear Algebra Functions
// 
///////////////////////////////////////////////////////////

namespace linalg {

    /*
     * Calculate the standard dot product of two vectors
     *
     *
     * @param   [IN]   a                            first vector
     * @param   [IN]   b                            second vector
     * @return                                      dot product
     */
    double dot(const Vector3d& p_a, const Vector3d& p_b) {
    
        double c = 0;
        for (size_t i=0; i<3; i++) {
            c += p_a[i] * p_b[i];
        }
        return c;
    }
    
    
    /*
     * Calculate the standard cross product of two vectors
     *
     *
     * @param   [IN]   a                            first vector
     * @param   [IN]   b                            second vector
     * @return                                      cross product
     */
    Vector3d cross(const Vector3d& p_a, const Vector3d& p_b) {
    
        Vector3d c = Vector3d(0, 0, 0);
    
        c[0] = p_a[1] * p_b[2] - p_a[2] * p_b[1];
        c[1] = p_a[2] * p_b[0] - p_a[0] * p_b[2];
        c[2] = p_a[0] * p_b[1] - p_a[1] * p_b[0];
    
        return c;
    }

    /*
     * Calculate the angle between two vectors.
     *
     *
     * @param   [IN]   a                            first vector
     * @param   [IN]   b                            second vector
     * @return                                      angle between them, in radians
     */
    double angleBetween(const Vector3d& p_a, const Vector3d& p_b) {
        // Angle between 2 vectors, between [0, pi] 
        return std::acos(linalg::dot(p_a, p_b) / (p_a.Magnitude() * p_b.Magnitude()));
    }
}

