#ifndef __vector3d_h__
#define __vector3d_h__

#include "constants.h"
#include <vector>
#include <array>

///////////////////////////////////////////////////////////
// 
// Vector3d object
// 
///////////////////////////////////////////////////////////

class Vector3d {

public:

    // constructors & initializers
    Vector3d();
    Vector3d(const double p_x,
             const double p_y,
             const double p_z);
    Vector3d(DBL_VECTOR p_v);

    // object identifiers - all classes have these 
    OBJECT_ID    ObjectId()    { return m_ObjectId; }           // object id for vectors - ordinal value from enum
    OBJECT_TYPE  ObjectType()  { return OBJECT_TYPE::NONE; }    // object type for vectors - always "NONE"
    STELLAR_TYPE StellarType() { return STELLAR_TYPE::NONE; }   // stellar type for vectors - always "NONE"

    // destructor
    virtual ~Vector3d() {}

    // getters 
    double      xValue() const       { return m_x; }
    double      yValue() const       { return m_y; }
    double      zValue() const       { return m_z; }
    double      Magnitude() const;

    DBL_VECTOR  asDBL_VECTOR();
    
    // member functions 
    Vector3d    RotateVector( const double p_ThetaE, 
                              const double p_PhiE, 
                              const double p_PsiE);
    Vector3d    UnitVector();

    ///////////////////////////////
    // Overload operators
    
    // Overload Indexing operator
    double& operator[] (size_t p_i);
    double operator[] (size_t p_i) const { return (*const_cast<Vector3d*>(this))[p_i]; }
    
    // Overload setting operator
    void operator =(Vector3d p_NewVec);

    // Overload += setting vector
    void operator +=(Vector3d p_Vec) { (*this) = (*this) + p_Vec; }

    // Overload vector addition operator
    Vector3d operator +(Vector3d p_Vec);
    
    // Overload vector subtraction operator
    Vector3d operator -(Vector3d p_Vec);
    
    // Overload scalar multiplication operator (only works as v*s, see below for s*v)
    Vector3d operator *(const double p_Scalar); 

    // Overload scalar division operator
    Vector3d operator /(const double p_Scalar);


    
    

protected:

    OBJECT_ID   m_ObjectId;

    // member variables
    double m_x;
    double m_y;
    double m_z;

    // member functions
    void UpdateVector(const double p_x, const double p_y, const double p_z);
};

// Free functions, 
// Overload * with reversed order of inputs
Vector3d operator *(double p_Scalar, Vector3d p_Vec);

// Overload output << operator for Vector3d
std::ostream &operator <<(std::ostream &os, Vector3d const p_Vec);



///////////////////////////////////////////////////////////
// 
// Linear Algebra Functions
// 
///////////////////////////////////////////////////////////

namespace linalg {

    double      dot(const Vector3d& p_a, const Vector3d& p_b);

    Vector3d    cross(const Vector3d& p_a, const Vector3d& p_b);

    double      angleBetween(const Vector3d& p_a, const Vector3d& p_b);
}

#endif // __vector3d_h__
