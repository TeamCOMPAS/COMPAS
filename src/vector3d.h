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
    Vector3d(const double x,
             const double y,
             const double z);
    Vector3d(DBL_VECTOR v);

    // object identifiers - all classes have these 
    OBJECT_ID    ObjectId()    { return m_ObjectId; }           // object id for vectors - ordinal value from enum
    OBJECT_TYPE  ObjectType()  { return OBJECT_TYPE::NONE; }    // object type for vectors - always "NONE"
    STELLAR_TYPE StellarType() { return STELLAR_TYPE::NONE; }   // stellar type for vectors - always "NONE"

    // destructor
    virtual ~Vector3d() {}

    // getters 
    double      xValue() const       { return m_0; }
    double      yValue() const       { return m_1; }
    double      zValue() const       { return m_2; }
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
    double& operator[] (size_t i);
    double operator[] (size_t i) const { return (*const_cast<Vector3d*>(this))[i]; }
    
    // Overload setting operator
    void operator =(Vector3d vNew);

    // Overload += setting vector
    void operator +=(Vector3d vec) { (*this) = (*this) + vec; }

    // Overload vector addition operator
    Vector3d operator +(Vector3d vec);
    
    // Overload vector subtraction operator
    Vector3d operator -(Vector3d vec);
    
    // Overload scalar multiplication operator (only works as v*s, see below for s*v)
    Vector3d operator *(const double scalar); 

    // Overload scalar division operator
    Vector3d operator /(const double scalar);


    
    

protected:

    OBJECT_ID   m_ObjectId;

    // member variables
    double m_0;
    double m_1;
    double m_2;

    // member functions
    void updateVector( const double x, const double y, const double z);
};

// Free functions, 
// Overload * with reversed order of inputs
Vector3d operator *(double scalar, Vector3d vec);

// Overload output << operator for Vector3d
std::ostream &operator <<(std::ostream &os, Vector3d const vec);



///////////////////////////////////////////////////////////
// 
// Linear Algebra Functions
// 
///////////////////////////////////////////////////////////

namespace linalg {

    double      dot(const Vector3d& a, const Vector3d& b);

    Vector3d    cross(const Vector3d& a, const Vector3d& b);

    double      angleBetween(const Vector3d& a, const Vector3d& b);
}

#endif // __vector3d_h__
