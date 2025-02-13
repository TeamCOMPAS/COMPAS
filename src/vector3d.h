#ifndef __vector3d_h__
#define __vector3d_h__

#include <vector>
#include <array>

#include "constants.h"
#include "Errors.h"


class Vector3d {

public:

    // constructors
    Vector3d();
    Vector3d(const double p_x, const double p_y, const double p_z);
    Vector3d(DBL_VECTOR p_Vec);

    // destructor
    virtual ~Vector3d() {}


    // object identifiers - all classes have these 
    OBJECT_ID          ObjectId()          { return m_ObjectId; }                                       // object id for vectors - ordinal value from enum
    OBJECT_TYPE        ObjectType()        { return OBJECT_TYPE::NONE; }                                // object type for vectors - always "NONE"
    OBJECT_PERSISTENCE ObjectPersistence() { return OBJECT_PERSISTENCE::PERMANENT; }                    // object persistence for vectors - always "PERMANENT"
    STELLAR_TYPE       StellarType()       { return STELLAR_TYPE::NONE; }                               // stellar type for vectors - always "NONE"


    // getters 
    DBL_VECTOR  asDBL_VECTOR() const { return { (*this)[0], (*this)[1] , (*this)[2]} ; }
    double      Magnitude() const    { return std::sqrt(Dot(*this, *this)); }
    double      xValue() const       { return m_x; }
    double      yValue() const       { return m_y; }
    double      zValue() const       { return m_z; }

    // Operator overloads
    double   operator [] (const size_t p_i) const { return (*const_cast<Vector3d*>(this))[p_i]; }
    double&  operator [] (const size_t p_i) {

        THROW_ERROR_IF(p_i < 0 || p_i > 2, ERROR::INDEX_OUT_OF_RANGE);                                  // this is a code defect
        
             if (p_i == 0) return m_x;
        else if (p_i == 1) return m_y;
        else               return m_z;
    }
   
    void     operator = (const Vector3d p_Vec) { UpdateVector(p_Vec[0], p_Vec[1], p_Vec[2]); }

    void     operator += (const Vector3d p_Vec) { *this = (*this) + p_Vec; }

    Vector3d operator + (const Vector3d p_Vec) {
        Vector3d vec;
        for (size_t i=0; i<3; i++) vec[i] = (*this)[i] + p_Vec[i];
        return vec;
    }
    
    Vector3d operator - (const Vector3d p_Vec) {
        Vector3d vec;
        for (size_t i = 0; i < 3; i++) vec[i] = (*this)[i] - p_Vec[i];
        return vec;
    }
    
    Vector3d operator * (const double p_Scalar) {                                                       // for Vector3d * scalar
        Vector3d vec;
        for (size_t i = 0; i < 3; i++) vec[i] = (*this)[i] * p_Scalar; 
        return vec;
    }
    friend Vector3d operator * (const double p_Scalar, Vector3d p_Vec) { return p_Vec * p_Scalar; };    // for scalar * Vector3d

    Vector3d operator / (const double p_Scalar) {                                                       // for Vector3d / scalar
        Vector3d vec;
        for (size_t i = 0; i < 3; i++) vec[i] = (*this)[i] / p_Scalar;
        return vec;
    }

    friend std::ostream &operator << (std::ostream &p_os, Vector3d const p_Vec) { return p_os << "{" << p_Vec[0] << ", " << p_Vec[1] << ", " << p_Vec[2] << "}"; }


    // member functions 
    static double   AngleBetween(const Vector3d& p_Vec1, const Vector3d& p_Vec2) { return std::acos(Dot(p_Vec1, p_Vec2) / (p_Vec1.Magnitude() * p_Vec2.Magnitude())); }

    Vector3d        ChangeBasis(const double p_ThetaE, const double p_PhiE, const double p_PsiE);

    static Vector3d Cross(const Vector3d& p_a, const Vector3d& p_b) {
    
        Vector3d result = Vector3d(0.0, 0.0, 0.0);
    
        result[0] = p_a[1] * p_b[2] - p_a[2] * p_b[1];
        result[1] = p_a[2] * p_b[0] - p_a[0] * p_b[2];
        result[2] = p_a[0] * p_b[1] - p_a[1] * p_b[0];
    
        return result;
    }

    static double   Dot(const Vector3d& p_Vec1, const Vector3d& p_Vec2) {
        double result = 0.0;
        for (size_t i = 0; i < 3; i++) result += p_Vec1[i] * p_Vec2[i];
        return result;
    }

    Vector3d        MatrixMult(const std::vector<DBL_VECTOR>& p_matrix, const Vector3d& p_vector);

    Vector3d        RotateVectorAboutX(const double p_Theta);
    Vector3d        RotateVectorAboutY(const double p_Theta);
    Vector3d        RotateVectorAboutZ(const double p_Theta);

    Vector3d        UnitVector() { return *this / (*this).Magnitude(); } 


protected:

    OBJECT_ID   m_ObjectId;

    // member variables
    double m_x;
    double m_y;
    double m_z;


    // member functions
    void UpdateVector(const double p_x, const double p_y, const double p_z) { m_x = p_x; m_y = p_y; m_z = p_z; }
};

#endif // __vector3d_h__
