#include <iostream>
#include "vector3d.h"

int main() {

    Vector3d v1(1,2,3);
    Vector3d v2(4,5,6);

    v1 += v2;

    std::cout << v1[0] << v1[1] << v1[2] << std::endl;

    //Vector3d v3 = Vector3d(0, 1, 2);

    //Vector3d v4(3,4,5);
    //v4 = v3;

    //Vector3d v5;

    //v5 = Vector3d();

    //std::cout << v3[0] << v3[1] << v3[2] << std::endl;
    //std::cout << v4[0] << v4[1] << v4[2] << std::endl;
    //std::cout << v5[0] << v5[1] << v5[2] << std::endl;

    //double DP = linalg::dot(v1, v2);
    //Vector3d CP = linalg::cross(v1, v2);

    //std::cout << DP << std::endl;
    //std::cout << v1[0] << std::endl;
    //std::cout << v1.Magnitude() << std::endl;
    //std::cout << CP << std::endl;
    //std::cout << CP[0] << std::endl;
    //std::cout << CP[1] << std::endl;
    //std::cout << CP[2] << std::endl;

    //double pi = 3.1415926;

    //double theta = .01;
    //double phi = .01;
    //double psi = .01;

    //Vector3d rotVec = v1.RotateVector(theta, phi, psi);
    ////std::cout << rotVec[0] << std::endl;
    ////std::cout << rotVec[1] << std::endl;
    ////std::cout << rotVec[2] << std::endl;

    //Vector3d sumVec1 = v1 + v2;
    //Vector3d sumVec2 = v2 + v1;
    //Vector3d diffVec1 = v1 - v2;
    //Vector3d diffVec2 = v2 - v1;
    //Vector3d prodVec1 = v1 * 3.0;
    //Vector3d prodVec2 = 3.0 * v1;

    ////std::cout << v1[0] << v1[1] << v1[2] << std::endl;
    ////std::cout << v2[0] << v2[1] << v2[2] << std::endl;

    ////std::cout << sumVec1[0] << sumVec1[1] << sumVec1[2] << std::endl;
    ////std::cout << sumVec2[0] << sumVec2[1] << sumVec2[2] << std::endl;

    ////std::cout << diffVec1[0] << diffVec1[1] << diffVec1[2] << std::endl;
    ////std::cout << diffVec2[0] << diffVec2[1] << diffVec2[2] << std::endl;

    ////std::cout << prodVec1[0] << " " << prodVec1[1] << " " << prodVec1[2] << " " << std::endl;
    ////std::cout << prodVec2[0] << " " << prodVec2[1] << " " << prodVec2[2] << " " << std::endl;

    //std::cout << rotVec[0] << " " << rotVec[1] << " " << rotVec[2] << std::endl;

    return 0;
}
