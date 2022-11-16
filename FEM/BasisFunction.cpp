#include "BasisFunction.h"

BasisFunction::BasisFunction()
{
}

BasisFunction::~BasisFunction()
{
}

BasisFunction_1D_Linear::BasisFunction_1D_Linear()
{
	mOrder = 1;
	mDim = 1;
}

BasisFunction_1D_Linear::~BasisFunction_1D_Linear()
{
}

Eigen::VectorXd BasisFunction_1D_Linear::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	Eigen::VectorXd N(2);
	N(0) = 0.5 * (1 - xi);
	N(1) = 0.5 * (1 + xi);
	return N;
}

Eigen::VectorXd BasisFunction_1D_Linear::getdN(Eigen::VectorXd coordinates)
{
	Eigen::VectorXd dN(2);
	dN(0) = -0.5;
	dN(1) = 0.5;
	return dN;
}

BasisFunction_1D_Quadratic::BasisFunction_1D_Quadratic()
{
	mOrder = 2;
	mDim = 1;
}

BasisFunction_1D_Quadratic::~BasisFunction_1D_Quadratic()
{
}

Eigen::VectorXd BasisFunction_1D_Quadratic::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	Eigen::VectorXd N(3);
	N(0) = 0.5 * xi * (xi - 1);
	N(1) = 1 - xi * xi;
	N(2) = 0.5 * xi * (xi + 1);
	return N;
}

Eigen::VectorXd BasisFunction_1D_Quadratic::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	Eigen::VectorXd dN(3);
	dN(0) = xi - 0.5;
	dN(1) = -2 * xi;
	dN(2) = xi + 0.5;
	return dN;
}

BasisFunction_2D_Triangle_Linear::BasisFunction_2D_Triangle_Linear()
{
	mOrder = 1;
	mDim = 2;
}

BasisFunction_2D_Triangle_Linear::~BasisFunction_2D_Triangle_Linear()
{
}

Eigen::VectorXd BasisFunction_2D_Triangle_Linear::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd N(3);
	N(0) = 1 - xi - eta;
	N(1) = xi;
	N(2) = eta;
	return N;
}

Eigen::VectorXd BasisFunction_2D_Triangle_Linear::getdN(Eigen::VectorXd coordinates)
{
	Eigen::VectorXd dN(6);
	dN(0) = -1;
	dN(1) = -1;
	dN(2) = 1;
	dN(3) = 0;
	dN(4) = 0;
	dN(5) = 1;
	return dN;
}

BasisFunction_2D_Triangle_Quadratic::BasisFunction_2D_Triangle_Quadratic()
{
	mOrder = 2;
	mDim = 2;
}

BasisFunction_2D_Triangle_Quadratic::~BasisFunction_2D_Triangle_Quadratic()
{
}

Eigen::VectorXd BasisFunction_2D_Triangle_Quadratic::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd N(6);
	N(0) = 2 * xi * eta;
	N(1) = 2 * xi * (1 - xi - eta);
	N(2) = 2 * eta * (1 - xi - eta);
	N(3) = 4 * xi * (1 - xi - eta);
	N(4) = 4 * xi * eta;
	N(5) = 4 * eta * (1 - xi - eta);
	return N;
}

Eigen::VectorXd BasisFunction_2D_Triangle_Quadratic::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd dN(12);
	dN(0) = 2 * eta;
	dN(1) = 2 * xi;
	dN(2) = 0;
	dN(3) = -4 * xi - 4 * eta + 4;
	dN(4) = 4 * xi;
	dN(5) = 4 * eta;
	dN(6) = 0;
	dN(7) = -4 * xi - 4 * eta + 4;
	dN(8) = 4 * xi;
	dN(9) = 4 * eta;
	dN(10) = 4 * xi;
	dN(11) = 4 * eta;
	return dN;
}

BasisFunction_2D_Quad_Linear::BasisFunction_2D_Quad_Linear()
{
	mOrder = 1;
	mDim = 2;
}

BasisFunction_2D_Quad_Linear::~BasisFunction_2D_Quad_Linear()
{
}

Eigen::VectorXd BasisFunction_2D_Quad_Linear::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd N(4);
	N(0) = 0.25 * (1 - xi) * (1 - eta);
	N(1) = 0.25 * (1 + xi) * (1 - eta);
	N(2) = 0.25 * (1 + xi) * (1 + eta);
	N(3) = 0.25 * (1 - xi) * (1 + eta);
	return N;
}

Eigen::VectorXd BasisFunction_2D_Quad_Linear::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd dN(8);
	dN(0) = -0.25 * (1 - eta);
	dN(1) = -0.25 * (1 - xi);
	dN(2) = 0.25 * (1 - eta);
	dN(3) = -0.25 * (1 + xi);
	dN(4) = 0.25 * (1 + eta);
	dN(5) = 0.25 * (1 + xi);
	dN(6) = -0.25 * (1 + eta);
	dN(7) = 0.25 * (1 - xi);
	return dN;
}

BasisFunction_2D_Quad_Quadratic::BasisFunction_2D_Quad_Quadratic()
{
	mOrder = 2;
	mDim = 2;
}

BasisFunction_2D_Quad_Quadratic::~BasisFunction_2D_Quad_Quadratic()
{
}

Eigen::VectorXd BasisFunction_2D_Quad_Quadratic::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd N(9);
	N(0) = 0.25 * (1 - xi) * (1 - eta) * (1 - xi - eta);
	N(1) = 0.25 * (1 + xi) * (1 - eta) * (1 - xi + eta);
	N(2) = 0.25 * (1 + xi) * (1 + eta) * (1 + xi + eta);
	N(3) = 0.25 * (1 - xi) * (1 + eta) * (1 + xi - eta);
	N(4) = 0.5 * (1 - xi * xi) * (1 - eta);
	N(5) = 0.5 * (1 + xi) * (1 - eta * eta);
	N(6) = 0.5 * (1 - xi * xi) * (1 + eta);
	N(7) = 0.5 * (1 - xi) * (1 - eta * eta);
	N(8) = (1 - xi * xi) * (1 - eta * eta);
	return N;
}

Eigen::VectorXd BasisFunction_2D_Quad_Quadratic::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	Eigen::VectorXd dN(18);
	dN(0) = -0.25 * (1 - eta) * (2 * xi + eta);
	dN(1) = -0.25 * (1 - xi) * (xi + 2 * eta);
	dN(2) = 0.25 * (1 - eta) * (2 * xi + eta);
	dN(3) = -0.25 * (1 + xi) * (xi - 2 * eta);
	dN(4) = 0.5 * (-2 * xi) * (1 - eta);
	dN(5) = -0.5 * (1 - eta * eta);
	dN(6) = 0.25 * (1 + eta) * (2 * xi + eta);
	dN(7) = 0.25 * (1 + xi) * (xi + 2 * eta);
	dN(8) = 0.5 * (1 - eta * eta);
	dN(9) = 0.25 * (1 + eta) * (2 * xi + eta);
	dN(10) = 0.25 * (1 + xi) * (xi + 2 * eta);
	dN(11) = 0.5 * (1 - eta * eta);
	dN(12) = -0.25 * (1 + eta) * (2 * xi - eta);
	dN(13) = 0.25 * (1 - xi) * (xi - 2 * eta);
	dN(14) = -0.5 * (1 - eta * eta);
	dN(15) = -0.25 * (1 - xi) * (xi + 2 * eta);
	dN(16) = -0.5 * (-2 * xi) * (1 + eta);
	dN(17) = 0.25 * (1 - eta) * (2 * xi - eta);
	return dN;
}

BasisFunction_3D_Tetrahedron_Linear::BasisFunction_3D_Tetrahedron_Linear()
{
	mOrder = 1;
	mDim = 3;
}

BasisFunction_3D_Tetrahedron_Linear::~BasisFunction_3D_Tetrahedron_Linear()
{
}

Eigen::VectorXd BasisFunction_3D_Tetrahedron_Linear::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	double zeta = coordinates(3);
	Eigen::VectorXd N(4);
	N(0) = 1 - xi - eta - zeta;
	N(1) = xi;
	N(2) = eta;
	N(3) = zeta;
	return N;
}

Eigen::VectorXd BasisFunction_3D_Tetrahedron_Linear::getdN(Eigen::VectorXd coordinates)
{
	Eigen::VectorXd dN(12);
	dN(0) = -1;
	dN(1) = -1;
	dN(2) = -1;
	dN(3) = 1;
	dN(4) = 0;
	dN(5) = 0;
	dN(6) = 0;
	dN(7) = 1;
	dN(8) = 0;
	dN(9) = 0;
	dN(10) = 0;
	dN(11) = 1;
	return dN;
}

BasisFunction_3D_Tetrahedron_Quadratic::BasisFunction_3D_Tetrahedron_Quadratic()
{
	mOrder = 2;
	mDim = 3;
}

BasisFunction_3D_Tetrahedron_Quadratic::~BasisFunction_3D_Tetrahedron_Quadratic()
{
}

Eigen::VectorXd BasisFunction_3D_Tetrahedron_Quadratic::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	double zeta = coordinates(3);
	Eigen::VectorXd N(10);
	N(0) = 1 - xi - eta - zeta;
	N(1) = xi;
	N(2) = eta;
	N(3) = zeta;
	N(4) = xi * (1 - xi - eta - zeta);
	N(5) = eta * (1 - xi - eta - zeta);
	N(6) = zeta * (1 - xi - eta - zeta);
	N(7) = xi * eta;
	N(8) = xi * zeta;
	N(9) = eta * zeta;
	return N;
}

Eigen::VectorXd BasisFunction_3D_Tetrahedron_Quadratic::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	double zeta = coordinates(3);
	Eigen::VectorXd dN(30);
	dN(0) = -1;
	dN(1) = -1;
	dN(2) = -1;
	dN(3) = 1;
	dN(4) = 0;
	dN(5) = 0;
	dN(6) = 0;
	dN(7) = 1;
	dN(8) = 0;
	dN(9) = 0;
	dN(10) = 0;
	dN(11) = 1;
	dN(12) = -1 + 2 * xi + eta + zeta;
	dN(13) = -xi;
	dN(14) = -xi;
	dN(15) = -xi;
	dN(16) = 1 - 2 * xi;
	dN(17) = 0;
	dN(18) = 0;
	dN(19) = 0;
	dN(20) = 0;
	dN(21) = 0;
	dN(22) = 0;
	dN(23) = 0;
	dN(24) = 0;
	dN(25) = 0;
	dN(26) = 0;
	dN(27) = 0;
	dN(28) = 0;
	dN(29) = 0;
	return dN;
}

BasisFunction_3D_Hexahedron_Linear::BasisFunction_3D_Hexahedron_Linear()
{
	mOrder = 1;
	mDim = 3;
}

BasisFunction_3D_Hexahedron_Linear::~BasisFunction_3D_Hexahedron_Linear()
{
}

Eigen::VectorXd BasisFunction_3D_Hexahedron_Linear::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	double zeta = coordinates(3);
	Eigen::VectorXd N(8);
	N(0) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
	N(1) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
	N(2) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
	N(3) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
	N(4) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
	N(5) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
	N(6) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
	N(7) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);
	return N;
}

Eigen::VectorXd BasisFunction_3D_Hexahedron_Linear::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double eta = coordinates(1);
	double zeta = coordinates(3);
	Eigen::VectorXd dN(24);
	//N0
	dN( 0) = -0.125 * (1 - eta) * (1 - zeta);
	dN( 1) = -0.125 * (1 -  xi) * (1 - zeta);
	dN( 2) = -0.125 * (1 -  xi) * (1 -  eta);
	//N1
	dN( 3) =  0.125 * (1 - eta) * (1 - zeta);
	dN( 4) = -0.125 * (1 +  xi) * (1 - zeta);
	dN( 5) = -0.125 * (1 +  xi) * (1 -  eta);
	//N2
	dN( 6) =  0.125 * (1 + eta) * (1 - zeta);
	dN( 7) =  0.125 * (1 +  xi) * (1 - zeta);
	dN( 8) = -0.125 * (1 +  xi) * (1 +  eta);
	//N3
	dN( 9) = -0.125 * (1 + eta) * (1 - zeta);
	dN(10) =  0.125 * (1 -  xi) * (1 - zeta);
	dN(11) = -0.125 * (1 -  xi) * (1 +  eta);
	//N4
	dN(12) = -0.125 * (1 - eta) * (1 + zeta);
	dN(13) = -0.125 * (1 -  xi) * (1 + zeta);
	dN(14) =  0.125 * (1 -  xi) * (1 -  eta);
	//N5
	dN(15) =  0.125 * (1 - eta) * (1 + zeta);
	dN(16) = -0.125 * (1 +  xi) * (1 + zeta);
	dN(17) =  0.125 * (1 +  xi) * (1 -  eta);
	//N6
	dN(18) =  0.125 * (1 + eta) * (1 + zeta);
	dN(19) =  0.125 * (1 +  xi) * (1 + zeta);
	dN(20) =  0.125 * (1 +  xi) * (1 +  eta);
	//N7
	dN(21) = -0.125 * (1 + eta) * (1 + zeta);
	dN(22) =  0.125 * (1 -  xi) * (1 + zeta);
	dN(23) =  0.125 * (1 -  xi) * (1 +  eta);
	return dN;
}

BasisFunction_1D_Hermite::BasisFunction_1D_Hermite()
{
	mOrder = 3;
	mDim = 1;
}

BasisFunction_1D_Hermite::~BasisFunction_1D_Hermite()
{
}

Eigen::VectorXd BasisFunction_1D_Hermite::getN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double h = coordinates(1);
	Eigen::VectorXd N(4);
	N(0) = 0.25 * (2. - 3. * xi + xi * xi * xi);
	N(1) = h / 8. * (1. - xi - xi * xi + xi * xi * xi);
	N(2) = 0.25 * (2. + 3. * xi - xi * xi * xi);
	N(3) = h / 8. * (-1. - xi + xi * xi + xi * xi * xi);
	return N;
}

Eigen::VectorXd BasisFunction_1D_Hermite::getdN(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double h = coordinates(1);
	Eigen::VectorXd dN(4);
	dN(0) = 0.75 * (-1. + xi * xi);
	dN(1) = h / 8. * (-1. - 2. * xi + 3. * xi * xi);
	dN(2) = 0.75 * (1. - xi * xi);
	dN(3) = h / 8. * (-1. + 2. * xi + 3. * xi * xi);
	return dN;
}

Eigen::VectorXd BasisFunction_1D_Hermite::getd2N(Eigen::VectorXd coordinates)
{
	double xi = coordinates(0);
	double h = coordinates(1);
	Eigen::VectorXd d2N(4);
	d2N(0) = 1.5 * xi;
	d2N(1) = h / 4. * (-1 + 3. * xi);
	d2N(2) = -1.5 * xi;
	d2N(3) = h / 4. * (1 + 3. * xi);
	return d2N;
}
