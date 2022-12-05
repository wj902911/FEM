#pragma once

#include <eigen/Eigen/Core>
#include <iostream>

struct Quadrature
{
	Eigen::MatrixXd mPoints;
	Eigen::VectorXd mWeight;
	//int mDim;

	Quadrature(int mDim, int mOrder)
	{
		if (mDim == 1)
		{
			if (mOrder == 1)
			{
				mPoints.resize(2, 1);
				mWeight.resize(2);
				mPoints << -sqrt(1. / 3.), sqrt(1. / 3.);
				mWeight << 1., 1.;
			}
			else if (mOrder == 2)
			{
				mPoints.resize(3, 1);
				mWeight.resize(3);
				mPoints << -sqrt(3. / 5.), 0, sqrt(3. / 5.);
				mWeight << 0.555555555555556, 0.888888888888889, 0.555555555555556;
			}
			else if (mOrder == 3)
			{
				mPoints.resize(4, 1);
				mWeight.resize(4);
				mPoints << -sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.)), -sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)), sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)), sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.));
				mWeight << (18. - sqrt(30.)) / 36., (18. + sqrt(30.)) / 36., (18. + sqrt(30.)) / 36., (18. - sqrt(30.)) / 36.;
			}
		}
		else if (mDim == 2)
		{
			mPoints.resize(4,2);
			mWeight.resize(4);
			mPoints << -1. / sqrt(3.), -1. / sqrt(3.),
				        1. / sqrt(3.), -1. / sqrt(3.),
				        1. / sqrt(3.),  1. / sqrt(3.), 
				       -1. / sqrt(3.),  1. / sqrt(3.);
			mWeight << 1., 1., 1., 1.;
			
		}
		else if (mDim == 3)
		{
			mPoints.resize(8, 3);
			mWeight.resize(8);
			mPoints << -1. / sqrt(3.), -1. / sqrt(3.),  1. / sqrt(3.),
				    1. / sqrt(3.), -1. / sqrt(3.),  1. / sqrt(3.),
				    1. / sqrt(3.),  1. / sqrt(3.),  1. / sqrt(3.),
				   -1. / sqrt(3.),  1. / sqrt(3.),  1. / sqrt(3.),
				   -1. / sqrt(3.), -1. / sqrt(3.), -1. / sqrt(3.),
				    1. / sqrt(3.), -1. / sqrt(3.), -1. / sqrt(3.),
				    1. / sqrt(3.),  1. / sqrt(3.), -1. / sqrt(3.),
				   -1. / sqrt(3.),  1. / sqrt(3.), -1. / sqrt(3.),
			mWeight << 1., 1., 1., 1., 1., 1., 1., 1.;
		}
		
	}
	
	~Quadrature()
	{
	}

	//Can only use in 1D
	Eigen::MatrixXd getXInPhysicalSpace(Eigen::MatrixXd nodalCoordinates)
	{
		Eigen::MatrixXd xiInPhysicalSpace;
		xiInPhysicalSpace.resize(mPoints.rows(), mPoints.cols());
		xiInPhysicalSpace = nodalCoordinates(0) + (nodalCoordinates(1) - nodalCoordinates(0)) / 2. * (mPoints.array() + 1.);
		return xiInPhysicalSpace;
	}
};

class BasisFunction
{
public:
	BasisFunction();
	~BasisFunction();
	virtual Eigen::VectorXd getN(Eigen::VectorXd coordinates) = 0;
	virtual Eigen::VectorXd getdN(Eigen::VectorXd coordinates) = 0;

	void SetNumDofEachNode(int numDofEachNode);

	int mDim;
	int mOrder;
	int mNumFunctions;
	int mNumNodes;
	Eigen::VectorXi mNumDofEachNode;
};

class BasisFunction_1D_Linear : public BasisFunction
{
public:
	BasisFunction_1D_Linear();
	~BasisFunction_1D_Linear();
	
	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_1D_Quadratic : public BasisFunction
{
public:
	BasisFunction_1D_Quadratic();
	~BasisFunction_1D_Quadratic();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_1D_Hermite : public BasisFunction
{
public:
	BasisFunction_1D_Hermite();
	~BasisFunction_1D_Hermite();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getd2N(Eigen::VectorXd coordinates);
};

class BasisFunction_2D_Triangle_Linear : public BasisFunction
{
public:
	BasisFunction_2D_Triangle_Linear();
	~BasisFunction_2D_Triangle_Linear();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_2D_Triangle_Quadratic : public BasisFunction
{
public:
	BasisFunction_2D_Triangle_Quadratic();
	~BasisFunction_2D_Triangle_Quadratic();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_2D_Quad_Linear : public BasisFunction
{
public:
	BasisFunction_2D_Quad_Linear();
	~BasisFunction_2D_Quad_Linear();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_2D_Quad_Quadratic : public BasisFunction
{
public:
	BasisFunction_2D_Quad_Quadratic();
	~BasisFunction_2D_Quad_Quadratic();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_3D_Tetrahedron_Linear : public BasisFunction
{
public:
	BasisFunction_3D_Tetrahedron_Linear();
	~BasisFunction_3D_Tetrahedron_Linear();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_3D_Tetrahedron_Quadratic : public BasisFunction
{
public:
	BasisFunction_3D_Tetrahedron_Quadratic();
	~BasisFunction_3D_Tetrahedron_Quadratic();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};

class BasisFunction_3D_Hexahedron_Linear : public BasisFunction
{
public:
	BasisFunction_3D_Hexahedron_Linear();
	~BasisFunction_3D_Hexahedron_Linear();

	Eigen::VectorXd getN(Eigen::VectorXd coordinates);
	Eigen::VectorXd getdN(Eigen::VectorXd coordinates);
};


class Function
{
public:
	Function() {}
	~Function() {}
	virtual double getValue(Eigen::VectorXd coordinates) = 0;
};

class BodyForceFunction : public Function
{
public:
	BodyForceFunction(Eigen::VectorXd value):mValue(value) { }
	~BodyForceFunction() { }
	
	double getValue(Eigen::VectorXd coordinates) { return mValue(0) * 2. * (2 - coordinates(0) / mValue(1)); }

	Eigen::VectorXd getElementLoadVector(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, Eigen::MatrixXd nodalCoordinates)
	{
		Eigen::VectorXd elementLoadVector;
		int nShapeF = basis->mNumFunctions;
		elementLoadVector.resize(nShapeF);
		elementLoadVector.setZero();
		for (int i = 0; i < quadrature->mPoints.rows(); i++)
		{
			Eigen::VectorXd CaussPt = quadrature->mPoints.row(i);
			Eigen::VectorXd xInPhysicalSpace = basis->getN(CaussPt).transpose() * nodalCoordinates;
			Eigen::Matrix2d J = (basis->getdN(CaussPt).reshaped(2, nShapeF) * nodalCoordinates).transpose();
			elementLoadVector += quadrature->mWeight[i] * getValue(xInPhysicalSpace) * basis->getN(CaussPt) * J.determinant();
		}
		//std::cout << elementLoadVector << std::endl;
		return elementLoadVector;
	}

	Eigen::VectorXd mValue;
};

class SextionAreaFunction : public Function
{
public:
	SextionAreaFunction(double A0, double L) :mA0(A0), mL(L) { }
	~SextionAreaFunction() { }

	double getValue(Eigen::VectorXd coordinates) 
	{ 
		double x = coordinates[0];
		if (x <= mL)
		{
			return mA0 * (1. - coordinates[0] / 2. / mL);
		}
		else
		{
			return mA0 / 2.;
		}
	}

	double mA0, mL;
};