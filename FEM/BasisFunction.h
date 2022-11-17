#pragma once

#include <eigen/Eigen/Core>

struct Quadrature
{
	Eigen::VectorXd mXi;
	Eigen::VectorXd mWeight;
	//int mDim;

	Quadrature(int mDim, int mOrder)
	{
		if (mDim == 1)
		{
			if (mOrder == 1)
			{
				mXi.resize(3);
				mWeight.resize(3);
				mXi << -sqrt(3. / 5.), 0, sqrt(3. / 5.);
				mWeight << 0.555555555555556, 0.888888888888889, 0.555555555555556;
			}
			else if (mOrder == 2)
			{
				mXi.resize(4);
				mWeight.resize(4);
				mXi << -sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.)), -sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)), sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)), sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.));
				mWeight << (18. - sqrt(30.)) / 36., (18. + sqrt(30.)) / 36., (18. + sqrt(30.)) / 36., (18. - sqrt(30.)) / 36.;
			}
			else if (mOrder == 3)
			{
				mXi.resize(5);
				mWeight.resize(5);
				mXi << -sqrt(5. + 2. * sqrt(10. / 7.)) / 3., -sqrt(5. - 2. * sqrt(10. / 7.)) / 3., 0, sqrt(5. - 2. * sqrt(10. / 7.)) / 3., sqrt(5. + 2. * sqrt(10. / 7.)) / 3.;
				mWeight << (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225., (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900.;
			}
		}
		else if (mDim == 2)
		{
			mXi.resize(8);
			mWeight.resize(4);
			mXi << -1. / sqrt(3.), -1. / sqrt(3.), 
				    1. / sqrt(3.), -1. / sqrt(3.),
				    1. / sqrt(3.),  1. / sqrt(3.), 
				   -1. / sqrt(3.),  1. / sqrt(3.);
			mWeight << 1., 1., 1., 1.;
			
		}
		else if (mDim == 3)
		{
			mXi.resize(24);
			mWeight.resize(8);
			mXi << -1. / sqrt(3.), -1. / sqrt(3.),  1. / sqrt(3.),
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

	Eigen::VectorXd getXInPhysicalSpace(Eigen::VectorXd nodalCoordinates)
	{
		Eigen::VectorXd xiInPhysicalSpace;
		xiInPhysicalSpace.resize(mXi.size());
		xiInPhysicalSpace = nodalCoordinates(0) + (nodalCoordinates(1) - nodalCoordinates(0)) / 2. * (mXi.array() + 1.);
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
	BodyForceFunction(double value):mValue(value) { }
	~BodyForceFunction() { }
	
	double getValue(Eigen::VectorXd coordinates) { return mValue; }

	Eigen::VectorXd getElementLoadVector(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, Eigen::VectorXd nodalCoordinates)
	{
		Eigen::VectorXd elementLoadVector;
		elementLoadVector.resize(basis->mOrder + 1);
		elementLoadVector.setZero();
		Eigen::VectorXd xInPhysicalSpace = quadrature->getXInPhysicalSpace(nodalCoordinates);
		for (int i = 0; i < quadrature->mXi.size(); i++)
		{
			Eigen::VectorXd x;
			x.resize(1);
			x << xInPhysicalSpace[i];
			Eigen::VectorXd xi;
			xi.resize(1);
			xi << quadrature->mXi[i];
			elementLoadVector += quadrature->mWeight[i] * getValue(x) * basis->getN(xi);
		}
		elementLoadVector *= (nodalCoordinates(1) - nodalCoordinates(0)) / 2.;
		return elementLoadVector;
	}

	double mValue;
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