#include "Element.h"
#include <iostream>

Element::Element(int index, int materialIndex, int basisFunctionIndex, int quadratureIndex)
{
	mIndex = index;
	mMaterialIndex = materialIndex;
	mBasisFunctionIndex = basisFunctionIndex;
	mQuadratureIndex = quadratureIndex;
}

Element::~Element()
{
}

void Element::setSize(int size)
{
	mKe.resize(size, size);
	mKe.setZero();
	mDofIndexes.resize(size);
	for (int i = 0; i < size; i++)
	{
		mDofIndexes(i) = -1;
	}
	mNodalU.resize(size);
	mNodalU.setZero();
}

Element_1D::Element_1D(int index, int materailIndex, double A, Eigen::VectorXi endNodeIndex, int basisFunctionIndex, int quadratureIndex) :Element(index, materailIndex, basisFunctionIndex, quadratureIndex)
{
	mSectionArea = A;
	mLength = 0;
	mEndNodeIndex = endNodeIndex;
}

Element_1D::~Element_1D()
{
}

void Element_1D::computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mXi.size();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd dN;
	Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(n + 1, n + 1);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mXi[i];
		dN = basis->getdN(xi);
		Ke += w(i) * dN * dN.transpose();
	}
	mKe = 2. * mat->GetE() * mSectionArea / mLength * Ke;
	//std::cout << mKe << std::endl;
	//std::cout << std::endl;
}

void Element_1D::computeInternalForce(std::shared_ptr<BasisFunction> basis, int outputNumber, std::shared_ptr<Material> mat)
{
	mStrain.resize(outputNumber);
	mStrain.setZero();
	mStress.resize(outputNumber);
	mStress.setZero();
	mInternalForce.resize(outputNumber);
	mInternalForce.setZero();
	mU.resize(outputNumber);
	double step = 2.0 / (outputNumber - 1);
	for (int i = 0; i < outputNumber; i++)
	{
		double x = -1.0 + i * step;
		Eigen::VectorXd xi(1);
		xi << x;
		Eigen::VectorXd phi = basis->getN(xi);
		mU(i) = phi.transpose() * mNodalU;
		Eigen::VectorXd dphi = basis->getdN(xi);
		mStrain(i) = 2.0 / mLength * dphi.transpose() * mNodalU;
		mStress(i) = mat->GetE() * mStrain(i);
		mInternalForce(i) = mStress(i) * mSectionArea;
	}
}

double Element_1D::getU(double xi, std::shared_ptr<BasisFunction> basis)
{
	Eigen::VectorXd x(1);
	x << xi;
	Eigen::VectorXd phi = basis->getN(x);
	return phi.transpose() * mNodalU;
}

double Element_1D::getdU(double xi, std::shared_ptr<BasisFunction> basis)
{
	Eigen::VectorXd x(1);
	x << xi;
	Eigen::VectorXd dphi = basis->getdN(x);
	return 2.0 / mLength * dphi.transpose() * mNodalU;
}

Element_1D_ununiformSection::Element_1D_ununiformSection(int index, int materailIndex, std::shared_ptr<SextionAreaFunction> A, Eigen::VectorXi endNodeIndex, Eigen::VectorXd nodalCoordinates, int basisFunctionIndex, int quadratureIndex) :Element(index, materailIndex, basisFunctionIndex, quadratureIndex)
{
	mSectionArea = A;
	mLength = 0;
	mEndNodeIndex = endNodeIndex;
	mNodalCoordinates = nodalCoordinates;
}

Element_1D_ununiformSection::~Element_1D_ununiformSection()
{
}

void Element_1D_ununiformSection::computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mXi.size();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd dN;
	Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(n + 1, n + 1);
	Eigen::MatrixXd x = quadrature->getXInPhysicalSpace(mNodalCoordinates);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mXi[i];
		dN = basis->getdN(xi);
		Eigen::VectorXd X(1);
		X << x(i);
		Ke += w(i) * mSectionArea->getValue(X) * dN * dN.transpose();
	}
	mKe = 2. * mat->GetE() / mLength * Ke;
}

void Element_1D_ununiformSection::computeInternalForce(std::shared_ptr<BasisFunction> basis, int outputNumber, std::shared_ptr<Material> mat)
{
	mStrain.resize(outputNumber);
	mStrain.setZero();
	mStress.resize(outputNumber);
	mStress.setZero();
	mInternalForce.resize(outputNumber);
	mInternalForce.setZero();
	mU.resize(outputNumber);
	double step = 2.0 / (outputNumber - 1);
	for (int i = 0; i < outputNumber; i++)
	{
		double x = -1.0 + i * step;
		Eigen::VectorXd xi(1);
		xi << x;
		Eigen::VectorXd phi = basis->getN(xi);
		mU(i) = phi.transpose() * mNodalU;
		Eigen::VectorXd dphi = basis->getdN(xi);
		mStrain(i) = 2.0 / mLength * dphi.transpose() * mNodalU;
		mStress(i) = mat->GetE() * mStrain(i);
		Eigen::VectorXd X(1);
		X << mNodalCoordinates(0) + mLength / 2. * (1. + x);
		mInternalForce(i) = mStress(i) * mSectionArea->getValue(X);
	}
}

double Element_1D_ununiformSection::getU(double xi, std::shared_ptr<BasisFunction> basis)
{
	Eigen::VectorXd x(1);
	x << xi;
	Eigen::VectorXd phi = basis->getN(x);
	return phi.transpose() * mNodalU;
}

double Element_1D_ununiformSection::getdU(double xi, std::shared_ptr<BasisFunction> basis)
{
	Eigen::VectorXd x(1);
	x << xi;
	Eigen::VectorXd dphi = basis->getdN(x);
	return 2.0 / mLength * dphi.transpose() * mNodalU;
}
