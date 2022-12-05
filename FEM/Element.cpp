#include "Element.h"
#include <iostream>
#include <eigen/Eigen/LU>

Element::Element(int index, int materialIndex, int basisFunctionIndex, int quadratureIndex)
{
	mIndex = index;
	mMaterialIndex = materialIndex;
	mBasisFunctionIndex = basisFunctionIndex;
	mQuadratureIndex = quadratureIndex;
	mLength = 0;
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
	mEndNodeIndex = endNodeIndex;
}

Element_1D::~Element_1D()
{
}

void Element_1D::evaluate(Eigen::MatrixXd cords,
	std::shared_ptr<BasisFunction> basis,
	std::shared_ptr<Quadrature> quadrature)
{
}

void Element_1D::computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mPoints.rows();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd dN;
	Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(n + 1, n + 1);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mPoints(i, 0);
		dN = basis->getdN(xi);
		Ke += w(i) * dN * dN.transpose();
	}
	mKe = 2. * mat->GetE() * mSectionArea / mLength * Ke;
	//std::cout << mKe << std::endl;
	//std::cout << std::endl;
}

void Element_1D::computeMassMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mPoints.rows();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd N;
	Eigen::MatrixXd Me = Eigen::MatrixXd::Zero(n + 1, n + 1);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mPoints(i, 0);
		N = basis->getN(xi);
		Me += w(i) * N * N.transpose();
	}
	mKe = mat->GetDensity() * mSectionArea * mLength / 2. * Me;
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

double Element_1D::getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	//Eigen::VectorXd x(1);
	//x << xi;
	Eigen::VectorXd phi = basis->getN(xi);
	return phi.transpose() * mNodalU;
}

double Element_1D::getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	//Eigen::VectorXd x(1);
	//x << xi;
	Eigen::VectorXd dphi = basis->getdN(xi);
	return 2.0 / mLength * dphi.transpose() * mNodalU;
}

Element_1D_ununiformSection::Element_1D_ununiformSection(int index, int materailIndex, std::shared_ptr<SextionAreaFunction> A, Eigen::VectorXi endNodeIndex, Eigen::VectorXd nodalCoordinates, int basisFunctionIndex, int quadratureIndex) :Element(index, materailIndex, basisFunctionIndex, quadratureIndex)
{
	mSectionArea = A;
	mEndNodeIndex = endNodeIndex;
	mNodalCoordinates = nodalCoordinates;
}

Element_1D_ununiformSection::~Element_1D_ununiformSection()
{
}

void Element_1D_ununiformSection::evaluate(Eigen::MatrixXd cords,
	std::shared_ptr<BasisFunction> basis,
	std::shared_ptr<Quadrature> quadrature)
{
}

void Element_1D_ununiformSection::computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mPoints.rows();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd dN;
	Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(n + 1, n + 1);
	Eigen::MatrixXd x = quadrature->getXInPhysicalSpace(mNodalCoordinates);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mPoints(i, 0);
		dN = basis->getdN(xi);
		Eigen::VectorXd X(1);
		X << x(i);
		Ke += w(i) * mSectionArea->getValue(X) * dN * dN.transpose();
	}
	mKe = 2. * mat->GetE() / mLength * Ke;
}

void Element_1D_ununiformSection::computeMassMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mPoints.rows();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd N;
	Eigen::MatrixXd Me = Eigen::MatrixXd::Zero(n + 1, n + 1);
	Eigen::MatrixXd x = quadrature->getXInPhysicalSpace(mNodalCoordinates);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mPoints(i, 0);
		N = basis->getN(xi);
		Eigen::VectorXd X(1);
		X << x(i);
		Me += w(i) * mSectionArea->getValue(X) * N * N.transpose();
	}
	mMe = mat->GetDensity() * mLength / 2. * Me;
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

double Element_1D_ununiformSection::getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	//Eigen::VectorXd x(1);
	//x << xi;
	Eigen::VectorXd phi = basis->getN(xi);
	return phi.transpose() * mNodalU;
}

double Element_1D_ununiformSection::getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	//Eigen::VectorXd x(1);
	//x << xi;
	Eigen::VectorXd dphi = basis->getdN(xi);
	return 2.0 / mLength * dphi.transpose() * mNodalU;
}

Element_1D_Beam::Element_1D_Beam(int index, int materailIndex, double A, double I, Eigen::VectorXi endNodeIndex, int basisFunctionIndex, int quadratureIndex) :Element(index, materailIndex, basisFunctionIndex, quadratureIndex)
{
	mSectionArea = A;
	mI = I;
	mEndNodeIndex = endNodeIndex;
}

Element_1D_Beam::~Element_1D_Beam()
{
}

void Element_1D_Beam::evaluate(Eigen::MatrixXd cords, std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature)
{
}

void Element_1D_Beam::computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int nDof = basis->mNumDofEachNode.sum();
	int nGauss = quadrature->mPoints.rows();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd d2N;
	Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(nDof, nDof);
	std::shared_ptr<BasisFunction_1D_Hermite> pBasis = std::dynamic_pointer_cast<BasisFunction_1D_Hermite>(basis);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(2);
		xi << quadrature->mPoints(i, 0), mLength;
		d2N = pBasis->getd2N(xi);
		Ke += w(i) * d2N * d2N.transpose();
	}
	mKe = 8. * mat->GetE() * mI / pow(mLength, 3) * Ke;
}

void Element_1D_Beam::computeMassMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int nDof = basis->mNumDofEachNode.sum();
	int nGauss = quadrature->mPoints.rows();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd N;
	Eigen::MatrixXd Me = Eigen::MatrixXd::Zero(nDof, nDof);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(2);
		xi << quadrature->mPoints(i, 0), mLength;
		N = basis->getN(xi);
		Me += w(i) * N * N.transpose();
	}
	mMe = mat->GetDensity() * mSectionArea * mLength / 2. * Me;
}

void Element_1D_Beam::computeInternalForce(std::shared_ptr<BasisFunction> basis, int outputNumber, std::shared_ptr<Material> mat)
{
}

double Element_1D_Beam::getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	return 0.0;
}

double Element_1D_Beam::getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	return 0.0;
}

Element_2D_Membrane::Element_2D_Membrane(int index, int materailIndex, Eigen::VectorXi endNodeIndex, int basisFunctionIndex, int quadratureIndex) :Element(index, materailIndex, basisFunctionIndex, quadratureIndex)
{
	mEndNodeIndex = endNodeIndex;
}

Element_2D_Membrane::~Element_2D_Membrane()
{
}

void Element_2D_Membrane::evaluate(Eigen::MatrixXd cords, std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature)
{
	mNodalCoordinates = cords;
	int nGauss = quadrature->mPoints.rows();
	int nShapeF = basis->mNumFunctions;
	mValues.resize(nShapeF, nGauss);
	mDerivatives.resize(nShapeF * 2, nGauss);
	mMp.resize(4, nGauss);
	for (int i = 0; i < nGauss; i++)
	{
		mValues.col(i) = basis->getN(quadrature->mPoints.row(i));
		mDerivatives.col(i) = basis->getdN(quadrature->mPoints.row(i));
		mMp.col(i) = (mDerivatives.col(i).reshaped(2, nShapeF) * cords).reshaped(4, 1);
		//std::cout << mMp.col(i) << std::endl;
		//std::cout << std::endl;
	}
}

void Element_2D_Membrane::computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
	int nGauss = quadrature->mPoints.rows();
	int nShapeF = basis->mNumFunctions;
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(nShapeF, nShapeF);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::Matrix2d J = mMp.col(i).reshaped(2, 2).transpose();
		Eigen::Matrix2d Jinv = J.inverse();
		Eigen::MatrixXd dN = mDerivatives.col(i).reshaped(2, nShapeF).transpose();
		//Eigen::MatrixXd temp = (dN.col(0) * Jinv(0, 0) + dN.col(1) * Jinv(1, 0)) * (dN.col(0).transpose() * Jinv(0, 0) + dN.col(1).transpose() * Jinv(1, 0))
		//	                 + (dN.col(0) * Jinv(0, 1) + dN.col(1) * Jinv(1, 1)) * (dN.col(0).transpose() * Jinv(0, 1) + dN.col(1).transpose() * Jinv(1, 1));
		Eigen::MatrixXd dNdx = dN * Jinv;
		//Ke += w(i) * mat->GetT() * temp * J.determinant();
		Ke += w(i) * dNdx * mat->GetT() * dNdx.transpose() * J.determinant();
	}
	mKe = Ke;
	//std::cout << Ke << std::endl;
}

void Element_2D_Membrane::computeMassMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat)
{
}

void Element_2D_Membrane::computeInternalForce(std::shared_ptr<BasisFunction> basis, int outputNumber, std::shared_ptr<Material> mat)
{
}

double Element_2D_Membrane::getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	return basis->getN(xi).transpose() * mNodalU;
}

double Element_2D_Membrane::getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis)
{
	return 0.0;
}
