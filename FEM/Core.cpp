#pragma once

#include "Core.h"
#include <iostream>

Core::Core()
{
	//mBodyForce = 0;
	mOutputNumber = 10;
	mNDof_unkonwn = 0;
	mNDof_konwn = 0;
	mNDof_all = 0;
}

Core::~Core()
{
}

void Core::setOutputNumber(int number)
{
	mOutputNumber = number;
}

void Core::assignDofIndexes()
{
	int n = 0;
	for (auto& pElement : mElements)
	{
		int nNodes = pElement->mEndNodeIndex.size();
		pElement->setSize(nNodes);
		for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
		{
			int Ni = pElement->mEndNodeIndex(i);
			mNodes[Ni]->mElementIndices.insert(pElement->mIndex);
		}
	}
	for (auto& node : mNodes)
	{
		if (!node->mIsFixed)
		{
			node->mDofIndex = n;
			for (auto& p : node->mElementIndices)
			{
				std::shared_ptr<Element> pElement = mElements[p];
				for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
				{
					if (pElement->mEndNodeIndex(i) == node->mIndex)
					{
						pElement->mDofIndexes(i) = n;
					}
				}
			}
			n++;
		}
	}
	for (auto& pBC : mBoundaryConditions)
	{
		if (pBC->mType == BoundaryConditionType::SPRING_CONCENTRATE)
		{
			int nodeIndex = pBC->mBearerIndex;
			std::shared_ptr<Node> pNode = mNodes[nodeIndex];
			pNode->mDofIndex = n;
			std::shared_ptr<ConcentrateSpring> pSpring = std::dynamic_pointer_cast<ConcentrateSpring>(pBC);
			pSpring->mDofIndex = n;
			for (auto& pElIdix : pNode->mElementIndices)
			{
				std::shared_ptr<Element> pElement = mElements[pElIdix];
				if (pElement->mEndNodeIndex(0) == nodeIndex)
				{
					pElement->mDofIndexes(0) = n;
				}
				else
				{
					pElement->mDofIndexes(1) = n;
				}
			}
			n++;
		}
	}
	int m = 0;
	vector<double> UKnown;
	for (auto& pBC : mBoundaryConditions)
	{
		if (pBC->mType == BoundaryConditionType::DISPLACEMENT)
		{
			int nodeIndex = pBC->mBearerIndex;
			std::shared_ptr<Node> pNode = mNodes[nodeIndex];
			pNode->mDofIndex = n + m;
			pNode->mU = pBC->mValue;
			UKnown.emplace_back(pBC->mValue);
			for (auto& pElIdix : pNode->mElementIndices)
			{
				std::shared_ptr<Element> pElement = mElements[pElIdix];
				//std::shared_ptr<Element_1D> pElement = std::dynamic_pointer_cast<Element_1D>(mElements[pElIdix]);
				if (pElement->mEndNodeIndex(0) == nodeIndex)
				{
					pElement->mDofIndexes(0) = n + m;
				}
				else
				{
					pElement->mDofIndexes(1) = n + m;
				}
			}
			m++;
		}
	}
	
	mUKnown.resize(UKnown.size());
	for (int i = 0; i < UKnown.size(); i++)
	{
		mUKnown(i) = UKnown[i];
	}
	mNDof_unkonwn = n;
	mNDof_konwn = m;
	mNDof_all = n + m;

	//for (auto& pEle : mElements)
	//{
	//	std::cout << pEle->mDofIndexes << std::endl;
	//	std::cout << std::endl;
	//}
}

TripletArray Core::computeStiffMatrix()
{
	for (auto& pElement : mElements)
	{
		int N1 = pElement->mEndNodeIndex(0);
		int N3 = pElement->mEndNodeIndex(pElement->mEndNodeIndex.size() - 1);

		double h = abs(mNodes[N3]->mX - mNodes[N1]->mX);
		//std::shared_ptr<Element_1D> pElement_1D = std::dynamic_pointer_cast<Element_1D>(pElement);
		pElement->mLength = h;
		
		pElement->computeStiffnessMatrix(mBasisFunctions[pElement->mBasisFunctionIndex], mQuadratures[pElement->mQuadratureIndex], mMaterials[pElement->mMaterialIndex]);
	}
	TripletArray Kijs;
	for (auto& pElement : mElements)
	{	
		for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
		{
			int m = pElement->mDofIndexes(i);
			if (m >= 0)
			{
				for (int j = 0; j < pElement->mEndNodeIndex.size(); j++)
				{
					int n = pElement->mDofIndexes(j);
					if (n >= 0)
					{
						Kijs.emplace_back(m, n, pElement->mKe(i, j));
					}
				}
			}
		}
	}
	for (auto& pBC : mBoundaryConditions)
	{
		if (pBC->mType == BoundaryConditionType::SPRING_CONCENTRATE)
		{
			std::shared_ptr<Spring> pSpring = std::dynamic_pointer_cast<Spring>(pBC);
			int i = pSpring->mDofIndex;
			Kijs.emplace_back(i, i, pSpring->mValue);
		}
		else if (pBC->mType == BoundaryConditionType::SPRING_DISTRIBUTED)
		{
			std::shared_ptr<EvenlyDistributedSpring> pSpring = std::dynamic_pointer_cast<EvenlyDistributedSpring>(pBC);
			int elementIndex = pSpring->mBearerIndex;
			std::shared_ptr<Element> pElement = mElements[elementIndex];
			//std::shared_ptr<Element_1D> pElement = std::dynamic_pointer_cast<Element_1D>(mElements[elementIndex]);
			Eigen::MatrixXd k = pSpring->getStiffnessMatrix(pElement->mLength);
			for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
			{
				int m = pElement->mDofIndexes(i);
				if (m >= 0)
				{
					for (int j = 0; j < pElement->mEndNodeIndex.size(); j++)
					{
						int n = pElement->mDofIndexes(j);
						if (n >= 0)
						{
							Kijs.emplace_back(m, n, k(i, j));
						}
					}
				}
			}
		}
		//other BCs
	}
	return Kijs;
}

void Core::computeLoadVector()
{
	mF.resize(mNDof_all);
	mF.setZero();
	for (auto& pBdForce : mBodyForceFunctions)
	{
		for (auto& pElement : mElements)
		{
			Eigen::Vector2d nodalCoordinates(mNodes[pElement->mEndNodeIndex[0]]->mX, 
				                             mNodes[pElement->mEndNodeIndex[pElement->mEndNodeIndex.size() - 1]]->mX);
			Eigen::VectorXd f = pBdForce->getElementLoadVector(mBasisFunctions[pElement->mBasisFunctionIndex],
				                                               mQuadratures[pElement->mQuadratureIndex],
				                                               nodalCoordinates);
			for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
			{
				int m = pElement->mDofIndexes(i);
				if (m >= 0)
				{
					mF(m) += f(i);
				}
			}
			
		}
	}
	
	for (auto& pBC : mBoundaryConditions)
	{
		if (pBC->mType == BoundaryConditionType::LOAD_CONCENTRATE_NODE)
		{
			int i = mNodes[pBC->mBearerIndex]->mDofIndex;
			mF(i) += pBC->mValue;
		}
		else if (pBC->mType == BoundaryConditionType::LOAD_DISTRIBUTED_LINEARLY)
		{
			std::shared_ptr<LinearlyDistributedLoad> pLoad = std::dynamic_pointer_cast<LinearlyDistributedLoad>(pBC);
			int elementIndex = pLoad->mBearerIndex;
			//std::shared_ptr<Element_1D> pElement = std::dynamic_pointer_cast<Element_1D>(mElements[elementIndex]);
			std::shared_ptr<Element> pElement = mElements[elementIndex];
			Eigen::VectorXd f = pLoad->getElementLoadVector(pElement->mLength);
			for (int i = 0; i < f.size(); i++)
			{
				int m = pElement->mDofIndexes(i);
				if (m >= 0)
				{
					mF(m) += f(i);
				}
			}
		}
		//other load
	}
}

void Core::solve()
{
	assignDofIndexes();
	TripletArray Kijs = computeStiffMatrix();
	mK.resize(mNDof_all, mNDof_all);
	mK.setZero();
	mK.setFromTriplets(Kijs.begin(), Kijs.end());
	computeLoadVector();
	solveForUnknownU();
	assignU();
	computeInternalForce();
}

void Core::solveForUnknownU()
{
	Eigen::MatrixXd K00 = mK.topLeftCorner(mNDof_unkonwn, mNDof_unkonwn);
	Eigen::MatrixXd K01 = mK.topRightCorner(mNDof_unkonwn, mNDof_konwn);
	Eigen::MatrixXd K10 = mK.bottomLeftCorner(mNDof_konwn, mNDof_unkonwn);
	Eigen::MatrixXd K11 = mK.bottomRightCorner(mNDof_konwn, mNDof_konwn);
	Eigen::VectorXd F0 = mF.head(mNDof_unkonwn);
	Eigen::VectorXd F1 = mF.tail(mNDof_konwn);
	Eigen::VectorXd b = F0 - K01 * mUKnown;
	mUUnknown = K00.colPivHouseholderQr().solve(b);
}

void Core::assignU()
{
	for (auto& pNode : mNodes)
	{
		if (pNode->mDofIndex >= 0)
		{
			pNode->mU = mUUnknown(pNode->mDofIndex);
		}
	}
	for (auto& pElement : mElements)
	{
		int n = pElement->mEndNodeIndex.size();
		pElement->mNodalU.resize(n);
		for (int i = 0; i < n; i++)
		{
			pElement->mNodalU(i) = mNodes[pElement->mEndNodeIndex(i)]->mU;
		}
	}
}

void Core::computeInternalForce()
{
	for (auto& pElement : mElements)
	{
		//std::shared_ptr<Element_1D> pElement_1D = std::dynamic_pointer_cast<Element_1D>(pElement);
		Eigen::VectorXd U;
		int numberOfNodes = pElement->mEndNodeIndex.size();
		U.resize(numberOfNodes);
		for (int i = 0; i < numberOfNodes; i++)
		{
			U(i)= mNodes[pElement->mEndNodeIndex(i)]->mU;
		}
		pElement->computeInternalForce(mBasisFunctions[pElement->mBasisFunctionIndex], mOutputNumber, mMaterials[pElement->mMaterialIndex]);
	}
}

double Core::getU(double x)
{
	double u = 0;
	for (auto& pElement : mElements)
	{
		int n = pElement->mEndNodeIndex.size();
		double xmin = mNodes[pElement->mEndNodeIndex[0]]->mX;
		double xmax = mNodes[pElement->mEndNodeIndex[n - 1]]->mX;
		if (x >= xmin && x <= xmax)
		{
			//std::shared_ptr<Element_1D> pElement_1D = std::dynamic_pointer_cast<Element_1D>(pElement);
			double xi = (2.0 * x - xmin - xmax) / pElement->mLength;
			u = pElement->getU(xi, mBasisFunctions[pElement->mBasisFunctionIndex]);
			break;
		}
	}
	return u;
}

double Core::getdU(double x)
{
	double du = 0;
	for (auto& pElement : mElements)
	{
		int n = pElement->mEndNodeIndex.size();
		double xmin = mNodes[pElement->mEndNodeIndex[0]]->mX;
		double xmax = mNodes[pElement->mEndNodeIndex[n - 1]]->mX;
		if (x >= xmin && x <= xmax)
		{
			//std::shared_ptr<Element_1D> pElement_1D = std::dynamic_pointer_cast<Element_1D>(pElement);
			double xi = (2.0 * x - xmin - xmax) / pElement->mLength;
			du = pElement->getdU(xi, mBasisFunctions[pElement->mBasisFunctionIndex]);
			break;
		}
	}
	return du;
}

void Core::setDimension(int dimension)
{
	mDim = dimension;
}
