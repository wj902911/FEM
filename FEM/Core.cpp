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
		int nDof = mBasisFunctions[pElement->mBasisFunctionIndex]->mNumDofEachNode.sum();
		pElement->setSize(nDof);
		for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
		{
			int Ni = pElement->mEndNodeIndex(i);
			mNodes[Ni]->mElementIndices.insert(pElement->mIndex);
		}
	}
	for (auto& node : mNodes)
	{
		for (int i = 0; i < node->mDofIndex.size(); i++)
		{
			if (!node->mIsFixed[i])
			{
				node->mDofIndex[i] = n;
				for (auto& p : node->mElementIndices)
				{
					std::shared_ptr<Element> pElement = mElements[p];
					Eigen::Index orderOfDof = 0;
					for (int j = 0; j < pElement->mEndNodeIndex.size(); j++)
					{
						if (pElement->mEndNodeIndex(j) == node->mIndex)
						{
							pElement->mDofIndexes(orderOfDof + i) = n;
						}
						else
						{
							orderOfDof += mBasisFunctions[pElement->mBasisFunctionIndex]->mNumDofEachNode[j];
						}
					}
				}
				n++;
			}
		}
		
	}
	for (auto& pBC : mBoundaryConditions)
	{
		if (pBC->mType == BoundaryConditionType::SPRING_CONCENTRATE)
		{
			int nodeIndex = pBC->mBearerIndex;
			std::shared_ptr<Node> pNode = mNodes[nodeIndex];
			//pNode->mDofIndex = n;
			std::shared_ptr<ConcentrateSpring> pSpring = std::dynamic_pointer_cast<ConcentrateSpring>(pBC);
			pSpring->mDofIndex = pNode->mDofIndex(pSpring->mDirection);
			/*
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
			*/
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
			std::shared_ptr<Displacement> pDisp = std::dynamic_pointer_cast<Displacement>(pBC);
			int dir = pDisp->mDirection;
			pNode->mDofIndex(dir) = n + m;
			pNode->mU(dir) = pBC->mValue;
			UKnown.emplace_back(pBC->mValue);
			for (auto& pElIdix : pNode->mElementIndices)
			{
				std::shared_ptr<Element> pElement = mElements[pElIdix];
				Eigen::Index orderOfDof = 0;
				//std::shared_ptr<Element_1D> pElement = std::dynamic_pointer_cast<Element_1D>(mElements[pElIdix]);
				for (int i = 0; i < pElement->mEndNodeIndex.size(); i++)
				{
					if (pElement->mEndNodeIndex(i) == nodeIndex)
						pElement->mDofIndexes(orderOfDof + dir) = n + m;
					else
						orderOfDof += mBasisFunctions[pElement->mBasisFunctionIndex]->mNumDofEachNode[i];
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

	/*
	for (auto& pEle : mElements)
	{
		std::cout << pEle->mDofIndexes << std::endl;
		std::cout << std::endl;
	}
	*/
}

void Core::computeStiffMatrix()
{
	for (auto& pElement : mElements)
	{
		Eigen::MatrixXd cords;
		int nNode = pElement->mEndNodeIndex.size();
		int geoDim = mNodes[pElement->mEndNodeIndex(0)]->mX.size();
		cords.resize(nNode, geoDim);
		for (int i = 0; i < nNode; i++)
		{
			int nodeIndex = pElement->mEndNodeIndex(i);
			cords.row(i) = mNodes[nodeIndex]->mX;
		}
		pElement->evaluate(cords, mBasisFunctions[pElement->mBasisFunctionIndex], mQuadratures[pElement->mQuadratureIndex]);
		pElement->computeStiffnessMatrix(mBasisFunctions[pElement->mBasisFunctionIndex], mQuadratures[pElement->mQuadratureIndex], mMaterials[pElement->mMaterialIndex]);
	}
	TripletArray Kijs;
	for (auto& pElement : mElements)
	{	
		int numDofs = pElement->mDofIndexes.size();
		for (int i = 0; i < numDofs; i++)
		{
			int m = pElement->mDofIndexes(i);
			if (m >= 0)
			{
				for (int j = 0; j < numDofs; j++)
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
			std::shared_ptr<EvenlyDistributedSpring_2D> pSpring = std::dynamic_pointer_cast<EvenlyDistributedSpring_2D>(pBC);
			int elementIndex = pSpring->mBearerIndex;
			std::shared_ptr<Element> pElement = mElements[elementIndex];
			//std::shared_ptr<Element_1D> pElement = std::dynamic_pointer_cast<Element_1D>(mElements[elementIndex]);
			std::shared_ptr<BasisFunction> pBasisFunction = mBasisFunctions[pElement->mBasisFunctionIndex];
			std::shared_ptr<Quadrature> pQuadrature = mQuadratures[pElement->mQuadratureIndex + 1];
			//double eleLength = pElement->mLength;
			Eigen::MatrixXd k = pSpring->getStiffnessMatrix(pBasisFunction, pQuadrature, pElement->mNodalCoordinates);
			int numDofs = pElement->mDofIndexes.size();
			for (int i = 0; i < numDofs; i++)
			{
				int m = pElement->mDofIndexes(i);
				if (m >= 0)
				{
					for (int j = 0; j < numDofs; j++)
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
	mK.resize(mNDof_all, mNDof_all);
	mK.setZero();
	mK.setFromTriplets(Kijs.begin(), Kijs.end());
}

void Core::computeMassMatrix()
{
	for (auto& pElement : mElements)
	{
		pElement->computeMassMatrix(mBasisFunctions[pElement->mBasisFunctionIndex], mQuadratures[pElement->mQuadratureIndex], mMaterials[pElement->mMaterialIndex]);
	}
	TripletArray Mijs;
	for (auto& pElement : mElements)
	{
		int numDofs = pElement->mDofIndexes.size();
		for (int i = 0; i < numDofs; i++)
		{
			int m = pElement->mDofIndexes(i);
			if (m >= 0)
			{
				for (int j = 0; j < numDofs; j++)
				{
					int n = pElement->mDofIndexes(j);
					if (n >= 0)
					{
						Mijs.emplace_back(m, n, pElement->mMe(i, j));
					}
				}
			}
		}
	}
	mM.resize(mNDof_all, mNDof_all);
	mM.setZero();
	mM.setFromTriplets(Mijs.begin(), Mijs.end());
}

void Core::computeLoadVector()
{
	mF.resize(mNDof_all);
	mF.setZero();
	for (auto& pBdForce : mBodyForceFunctions)
	{
		for (auto& pElement : mElements)
		{
			Eigen::VectorXd f = pBdForce->getElementLoadVector(mBasisFunctions[pElement->mBasisFunctionIndex],
				                                               mQuadratures[pElement->mQuadratureIndex],
															   pElement->mNodalCoordinates);
			for (int i = 0; i < pElement->mDofIndexes.size(); i++)
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
			std::shared_ptr<NodeLoad> pLoad = std::dynamic_pointer_cast<NodeLoad>(pBC);
			int i = mNodes[pBC->mBearerIndex]->mDofIndex(pLoad->mDirection);
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
	computeStiffMatrix();
	computeLoadVector();
	solveForUnknownU();
	assignU();
	computeInternalForce();
}

void Core::solveForFrequency()
{
	assignDofIndexes();
	computeStiffMatrix();
	computeMassMatrix();
	computeNaturalFrequency();
}

void Core::computeNaturalFrequency()
{
	Eigen::MatrixXd dmK = Eigen::MatrixXd(mK);
	Eigen::MatrixXd dmM = Eigen::MatrixXd(mM);
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dmK, dmM);
	mNaturalFrequencies = sqrt(es.eigenvalues().array());
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
		for (int i = 0; i < pNode->mDofIndex.size(); i++)
		{
			int dofIndex = pNode->mDofIndex(i);
			if (dofIndex >= 0 && dofIndex < mUUnknown.size())
			{
				pNode->mU(i) = mUUnknown(dofIndex);
			}
		}
	}
	for (auto& pElement : mElements)
	{
		int n = pElement->mEndNodeIndex.size();
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < mNodes[pElement->mEndNodeIndex(i)]->mU.size(); j++)
			{
				pElement->mNodalU(i) = mNodes[pElement->mEndNodeIndex(i)]->mU(j);
			}
		}
	}
}

void Core::computeInternalForce()
{
	for (auto& pElement : mElements)
	{
		//std::shared_ptr<Element_1D> pElement_1D = std::dynamic_pointer_cast<Element_1D>(pElement);
		//Eigen::VectorXd U;
		//int numberOfNodes = pElement->mEndNodeIndex.size();
		//U.resize(numberOfNodes);
		//for (int i = 0; i < numberOfNodes; i++)
		//{
			//U(i)= mNodes[pElement->mEndNodeIndex(i)]->mU;
		//}
		pElement->computeInternalForce(mBasisFunctions[pElement->mBasisFunctionIndex], mOutputNumber, mMaterials[pElement->mMaterialIndex]);
	}
}

double Core::getU(Eigen::VectorXd x)
{
	double u = 0;
	for (auto& pElement : mElements)
	{
		int n = pElement->mEndNodeIndex.size();
		Eigen::VectorXd min = pElement->mNodalCoordinates.colwise().minCoeff();
		Eigen::VectorXd max = pElement->mNodalCoordinates.colwise().maxCoeff();
		if (x(0) >= min(0) && x(0) <= max(0))
		{
			if (x(1) >= min(1) && x(1) <= max(1))
			{
				//double xi = (2.0 * x - xmin - xmax) / pElement->mLength;
				//u = pElement->getU(xi, mBasisFunctions[pElement->mBasisFunctionIndex]);
				break;
			}
		}
	}
	return u;
}

double Core::getdU(Eigen::VectorXd x)
{
	double du = 0;
	for (auto& pElement : mElements)
	{
		int n = pElement->mEndNodeIndex.size();
		double xmin = mNodes[pElement->mEndNodeIndex[0]]->mX(0);
		double xmax = mNodes[pElement->mEndNodeIndex[n - 1]]->mX(0);
		if (x(0) >= xmin && x(0) <= xmax)
		{
			//std::shared_ptr<Element_1D> pElement_1D = std::dynamic_pointer_cast<Element_1D>(pElement);
			//double xi = (2.0 * x - xmin - xmax) / pElement->mLength;
			//du = pElement->getdU(xi, mBasisFunctions[pElement->mBasisFunctionIndex]);
			break;
		}
	}
	return du;
}

void Core::setDimension(int dimension)
{
	//mDim = dimension;
}
