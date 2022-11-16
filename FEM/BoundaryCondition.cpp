#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition(int index, BoundaryConditionType type, int BearerIndex, double value)
{
	mIndex = index;
	mType = type;
	mBearerIndex = BearerIndex;
	mValue = value;
}

BoundaryCondition::~BoundaryCondition()
{
}

NodeLoad::NodeLoad(int index, int nodeIndex, double value):BoundaryCondition(index, BoundaryConditionType::LOAD_CONCENTRATE_NODE, nodeIndex, value)
{
}

NodeLoad::~NodeLoad()
{
}

Displacement::Displacement(int index, int nodeIndex, double value):BoundaryCondition(index, BoundaryConditionType::DISPLACEMENT, nodeIndex, value)
{
}

Displacement::~Displacement()
{
}

Spring::Spring(int index, BoundaryConditionType type, int nodeIndex, double k):BoundaryCondition(index, type, nodeIndex, k)
{
	mDofIndex = -1;
}

Spring::~Spring()
{
}

ElementLoad::ElementLoad(int index, int elementIndex, double value):BoundaryCondition(index, BoundaryConditionType::LOAD_CONCENTRATE_ELEMENT, elementIndex, value)
{
}

ElementLoad::~ElementLoad()
{
}

DistributedLoad::DistributedLoad(int index, int elementIndex, double value):BoundaryCondition(index, BoundaryConditionType::LOAD_DISTRIBUTED, elementIndex, value)
{
}

DistributedLoad::~DistributedLoad()
{
}

/*
EvenlyDistributedSpring_1st::EvenlyDistributedSpring_1st(int index, int elementIndex, double k) :EvenlyDistributedSpring(index, elementIndex, k)
{
}

EvenlyDistributedSpring_1st::~EvenlyDistributedSpring_1st()
{
}

Eigen::MatrixXd EvenlyDistributedSpring_1st::getStiffnessMatrix(double length)
{
	Eigen::Matrix2d stiffnessMatrix;
	stiffnessMatrix << mValue * length / 3.0, mValue * length / 6.0, mValue * length / 6.0, mValue * length / 3.0;
	return stiffnessMatrix;
}
*/

LinearlyDistributedLoad::LinearlyDistributedLoad(int index, int elementIndex, double startValue, double endValue) :BoundaryCondition(index, BoundaryConditionType::LOAD_DISTRIBUTED_LINEARLY, elementIndex, startValue)
{
	mEndValue = endValue;
}


LinearlyDistributedLoad::~LinearlyDistributedLoad()
{
}

LinearlyDistributedLoad_1st::LinearlyDistributedLoad_1st(int index, int elementIndex, double startValue, double endValue):LinearlyDistributedLoad(index, elementIndex, startValue, endValue)
{
}

LinearlyDistributedLoad_1st::~LinearlyDistributedLoad_1st()
{
}

Eigen::VectorXd LinearlyDistributedLoad_1st::getElementLoadVector(double length)
{
	Eigen::Vector2d loadVector;
	loadVector << (mValue / 3.0 + mEndValue / 6.0) * length, (mValue / 6.0 + mEndValue / 3.0)* length;
	return loadVector;
}

/*
EvenlyDistributedSpring_2nd::EvenlyDistributedSpring_2nd(int index, int elementIndex, double k) :EvenlyDistributedSpring(index, elementIndex, k)
{
}

EvenlyDistributedSpring_2nd::~EvenlyDistributedSpring_2nd()
{
}

Eigen::MatrixXd EvenlyDistributedSpring_2nd::getStiffnessMatrix(double length)
{
	//Eigen::MatrixXd stiffnessMatrix;
	//stiffnessMatrix.resize(3,3);
	Eigen::Matrix3d stiffnessMatrix;
	stiffnessMatrix << 4.0,  2.0,  -1,
		               2.0, 16.0, 2.0,
		                -1,  2.0, 4.0;
	stiffnessMatrix *= mValue * length / 30.0;
	return stiffnessMatrix;
}
*/

LinearlyDistributedLoad_2nd::LinearlyDistributedLoad_2nd(int index, int elementIndex, double startValue, double endValue) :LinearlyDistributedLoad(index, elementIndex, startValue, endValue)
{
}

LinearlyDistributedLoad_2nd::~LinearlyDistributedLoad_2nd()
{
}

Eigen::VectorXd LinearlyDistributedLoad_2nd::getElementLoadVector(double length)
{
	Eigen::Vector3d loadVector;
	loadVector << mValue, 2.0 * (mValue + mEndValue), mEndValue;
	loadVector *= length / 6.0;
	return loadVector;
}

ConcentrateSpring::ConcentrateSpring(int index, int nodeIndex, double k):Spring(index, BoundaryConditionType::SPRING_CONCENTRATE, nodeIndex, k)
{
}

ConcentrateSpring::~ConcentrateSpring()
{
}

EvenlyDistributedSpring::EvenlyDistributedSpring(int index, int elementIndex, double k) :Spring(index, BoundaryConditionType::SPRING_DISTRIBUTED, elementIndex, k)
{
}

EvenlyDistributedSpring::~EvenlyDistributedSpring()
{
}

Eigen::MatrixXd EvenlyDistributedSpring::getStiffnessMatrix(std::shared_ptr<BasisFunction> basis, std::shared_ptr<Quadrature> quadrature, std::shared_ptr<Material> mat, double elementLength)
{
	int n = basis->mOrder;
	int nGauss = quadrature->mXi.size();
	Eigen::VectorXd w = quadrature->mWeight;
	Eigen::VectorXd N;
	Eigen::MatrixXd Kbc = Eigen::MatrixXd::Zero(n + 1, n + 1);
	for (int i = 0; i < nGauss; i++)
	{
		Eigen::VectorXd xi(1);
		xi << quadrature->mXi[i];
		N = basis->getN(xi);
		Kbc += w(i) * N * N.transpose();
	}
	return Kbc * mValue * elementLength / 2.;
}
