#pragma once

#include <eigen/Eigen/Core>
#include "BasisFunction.h"
#include "Material.h"

class BasisFunction;
class Material;
class SextionAreaFunction;

enum BoundaryConditionType
{
	LOAD_CONCENTRATE_NODE,
	LOAD_CONCENTRATE_ELEMENT,
	LOAD_DISTRIBUTED,
	LOAD_DISTRIBUTED_LINEARLY,
	DISPLACEMENT,
	//SPRING,
	SPRING_DISTRIBUTED,
	SPRING_CONCENTRATE,
};

class BoundaryCondition
{
public:
	BoundaryCondition(int index, BoundaryConditionType type, int mBearerIndex, double value);
	virtual ~BoundaryCondition();

	//virtual Eigen::MatrixXd getStiffnessMatrix(double length) = 0;

	int mIndex;
	BoundaryConditionType mType;
	int mBearerIndex;
	double mValue;
};

class NodeLoad : public BoundaryCondition
{
public:
	NodeLoad(int index, int nodeIndex, double value, int direction );
	~NodeLoad();

	int mDirection;
};

class ElementLoad :public BoundaryCondition
{
public:
	ElementLoad(int index, int elementIndex, double value);
	~ElementLoad();
};

class DistributedLoad :public BoundaryCondition
{
public:
	DistributedLoad(int index, int elementIndex, double value);
	~DistributedLoad();
};

class LinearlyDistributedLoad :public BoundaryCondition
{
public:
	LinearlyDistributedLoad(int index, int elementIndex, double startValue, double endValue);
	~LinearlyDistributedLoad();

	virtual Eigen::VectorXd getElementLoadVector(double length) = 0;

	double mEndValue;
};

class LinearlyDistributedLoad_1st :public LinearlyDistributedLoad
{
public:
	LinearlyDistributedLoad_1st(int index, int elementIndex, double startValue, double endValue);
	~LinearlyDistributedLoad_1st();

	Eigen::VectorXd getElementLoadVector(double length);
};

class LinearlyDistributedLoad_2nd :public LinearlyDistributedLoad
{
public:
	LinearlyDistributedLoad_2nd(int index, int elementIndex, double startValue, double endValue);
	~LinearlyDistributedLoad_2nd();

	Eigen::VectorXd getElementLoadVector(double length);
};

class Displacement : public BoundaryCondition
{
public:
	Displacement(int index, int nodeIndex, double value, int direction);
	~Displacement();
	
	int mDirection;
};

class Spring : public BoundaryCondition
{
public:
	Spring(int index, BoundaryConditionType type, int nodeIndex, double k);
	~Spring();
	

	int mDofIndex;
};
class ConcentrateSpring : public Spring
{
public:
	ConcentrateSpring(int index, int nodeIndex, double k, int direction);
	~ConcentrateSpring();
	
	int mDirection;
};

class EvenlyDistributedSpring : public Spring
{
public:
	EvenlyDistributedSpring(int index, int elementIndex, double k);
	~EvenlyDistributedSpring();
	
	virtual Eigen::MatrixXd getStiffnessMatrix(std::shared_ptr<BasisFunction> basis,
		                                       std::shared_ptr<Quadrature> quadrature,
	                                           double elementLength);
};

class EvenlyDistributedSpring_2D : public EvenlyDistributedSpring
{
public:
	EvenlyDistributedSpring_2D(int index, int elementIndex, double k, int dir, int orentationSign, double valueOfFixDir);
	~EvenlyDistributedSpring_2D();

	virtual Eigen::MatrixXd getStiffnessMatrix(std::shared_ptr<BasisFunction> basis,
		std::shared_ptr<Quadrature> quadrature,
		Eigen::MatrixXd cords);

	int mDir, mOrentationSign;
	double mValueOfFixDir;
};



/*
class EvenlyDistributedSpring_1st : public EvenlyDistributedSpring
{
public:
	EvenlyDistributedSpring_1st(int index, int elementIndex, double k);
	~EvenlyDistributedSpring_1st();

	Eigen::MatrixXd getStiffnessMatrix(double length);
};

class EvenlyDistributedSpring_2nd : public EvenlyDistributedSpring
{
public:
	EvenlyDistributedSpring_2nd(int index, int elementIndex, double k);
	~EvenlyDistributedSpring_2nd();

	Eigen::MatrixXd getStiffnessMatrix(double length);
};
*/


