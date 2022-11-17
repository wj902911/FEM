#include "Node.h"

Node::Node(int index, int dim, Eigen::VectorXd x, Eigen::VectorXi isFixed, bool hasRotation)
{
	mIndex = index;
	mDim = dim;
	mX = x;
	mIsFixed = isFixed;
	mHasRotation = hasRotation;
	int numDof = mDim;
	if (hasRotation)
	{
		if (mDim < 3)
			numDof += 1;
		else
			numDof += 3;
	}
	mDofIndex = Eigen::VectorXi::Constant(numDof, -1);
	mU.setZero(numDof);
}

Node::~Node()
{
}

void Node::fixDof()
{
	mIsFixed = Eigen::VectorXi::Constant(mIsFixed.size(), 1);
}
