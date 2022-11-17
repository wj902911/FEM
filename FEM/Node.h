#pragma once

#include <eigen/Eigen/Core>
#include <unordered_set>

class Node
{
public:
	Node(int index, int dim, Eigen::VectorXd x, Eigen::VectorXi isFixed, bool hasRotation);
	~Node();

	void fixDof();

	int mIndex, mDim;
	Eigen::VectorXd mX, mU;
	std::unordered_set<int> mElementIndices;
	bool mHasRotation;
	Eigen::VectorXi mIsFixed, mDofIndex;
};

