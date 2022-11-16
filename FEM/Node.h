#pragma once

#include <eigen/Eigen/Core>
#include <unordered_set>

class Node
{
public:
	Node(int index, int dim, Eigen::VectorXd x, bool isFixed, bool hasRotation);
	~Node();

	void fixDof();

	int mIndex, mDim;
	Eigen::VectorXd mX, mU;
	std::unordered_set<int> mElementIndices;
	bool mIsFixed, mHasRotation;
	Eigen::VectorXi mDofIndex;
};

