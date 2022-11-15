#pragma once

#include <eigen/Eigen/Core>
#include <unordered_set>

class Node
{
public:
	Node(int index, double x, bool isFixed);
	~Node();

	void fixDof();

	int mIndex;
	double mX;
	std::unordered_set<int> mElementIndices;
	bool mIsFixed;
	int mDofIndex;
	double mU;
};

