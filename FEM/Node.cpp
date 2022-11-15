#include "Node.h"

Node::Node(int index, double x, bool isFixed)
{
	mIndex = index;
	mX = x;
	mIsFixed = isFixed;
	mDofIndex = -1;
	mU = 0;
}

Node::~Node()
{
}

void Node::fixDof()
{
	mIsFixed = true;
}
