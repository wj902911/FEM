#pragma once

#include <vector>
#include <tuple>
#include <eigen/Eigen/Sparsecore>
#include <eigen/Eigen/Dense>
#include "Node.h"
#include "Element.h"
#include "Material.h"
#include "BoundaryCondition.h"

using namespace std;
using TripletArray = std::vector<Eigen::Triplet<double>>;

class Core
{
public:
	Core();
	~Core();

	void setOutputNumber(int number);

	void assignDofIndexes();
	void computeStiffMatrix();
	void computeMassMatrix();
	void computeLoadVector();
	void computeNaturalFrequency();

	void solve();
	void solveForFrequency();

	void solveForUnknownU();
	void assignU();
	void computeInternalForce();
	double getU(double x);
	double getdU(double x);
	void setDimension(int dimension);
	
	vector<shared_ptr<Node>> mNodes;
	vector<shared_ptr<Element>> mElements;
	vector<shared_ptr<BasisFunction>> mBasisFunctions;
	vector<shared_ptr<Quadrature>> mQuadratures;
	vector<shared_ptr<Material>> mMaterials;
	vector<shared_ptr<BoundaryCondition>> mBoundaryConditions;
	vector<shared_ptr<BodyForceFunction>> mBodyForceFunctions;
	Eigen::SparseMatrix<double> mK, mM;
	Eigen::VectorXd mF;
	Eigen::VectorXd mUKnown;
	Eigen::VectorXd mUUnknown;
	Eigen::VectorXd mNaturalFrequencies;
	int mNDof_unkonwn, mNDof_konwn, mNDof_all;
	int mOutputNumber;
	int mDim;
};

