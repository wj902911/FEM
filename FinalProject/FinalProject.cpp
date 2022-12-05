#include <iostream>
#include <fstream>
#include <eigen/Eigen/Dense>
#include "Core.h"

using namespace std;

int main()
{
	int nElementsOnEachDir = 8;
	double T = 1.;
	double L = 1.;
	double p0 = 1.;
	int outputNumber = 200;
	
	double h = L / nElementsOnEachDir;
	double k = T / L;

	Core core;
	Eigen::VectorXd p(2);
	p << p0, L;
	BodyForceFunction bodyForce(p);
	core.mBodyForceFunctions.emplace_back(make_shared<BodyForceFunction>(bodyForce));

	Material material;
	material.SetT(T);
	core.mMaterials.emplace_back(make_shared<Material>(material));
	
	Quadrature quadrature_2D(2, 1);
	core.mQuadratures.emplace_back(make_shared<Quadrature>(quadrature_2D));
	Quadrature quadrature_1D(1, 1);
	core.mQuadratures.emplace_back(make_shared<Quadrature>(quadrature_1D));

	BasisFunction_2D_Quad_Linear basisFunction;
	basisFunction.SetNumDofEachNode(1);
	core.mBasisFunctions.emplace_back(make_shared<BasisFunction_2D_Quad_Linear>(basisFunction));

	int dispIndex = 0;
	for (int i = 0; i <= nElementsOnEachDir; i++)
	{
		for (int j = 0; j <= nElementsOnEachDir; j++)
		{
			int nodeIndex = i * (nElementsOnEachDir + 1) + j;
			Eigen::Vector2d coord(j * h, i * h);
			Eigen::VectorXi isFixed(1);
			if (i == 0 || i == nElementsOnEachDir || j == 0)
			{
				isFixed << 1;
			}
			else
			{
				isFixed << 0;
			}
			if (j == 0 && i != 0 && i != nElementsOnEachDir)
			{
				double disp = 2.0 * p0 * L * L / T * coord(1) / L * (1 - coord(1) / L);
				Displacement displacement(dispIndex, nodeIndex, disp, 0);
				core.mBoundaryConditions.emplace_back(make_shared<Displacement>(displacement));
				dispIndex++;
			}
			Node node(nodeIndex, 1, coord, isFixed, 0);
			core.mNodes.emplace_back(make_shared<Node>(node));
		}
	}
	int springIndex = 0;
	for (int i = 0; i < nElementsOnEachDir; i++)
	{
		for (int j = 0; j < nElementsOnEachDir; j++)
		{
			int eleIndex = i * nElementsOnEachDir + j;
			Eigen::VectorXi nodeIds(4);
			nodeIds << i * (nElementsOnEachDir + 1) + j,
				i* (nElementsOnEachDir + 1) + j + 1,
				(i + 1)* (nElementsOnEachDir + 1) + j + 1,
				(i + 1)* (nElementsOnEachDir + 1) + j;
			Element_2D_Membrane element(eleIndex, 0, nodeIds);
			core.mElements.emplace_back(make_shared<Element_2D_Membrane>(element));
			if (j == nElementsOnEachDir - 1)
			{
				EvenlyDistributedSpring_2D Spring(springIndex, eleIndex, k, 1, 1, 1);
				core.mBoundaryConditions.emplace_back(make_shared<EvenlyDistributedSpring_2D>(Spring));
				springIndex++;
			}
		}
	}

	core.solve();

	int nNode = core.mNodes.size();
	Eigen::MatrixXd coords(nNode, 2);
	for (int i = 0; i < nNode; i++)
	{
		cout << core.mNodes[i]->mU << endl;
		coords.row(i) = core.mNodes[i]->mX;
	}

	Eigen::VectorXd uex = p0 * L * L / T * (2. - coords.col(0).array() / L) * coords.col(1).array() / L * (1. - coords.col(1).array() / L);
	cout << uex << endl;

	double s = 2. / outputNumber;
	double norm2 = 0;
	for (int i = 0; i < core.mElements.size(); i++)
	{
		for (int j = 0; j < outputNumber; j++)
		{
			for (int k = 0; k < outputNumber; k++)
			{
				shared_ptr<Element> pEle = core.mElements[i];
				shared_ptr<BasisFunction> pBasis = core.mBasisFunctions[pEle->mBasisFunctionIndex];
				Eigen::MatrixXd coords = pEle->mNodalCoordinates;
				Eigen::Vector2d xi_eta(-1 + k * s, -1 + j * s);
				double u = core.mElements[i]->getU(xi_eta, pBasis);
				Eigen::Vector2d xy = pBasis->getN(xi_eta).transpose() * coords;
				double ue = p0 * L * L / T * (2. - xy(0) / L) * xy(1) / L * (1. - xy(1) / L);
				//cout << ue-u << endl;
				norm2 += (ue - u) * (ue - u) * s * s;
			}
		}
	}

	cout << sqrt(norm2) << endl;
}