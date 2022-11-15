#include <iostream>
#include <fstream>
#include <eigen/Eigen/Dense>
#include "Core.h"

#if 0 //HW3
double computeUex(double x, double p0, double A, double E, double L)
{
	return p0 * E * A / L / L * (x / L - sinh(x / L) / sinh(1.0));
}

double computeDudxex(double x, double p0, double A, double E, double L)
{	
	return p0 * E * A / L / L * (1.0 / L - cosh(x / L) / sinh(1.0));
}
#endif

#if 0 //HW4
double computeUex(double x, double alpha)
{
	return 1. / (1. + alpha) * pow(x, 1. + alpha);
}

double computeDudxex(double x, double alpha)
{
	return pow(x, alpha);
}

double computeSexact(double alpha)
{
	return 1. / (2. * alpha + 1.);
}
#endif

#if 1 //EXAM1
double computeUex(double x, double p0, double A0, double E, double L)
{
	if (x <= L)
		return p0 * L * L / E / A0 * (4. * (x / 2. / L + log(1. - x / 2. / L)) - (4. * log(1. - x / 2. / L) - 1.) * (3. + 4. * log(2.)) / (7. + 4. * log(2.)));
	else
		return p0 * L * L / E / A0 * (-(x / L) * (x / L) + 4. * x / L * (3. + 4. * log(2.)) / (7. + 4. * log(2.)) - 4. * (-3. + 4. * log(2.)) / (7. + 4. * log(2.)));
}

double computeDudxex(double x, double p0, double A0, double E, double L)
{
	if (x <= L)
		return p0 * L * L / E / A0 * (4. * (1. / 2. / L - 1. / (2. * L - x)) + (4. / (2. * L - x)) * (3. + 4. * log(2.)) / (7. + 4. * log(2.)));
	else
		return p0 * L * L / E / A0 * (-2. / L * (x / L) + 4. / L * (3. + 4. * log(2.)) / (7. + 4. * log(2.)));

}
#endif

double Simpson(Eigen::VectorXd as, Eigen::VectorXd bs, Eigen::VectorXd fas, Eigen::VectorXd fab2s, Eigen::VectorXd fbs)
{
	double sum = 0;
	for (int i = 0; i < as.size(); i++)
	{
		sum += (bs[i] - as[i]) / 6.0 * (fas[i] + 4 * fab2s[i] + fbs[i]);
	}
	return sum;
}

int main()
{
	int nElements = 8;
	int elementOrder = 1;
	int outputNumber = 20;
	outputNumber += 1;
	double E = 1.;
	double A0 = 1.;
	double L = 1.;
	double p0 = 1.;
	double k2 = E * A0 / L;
	double k1 = 2. * k2;
	double h = L / nElements;

	Core core;
	core.setOutputNumber(outputNumber);
	BodyForceFunction bodyForce(p0);
	core.mBodyForceFunctions.emplace_back(make_shared<BodyForceFunction>(bodyForce));
	
	Material material;
	material.SetE(E);
	core.mMaterials.emplace_back(make_shared<Material>(material));

	Quadrature quadrature_1D(1, elementOrder);
	core.mQuadratures.emplace_back(make_shared<Quadrature>(quadrature_1D));
	Quadrature quadrature_1D_raiseOrder1(1, elementOrder + 1);
	core.mQuadratures.emplace_back(make_shared<Quadrature>(quadrature_1D_raiseOrder1));
	
	if (elementOrder == 1)
	{
		BasisFunction_1D_Linear basisFunction;
		core.mBasisFunctions.emplace_back(make_shared<BasisFunction_1D_Linear>(basisFunction));
		for (int i = 0; i < 2 * nElements; i++)
		{
			Node node(i, i * h, false);
			core.mNodes.emplace_back(make_shared<Node>(node));
		}
		core.mNodes[0]->fixDof();
		Node node(2 * nElements, 2. * L, true);
		core.mNodes.emplace_back(make_shared<Node>(node));
	}
	else if (elementOrder == 2)
	{
		BasisFunction_1D_Quadratic basisFunction;
		core.mBasisFunctions.emplace_back(make_shared<BasisFunction_1D_Quadratic>(basisFunction));
		for (int i = 0; i < 4 * nElements; i++)
		{
			Node node(i, i* h / 2, false);
			core.mNodes.emplace_back(make_shared<Node>(node));
		}
		core.mNodes[0]->fixDof();
		Node node(4 * nElements, 2. * L, true);
		core.mNodes.emplace_back(make_shared<Node>(node));
	}
	else
	{
		std::cout << "Element order not supported!" << std::endl;
		return 0;
	}
	
	SextionAreaFunction A(A0, L);

	for (int i = 0; i < 2 * nElements; i++)
	{
		if (i < nElements)
		{
			if (elementOrder == 1)
			{
				Element_1D_ununiformSection element(i, 
					                                0, 
					                                make_shared<SextionAreaFunction>(A), 
					                                Eigen::Vector2i(i, i + 1), 
					                                Eigen::Vector2d(core.mNodes[i]->mX, core.mNodes[i + 1]->mX),
				                                    0,
					                                1);
				
				core.mElements.emplace_back(make_shared<Element_1D_ununiformSection>(element));
			}
			else if (elementOrder == 2)
			{
				Element_1D_ununiformSection element(i, 
													0, 
					                                make_shared<SextionAreaFunction>(A),
													Eigen::Vector3i(2 * i, 2 * i + 1, 2 * i + 2),
					                                Eigen::Vector2d(core.mNodes[2 * i]->mX, core.mNodes[2 * i + 2]->mX),
					                                0,
					                                1);
				
				core.mElements.emplace_back(make_shared<Element_1D_ununiformSection>(element));
			}
		}
		else
		{
			if (elementOrder == 1)
			{
				Element_1D element(i, 0, A0 * 0.5, Eigen::Vector2i(i, i + 1), 0, 1);
				core.mElements.emplace_back(make_shared<Element_1D>(element));
			}
			else if (elementOrder == 2)
			{
				Element_1D element(i, 0, A0 * 0.5, Eigen::Vector3i(2 * i, 2 * i + 1, 2 * i + 2), 0, 1);
				core.mElements.emplace_back(make_shared<Element_1D>(element));
			}
		}
		
	}
	
	ConcentrateSpring Spring1(0, 0, k1);
	core.mBoundaryConditions.emplace_back(make_shared<ConcentrateSpring>(Spring1));
	ConcentrateSpring Spring2(1, core.mNodes.size() - 1, k2);
	core.mBoundaryConditions.emplace_back(make_shared<ConcentrateSpring>(Spring2));
	
	core.solve();
	
	cout << core.mK << endl;
	cout << endl;
	cout << core.mF << endl;
	cout << endl;
	//Eigen::VectorXd F = core.mF;
	//F(nElements) += core.mUUnknown(nElements) * k;
	//cout << "Sapp = " << 0.5 * core.mUUnknown.dot(F) << endl;
	//cout << endl;
	
	int nSegments = 40;
	Eigen::VectorXd as(nSegments);
	Eigen::VectorXd bs(nSegments);
	Eigen::VectorXd fas(nSegments);
	Eigen::VectorXd fab2s(nSegments);
	Eigen::VectorXd fbs(nSegments);
	double dx = 2. * L / nSegments;
	for (int i = 0; i < nSegments; i++)
	{
		double a = i * dx;
		as(i) = a;
		double temp = computeUex(a, p0, A0, E, L) - core.getU(a);
		fas(i) = temp * temp;
	}
	for (int i = 0; i < nSegments - 1; i++)
	{
		bs(i) = as(i + 1);
		fbs(i) = fas(i + 1);
	}
	bs(nSegments - 1) = 2. * L;
	double temp = computeUex(2. * L, p0, A0, E, L) - core.mNodes[core.mNodes.size() - 1]->mU;
	fbs(nSegments - 1) = temp * temp;
	for (int i = 0; i < nSegments; i++)
	{
		double ab2 = (as(i) + bs(i)) / 2.0;
		double temp = computeUex(ab2, p0, A0, E, L) - core.getU(ab2);
		fab2s(i) = temp * temp;
	}
	double L2Norm = log(sqrt(Simpson(as, bs, fas, fab2s, fbs)));
	cout << "ln(L2Norm) = " << L2Norm << endl;
	cout << endl;
	
	/*
	for (int i = 0; i < nSegments; i++)
	{
		double temp = computeDudxex(as(i), alpha) - core.getdU(as(i));
		//fas(i) *= k;
		//fas(i) += E * A * temp1 * temp1;
		fas(i) = temp * temp;
	}
	for (int i = 0; i < nSegments - 1; i++)
	{
		fbs(i) = fas(i + 1);
	}
	temp = computeDudxex(L, alpha) - core.getdU(L);
	//fbs(nSegments - 1) *= k;
	//fbs(nSegments - 1) += E * A * temp * temp;
	fbs(nSegments - 1) = temp * temp;
	for (int i = 0; i < nSegments; i++)
	{
		double ab2 = (as(i) + bs(i)) / 2.0;
		double temp = computeDudxex(ab2, alpha) - core.getdU(ab2);
		//fab2s(i) *= k;
		//fab2s(i) += E * A * temp * temp;
		fab2s(i) = temp * temp;
	}
	*/
	
	//double energyNorm = log(sqrt(Simpson(as, bs, fas, fab2s, fbs)));
	/*
	if (alpha > -0.5)
	{
		double energyNorm = log(sqrt(computeSexact(alpha) - core.mUUnknown.dot(F)));
		cout << "ln(energyNorm) = " << energyNorm << endl;
		cout << endl;
	}
	*/
	
#if 1
	std::ofstream uOut("u.txt");
	std::ofstream uexOut("uex.txt");
	std::ofstream fOut("f.txt");
	std::ofstream fexOut("fex.txt");
	std::ofstream xOut("x.txt");
	std::ofstream xexOut("xex.txt");

	int nInElement = outputNumber - 1;
	double step = h / nInElement;

	for (int i = 0; i < 2 * nElements; i++)
	{
		if (elementOrder == 1)
			std::cout << "u" << i << " = " << core.mNodes[i]->mU << std::endl;
		else
		{
			std::cout << "u" << 2 * i << " = " << core.mNodes[2 * i]->mU << std::endl;
			std::cout << "u" << 2 * i + 1 << " = " << core.mNodes[2 * i + 1]->mU << std::endl;
		}
		std::shared_ptr<Element> pElement = core.mElements[i];
		for (int j = 0; j < outputNumber; j++)
		{
			double x = i * h + j * step;
			xOut << x << std::endl;
			double u = pElement->mU[j];
			uOut << u << std::endl;
			double f = pElement->mInternalForce[j];
			fOut << f << std::endl;
		}
	}
	if (elementOrder == 1)
		std::cout << "u" << 2 * nElements << " = " << core.mNodes[2 * nElements]->mU << std::endl;
	else
		std::cout << "u" << 4 * nElements << " = " << core.mNodes[4 * nElements]->mU << std::endl;
	
	uOut.close();
	fOut.close();
	xOut.close();
	
	int n = 100;
	for (int i = 0; i <= n; i++)
	{
		double x = i * 2. * L / n;
		xexOut << x << std::endl;
		//double Lambda = sqrt(k / E / A);
		//double uex = p0 / k * (x / L - sinh(Lambda * x) / sinh(Lambda * L));
		//double uex = p0 * L * L / (4.0 * E * A) * (3.0 - 2.0 * (x / L) * (x / L));
		double uex = computeUex(x, p0, A0, E, L);
		uexOut << uex << std::endl;
		Eigen::VectorXd X(1);
		X << x;
		double fex = E * A.getValue(X) * computeDudxex(x, p0, A0, E, L);
		fexOut << fex << std::endl;
	}
	uexOut.close();
	xexOut.close();
	fexOut.close();
#endif
	
}