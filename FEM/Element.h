#pragma once

#include <eigen/Eigen/Core>
#include <vector>
#include "BasisFunction.h"
#include "Material.h"

class BasisFunction;
class Material;
class SextionAreaFunction;
//class Node;

class Element
{
public:
	Element(int index, int materailIndex, int basisFunctionIndex, int quadratureIndex);
	virtual ~Element();
	
	virtual void evaluate(Eigen::MatrixXd cords,
		std::shared_ptr<BasisFunction> basis,
		std::shared_ptr<Quadrature> quadrature) = 0;
	
	//virtual void getJacobian();
	virtual void computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, 
		                                std::shared_ptr<Quadrature> quadrature, 
										std::shared_ptr<Material> mat) = 0;

	virtual void computeMassMatrix(std::shared_ptr<BasisFunction> basis,
		                                std::shared_ptr<Quadrature> quadrature,
		                                std::shared_ptr<Material> mat) = 0;
	
	virtual void computeInternalForce(std::shared_ptr<BasisFunction> basis, 
		                              int outputNumber, 
		                              std::shared_ptr<Material> mat) = 0;
	
	virtual double getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis) = 0;
	virtual double getdU(Eigen::VectorXd, std::shared_ptr<BasisFunction> basis) = 0;

	void setSize(int size);

	int mIndex;
	//double mE;
	int mMaterialIndex;
	int mBasisFunctionIndex;
	int mQuadratureIndex;

	Eigen::VectorXi mEndNodeIndex;
	Eigen::MatrixXd mKe, mMe, mNodalCoordinates;
	Eigen::VectorXi mDofIndexes;
	Eigen::VectorXd mNodalU;

	Eigen::VectorXd mInternalForce;
	Eigen::VectorXd mStrain;
	Eigen::VectorXd mStress;
	Eigen::VectorXd mU;

	double mLength;
};

class Element_1D :public Element
{
public:
	Element_1D(int index, 
		       int materailIndex, 
			   double A, 
		       Eigen::VectorXi endNodeIndex, 
		       int basisFunctionIndex = 0, 
		       int quadratureIndex = 0);
	
	~Element_1D();
	
	void evaluate(Eigen::MatrixXd cords,
		std::shared_ptr<BasisFunction> basis,
		std::shared_ptr<Quadrature> quadrature);
	
	void computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, 
		                        std::shared_ptr<Quadrature> quadrature, 
								std::shared_ptr<Material> mat);
	
	void computeMassMatrix(std::shared_ptr<BasisFunction> basis,
		                   std::shared_ptr<Quadrature> quadrature,
		                   std::shared_ptr<Material> mat);
	
	void computeInternalForce(std::shared_ptr<BasisFunction> basis, 
							  int outputNumber, 
							  std::shared_ptr<Material> mat);
	
	double getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);
	double getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);
	
	double mSectionArea;
	//double mLength;
};

class Element_1D_ununiformSection :public Element
{
public:
	Element_1D_ununiformSection(int index, 
		                        int materailIndex, 
								std::shared_ptr<SextionAreaFunction> A, 
								Eigen::VectorXi endNodeIndex, 
								Eigen::VectorXd nodalCoordinates, 
								int basisFunctionIndex = 0, 
								int quadratureIndex = 0);
	
	~Element_1D_ununiformSection();

	void evaluate(Eigen::MatrixXd cords,
		std::shared_ptr<BasisFunction> basis,
		std::shared_ptr<Quadrature> quadrature);
	
	void computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis, 
								std::shared_ptr<Quadrature> quadrature, 
								std::shared_ptr<Material> mat);

	void computeMassMatrix(std::shared_ptr<BasisFunction> basis,
		                   std::shared_ptr<Quadrature> quadrature,
		                   std::shared_ptr<Material> mat);
	
	void computeInternalForce(std::shared_ptr<BasisFunction> basis, 
							  int outputNumber, 
							  std::shared_ptr<Material> mat);
	
	double getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);
	double getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);

	std::shared_ptr<SextionAreaFunction> mSectionArea;
	//double mLength;
	Eigen::VectorXd mNodalCoordinates;
};

class Element_1D_Beam :public Element
{
public:
	Element_1D_Beam(int index,
		            int materailIndex,
		            double A,
		            double I,
		            Eigen::VectorXi endNodeIndex,
		            int basisFunctionIndex = 0,
		            int quadratureIndex = 0);

	~Element_1D_Beam();

	void evaluate(Eigen::MatrixXd cords, 
		std::shared_ptr<BasisFunction> basis, 
		std::shared_ptr<Quadrature> quadrature);
	
	void computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis,
		                        std::shared_ptr<Quadrature> quadrature,
		                        std::shared_ptr<Material> mat);

	void computeMassMatrix(std::shared_ptr<BasisFunction> basis,
		                   std::shared_ptr<Quadrature> quadrature,
		                   std::shared_ptr<Material> mat);

	void computeInternalForce(std::shared_ptr<BasisFunction> basis,
		                      int outputNumber,
		                      std::shared_ptr<Material> mat);

	double getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);
	double getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);

	double mSectionArea, mI;
};

class Element_2D_Membrane :public Element
{
public:
	Element_2D_Membrane(int index,
		int materailIndex,
		Eigen::VectorXi endNodeIndex,
		int basisFunctionIndex = 0,
		int quadratureIndex = 0);

	~Element_2D_Membrane();

	void evaluate(Eigen::MatrixXd cords, 
		std::shared_ptr<BasisFunction> basis, 
		std::shared_ptr<Quadrature> quadrature);

	void computeStiffnessMatrix(std::shared_ptr<BasisFunction> basis,
		std::shared_ptr<Quadrature> quadrature,
		std::shared_ptr<Material> mat);

	void computeMassMatrix(std::shared_ptr<BasisFunction> basis,
		std::shared_ptr<Quadrature> quadrature,
		std::shared_ptr<Material> mat);

	void computeInternalForce(std::shared_ptr<BasisFunction> basis,
		int outputNumber,
		std::shared_ptr<Material> mat);

	double getU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);
	double getdU(Eigen::VectorXd xi, std::shared_ptr<BasisFunction> basis);

	Eigen::MatrixXd mValues, mDerivatives, mMp;
};