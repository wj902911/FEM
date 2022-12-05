#pragma once
class Material
{
public:
	Material();
	~Material();
	void SetE(double E);
	void SetNu(double Nu);
	void SetDensity(double Density);
	void SetT(double T);
	double GetE();
	double GetNu();
	double GetDensity();
	double GetT();
	
private:
	double mE;
	double mNu;
	double mDensity;
	double mT;
};

