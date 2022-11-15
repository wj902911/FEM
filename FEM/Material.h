#pragma once
class Material
{
public:
	Material();
	~Material();
	void SetE(double E);
	void SetNu(double Nu);
	void SetDensity(double Density);
	double GetE();
	double GetNu();
	double GetDensity();
	
private:
	double mE;
	double mNu;
	double mDensity;
};

