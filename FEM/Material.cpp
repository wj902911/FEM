#include "Material.h"

Material::Material()
{
	mE = 0;
	mNu = 0;
	mDensity = 0;
	mT = 0;
}

Material::~Material()
{
}

void Material::SetE(double E)
{
	mE = E;
}

void Material::SetNu(double Nu)
{
	mNu = Nu;
}

void Material::SetDensity(double Density)
{
	mDensity = Density;
}

void Material::SetT(double T)
{
	mT = T;
}

double Material::GetE()
{
	return mE;
}

double Material::GetNu()
{
	return mNu;
}

double Material::GetDensity()
{
	return mDensity;
}

double Material::GetT()
{
	return mT;
}
