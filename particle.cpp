#include <cmath>
#include <iostream>
#include <cassert>
#include "particle.hpp"
#include <numeric>
#include <gsl/gsl_util>

float q_e = -1.602176487e-19;
float c   = 299792458;
float m_e = 9.10938215e-31;
float PI  = 3.141592654;
//constructor

Set::Set(int size)
{
	assert(size > 0);
	mSize = size;
	for (int i = 0; i < mSize; ++i) mData.push_back(Particle());
}

Set& Set::operator=(const Set& s)
{
	assert(s.mSize == mSize);
	for (int i = 0; i < mSize; ++i) mData[i] = s.mData[i];
	return *this;
}
//overloading square brackets, zerobased indexing

Particle& Set::operator[](int i)
{
	assert(i > -1);
	assert(i < mSize);
	return mData[i];
}

//overloading curved brackets, one based indexing

Particle& Set::operator()(int i)
{
	assert(i > 0);
	assert(i < mSize + 1);
	return mData[static_cast<_int64>(i) - 1];
};

Set Set::operator+(const Set& s1) const
{
	assert(mSize == s1.mSize);
	Set s(mSize);
	for (int i = 0; i < mSize; ++i) {
		s[i].x = mData[i].x + s1.mData[i].x;
		s[i].GB_x = mData[i].GB_x + s1.mData[i].GB_x;
		s[i].y = mData[i].y + s1.mData[i].y;
		s[i].GB_y = mData[i].GB_y + s1.mData[i].GB_y;
		s[i].z = mData[i].z + s1.mData[i].z;
		s[i].GB_z = mData[i].GB_z + s1.mData[i].GB_z;
	}
	return s;
}

Set Set::operator-(const Set& s1) const
{
	assert(mSize == s1.mSize);
	Set s(mSize);
	for (int i = 0; i < mSize; ++i) {
		s[i].x = mData[i].x - s1.mData[i].x;
		s[i].GB_x = mData[i].GB_x - s1.mData[i].GB_x;
		s[i].y = mData[i].y - s1.mData[i].y;
		s[i].GB_y = mData[i].GB_y - s1.mData[i].GB_y;
		s[i].z = mData[i].z - s1.mData[i].z;
		s[i].GB_z = mData[i].GB_z - s1.mData[i].GB_z;
	}
	return s;
}
void Set::add(const Set& s1)
{
	assert(mSize == s1.mSize);
	for (int i = 0; i < mSize; ++i) {
		mData[i].x += s1.mData[i].x;
		mData[i].GB_x += s1.mData[i].GB_x;
		mData[i].y += s1.mData[i].y;
		mData[i].GB_y += s1.mData[i].GB_y;
		mData[i].z += s1.mData[i].z;
		mData[i].GB_z += s1.mData[i].GB_z;
	}
}
void Set::subtract(const Set& s1)
{
	assert(mSize == s1.mSize);
	for (int i = 0; i < mSize; ++i) {
		mData[i].x -= s1.mData[i].x;
		mData[i].GB_x -= s1.mData[i].GB_x;
		mData[i].y -= s1.mData[i].y;
		mData[i].GB_y -= s1.mData[i].GB_y;
		mData[i].z -= s1.mData[i].z;
		mData[i].GB_z -= s1.mData[i].GB_z;
	}
}

void Set::SetScale(float a1, float a2, float a3, float b1, float b2, float b3)
{
	for (int i = 0; i < mSize; ++i) {
		mData[i].x *= a1;
		mData[i].GB_x *= b1;
		mData[i].y *= a2;
		mData[i].GB_y *= b2;
		mData[i].z *= a3;
		mData[i].GB_z *= b3;
	}
}

void Set::SetRXYPhidist(float mu, float sigma, int sout)
{
	unsigned seed1 = (unsigned int) std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine e(seed1);
	std::normal_distribution<float> distN(mu, sigma);

	unsigned seed2 = (unsigned int) std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine ee(seed2);
	std::uniform_real_distribution<float> distU(0.0, 1.0);
	
	for (int i = 0; i < mSize; ++i) {
		float r = fabs(distN(e));
		float phi = 2.0*PI*fabs(distU(ee));
		if (r > mu + sout * sigma) r -= sout * sigma;
		if (r + sout * sigma < mu) r += sout * sigma;
		mData[i].x = r * cos(phi);
		mData[i].y = r * sin(phi);
	}
}

void Set::SetZdist(float mu, float sigma, float sleft, float sright)
{
	unsigned seed = (unsigned int) std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine e(seed);
	std::normal_distribution<float> distN(mu, sigma);
	for (int i = 0; i < mSize; ++i) {
		float z = distN(e);
		if (z > mu + sright * sigma) z -= sright * sigma;
		if (z + sleft * sigma < mu) z += sleft * sigma;
		mData[i].z = z;
	}
}

void Set::SetGBRXYPhidist(float mu, float sigma, int sout)
{
	unsigned seed1 = (unsigned int) std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine e(seed1);
	std::normal_distribution<float> distN(mu, sigma);

	unsigned seed2 = (unsigned int) std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine ee(seed2);
	std::uniform_real_distribution<float> distU(0.0, 1.0);

	for (int i = 0; i < mSize; ++i) {
		float GBr = fabs(distN(e));
		float phi = 2.0 * PI * fabs(distU(ee));
		if (GBr > mu + sout * sigma) GBr -= sout * sigma;
		if (GBr + sout * sigma < mu) GBr += sout * sigma;

		mData[i].GB_x = GBr * cos(phi);
		mData[i].GB_y = GBr * sin(phi);
	}
}

void Set::SetGBZdist(float center, float width)
{
	unsigned seed = (unsigned int) std::chrono::steady_clock::now().time_since_epoch().count();
	std::default_random_engine e(seed);
	std::uniform_real_distribution<float> distU(center-width/2, center+width/2);
	for (int i = 0; i < mSize; ++i) {
		mData[i].GB_z = distU(e);
	}
}

void Set::AddXdiv(float center, float div)
{
	for (int i = 0; i < mSize; ++i) {
		float Gold = sqrt(1 + mData[i].GB_x * mData[i].GB_x + mData[i].GB_y * mData[i].GB_y + mData[i].GB_z * mData[i].GB_z);
		mData[i].GB_x += div * (mData[i].x - center);
		float GBz2 = Gold * Gold - mData[i].GB_x * mData[i].GB_x - mData[i].GB_y * mData[i].GB_y - 1;
		assert(GBz2 > 0);
		mData[i].GB_z = sqrt(GBz2);
	}
}

void Set::AddYdiv(float center, float div)
{
	for (int i = 0; i < mSize; ++i) {
		float Gold = sqrt(1 + mData[i].GB_x * mData[i].GB_x + mData[i].GB_y * mData[i].GB_y + mData[i].GB_z * mData[i].GB_z);
		mData[i].GB_y += div * (mData[i].y - center);
		float GBz2 = Gold * Gold - mData[i].GB_x * mData[i].GB_x - mData[i].GB_y * mData[i].GB_y - 1;
		assert(GBz2 > 0);
		mData[i].GB_z = sqrt(GBz2);
	}
}

void Set::AddZdiv(float center, float div)
{
	for (int i = 0; i < mSize; ++i) {
		float Gset = sqrt(1 + mData[i].GB_x * mData[i].GB_x + mData[i].GB_y * mData[i].GB_y + mData[i].GB_z * mData[i].GB_z);
		Gset += div * (mData[i].z - center) / (m_e * c * c / -q_e);
		assert(Gset > 1);
		float GBz2 = Gset * Gset - mData[i].GB_x * mData[i].GB_x - mData[i].GB_y * mData[i].GB_y - 1;
		mData[i].GB_z = sqrt(GBz2);
	}
}

void Set::SetGBXEmittance(float nemix)
{
	float xx, xpxp, xxp;
	//xc = (float*)malloc(mSize * sizeof(float));
	//xpc = (float*)malloc(mSize * sizeof(float));

	float* xc{ new float[mSize] {} };
	float* xpc{ new float[mSize] {} };
	//auto xc = std::make_unique<float[]>(mSize);
	//auto xpc = std::make_unique<float[]>(mSize);

	for (int i = 0; i < mSize; ++i) {
		xpc[i] = mData[i].GB_x;
		xc[i] = mData[i].x;
	}
	float mean_xpc = std::accumulate(xpc, xpc+mSize, 0.0) / mSize;
	float mean_xc = std::accumulate(xc, xc+mSize, 0.0) / mSize;

	std::transform(xpc, xpc + mSize, xpc, [mean_xpc](float x) {return x - mean_xpc; });
	std::transform(xc, xc + mSize, xc, [mean_xc](float x) {return x - mean_xc; });

	xx   = std::inner_product(xc, xc+mSize, xc, 0.0) / mSize;
	xpxp = std::inner_product(xpc, xpc+mSize, xpc, 0.0) / mSize;
	xxp  = std::inner_product(xc,xc+mSize, xpc, 0.0) / mSize;

	float oldnemix = std::sqrt(xx * xpxp - xxp * xxp);

	delete[] xc;
	delete[] xpc;

	assert(oldnemix > 0);
	float scale = nemix / oldnemix;
	int mags = (int)fabs(log10(scale));
	for (int i = 0; i < mSize; ++i) mData[i].GB_x *= scale;
}

void Set::SetGBYEmittance(float nemiy)
{
	float yy, ypyp, yyp;
	float* yc{ new float[mSize] {} };
	float* ypc{ new float[mSize] {} };

	for (int i = 0; i < mSize; ++i) {
		ypc[i] = mData[i].GB_y;
		yc[i]  = mData[i].y;
	}
	float mean_ypc = std::accumulate(ypc, ypc + mSize, 0.0) / mSize;
	float mean_yc = std::accumulate(yc, yc + mSize, 0.0) / mSize;

	std::transform(ypc, ypc + mSize, ypc, [mean_ypc](float x) {return x - mean_ypc; });
	std::transform(yc, yc + mSize, yc, [mean_yc](float x) {return x - mean_yc; });

	yy = std::inner_product(yc, yc + mSize, yc, 0.0) / mSize;
	ypyp = std::inner_product(ypc, ypc + mSize, ypc, 0.0) / mSize;
	yyp = std::inner_product(yc, yc + mSize, ypc, 0.0) / mSize;
	float oldnemiy = std::sqrt(yy * ypyp - yyp * yyp);
	delete[] ypc;
	delete[] yc;
	assert(oldnemiy > 0);
	float scale = nemiy / oldnemiy;
	int mags = (int)fabs(log10(scale));
	for (int i = 0; i < mSize; ++i) mData[i].GB_y *= scale;
}
int Set::size() const
{
	return mSize;
}