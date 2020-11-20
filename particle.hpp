#pragma once
#include<iostream>
#include<chrono>
#include<random>
#include <vector>
#include <string>



struct Particle {
	float x;
	float GB_x;
	float y;
	float GB_y;
	float z;
	float GB_z;
};
#ifndef SETHEADERDEF
#define SETHEADERDEF

class Set 
{
private:
	std::vector<Particle> mData;
	int mSize;
public:
	Set(int size);
	Set& operator=(const Set& s);
	Particle& operator[](int i); //zero based particle indexing
	Particle& operator()(int i); //one based particle indexing
	Set operator+(const Set& s1) const; //binary +
	Set operator-(const Set& s1) const; //binary -
	void add(const Set& s1); //change local memory, instead of generating new 'set' object
	void subtract(const Set& s1);

	void SetRXYPhidist(float mu, float sigma, int sout); //gaussian in R, uniform in phi
	void SetZdist(float mu, float sigma, float sleft, float sright);
	//void SetPhidist(float center, float width);
	void SetScale(float a1, float a2, float a3, float b1, float b2, float b3);

	void SetGBRXYPhidist(float mu, float sigma, int sout);
	void SetGBZdist(float center, float width);

	void AddXdiv(float center, float div);
	void AddYdiv(float center, float div);
	void AddZdiv(float center, float div);

	void SetGBXEmittance(float emit);
	void SetGBYEmittance(float emit);
	int size() const;
};
#endif