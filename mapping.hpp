#include <math.h>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "particle.hpp"
#include <map>
#include <string>
#include <boost/function.hpp>

//typedef float (*FnPtr)(Set&, float, char);
//typedef boost::function<float (Set& s, float time, char t, float& result) noexcept> fun_t;
//typedef std::map<std::string, fun_t> funs_t;


#ifndef MAPHEADERDEF
#define MAPVDEADERDEF
class Map
{
private:
	//fun_t funcArray[36];
	//float (funcArray[36])(Set& s, float time, char t, float& result) noexcept;
	float* mData;
	float* mData_scaled;
	int mSize;
	char type;
	float L;
	float phi0;
	float omega;
	float B0;
	float v0;
	unsigned int nP;
	//std::map<std::string, float(*)(Set& s, float time, char t, float& result) noexcept> mMap;
	//funs_t mMap;
	typedef void (Map::* method)(Set&, float, char, float&);
	method funcArray[36] = {	&Map::M_00, &Map::M_01, &Map::M_02, &Map::M_03, &Map::M_04, &Map::M_05,
								&Map::M_10, &Map::M_11, &Map::M_12, &Map::M_13, &Map::M_14, &Map::M_15,
								&Map::M_20, &Map::M_21, &Map::M_22, &Map::M_23, &Map::M_24, &Map::M_25,
								&Map::M_30, &Map::M_31, &Map::M_32, &Map::M_33, &Map::M_34, &Map::M_35,
								&Map::M_40, &Map::M_41, &Map::M_42, &Map::M_43, &Map::M_44, &Map::M_45,
								&Map::M_50, &Map::M_51, &Map::M_52, &Map::M_53, &Map::M_54, &Map::M_55 };

public:
	Map(int size, char t, int numberParticles); //constructor 	float* yc{ new float[mSize] {} };
	~Map(); //destructor to (use delete[] mData)
	//define matrix elements
	void M_00(Set& s, float time, char t, float& result) noexcept; void M_01(Set& s, float time, char t, float& result) noexcept; void M_02(Set& s, float time, char t, float& result) noexcept; void M_03(Set& s, float time, char t, float& result) noexcept; void M_04(Set& s, float time, char t, float& result) noexcept; void M_05(Set& s, float time, char t, float& result) noexcept; 
	void M_10(Set& s, float time, char t, float& result) noexcept; void M_11(Set& s, float time, char t, float& result) noexcept; void M_12(Set& s, float time, char t, float& result) noexcept; void M_13(Set& s, float time, char t, float& result) noexcept; void M_14(Set& s, float time, char t, float& result) noexcept; void M_15(Set& s, float time, char t, float& result) noexcept; 
	void M_20(Set& s, float time, char t, float& result) noexcept; void M_21(Set& s, float time, char t, float& result) noexcept; void M_22(Set& s, float time, char t, float& result) noexcept; void M_23(Set& s, float time, char t, float& result) noexcept; void M_24(Set& s, float time, char t, float& result) noexcept; void M_25(Set& s, float time, char t, float& result) noexcept; 
	void M_30(Set& s, float time, char t, float& result) noexcept; void M_31(Set& s, float time, char t, float& result) noexcept; void M_32(Set& s, float time, char t, float& result) noexcept; void M_33(Set& s, float time, char t, float& result) noexcept; void M_34(Set& s, float time, char t, float& result) noexcept; void M_35(Set& s, float time, char t, float& result) noexcept; 
	void M_40(Set& s, float time, char t, float& result) noexcept; void M_41(Set& s, float time, char t, float& result) noexcept; void M_42(Set& s, float time, char t, float& result) noexcept; void M_43(Set& s, float time, char t, float& result) noexcept; void M_44(Set& s, float time, char t, float& result) noexcept; void M_45(Set& s, float time, char t, float& result) noexcept; 
	void M_50(Set& s, float time, char t, float& result) noexcept; void M_51(Set& s, float time, char t, float& result) noexcept; void M_52(Set& s, float time, char t, float& result) noexcept; void M_53(Set& s, float time, char t, float& result) noexcept; void M_54(Set& s, float time, char t, float& result) noexcept; void M_55(Set& s, float time, char t, float& result) noexcept; 
	//update mData elements for given time and 6-vectors
	void Update(Set& s, float time, char t);
	void ScaleMap();
	void Set_L(float length) noexcept;
	void Set_Phi0(float phi0) noexcept;
	void Set_Omega(float omega) noexcept;
	void Set_B0(float b) noexcept;
	void Set_v0(float v) noexcept;
	float* GetScaled();
	float* GetData();

	unsigned int GetNp();
};
#endif
