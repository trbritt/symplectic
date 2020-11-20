#include "mapping.hpp"
#include "particle.hpp"
#include <map>
#include <string>
#include <omp.h>
#include <boost/lexical_cast.hpp>

extern float q_e, c, m_e, PI;
Map::Map(int size, char t, int numberParticles)
{
	if ((t != 'c') ||( t != 's')) type = 'c';
	else type = t;
	mSize = size;
	nP = numberParticles;
	mData = new float[size] {0.0f};
	mData_scaled = new float[size * numberParticles * numberParticles]{ 0.0f };

	/*
	mMap["M_00") = &Map::M_00;  mMap["M_01") = &Map::M_01; mMap["M_02") = &Map::M_02;	 mMap["M_03") = &Map::M_03; mMap["M_04") = &Map::M_04; mMap["M_05") = &Map::M_05;
	mMap["M_10") = &Map::M_10;  mMap["M_11") = &Map::M_11; mMap["M_12") = &Map::M_12;	 mMap["M_13") = &Map::M_13; mMap["M_14") = &Map::M_14; mMap["M_15") = &Map::M_15;
	mMap["M_20") = &Map::M_20;  mMap["M_21") = &Map::M_21; mMap["M_22") = &Map::M_22;	 mMap["M_23") = &Map::M_23; mMap["M_24") = &Map::M_24; mMap["M_25") = &Map::M_25;
	mMap["M_30") = &Map::M_30;  mMap["M_31") = &Map::M_31; mMap["M_32") = &Map::M_32;	 mMap["M_33") = &Map::M_33; mMap["M_34") = &Map::M_34; mMap["M_35") = &Map::M_35;
	mMap["M_40") = &Map::M_40;  mMap["M_41") = &Map::M_41; mMap["M_42") = &Map::M_42;	 mMap["M_43") = &Map::M_43; mMap["M_44") = &Map::M_44; mMap["M_45") = &Map::M_45;
	mMap["M_50") = &Map::M_50;  mMap["M_51") = &Map::M_51; mMap["M_52") = &Map::M_52;	 mMap["M_53") = &Map::M_53; mMap["M_54") = &Map::M_54; mMap["M_55") = &Map::M_55;*/

	/*funcArray[0]  = &Map::M_00;  funcArray[1]  = &Map::M_01; funcArray[2]  = &Map::M_02; funcArray[3]  = &Map::M_03; funcArray[4]  = &Map::M_04; funcArray[5]  = &Map::M_05;
	funcArray[6]  = &Map::M_10;  funcArray[7]  = &Map::M_11; funcArray[8]  = &Map::M_12; funcArray[9]  = &Map::M_13; funcArray[10) = &Map::M_14; funcArray[11) = &Map::M_15;
	funcArray[12) = &Map::M_20;  funcArray[13) = &Map::M_21; funcArray[14) = &Map::M_22; funcArray[15) = &Map::M_23; funcArray[16) = &Map::M_24; funcArray[17) = &Map::M_25;
	funcArray[18) = &Map::M_30;  funcArray[19) = &Map::M_31; funcArray[20) = &Map::M_32; funcArray[21) = &Map::M_33; funcArray[22) = &Map::M_34; funcArray[23) = &Map::M_35;
	funcArray[24) = &Map::M_40;  funcArray[25) = &Map::M_41; funcArray[26) = &Map::M_42; funcArray[27) = &Map::M_43; funcArray[28) = &Map::M_44; funcArray[29) = &Map::M_45;
	funcArray[30) = &Map::M_50;  funcArray[31) = &Map::M_51; funcArray[32) = &Map::M_52; funcArray[33) = &Map::M_53; funcArray[34) = &Map::M_54; funcArray[35) = &Map::M_55;*/


}
Map::~Map()
{
	delete [] mData;
	if (*mData_scaled != NULL) delete [] mData_scaled;
}
void Map::Set_L(float length) noexcept 
{
	L = length;
}
void Map::Set_Phi0(float p) noexcept
{
	phi0 = p;
}
void Map::Set_Omega(float o) noexcept
{
	omega = o;
}
void Map::Set_B0(float b) noexcept
{
	B0 = b;
}
void Map::Set_v0(float v) noexcept
{
	v0 = v;
}
//Below are the massive 36 functions defining the elements of the maps
//throw std::invalid_argument(boost::lexical_cast<std::string>(EXIT_FAILURE));

//first row
void Map::M_00(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		result = 1.0f;
		*(mData + 6 * 0 + 0) = result; break; }
		
	case 's': {
		result = 0;
		*(mData + 6 * 0 + 0) = result; break; }

	default:
	result=0.0f;
	}
}
void Map::M_01(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		result = L;
		*(mData + 6 * 0 + 1) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 0 + 1) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_02(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		result = 0.0f;
		*(mData + 6 * 0 + 2) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 0 + 2) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_03(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		result = 0.0f;
		*(mData + 6 * 0 + 3) = result; break; }
	case 's': {
		result = 0.0f;
		*(mData + 6 * 0 + 3) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_04(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		result = -(omega_c / omega) * (LAM * cos(phi0) + sin(phi0) - sin(LAM + phi0));
		*(mData + 6 * 0 + 4) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 0 + 4) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_05(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		result = (L / LAM) * (omega_c / omega) * (cos(LAM + phi0) - cos(phi0) + LAM * sin(LAM + phi0));
		*(mData + 6 * 0 + 5) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 0 + 5) = result; break; }
	default:
	result=0.0f;
	}
}
//second row
void Map::M_10(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		result = 0.0f;
		*(mData + 6 * 1 + 0) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 1 + 0) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_11(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		result = 1.0f;
		*(mData + 6 * 1 + 1) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 1 + 1) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_12(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 1 + 2) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 1 + 2) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_13(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 1 + 3) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 1 + 3) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_14(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		result= -(LAM / L) * (omega_c / omega) * (cos(phi0) - cos(LAM + phi0)); 
		*(mData + 6 * 1 + 4) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 1 + 4) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_15(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		result= (omega_c / omega) * (LAM * cos(LAM + phi0) + sin(phi0) - sin(LAM + phi0)); 
		*(mData + 6 * 1 + 5) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 1 + 5) = result; break; }
	default:
	result=0.0f;
	}
}
//third row
void Map::M_20(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 2 + 0) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 2 + 0) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_21(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 2 + 1) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 2 + 1) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_22(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 1.0f;
		*(mData + 6 * 2 + 2) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 2 + 2) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_23(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= L;
		*(mData + 6 * 2 + 3) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 2 + 3) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_24(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 2 + 4) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 2 + 4) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_25(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 2 + 5) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 2 + 5) = result; break; }
	default:
	result=0.0f;
	}
}
//fourth row
void Map::M_30(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 3 + 0) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 3 + 0) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_31(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 3 + 1) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 3 + 1) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_32(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 3 + 2) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 3 + 2) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_33(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 1.0f;
		*(mData + 6 * 3 + 3) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 3 + 3) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_34(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 3 + 4) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 3 + 4) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_35(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 3 + 5) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 3 + 5) = result; break; }
	default:
	result=0.0f;
	}
}
//fifth row
void Map::M_40(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		result= -(omega_c / omega) * (LAM * cos(phi0) + sin(phi0) - sin(LAM + phi0)) / (g0 * g0); 
		*(mData + 6 * 4 + 0) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 4 + 0) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_41(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= (L / LAM) * (omega_c / omega) * ((1 - 2 * beta0 * beta0) * (cos(LAM + phi0)
			- cos(phi0)) - beta0 * beta0 * LAM * sin(phi0)
			+ (LAM / (g0 * g0)) * sin(LAM + phi0)); 
		*(mData + 6 * 4 + 1) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 4 + 1) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_42(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 4 + 2) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 4 + 2) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_43(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 4 + 3) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 4 + 3) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_44(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= 1 +
			(omega_c / omega) * (omega_c / omega) * (
				(0.5f - 1.25f * beta0 * beta0) * cos(2 * phi0)
				+ (-0.5f + 0.25f * beta0 * beta0) * cos(2 * (LAM + phi0))
				+ beta0 * beta0 * (cos(LAM + 2 * phi0) + (LAM / 2) * sin(2 * phi0))
				- (LAM / (g0 * g0)) * sin(LAM + 2 * phi0)
				); 
		*(mData + 6 * 4 + 4) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 4 + 4) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_45(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= L +
			L * (omega_c / omega) * (omega_c / omega) *
			(
				0.5f * (1.0f - 4.0f * beta0 * beta0) * cos(LAM)
				+ 0.25f * beta0 * beta0 * cos(2 * phi0)
				+ 0.5f * (1.0f - 0.5f * beta0 * beta0) * cos(2 * (LAM + phi0))
				+ 0.5f * (1.0f + beta0 * beta0) * (cos(LAM + 2 * phi0) - 1)
				+ (3.0f * beta0 * beta0 / LAM + LAM / (2 * g0 * g0)) * sin(LAM)
				+ 5.0f * beta0 * beta0 / (4 * LAM) * sin(2 * phi0)
				- beta0 * beta0 / (4 * LAM) * sin(2 * (LAM + phi0))
				); 
		*(mData + 6 * 4 + 5) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 4 + 5) = result; break; }
	default:
	result=0.0f;
	}
}
//sixth row
void Map::M_50(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= -(LAM / L) * (omega_c / omega) * (cos(phi0) - cos(LAM + phi0)) / (g0 * g0); 
		*(mData + 6 * 5 + 0) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 5 + 0) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_51(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= (omega_c / omega) * (
			(LAM / (g0 * g0)) * cos(LAM + phi0)
			+ beta0 * beta0 * (sin(LAM + phi0) - sin(phi0))
			); 
		*(mData + 6 * 5 + 1) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 5 + 1) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_52(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 5 + 2) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 5 + 2) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_53(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c':{
		result= 0.0f;
		*(mData + 6 * 5 + 3) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 5 + 3) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_54(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= (omega_c / omega) * (omega_c / omega) * (LAM / L) *
			(
				sin(2 * (LAM + phi0))
				- sin(LAM + 2 * phi0)
				- ((LAM / (g0 * g0)) + beta0 * beta0 * sin(LAM)) * cos(LAM + 2 * phi0)
				); 
		*(mData + 6 * 5 + 4) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 5 + 4) = result; break; }
	default:
	result=0.0f;
	}
}
void Map::M_55(Set& s, float time, char t, float& result) noexcept
{
	switch (t) {
	case 'c': {
		float g0 = 1 / sqrt(1 - v0 * v0 / (c * c));
		float omega_c = q_e * B0 / (g0 * m_e);
		float LAM = omega * L / v0;
		float beta0 = v0 / c;
		result= 1 + (omega_c / omega) * (omega_c / omega) *
			(
				LAM *
				(
					sin(2 * (LAM + phi0)) - sin(LAM + 2 * phi0) + LAM * sin(phi0) * sin(LAM + phi0)
					)
				- beta0 * beta0 *
				(
					cos(LAM + phi0) - cos(phi0) + LAM * sin(phi0)
					) * LAM * sin(LAM + phi0)
				- 2 * beta0 * beta0 *
				(
					1.0f - cos(LAM) + LAM * cos(LAM + phi0) * sin(phi0)
					+ 0.25f *
					(
						cos(2 * (LAM + phi0)) - cos(2 * phi0)
						)
					)
				); 
		*(mData + 6 * 5 + 5) = result; break; }
	case 's': {
		result = 0;
		*(mData + 6 * 5 + 5) = result; break; }
	default:
	result=0.0f;
	}
}


void Map::Update(Set& s, float time, char t)
{
	#pragma omp parallel for collapse(2) num_threads(11)
	for (unsigned int i = 0; i < 6; ++i) {
		for (unsigned int j = 0; j < 6; ++j) {
			//*(mData + 6 * i + j) = reinterpret_cast<float(*)()>(mMap[name])(s, time, t);
			const int idx = 6 * i + j;
			const method func = funcArray[idx];
			assert(func != NULL);
			float result;
			(this->*func)(s, time, t, result);
		}
	}
}
void Map::ScaleMap()
{
#pragma omp parallel for collapse(3) num_threads(11)
	for (unsigned int n = 0; n < nP; ++n) {
		for (unsigned int i = 0; i < 6; ++i) {
			for (unsigned int j = 0; j < 6; ++j) {
				//m[w*i+j) = m[w*i+j + n(w^2+h)]
				*(mData_scaled + static_cast<_int64>(nP*6*(i+n*6)+j+n*6)) = *(mData + static_cast<_int64>(6 * i + j));
			}
		}
	}
}
float* Map::GetScaled()
{
	return mData_scaled;
}
float* Map::GetData()
{
	return mData;
}
unsigned int Map::GetNp()
{
	return nP;
}