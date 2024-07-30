/*****************************************************************/
// Refs:
// a. 
// <An advanced time approach for acoustic analogy predictions> 
// Casalino, 2003
// b. 
// <Derivation of Formulations 1 and 1A of Farassat>
// Farassat, 2007
// 
// This version assumes: 
// 1. a static integral penetrable surface
// 2. freestream flow is along x+ direction
// 3. The cell number/area of the integral surface is constant
// 4. The viscous shear force over the data surface is neglected
// 
// TODO: 
//		Moving surface, even only on-body integration?
// 
//	MxInf, rhoInf, cInf, isPenetrable
/*****************************************************************/

#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <limits.h>
#include <complex>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>
#include "fftwLib/fftw3.h"
#include "vec3.h"
#include <cmath>

// for debug only
// show function names before printing information to terminal
#define fcout cout << __func__ << ":: "

using namespace std;

class Farrasat1A
{
public:
	Farrasat1A(double MxIn, double C0In, double rhoIn, double pIn, bool isPenetrableIn, int rankIn);
	~Farrasat1A();

	void Initialize(void)
	{
		if (isMaster) fcout << "Mach number x+ " << M0 << endl;
		if (isMaster) fcout << "Speed of sound " << C0 << endl;
		if (isMaster) fcout << "Density " << rho0 << endl;
		if (isMaster) fcout << "Pressure " << p0 << endl;
	}

	// Member functions declaration
	void ReadObserverLocation(vector<vec3> obsLocationIn) {
		obsLocation = obsLocationIn;
	}

	// Print all observers to terminal
	void PrintObserverLocations(void) {
		for (auto eachObsever : obsLocation)
		{
			if (isMaster) fcout << eachObsever << endl;
		}
	}

	// Get the number of observers
	size_t GetObserverNum(void) {
		return obsLocation.size();
	}

	// Define the time values for observers
	void InitializeObserverTime(double tMin, double tMax, size_t tNum);

	// Get observer-time values
	double GetObserverTimeValue(size_t iStep);

	// Get maximum observer-time value
	double GetMaximumObserverTime(void);

	// Get minimum observer-time value
	double GetMinimumObserverTime(void);

	// Get the total number of all observers
	size_t GetNumberOfObserverTime(void);

	// This is for debug and check if the surface
	//	is properly defined by looking at its averaged
	// center, normal vector and total area
	void CheckIntegralSurface(void);

	// Read integral surface, s, x, n
	// TODO:
	//		For parallel computation, this should be 
	//		isolated. Read outside this class
	//		and transfer partitioned data only
	void ReadSurface(
		vector<double>		  srcAreaIn,
		vector<vector<vec3>>  srcLocationIn,
		vector<vector<vec3>>  srcNormVecIn,
		bool isMovingIn);


	// Read 2D integral surface, s, x, n
	void ReadSurface2D(
		double AreaIn2D,
		vector<vector<vec3>>  LocationIn2D,
		vector<vector<vec3>>  NormVecIn2D,
		bool isMovingIn2D);


	// Read flow variables, rho, p, u from a file
	// TODO:
	//		For parallel computation, this should be 
	//		isolated. Read outside this class
	//		and transfer partitioned data only
	void SetFlowData(vector<double> srcTimeIn,
		vector<vector<double>>	srcDensityIn,
		vector<vector<double>>	srcPressureIn,
		vector<vector<vec3>>	srcVelocityIn
	);
	
	// Read 2D flow variables, rho, p, u from a file
	void SetflowData2D(vector<double> flowTimeIn2D,
		vector<vector<double>> flowDensityIn2D,
		vector<vector<double>> flowPressureIn2D,
		vector<vector<vec3>> flowVelocityIn2D
	);
	// Calculate source-observer radiation relations
	// this is useful when deciding the output time series
	// !!! Please note when the integral surface is moving,
	//	   make sure it is excuted every source time-step
	void CalTimeDelayAtStep(size_t it);

	// Calculate the minimum and maximum time-delay values.
	// It seems this function is needed only for GetTimeDelayMin()
	void CalTimeDelayMinMax(void);

	// Get the minimum time delay for sound radiation.
	double GetTimeDelayMin(void) {
		return timeDelayMin;
	}

	// Get the maximum time delay for sound radiation.
	double GetTimeDelayMax(void) {
		return timeDelayMax;
	}

	// Get how many time steps at sources
	// this is automatically defined when reading flow data.
	size_t GetNumberOfSourceTime(void) {
		return srcTime.size();
	}

	// Get how many time steps at sources
// this is automatically defined when reading flow data.
	size_t GetNumberOfSourceTime2D(void) {
		return srcTime2D.size();
	}

	// Get how many time steps at sources
	// this is automatically defined when reading flow data.
	size_t GetNumberOfSourceFrequency(void) {
		return srcFrequency.size();
	}

	size_t GetNumberOfSourceFrequency2D(void) {
		return srcFrequency2D.size();
	}


	// Calculate source terms Qn and Li
	void CalSourceTerms(void);

	// Calculate 2D source terms Qn and Li
	void CalSourceTerms2D(void);


	// Calculate source time-derivatives QnDot, LiDot
	// via interpolation
	// this may cause phase or amplitude error
	void InterpTimeDerivative(void);

	// do forward-time computation for pressure at observers
	void CalTimeSignal(void);

	// Solve FW-H equations using frequency-domain formulation
	void CalFreqSpectra(void);

	// Solve 2D FW-H  frequency-domain formulation
	void CalFreqSpectra2D(void);

	//2D green function
	complex<double> Green(double xpos1, double xpos2, double ypos1, double ypos2, double omega1);

	// Recover pressure signals via inverse Fourier transform
	// Note it is aligned with the source
	void RecoverSignals(void);
	// Recover pressure signals via inverse Fourier transform 2D
	void RecoverSignals2D(void);


	// get the minimum source time value
	double GetSourceTimeMin(void) {
		return *min_element(srcTime.begin(), srcTime.end());
	}

	// Get the maximum source time value
	double GetSourceTimeMax(void) {
		return *max_element(srcTime.begin(), srcTime.end());
	}

	// Write far-field pressure signals into a given file
	//	pFile defines the file name
	void SaveTimeSignals(string pFile);

	void SaveTimeSignals2D(string pFile);

	// Write far-field complex pressure amplitudes into a given file
	//	pFile defines the file name
	void SaveSpectrumMag(string pFile);

	void SaveSpectrumMag2D(string pFile);

	// This function currently saves complex pressure amplitudes to 
	// two files. It is expected to serve for BEM solvers.
	// Note: 
	//	data format is to be decided, suggest using HDF5 Binary files
	void SaveSpectrumCplx(string pRealFile, string pImagFile);

	void SaveSpectrumCplx2D(string pRealFile, string pImagFile);


	// output sound pressure signals at all observers
	vector<vector<double>>	pPrime;

	// complex pressure amplitudes at all observers
	vector<vector<complex<double>>> pPrimeCplx;

	vector<vector<complex<double>>> pPrimeCplx2D;
	
	vector<vector<double>>	pPrime2D;

private:

	/**************** define freestream conditions *******************/
	// freestream Mach number, [Mx, 0, 0]
	double		M0;

	// freestream density 
	double		rho0;

	// freestream speed of sound
	double		C0;

	// freestream pressure
	double		p0;
	/*****************************************************************/

	//
	bool isPenetrable;
	bool isMoving = false;

	/**************** define observer time ***************************/
	vector<double>		obsTime;		// observer-obsTimeime series 
	double		obsTimeMin;
	double		obsTimeMax;
	size_t		obsTimeNum;	// number of time steps at observers
	/*****************************************************************/

	/**************** define source time *****************************/
	vector<double>		srcTime;	// source-time series 
	vector<double> srcFrequency;	// source-frequency series 
	/*****************************************************************/

	/**************** define observers *******************************/
	// Cartesian coordinates, x, static
	vector<vec3> obsLocation;
	/*****************************************************************/

	/**************** define sources *********************************/
	// number of cells on the integral surface
	size_t		srcNum;

	// area of local surface, cell size
	vector<double>	srcArea;

	// moving integral surface with unchanged shapes
	// source points in Cartesian coordinates, y
	vector<vector<vec3>>	srcLocation;

	// surface normal vectors, n
	vector<vector<vec3>>	srcNormVec;

	// flow density, rho
	vector<vector<double>>	srcDensity;

	// flow pressure, p
	vector<vector<double>>	srcPressure;

	// flow velocity, U
	vector<vector<vec3>>	srcVelocity;

	/**************** define 2D source *****************************/

	vector<double> srcTime2D;
	
	vector<vector<double>> srcDensity2D;
	
	vector<vector<double>> srcPressure2D;
	
	vector<vector<vec3>> srcVelocity2D;
	
	double srcArea2D;
	
	vector<vector<vec3>>  srcLocation2D;
	
	vector<vector<vec3>>  srcNormVec2D;
	
	bool isMoving2D;
	
	// 2D sound sources 
	vector<vector<double>>	srcQn2D;
	vector<vector<vec3>>	srcLi2D;
	
	vector<double> srcFrequency2D;	// source-frequency series 



	// sound sources and their time-derivatives
	vector<vector<double>>	srcQn;
	vector<vector<double>>	srcQnDot;
	vector<vector<vec3>>	srcLi;
	vector<vector<vec3>>	srcLiDot;

	/**************** define sources-observer coupled ****************/
	vector<vector<double>>	timeDelay;
	double timeDelayMax;
	double timeDelayMin;

	// The acoustic phase distance
	vector<vector<double>>	R;

	// The acoustic amplitude distance
	vector<vector<double>>	Ra;

	// The unit radiation vector
	vector<vector<vec3>>	Rhat;
	/*****************************************************************/

	// if this rank is master processor, printing is activated
	bool isMaster;

};
