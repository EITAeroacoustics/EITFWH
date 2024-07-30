#include "Farrasat1A.h"

Farrasat1A::Farrasat1A(
	double MxIn,
	double C0In,
	double rhoIn,
	double pIn,
	bool isPenetrableIn,
	int rankIn) :
	M0(MxIn),
	C0(C0In),
	rho0(rhoIn),
	p0(pIn),
	isPenetrable(isPenetrableIn),
	isMaster(!rankIn)
{
	// Initialize();
};

Farrasat1A::~Farrasat1A()
{
}

// Define the time values for observers£»
void Farrasat1A::InitializeObserverTime(double tMin, double tMax, size_t tNum)
{
	obsTimeMin = tMin;
	obsTimeMax = tMax;
	obsTimeNum = tNum;
	double dt = (obsTimeMax - obsTimeMin) / (obsTimeNum - 1.0);
	for (size_t i = 0; i < obsTimeNum; i++)
		obsTime.push_back(double(i) * dt + obsTimeMin);
}

double Farrasat1A::GetObserverTimeValue(size_t iStep) {
	return obsTime[iStep];
}

double Farrasat1A::GetMaximumObserverTime() {
	return obsTimeMax;
}

double Farrasat1A::GetMinimumObserverTime() {
	return obsTimeMin;
}

//*********************************************
size_t Farrasat1A::GetNumberOfObserverTime() {
	return obsTimeNum;
}


void Farrasat1A::CheckIntegralSurface() {
	//
	vec3 avgNormal, avgCenter;
	// 
	double minArea = *min_element(srcArea.begin(), srcArea.end());
	double maxArea = *max_element(srcArea.begin(), srcArea.end());
	double allArea = accumulate(srcArea.begin(), srcArea.end(), 0.0);
	for (size_t ip = 0; ip < srcNum; ip++)
	{
		allArea += srcArea[ip];
		avgNormal += srcNormVec[0][ip] * srcArea[ip];
		avgCenter += srcLocation[0][ip] * srcArea[ip];
	}
	avgNormal /= allArea;
	avgCenter /= allArea;

	if (isMaster) fcout << "Number of surface panels " << srcNum << endl;
	if (isMaster) fcout << "Averaged normal vector " << avgNormal << endl;
	if (isMaster) fcout << "Averaged panel  center " << avgCenter << endl;
	if (isMaster) fcout << "Min / max panel area " << minArea << " / " << maxArea << endl;
	if (isMaster) fcout << "Total aera " << allArea << endl;
	if (isMaster) fcout << "Minumum / maximum length " << sqrt(minArea) << " / " << sqrt(maxArea) << endl;
}

void Farrasat1A::ReadSurface(
	vector<double>srcAreaIn,
	vector<vector<vec3>>  srcLocationIn,
	vector<vector<vec3>>  srcNormVecIn,
	bool isMovingIn) {
	srcArea = srcAreaIn;
	srcLocation = srcLocationIn;
	srcNormVec = srcNormVecIn;
	isMoving = isMovingIn;
	srcNum = srcArea.size();
}

void Farrasat1A::ReadSurface2D(
	double AreaIn2D,
	vector<vector<vec3>>  LocationIn2D,
	vector<vector<vec3>>  NormVecIn2D,
	bool isMovingIn2D) 
{
	
	srcArea2D = AreaIn2D;
	srcLocation2D = LocationIn2D;
	srcNormVec2D = NormVecIn2D;
	isMoving2D = isMovingIn2D;
	
}

void Farrasat1A::SetFlowData(vector<double> srcTimeIn,
	vector<vector<double>> srcDensityIn,
	vector<vector<double>> srcPressureIn,
	vector<vector<vec3>> srcVelocityIn
) {

	srcTime = srcTimeIn;
	srcDensity = srcDensityIn;
	srcPressure = srcPressureIn;
	srcVelocity = srcVelocityIn;
	/*****************************************************************/
}

void Farrasat1A::SetflowData2D(vector<double> flowTimeIn2D,
	vector<vector<double>> flowDensityIn2D,
	vector<vector<double>> flowPressureIn2D,
	vector<vector<vec3>> flowVelocityIn2D
) {

	srcTime2D = flowTimeIn2D;
	srcDensity2D = flowDensityIn2D;
	srcPressure2D = flowPressureIn2D;
	srcVelocity2D = flowVelocityIn2D;
	
}

void Farrasat1A::CalTimeDelayAtStep(size_t it) {

	// TODO: to put it into a initialize function
	// allocate memory space for radiation radius/vector
	timeDelay.resize(srcNum); 
	R.resize(srcNum);
	Rhat.resize(srcNum);
	for (size_t i = 0; i < srcNum; i++)
	{
		timeDelay[i].resize(obsLocation.size()); 
		R[i].resize(obsLocation.size());
		Rhat[i].resize(obsLocation.size());
	}

	const vec3 Mo = { -M0, 0., 0. };
	const double betaSqr = 1.0 - M0 * M0;

	for (size_t i = 0; i < srcNum; i++) 
	{
		for (size_t j = 0; j < obsLocation.size(); j++) 
		{
			const vec3& xPos = obsLocation[j];
			const vec3& yPos = srcLocation[it][i];
			vec3 rHat;
			rHat = xPos - yPos;
			double r = rHat.mag();
			double Mor = dot(Mo, rHat) / r;
			double td = r / C0 * (Mor + sqrt(Mor * Mor + betaSqr)) / betaSqr;

			// Ref: REF D. Casalino 2003 JSV Eq.(39)
			//  from the source location at emission time to the observer location at the advanced time 
			// 
			//	Below formulation assumes a static observer, but it could be easily extened to consider
			//	any subsonic motion
			rHat += Mo * C0 * td;
			r = rHat.mag();

			R[i][j] = r;
			Rhat[i][j] = rHat / r;
			timeDelay[i][j] = r / C0;
		}
	}

}

// This function computes the maximum/minimum time delay, and it
// works for sources in subsonic motion.
void Farrasat1A::CalTimeDelayMinMax() {

	CalTimeDelayAtStep((size_t)0); 

	timeDelayMin = DBL_MAX;
	timeDelayMax = 0.0;
	for (size_t i = 0; i < srcNum; i++)
	{
		timeDelayMin = min(timeDelayMin, *min_element(timeDelay[i].begin(), timeDelay[i].end()));
		timeDelayMax = max(timeDelayMax, *max_element(timeDelay[i].begin(), timeDelay[i].end()));
	}
}

void Farrasat1A::CalSourceTerms() {
	const double deltaT = obsTime[1] - obsTime[0];

	srcQn.resize(srcTime.size());
	srcLi.resize(srcTime.size());
	for (size_t i = 0; i < srcTime.size(); i++)
	{
		srcQn[i].resize(srcNum);
		srcLi[i].resize(srcNum);
	}

	for (size_t i = 0; i < srcTime.size(); i++)
	{
		if(isMaster)   cout << "\r";
		if(isMaster)  fcout << "Progress: " << floor(100 * (i + 1) / srcTime.size()) << " % " << flush;
		for (size_t j = 0; j < srcNum; j++)
		{
			const double& rho = srcDensity[i][j];
			const vec3& n = srcNormVec[i][j];
			const vec3& u = srcVelocity[i][j]+ vec3(-C0 * M0, 0.0, 0.0);
			double p = srcPressure[i][j] - p0;
			// relative to a static medium
			// TODO: consider on-body fw-h: v=u+vec3(-C0 * M0, 0.0, 0.0);
			vec3 v = vec3(-C0 * M0, 0.0, 0.0);
			if (!isPenetrable) v = u;

			// compute source terms
			vec3 Qi = rho0 * v + rho * (u - v);
			srcQn[i][j] = dot(Qi, n);

			// We have neglected the viscous shear force over the data surface 
			// acting on the flId exterior to the surface
			// TODO:
			// The viscous term, I*uj, needs to be included.
			srcLi[i][j] = p * n + rho * u * dot(u - v, n);
		}
	}
	if (isMaster)   cout << endl;
}

void Farrasat1A::CalSourceTerms2D(void)
{
	srcQn2D.resize(srcTime2D.size());
	srcLi2D.resize(srcTime2D.size());
	
	size_t srcNum2D = srcLocation2D[0].size();
	
	for (size_t i = 0; i < srcTime2D.size(); i++)
	{
		srcQn2D[i].resize(srcNum2D);
		srcLi2D[i].resize(srcNum2D);
	}

	for (size_t i = 0; i < srcTime2D.size(); i++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (i + 1) / srcTime2D.size()) << " % " << flush;
		for (size_t j = 0; j < srcNum2D; j++)
		{
			
			const double& rho = srcDensity2D[i][j];
			const vec3& n = srcNormVec2D[i][j];
			const vec3& u = srcVelocity2D[i][j];
			double p = srcPressure2D[i][j];   
			
			// relative to a uniform x0 medium
			vec3 U0= vec3(C0 * M0, 0.0, 0.0);
			
			// compute source terms
			vec3 Qi = rho * u - rho0 * U0;
			srcQn2D[i][j] = dot(Qi, n);

			// We have neglected the viscous shear force over the data surface 
			
			srcLi2D[i][j] = (p-p0) * n + rho * (u - 2 * U0) * dot(u, n)+ rho0 * U0 * dot(U0, n);

		
		}

	}
	if (isMaster)   cout << endl;

}

void Farrasat1A::InterpTimeDerivative() {
	// 1. this assumes a uniform sampling from cfd, i.e., constant time-step
	// 2. this will produce numerical errors, epsecially at high-frequencies
	const double deltaT = obsTime[1] - obsTime[0];
	srcQnDot.resize(srcTime.size());
	srcLiDot.resize(srcTime.size());
	for (size_t i = 0; i < srcTime.size(); i++)
	{
		srcQnDot[i].resize(srcNum);
		srcLiDot[i].resize(srcNum);
	}

	for (size_t i = 0; i < srcTime.size(); i++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (i + 1) / srcTime.size()) << " % " << flush;
		if (i == 0) // first step
		{
			// first-order forward difference
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] = 0.0;
				srcLiDot[i][j] = vec3(0, 0., 0.);
			}
		}
		else if (i == 1) // second step
		{
			// first-order backward difference
			double oneByDt = 1.0 / deltaT;
			const double bd1w0 = 1.0 * oneByDt;
			const double bd1w1 = -1.0 * oneByDt;
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] = bd1w0 * srcQn[i][j] + bd1w1 * srcQn[i - 1][j];
				srcLiDot[i][j] = bd1w0 * srcLi[i][j] + bd1w1 * srcLi[i - 1][j];
			}
		}
		else if (i == 2) // third step
		{
			// second-order backward difference
			double oneByDt = 1.0 / deltaT;
			const double bd2w0 = 1.5 * oneByDt;
			const double bd2w1 = -2.0 * oneByDt;
			const double bd2w2 = 0.5 * oneByDt;
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] =
					bd2w0 * srcQn[i][j]
					+ bd2w1 * srcQn[i - 1][j]
					+ bd2w2 * srcQn[i - 2][j];
				srcLiDot[i][j] =
					bd2w0 * srcLi[i][j]
					+ bd2w1 * srcLi[i - 1][j]
					+ bd2w2 * srcLi[i - 2][j];
			}
		}
		else if (i == 3) // forth step
		{
			// third-order backward difference
			double oneByDt = 1.0 / deltaT;
			const double bd3w0 = 11. / 6. * oneByDt;
			const double bd3w1 = -3.0 * oneByDt;
			const double bd3w2 = 1.5 * oneByDt;
			const double bd3w3 = -1. / 3. * oneByDt;
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] =
					bd3w0 * srcQn[i][j] + bd3w1 * srcQn[i - 1][j]
					+ bd3w2 * srcQn[i - 2][j]
					+ bd3w3 * srcQn[i - 3][j];
				srcLiDot[i][j] =
					bd3w0 * srcLi[i][j]
					+ bd3w1 * srcLi[i - 1][j]
					+ bd3w2 * srcLi[i - 2][j]
					+ bd3w3 * srcLi[i - 3][j];
			}
		}
		else if (i == 4) // fifth step
		{
			// forth-order backward difference
			double oneByDt = 1.0 / deltaT;
			const double bd4w0 = 25. / 12. * oneByDt;
			const double bd4w1 = -4.0 * oneByDt;
			const double bd4w2 = 3.0 * oneByDt;
			const double bd4w3 = -4. / 3. * oneByDt;
			const double bd4w4 = 1. / 4. * oneByDt;
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] =
					bd4w0 * srcQn[i][j]
					+ bd4w1 * srcQn[i - 1][j]
					+ bd4w2 * srcQn[i - 2][j]
					+ bd4w3 * srcQn[i - 3][j]
					+ bd4w4 * srcQn[i - 4][j];
				srcLiDot[i][j] =
					bd4w0 * srcLi[i][j]
					+ bd4w1 * srcLi[i - 1][j]
					+ bd4w2 * srcLi[i - 2][j]
					+ bd4w3 * srcLi[i - 3][j]
					+ bd4w4 * srcLi[i - 4][j];
			}
		}
		else if (i == 5) // sixth step
		{
			// fifth-order backward difference
			double oneByDt = 1.0 / deltaT;
			const double bd5w0 = 137. / 60. * oneByDt;
			const double bd5w1 = -5.0 * oneByDt;
			const double bd5w2 = 5.0 * oneByDt;
			const double bd5w3 = -10. / 3. * oneByDt;
			const double bd5w4 = 5. / 4. * oneByDt;
			const double bd5w5 = -1. / 5. * oneByDt;
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] =
					bd5w0 * srcQn[i][j]
					+ bd5w1 * srcQn[i - 1][j]
					+ bd5w2 * srcQn[i - 2][j]
					+ bd5w3 * srcQn[i - 3][j]
					+ bd5w4 * srcQn[i - 4][j]
					+ bd5w5 * srcQn[i - 5][j];
				srcLiDot[i][j] =
					bd5w0 * srcLi[i][j]
					+ bd5w1 * srcLi[i - 1][j]
					+ bd5w2 * srcLi[i - 2][j]
					+ bd5w3 * srcLi[i - 3][j]
					+ bd5w4 * srcLi[i - 4][j]
					+ bd5w5 * srcLi[i - 5][j];
			}
		}
		else
		{
			// sixth-order backward difference
			double oneByDt = 1.0 / deltaT;
			const double bd6w0 = 49. / 20. * oneByDt;
			const double bd6w1 = -6. * oneByDt;
			const double bd6w2 = 15. / 2. * oneByDt;
			const double bd6w3 = -20. / 3. * oneByDt;
			const double bd6w4 = 15. / 4. * oneByDt;
			const double bd6w5 = -6. / 5. * oneByDt;
			const double bd6w6 = 1. / 6. * oneByDt;
			for (size_t j = 0; j < srcNum; j++)
			{
				srcQnDot[i][j] =
					bd6w0 * srcQn[i][j]
					+ bd6w1 * srcQn[i - 1][j]
					+ bd6w2 * srcQn[i - 2][j]
					+ bd6w3 * srcQn[i - 3][j]
					+ bd6w4 * srcQn[i - 4][j]
					+ bd6w5 * srcQn[i - 5][j]
					+ bd6w6 * srcQn[i - 6][j];
				srcLiDot[i][j] =
					bd6w0 * srcLi[i][j]
					+ bd6w1 * srcLi[i - 1][j]
					+ bd6w2 * srcLi[i - 2][j]
					+ bd6w3 * srcLi[i - 3][j]
					+ bd6w4 * srcLi[i - 4][j]
					+ bd6w5 * srcLi[i - 5][j]
					+ bd6w6 * srcLi[i - 6][j];
			}
		}
	}
	if (isMaster)   cout << endl;
}

void Farrasat1A::CalTimeSignal() {
	CalSourceTerms();
	InterpTimeDerivative();
	// allocate memory space for pressure at observers
	pPrime.resize(obsTimeNum, vector<double>(obsLocation.size(), 0.0));

	// this implies a uniform pressure sampling at observers
	// it is ok, for the convenience of signal processing, eg. fft.
	const double deltaT = obsTime[1] - obsTime[0];

	vec3 M(-M0, 0, 0);
	double MSqr = M.magSqr();

	for (size_t i = 0; i < srcTime.size(); i++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (i + 1) / srcTime.size()) << " % " << flush;
		// TODO: for a static surface CalTimeDelay only runs once. 
		if (isMoving || i == 0)  CalTimeDelayAtStep(i);

		for (size_t j = 0; j < srcNum; j++)
		{
			const double& dS = srcArea[j];
			const double& Qn = srcQn[i][j];
			const double& QnDot = srcQnDot[i][j];
			const vec3& Li = srcLi[i][j];
			const vec3& LiDot = srcLiDot[i][j];
			double LM = dot(M, Li);//

			for (size_t k = 0; k < obsLocation.size(); k++)
			{
				const double& rMag = R[j][k];	// TODO: may be better to use [j][k]
				const vec3& rHat = Rhat[j][k];

				double Mr = dot(M, rHat);
				double LrDot = dot(LiDot, rHat);
				double Lr = dot(Li, rHat);
				// these did not include moving surface terms, e.g. MDot?
				double S1 = QnDot / rMag / (1 - Mr) / (1 - Mr);
				double S2 = Qn * C0 * (Mr - MSqr) / pow(rMag, 2) / pow(1 - Mr, 3);
				double S3 = LrDot / C0 / rMag / pow(1 - Mr, 2);
				double S4 = (Lr - LM) / pow(rMag, 2) / pow(1 - Mr, 2);
				double S5 = Lr * (Mr - MSqr) / pow(rMag, 2) / pow(1 - Mr, 3);

				const double pi = 2. * acos(0.0);

				double pTotal = (S1 + S2 + S3 + S4 + S5) * dS / 4. / pi;

				// advancing time at observers
				const double advTime = timeDelay[j][k] + srcTime[i];
				size_t advTimeLeftIndex = (size_t)floor((advTime - obsTimeMin) / deltaT);

				if (advTimeLeftIndex >= obsTimeNum)
					continue;
				double weight = 1. - (advTime - obsTime[advTimeLeftIndex]) / deltaT;
				pPrime[advTimeLeftIndex][k] += weight * pTotal;

				// write protection
				size_t advTimeRightIndex = advTimeLeftIndex + 1;
				if (advTimeRightIndex < obsTimeNum)
					pPrime[advTimeRightIndex][k] += (1. - weight) * pTotal;
			}
		}
	}
	if (isMaster) cout << endl;
}



void Farrasat1A::CalFreqSpectra2D()
{
	CalSourceTerms2D();
	size_t srcTNum = srcTime2D.size();
	size_t srcFrequencyNum = (size_t)floor(srcTNum / 2) + 1;
	
	size_t srcNum2D = srcLocation2D[0].size();
	double* in_real;
	fftw_complex* dft_cplx;
	
	in_real = (double*)fftw_malloc(sizeof(double) * srcTNum);
	dft_cplx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * srcFrequencyNum);
	
	for (int i = 0; i < srcTNum; i++)
		in_real[i] = 0.;

	for (int i = 0; i < srcFrequencyNum; i++)
	{
		dft_cplx[i][0] = 0.;
		dft_cplx[i][1] = 0.;
	}
	fftw_plan plan = fftw_plan_dft_r2c_1d((int)srcTNum, in_real, dft_cplx, FFTW_ESTIMATE);
	fftw_execute(plan); 

	// allocate memory space for pressure at observers
	vector<vector<complex<double>>> QnCplx2D(srcFrequencyNum, vector<complex<double>>(srcNum2D, complex<double>(0., 0.)));
	vector<vector<complex<double>>> Li0Cplx2D(srcFrequencyNum, vector<complex<double>>(srcNum2D, complex<double>(0., 0.)));
	vector<vector<complex<double>>> Li1Cplx2D(srcFrequencyNum, vector<complex<double>>(srcNum2D, complex<double>(0., 0.)));
	vector<vector<complex<double>>> Li2Cplx2D(srcFrequencyNum, vector<complex<double>>(srcNum2D, complex<double>(0., 0.)));

	double scaleFactor = 2.0 / srcTNum;
	if (isMaster)  fcout << "Calculate 2D Fourier transform " << endl;
	for (int ip = 0; ip < srcNum2D; ip++) 
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (ip + 1) / srcNum2D) << " % " << flush;

		// DFT for Qn
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcQn2D[it][ip];
		}
		fftw_execute(plan); 
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			QnCplx2D[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}

		// DFT for Li0
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcLi2D[it][ip][0];
		}
		fftw_execute(plan); 
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			Li0Cplx2D[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}

		// DFT for Li1
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcLi2D[it][ip][1];
		}
		fftw_execute(plan); 
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			Li1Cplx2D[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}

		// DFT for Li2
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcLi2D[it][ip][2];
		}
		fftw_execute(plan); 
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			Li2Cplx2D[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}
	}
	fftw_destroy_plan(plan);
	if (isMaster) cout << endl;
	//*******************************************//
	//*******************************************//
	const double pi = 2. * acos(0.0);
	//
	const double deltaT = srcTime2D[1] - srcTime2D[0];
	const double samplingFrequency = 1.0 / deltaT;
	srcFrequency2D.resize(srcFrequencyNum, 0.0);
	double deltaF = samplingFrequency / srcTNum;
	for (int i = 0; i < srcFrequencyNum; i++) {
		srcFrequency2D[i] = double(i) * deltaF;
	}
	srcFrequency2D[0] = 1E-8;

	pPrimeCplx2D.resize(srcFrequencyNum, vector<complex<double>>(obsLocation.size(), complex<double>(0., 0.)));
	const complex<double> I(0., 1.);
	if (isMaster)  fcout << "Calculate 2D pressure spectrum" << endl;

	for (int ik = 0; ik < srcFrequencyNum; ik++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (ik + 1) / srcFrequencyNum) << " % " << flush;
		double omega = 2 * pi * srcFrequency2D[ik];

		for (int ip = 0; ip < srcNum2D; ip++)
		{
			const double& dS = srcArea2D;
			const complex<double>& Qn = QnCplx2D[ik][ip];
			const complex<double>& Li0 = Li0Cplx2D[ik][ip];
			const complex<double>& Li1 = Li1Cplx2D[ik][ip];
			const complex<double>& Li2 = Li2Cplx2D[ik][ip];		

			for (int io = 0; io < obsLocation.size(); io++)
			{
				const vec3& xpos = obsLocation[io];
				const vec3& ypos = srcLocation2D[0][ip];
				double h = 1E-6;
				complex<double> G1 = Green(xpos[0],xpos[1],ypos[0],ypos[1], omega);
				complex<double> dG1 = Green(xpos[0], xpos[1], ypos[0] + h, ypos[1], omega);
				complex<double> dG2 = Green(xpos[0], xpos[1], ypos[0] - h, ypos[1], omega);
				complex<double> dGx = (dG1 - dG2) / (2.0 * h);
				complex<double> dG3 = Green(xpos[0], xpos[1], ypos[0], ypos[1] + h, omega);
				complex<double> dG4 = Green(xpos[0], xpos[1], ypos[0], ypos[1] - h, omega);
				complex<double> dGy = (dG3 - dG4) / (2.0 * h);
				//
				complex<double> S1 = I * omega * Qn * G1;
				complex<double> S2 = Li0 * dGx;
				complex<double> S3 = Li1 * dGy;
				complex<double> pTotal = - (S1 + S2 + S3) * dS;
				//
				pPrimeCplx2D[ik][io] += pTotal;
			}
		}
	}
	if (isMaster)   cout << endl;


}

complex<double> Farrasat1A::Green(double xpos1, double xpos2, double ypos1, double ypos2,double omega1)
{
	
	vec3 M2D(M0, 0, 0);
	double M = M2D.mag();
	double betasqr = 1.0 - M * M; 
	const complex<double> I(0., 1.); 
	
	double r1 = xpos1 - ypos1;
	double r2 = xpos2 - ypos2;
	double R2D = sqrt(pow(r1, 2) + betasqr * pow(r2, 2));
	
	double Hk = omega1 / C0 / betasqr * R2D;
	double H1 = cyl_bessel_j(0, Hk);
	double H2 = cyl_neumann(0, Hk);
	complex<double> Hkel = H1 - I * H2;
	
	complex<double> Gn = I / 4.0 / sqrt(betasqr) * exp(I * M * omega1 / C0 / betasqr * r1) * Hkel;
	return Gn;

}

void Farrasat1A::RecoverSignals2D()
{
	size_t srcTNum = srcTime2D.size();
	size_t srcFrequencyNum = (size_t)floor(srcTNum / 2) + 1;

	// allocate memory space for pressure at observers
	if (isMaster)  fcout << "Calculate pressure signal" << endl;

	// 	declare real number arrays for 1d-fft real input / complex output
	double* out_real;
	fftw_complex* ift_cplx;
	// allocate memory
	out_real = (double*)fftw_malloc(sizeof(double) * srcTNum);
	ift_cplx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * srcFrequencyNum);

	for (int i = 0; i < srcFrequencyNum; i++)
	{
		ift_cplx[i][0] = 0.;
		ift_cplx[i][1] = 0.;
	}
	fftw_plan ift_plan = fftw_plan_dft_c2r_1d((int)srcTNum, ift_cplx, out_real, FFTW_ESTIMATE);
	fftw_execute(ift_plan); // initialize the ift plan

	pPrime2D.resize(srcTNum, vector<double>(obsLocation.size(), 0.0));
	for (int io = 0; io < obsLocation.size(); io++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (io + 1) / obsLocation.size()) << " % " << flush;

		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{ 
			ift_cplx[ik][0] = pPrimeCplx2D[ik][io].real();
			ift_cplx[ik][1] = pPrimeCplx2D[ik][io].imag();
		}
		fftw_execute(ift_plan);
		for (int it = 0; it < srcTNum; it++)
		{
			pPrime2D[it][io] = out_real[it] / 2.0 ;
		}
	}
	fftw_destroy_plan(ift_plan);
	if (isMaster)   cout << endl;
}



void Farrasat1A::CalFreqSpectra() {
	CalSourceTerms();

	size_t srcTNum = srcTime.size();
	size_t srcFrequencyNum = (size_t) floor(srcTNum / 2) + 1;

	// declare real number arrays for 1d-fft real input / complex output
	double* in_real;
	fftw_complex* dft_cplx;
	// allocate memory 
	in_real = (double*)fftw_malloc(sizeof(double) * srcTNum);
	dft_cplx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * srcFrequencyNum);

	for (int i = 0; i < srcTNum; i++)
		in_real[i] = 0.;

	for (int i = 0; i < srcFrequencyNum; i++)
	{
		dft_cplx[i][0] = 0.;
		dft_cplx[i][1] = 0.;
	}
	fftw_plan plan = fftw_plan_dft_r2c_1d((int) srcTNum, in_real, dft_cplx, FFTW_ESTIMATE);
	fftw_execute(plan); 

	// allocate memory space for pressure at observers
	vector<vector<complex<double>>> QnCplx(srcFrequencyNum, vector<complex<double>>(srcNum, complex<double>(0., 0.)));
	vector<vector<complex<double>>> Li0Cplx(srcFrequencyNum, vector<complex<double>>(srcNum, complex<double>(0., 0.)));
	vector<vector<complex<double>>> Li1Cplx(srcFrequencyNum, vector<complex<double>>(srcNum, complex<double>(0., 0.)));
	vector<vector<complex<double>>> Li2Cplx(srcFrequencyNum, vector<complex<double>>(srcNum, complex<double>(0., 0.)));

	double scaleFactor = 2.0 / srcTNum;
	if (isMaster)  fcout << "Calculate Fourier transform " << endl;
	for (int ip = 0; ip < srcNum; ip++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (ip + 1) / srcNum) << " % " << flush;
		
		// DFT for Qn
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcQn[it][ip];
		}
		fftw_execute(plan); // run dft
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			QnCplx[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}

		// DFT for Li0
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcLi[it][ip][0];
		}
		fftw_execute(plan); // run dft
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			Li0Cplx[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}

		// DFT for Li1
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcLi[it][ip][1];
		}
		fftw_execute(plan); // run dft
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			Li1Cplx[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}

		// DFT for Li2
		for (int it = 0; it < srcTNum; it++)
		{
			in_real[it] = srcLi[it][ip][2];
		}
		fftw_execute(plan); // run dft
		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			Li2Cplx[ik][ip] = complex<double>(dft_cplx[ik][0], dft_cplx[ik][1]) * scaleFactor;
		}
	}
	fftw_destroy_plan(plan);
	if (isMaster) cout << endl;

	// this implies a uniform pressure sampling at observers
	// it is ok, for the convenience of signal processing, eg. fft.
	const double pi = 2. * acos(0.0);
	//
	const double deltaT = srcTime[1] - srcTime[0];
	const double samplingFrequency = 1.0 / deltaT;
	srcFrequency.resize(srcFrequencyNum, 0.0);
	double deltaF = samplingFrequency / srcTNum;
	for (int i = 0; i < srcFrequencyNum; i++) {
		srcFrequency[i] = double(i) * deltaF;
	}

	vec3 M(-M0, 0, 0);
	double MSqr = M.magSqr();
	CalTimeDelayAtStep(0);

	pPrimeCplx.resize(srcFrequencyNum, vector<complex<double>>(obsLocation.size(), complex<double>(0., 0.)));
	const complex<double> I(0., 1.);
	if (isMaster)  fcout << "Calculate pressure spectrum" << endl;
	for (int ik = 0; ik < srcFrequencyNum; ik++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (ik + 1) / srcFrequencyNum) << " % " << flush;
		double omega = 2 * pi * srcFrequency[ik];

		for (int ip = 0; ip < srcNum; ip++)
		{
			const double& dS = srcArea[ip];
			const complex<double>& Qn = QnCplx[ik][ip];
			const complex<double>& Li0 = Li0Cplx[ik][ip];
			const complex<double>& Li1 = Li1Cplx[ik][ip];
			const complex<double>& Li2 = Li2Cplx[ik][ip];

			complex<double> QnDot = I * omega * Qn;
			complex<double> Li0Dot = I * omega * Li0;
			complex<double> Li1Dot = I * omega * Li1;
			complex<double> Li2Dot = I * omega * Li2;

			complex<double> LM = M[0] * Li0 + M[1] * Li1 + M[2] * Li2;

			for (int io = 0; io < obsLocation.size(); io++)
			{
				// NOTE:
				// The radiation relation arrays can be computed at run-time,
				//	if their memory space causes issues.
				const double& rMag = R[ip][io];
				const vec3& rHat = Rhat[ip][io];
				const double& td = timeDelay[ip][io];
				double Mr = dot(M, rHat);

				// source-time, this version currently assumes a static integral surface
				complex<double> LrDot = Li0Dot * rHat[0] + Li1Dot * rHat[1] + Li2Dot * rHat[2];
				complex<double> Lr = Li0 * rHat[0] + Li1 * rHat[1] + Li2 * rHat[2];

				// these did not include moving surface terms, e.g. MDot, rDot, nDot, vDot?
				complex<double> S1 = QnDot / rMag / (1 - Mr) / (1 - Mr);
				complex<double> S2 = Qn * C0 * (Mr - MSqr) / pow(rMag, 2) / pow(1 - Mr, 3);
				complex<double> S3 = LrDot / C0 / rMag / pow(1 - Mr, 2);
				complex<double> S4 = (Lr - LM) / pow(rMag, 2) / pow(1 - Mr, 2);
				complex<double> S5 = Lr * (Mr - MSqr) / pow(rMag, 2) / pow(1 - Mr, 3);
				complex<double> pTotal = (S1 + S2 + S3 + S4 + S5) * dS / 4. / pi;

				// 
				pPrimeCplx[ik][io] += pTotal * exp(-I * omega * td);
			}
		}
	}
	if (isMaster)   cout << endl;
}

void Farrasat1A::RecoverSignals()
{
	size_t srcTNum = srcTime.size();
	size_t srcFrequencyNum = (size_t) floor(srcTNum / 2) + 1;

	// allocate memory space for pressure at observers
	if (isMaster)  fcout << "Calculate pressure signal" << endl;

	// 	declare real number arrays for 1d-fft real input / complex output
	double* out_real;
	fftw_complex* ift_cplx;
	// allocate memory
	out_real = (double*)fftw_malloc(sizeof(double) * srcTNum);
	ift_cplx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * srcFrequencyNum);

	for (int i = 0; i < srcFrequencyNum; i++)
	{
		ift_cplx[i][0] = 0.;
		ift_cplx[i][1] = 0.;
	}
	fftw_plan ift_plan = fftw_plan_dft_c2r_1d((int) srcTNum, ift_cplx, out_real, FFTW_ESTIMATE);
	fftw_execute(ift_plan); // initialize the ift plan

	pPrime.resize(obsTimeNum, vector<double>(obsLocation.size(), 0.0));
	for (int io = 0; io < obsLocation.size(); io++)
	{
		if (isMaster)   cout << "\r";
		if (isMaster)  fcout << "Progress: " << floor(100 * (io + 1) / obsLocation.size()) << " % " << flush;

		for (int ik = 0; ik < srcFrequencyNum; ik++)
		{
			ift_cplx[ik][0] = pPrimeCplx[ik][io].real();
			ift_cplx[ik][1] = pPrimeCplx[ik][io].imag();
		}
		fftw_execute(ift_plan);
		for (int it = 0; it < obsTime.size(); it++)
		{
			pPrime[it][io] = out_real[it] / 2.0;
		}
	}
	fftw_destroy_plan(ift_plan);
	if (isMaster)   cout << endl;
}


void Farrasat1A::SaveSpectrumMag(string pFile) {
	
	if (isMaster)   fcout << pFile << endl;

	ofstream spectrumFile_mag(pFile.c_str());
	for (int ik = 0; ik < srcFrequency.size(); ik++) {
		double freq = srcFrequency[ik];
		spectrumFile_mag << freq << " ";
		for (int io = 0; io < obsLocation.size(); io++)
		{
			spectrumFile_mag << abs(pPrimeCplx[ik][io]) << " ";
		}
		spectrumFile_mag << endl;
	}
}

void Farrasat1A::SaveSpectrumMag2D(string pFile) {

	if (isMaster)   fcout << pFile << endl;

	ofstream spectrumFile_mag(pFile.c_str());
	for (int ik = 0; ik < srcFrequency2D.size(); ik++) {
		double freq = srcFrequency2D[ik];
		spectrumFile_mag << freq << " ";
		for (int io = 0; io < obsLocation.size(); io++)
		{
			spectrumFile_mag << abs(pPrimeCplx2D[ik][io]) << " ";
		}
		spectrumFile_mag << endl;
	}
	spectrumFile_mag.close();
}



void Farrasat1A::SaveSpectrumCplx(string pRealFile, string pImagFile) {

	if (isMaster)   fcout << pRealFile << endl;

	ofstream spectrumRealFile(pRealFile.c_str());
	for (int ik = 0; ik < srcFrequency.size(); ik++) {
		double freq = srcFrequency[ik];
		spectrumRealFile << freq << " ";
		for (int io = 0; io < obsLocation.size(); io++)
		{
			spectrumRealFile << real(pPrimeCplx[ik][io]) << " ";
		}
		spectrumRealFile << endl;
	}

	if (isMaster)   fcout << pImagFile << endl;
	ofstream spectrumImagFile(pImagFile.c_str());
	for (int ik = 0; ik < srcFrequency.size(); ik++) {
		double freq = srcFrequency[ik];
		spectrumImagFile << freq << " ";
		for (int io = 0; io < obsLocation.size(); io++)
		{
			spectrumImagFile << imag(pPrimeCplx[ik][io]) << " ";
		}
		spectrumImagFile << endl;
	}
}


void Farrasat1A::SaveSpectrumCplx2D(string pRealFile, string pImagFile) {

	if (isMaster)   fcout << pRealFile << endl;

	ofstream spectrumRealFile(pRealFile.c_str());
	for (int ik = 0; ik < srcFrequency2D.size(); ik++) {
		double freq = srcFrequency2D[ik];
		spectrumRealFile << freq << " ";
		for (int io = 0; io < obsLocation.size(); io++)
		{
			spectrumRealFile << real(pPrimeCplx2D[ik][io]) << " ";
		}
		spectrumRealFile << endl;
	}

	if (isMaster)   fcout << pImagFile << endl;
	ofstream spectrumImagFile(pImagFile.c_str());
	for (int ik = 0; ik < srcFrequency2D.size(); ik++) {
		double freq = srcFrequency2D[ik];
		spectrumImagFile << freq << " ";
		for (int io = 0; io < obsLocation.size(); io++)
		{
			spectrumImagFile << imag(pPrimeCplx2D[ik][io]) << " ";
		}
		spectrumImagFile << endl;
	}
}




void Farrasat1A::SaveTimeSignals(string pFile) {
	if (isMaster)   fcout << pFile << endl;
	// TDODO: deine output path and file name
	ofstream outfile(pFile.c_str(), ios::out);
	for (size_t i = 0; i < obsTime.size(); i++) 
	{
		outfile << obsTime[i] << " ";
		for (size_t j = 0; j < obsLocation.size(); j++) {
			outfile << pPrime[i][j] << " ";
		}
		outfile << "\n";
	}
	outfile.close();
}



void Farrasat1A::SaveTimeSignals2D(string pFile) {
	if (isMaster)   fcout << pFile << endl;
	// TDODO: deine output path and file name
	ofstream outfile2D(pFile.c_str(), ios::out);
	for (size_t i = 0; i < srcTime2D.size(); i++)
	{
		outfile2D << srcTime2D[i]+20.0 << " ";
		for (size_t j = 0; j < obsLocation.size(); j++) {
			outfile2D << pPrime2D[i][j] << " ";
		}
		outfile2D << "\n";
	}
	outfile2D.close();
}
