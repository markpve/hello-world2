//
//


#include <stdexcept>
#include<ppl.h>

#include "emfield.h"


using namespace std;
using namespace concurrency;


namespace FDTD
{

	//*****************************************************************************************************************************************
	//*****************************************************************************************************************************************
	//Functions for EMField3D

	EMField3D::EMField3D(double(*IncidentEzField)(double t, double x), double(*IncidentEyField)(double t, double x), Structure3DWithGrid& aStr, int anumkVals, int amb) :
		
		Ex(aStr.sizeX() - 1, vector<vector<double>>(aStr.sizeY(), vector<double>(aStr.sizeZ(), 0.0))),
		Ey(aStr.sizeX(), vector<vector<double>>(aStr.sizeY() - 1, vector<double>(aStr.sizeZ(), 0.0))),
		Ez(aStr.sizeX(), vector<vector<double>>(aStr.sizeY(), vector<double>(aStr.sizeZ() - 1, 0.0))),
		Hx(aStr.sizeX(), vector<vector<double>>(aStr.sizeY() - 1, vector<double>(aStr.sizeZ() - 1, 0.0))),
		Hy(aStr.sizeX() - 1, vector<vector<double>>(aStr.sizeY(), vector<double>(aStr.sizeZ() - 1, 0.0))),
		Hz(aStr.sizeX() - 1, vector<vector<double>>(aStr.sizeY() - 1, vector<double>(aStr.sizeZ(), 0.0))),
		Number_of_power_monitors(aStr.numLx() + 2), 
		Pflow(anumkVals, Number_of_power_monitors)

	{
		int i = 0; int n = 0; int p = 0;
		int displayCode = 0;

		Location_index_for_incident_power = Number_of_power_monitors - 1;
		Location_index_for_incident_power = aStr.numLx();
		Location_index_for_reflected_power = 0;


		M = aStr.sizeX();
		N = aStr.sizeY();
		P = aStr.sizeZ();

		theGrid = *(aStr.grid());

		if(displayCode == 0)
			mb = min(amb, int(theGrid.nStepsXL(0) / 2 + 0.5));	//Note that nNds(0) must be at least two in order for a SF/TF boundary to be applied.
		else
			mb = int(theGrid.nStepsXL(0) / 2 + 0.5);	//Note that nNds(0) must be at least two in order for a SF/TF boundary to be applied.
		//End if

		mD = theGrid.nStepsXL(0) - 1;				//Should just assign value of mD from a Grid3D accessor. Must have mD >= mb + 2

		Eincz = IncidentEzField;
		Eincy = IncidentEyField;

		q = 0;
		delt = theGrid.minDelta()*Sc / c0;	//delt is in microseconds and delx in micrometers.

		ScExterior = c0*delt / theGrid.deltaxL(0);

		numkVals = anumkVals;



		/*
		//Initialize plane  of Ez and Ey nodes with x-index mb to incident field values
		for (n = 0; n < N; n++)
			for (p = 0; p < P - 1; p++)
			{
				Ez[mb][n][p] = Eincz(0, mb*theGrid.deltaxL(0));
			}//End for


		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P; p++)
			{
				Ey[mb][n][p] = Eincy(0, mb*theGrid.deltaxL(0));
			}//End for
		*/

	}//End Constructor EMField3D



	int EMField3D::propagateIn(Structure3DWithGrid& str, double aTinc, int maxTime, int ak0, int adeltak,
								int(*displayFunctionToCall)(const vector<vector<vector<double>>>& F, int q),
								int displayTime, int(*displayFunctionToCall_for_1D)(const EMField1D_zPolarized*, int timeIndex))
	{
		int i;
		int j;
		int k;
		int m;
		int n;
		int p;
		int pFlowIndex;
		int NT = maxTime;
		int Npulse;

		//**************************************************
		//objects needed for time-average power calculations

		k0 = ak0;
		deltak = adeltak;

		Npulse = int(aTinc / delt + 0.5);			//duration of incident pulse in time steps

		//const int num_xPlanes = 2;		//Number of x = const planes for which time-average power calculated in addition to incident power

		Number_of_power_monitors = str.numLx() + 2;

		Location_index_for_incident_power = Number_of_power_monitors - 1;

		Location_index_for_transmitted_power = str.numLx();

		Location_index_for_reflected_power = 0;

		vector <int> xIndexes_for_power_calculation(str.numLx() + 2, 0);		//Vector to store locations of x=const planes for which incident, 
																						//reflected power and transmitted power are to be calculated.

		xIndexes_for_power_calculation[Location_index_for_reflected_power] = mb - 2;		//Node plane index at which reflected power calculated

		xIndexes_for_power_calculation[1] = str.numSteps_xDir_throuh(0) - 1;		//Node plane index at which power flow just above structure is calculated

		for (i = 1; i < str.numLx(); i++)
		{
			xIndexes_for_power_calculation[i + 1] = str.numSteps_xDir_throuh(i - 1) + 2;		//Node plane index two steps into grid layer with index i
		}//End for

		xIndexes_for_power_calculation[Location_index_for_incident_power] = mb;		//Node plane index just beyond TFSF boundary


		//Arrays for real and imaginary parts of component phasors needed for time-average power calculations

		Array1D_Vector3D Re_Eys(numkVals, Number_of_power_monitors, N - 1, P);
		Array1D_Vector3D Im_Eys(numkVals, Number_of_power_monitors, N - 1, P);
		Array1D_Vector3D Re_C_Hzs(numkVals, Number_of_power_monitors, N - 1, P);
		Array1D_Vector3D Im_C_Hzs(numkVals, Number_of_power_monitors, N - 1, P);

		Array1D_Vector3D Re_Ezs(numkVals, Number_of_power_monitors, N, P - 1);
		Array1D_Vector3D Im_Ezs(numkVals, Number_of_power_monitors, N, P - 1);
		Array1D_Vector3D Re_C_Hys(numkVals, Number_of_power_monitors, N, P - 1);
		Array1D_Vector3D Im_C_Hys(numkVals, Number_of_power_monitors, N, P - 1);

		Array1D_Vector3D Pave_x(numkVals, Number_of_power_monitors, N, P);

		Array1D_Vector1D Pflow_freq(numkVals, Number_of_power_monitors);


		//*************************************************************************
		//Objects needed for analytical absorbing boundary conditions at front and rear x = const boundaries

		double C0;		//Coefficients for analytical absorbing boundary conditions
		double C1;
		double C2;
		double d;

		//****************************************************************************************
		//Vectors for field arrays for abc

		vector <double> Ez1D(P, 0);
		vector<vector<vector<double>>> EzOldLeft0(3, vector <vector<double>>(N, Ez1D));
		vector<vector<vector<double>>> EzOldLeft1(3, vector <vector<double>>(N, Ez1D));
		vector<vector<vector<double>>> EzOldRight0(3, vector <vector<double>>(N, Ez1D));
		vector<vector<vector<double>>> EzOldRight1(3, vector <vector<double>>(N, Ez1D));

		vector <double> Ey1D(P, 0);
		vector<vector<vector<double>>> EyOldLeft0(3, vector <vector<double>>(N, Ey1D));
		vector<vector<vector<double>>> EyOldLeft1(3, vector <vector<double>>(N, Ey1D));
		vector<vector<vector<double>>> EyOldRight0(3, vector <vector<double>>(N, Ey1D));
		vector<vector<vector<double>>> EyOldRight1(3, vector <vector<double>>(N, Ey1D));

		//Compute coefficients for second-order abc
		computeCoefficientsForSecondOrderABCs(d, C0, C1, C2);


		//*****************************************************************************************
		//objects needed for auxiliary one-dimensional incident field simulation to be used for corrections along the x=const TFSF boundary.
		//One-dimensional y-polarized and z-polarized incident fields will be simulated separately for the region (vacuum) above the structure. 
		//These simulated fields will be used to make corrections along the TFSF boundary in the three-dimensional simulation.

		//Instantiate a two-layer Grid object. The first layer is for a vacuum and the second is for a lossy layer for grid termination.
		Grid1D auxGrid(mD + int((M - mD)/2.0), theGrid.deltaxL(0));

		//To simulate the incident field, instantiate an EMField1DzPolar using the auxiliary two-layer Grid1D .
		EMField1D_zPolarized IFieldz(Eincz, &auxGrid, delt);

		//To simulate the incident field, instantiate an EMField1DyPolar field using the auxiliary two-layer Grid1D .
		EMField1D_yPolarized IFieldy(Eincy, &auxGrid, delt);

		//Instantiate a medium object representing vaccuum for the auxiliary 1D fields
		Medium aM(0, 1, 1);

		//1D coefficients
		vector <double> CEzE1D(auxGrid.size(), 0);
		vector <double> CEzH1D(auxGrid.size(), 0);
		vector <double> CHyE1D(auxGrid.size(), 0);
		vector <double> CHyH1D(auxGrid.size(), 0);

		vector <double> CEyE1D(auxGrid.size(), 0);
		vector <double> CEyH1D(auxGrid.size(), 0);
		vector <double> CHzE1D(auxGrid.size(), 0);
		vector <double> CHzH1D(auxGrid.size(), 0);

		IFieldz.compute_coefficients_for_medium_with_lossy_layer_underneath(aM, CEzE1D, CEzH1D, CHyE1D, CHyH1D);
		IFieldy.compute_coefficients_for_medium_with_lossy_layer_underneath(aM, CEyE1D, CEyH1D, CHzE1D, CHzH1D);


		//*******************************************************************************************
		//Objects needed for boundary conditions for periodicity in the y and z directions

		//Allocate memory for magnetic field components just outside the high y-boundary and just outside the high z-boundary
		//Used only if the structure is the unit cell of a larger periodic structure

		vector <vector<double>> HxOutsideY(M, vector<double>(P - 1, 0.0));			//Hx at nodes jsut outside high y-boundary. 
													//Used only if the structure is the unit cell of a larger structure periodic in the y-direction 
		vector <vector<double>> HzOutsideY(M - 1, vector<double>(P, 0.0));			//Hz at nodes jsut outside high y-boundary. 
													//Used only if the structure is the unit cell of a larger structure periodic in the y-direction 

		vector <vector<double>> HxOutsideZ(M, vector<double>(N - 1, 0.0));			//Hx at nodes jsut outside high z-boundary. 
													//Used only if the structure is the unit cell of a larger structure periodic in the z-direction 

		vector <vector<double>> HyOutsideZ(M - 1, vector<double>(N, 0.0));			//Hy at nodes jsut outside high z-boundary. 
													//Used only if the structure is the unit cell of a larger structure periodic in the z-direction



		//***************************************************************************************************************************
		//Time stepping

		for (q = 0; q < maxTime; q++)
		{

			updateMagneticField(str, IFieldz, IFieldy);

			correctMagneticField_alongTFSFxFrontBoundary(str, IFieldz, IFieldy);

			//for normal incidence
			if (str.Is_a_unitCell())
				updateMagneticField_outsideBoundaries_ofGridUnitCell(str, HxOutsideY, HzOutsideY, HxOutsideZ, HyOutsideZ);	//If the structure is not periodic in both y and z directions, need more code.
			else
			{
				//Need more code if structure not a unit cell or error condition
			}//End if

																															//Update auxilliary one-dimensional incident fields to use for correction to nodes adjacent to TF/SF boundary. 
			IFieldz.Update(CEzE1D, CEzH1D, CHyE1D, CHyH1D);
			IFieldy.Update(CEyE1D, CEyH1D, CHzE1D, CHzH1D);

			if (q % displayTime == 0)
			{
				displayFunctionToCall_for_1D(&IFieldz, q);
			}//End if


			updateElectricField(str, IFieldz, IFieldy);

			//For normal incidence
			if (str.Is_a_unitCell())
				updateElectricField_onBoundaries_ofGridUnitCell(str, HxOutsideY, HzOutsideY, HxOutsideZ, HyOutsideZ);		//If the structure is not periodic in both y and z directions, need more code.
			else
			{
				//Need more code if structure not a unit cell or error condition
			}//End if

			correctElectricField_alongTFSFxFrontBoundary(str, IFieldz, IFieldy);

			applySecondOrderABC_on_xBoundaries(EzOldLeft0, EzOldLeft1, EzOldRight0, EzOldRight1, EyOldLeft0, EyOldLeft1, EyOldRight0, EyOldRight1, d, C0, C1, C2);

			runningSum_for_Component_Phasors(Re_Eys, Im_Eys, Re_C_Hzs, Im_C_Hzs,
				Re_Ezs, Im_Ezs, Re_C_Hys, Im_C_Hys,
				IFieldz, IFieldy,
				xIndexes_for_power_calculation, Npulse);

			if (q % displayTime == 0)
			{
				//EMField3D* PointerToMe = this;
				displayFunctionToCall(Ez, q);
			}

		}//End Time stepping


		//**********************************************************************************************
		//Compute spectral time-average power flows

		adjust_component_phasors_to_Yee_cell_face_centers(Re_Eys, Im_Eys, Re_C_Hzs, Im_C_Hzs,
			Re_Ezs, Im_Ezs, Re_C_Hys, Im_C_Hys,
			xIndexes_for_power_calculation);		//component phasors adjusted to Yee cell face centers on each xp plane

		compute_spectral_power_densities(Re_Eys, Im_Eys, Re_C_Hzs, Im_C_Hzs,
			Re_Ezs, Im_Ezs, Re_C_Hys, Im_C_Hys, Pave_x,
			xIndexes_for_power_calculation);	//Time average x-component spectral power density calculated at Yee cell face centers on each xp plane


		compute_spectral_power_flows(Pave_x, str, xIndexes_for_power_calculation);


		return 0;
	}//End function propagateIn of class EMField3D




	void EMField3D::computeCoefficientsForSecondOrderABCs(double& d, double& C0, double& C1, double& C2)
	{
		d = ScExterior + 2.0 + 1 / ScExterior;

		C0 = -(ScExterior - 2.0 + 1 / ScExterior) / d;
		C1 = -2.0*(ScExterior - 1 / ScExterior) / d;
		C2 = 4.0*(ScExterior + 1 / ScExterior) / d;
	}//End function computeCoefficientsForSecondOrderABCs



	int EMField3D::updateMagneticField(Structure3DWithGrid& s, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy)
	{
		int m = 0;
		int n = 0;
		int p = 0;


		//Update Hx

		/************************************************************************************************************************************************************
		for (m = 0; m < M; m++)
			for (n = 0; n < N - 1; n++)
				for (p = 0; p < P - 1; p++)
					Hx[m][n][p] = s.CHxH(m, n, p)*Hx[m][n][p] - s.CHxEz(m, n, p)*(Ez[m][n + 1][p] - Ez[m][n][p]) + s.CHxEy(m, n, p)*(Ey[m][n][p + 1] - Ey[m][n][p]);

		*/

		parallel_for(0, M, [this, &s](size_t m)
							{parallel_for(0, N - 1, [this, &s, m](size_t n)
													{parallel_for(0, P - 1, [this, &s, m, n](size_t p)
																			{Hx[m][n][p] = s.CHxH(m, n, p)*Hx[m][n][p] - s.CHxEz(m, n, p)*(Ez[m][n + 1][p] - Ez[m][n][p])
																							+ s.CHxEy(m, n, p)*(Ey[m][n][p + 1] - Ey[m][n][p]); }); }); });
		

		//Update Hy

		/**

		for (m = 0; m < M - 1; m++)
			for (n = 0; n < N; n++)
				for (p = 0; p < P - 1; p++)
					Hy[m][n][p] = s.CHyH(m, n, p)*Hy[m][n][p] - s.CHyEx(m, n, p)*(Ex[m][n][p + 1] - Ex[m][n][p]) + s.CHyEz(m, n, p)*(Ez[m + 1][n][p] - Ez[m][n][p]);
		*/
	

		parallel_for(0, M - 1, [this, &s](size_t m)
								{parallel_for(0, N, [this, &s, m](size_t n)
													{parallel_for(0, P - 1, [this, &s, m, n](size_t p)
																			{Hy[m][n][p] = s.CHyH(m, n, p)*Hy[m][n][p] - s.CHyEx(m, n, p)*(Ex[m][n][p + 1] - Ex[m][n][p])
																							+ s.CHyEz(m, n, p)*(Ez[m + 1][n][p] - Ez[m][n][p]); }); }); });


		//Update Hz

		/**
		for (m = 0; m < M - 1; m++)
			for (n = 0; n < N - 1; n++)
				for (p = 0; p < P; p++)
					Hz[m][n][p] = s.CHzH(m, n, p) * Hz[m][n][p] - s.CHzEy(m, n, p) * (Ey[m + 1][n][p] - Ey[m][n][p]) + s.CHzEx(m, n, p)*(Ex[m][n + 1][p] - Ex[m][n][p]);

		*/

		parallel_for(0, M - 1, [this, &s](size_t m)
								{parallel_for(0, N - 1, [this, &s, m](size_t n)
														{parallel_for(0, P, [this, &s, m, n](size_t p)
																			{Hz[m][n][p] = s.CHzH(m, n, p) * Hz[m][n][p] - s.CHzEy(m, n, p) * (Ey[m + 1][n][p] - Ey[m][n][p]) + 
																							s.CHzEx(m, n, p)*(Ex[m][n + 1][p] - Ex[m][n][p]); }); }); });


		return 0;
	}//End function updateMagneticField



	void EMField3D::correctMagneticField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy)
	{
		int n, p;

		//Correct Hy along TFSF boundary
		/****************************************************************************************************************************************
		for (n = 0; n < N; n++)
			for (p = 0; p < P - 1; p++)
				Hy[mb - 1][n][p] -= s.CHyEz(mb - 1, n, p)*IFieldz.Ezdir(mb);
		*/

		parallel_for(0, N, [this, &s, &IFieldz](size_t n) 
							{parallel_for(0, P - 1, [this, &s, &IFieldz, n](size_t p)
													{Hy[mb - 1][n][p] -= s.CHyEz(mb - 1, n, p)*IFieldz.Ezdir(mb); }); });


		//Correct Hz along TFSF boundary
		/******************************************************************************************************************************************
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P; p++)
				Hz[mb - 1][n][p] += s.CHzEy(mb - 1, n, p)*IFieldy.Eydir(mb);
		*/


		parallel_for(0, N - 1, [this, &s, &IFieldy](size_t n)
							{parallel_for(0, P, [this, &s, &IFieldy, n](size_t p)
													{Hz[mb - 1][n][p] += s.CHzEy(mb - 1, n, p)*IFieldy.Eydir(mb); }); });

	}//End function correctMagneticField_alongTFSFxFrontBoundary



	int EMField3D::updateMagneticField_outsideBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
		vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ)
	{
		int m, n, p;

		if (s.Is_periodic_in_yDirection())
		{
			/***********************************************************************************************************************************************
			for (m = 0; m < M; m++)
				for (p = 0; p < P - 1; p++)
					HxOutsideY[m][p] = Hx[m][0][p];
			//End for
			//End for
			*/

			parallel_for(0, M, [this, &HxOutsideY](size_t m) {parallel_for(0, P - 1, [this, &HxOutsideY, m](size_t p) {HxOutsideY[m][p] = Hx[m][0][p]; }); });


			/************************************************************************************************************************************************
			for (m = 0; m < M - 1; m++)
				for (p = 0; p < P; p++)
					HzOutsideY[m][p] = Hz[m][0][p];
			//End for
			//End for
			*/

			parallel_for(0, M - 1, [this, &HzOutsideY](size_t m) {parallel_for(0, P, [this, &HzOutsideY, m](size_t p) {HzOutsideY[m][p] = Hz[m][0][p]; }); });

		}//End if


		if (s.Is_periodic_in_zDirection())
		{
			/**************************************************************************************************************************************************
			for (m = 0; m < M; m++)
				for (n = 0; n < N - 1; n++)
					HxOutsideZ[m][n] = Hx[m][n][0];
			//End for
			//End for
			*/

			parallel_for(0, M, [this, &HxOutsideZ](size_t m) {parallel_for(0, N - 1, [this, &HxOutsideZ, m](size_t n) {HxOutsideZ[m][n] = Hx[m][n][0]; }); });


			/**************************************************************************************************************************************************
			for (m = 0; m < M - 1; m++)
				for (n = 0; n < N; n++)
					HyOutsideZ[m][n] = Hy[m][n][0];
			//End for
			//End for
			*/

			parallel_for(0, M - 1, [this, &HyOutsideZ](size_t m) {parallel_for(0, N, [this, &HyOutsideZ, m](size_t n) {HyOutsideZ[m][n] = Hy[m][n][0]; }); });

		}//End if

		return 0;
	}//End functionupdateMagneticField_outsideBoundaries_ofGridUnitCell



	int EMField3D::updateElectricField(Structure3DWithGrid& s, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy)
	{
		int m;
		int n;
		int p;

		//Update Ez

		/*******************************************************************************************************************************************************************
		for (m = 1; m < M - 1; m++)
			for (n = 1; n < N - 1; n++)
				for (p = 0; p < P - 1; p++)
					Ez[m][n][p] = s.CEzE(m, n, p) * Ez[m][n][p] + s.CEzHy(m, n, p) * (Hy[m][n][p] - Hy[m - 1][n][p]) - s.CEzHx(m, n, p)*(Hx[m][n][p] - Hx[m][n - 1][p]);
		*/

		parallel_for(1, M - 1, [this, &s](size_t m) 
								{parallel_for(1, N - 1, [this, &s, m] (size_t n)
														{parallel_for(0, P - 1, [this, &s, m, n] (size_t p)
																				{Ez[m][n][p] = s.CEzE(m, n, p) * Ez[m][n][p] + s.CEzHy(m, n, p) * (Hy[m][n][p] - Hy[m - 1][n][p]) 
																								- s.CEzHx(m, n, p)*(Hx[m][n][p] - Hx[m][n - 1][p]);} );} );} );

		//Update Ey

		/********************************************************************************************************************************************************************
		for (m = 1; m < M - 1; m++)
			for (n = 0; n < N - 1; n++)
				for (p = 1; p < P - 1; p++)
					Ey[m][n][p] = s.CEyE(m, n, p) * Ey[m][n][p] + s.CEyHx(m, n, p) * (Hx[m][n][p] - Hx[m][n][p - 1]) - s.CEyHz(m, n, p)*(Hz[m][n][p] - Hz[m - 1][n][p]);

		*/

		parallel_for(1, M - 1, [this, &s](size_t m)
								{parallel_for(0, N - 1, [this, &s, m](size_t n)
														{parallel_for(1, P - 1, [this, &s, m, n](size_t p)
																				{Ey[m][n][p] = s.CEyE(m, n, p) * Ey[m][n][p] + s.CEyHx(m, n, p) * (Hx[m][n][p] - Hx[m][n][p - 1])
																								- s.CEyHz(m, n, p)*(Hz[m][n][p] - Hz[m - 1][n][p]); }); }); });


		//Update Ex

		/********************************************************************************************************************************************************************
		for (m = 0; m < M - 1; m++)
			for (n = 1; n < N - 1; n++)
				for (p = 1; p < P - 1; p++)
					Ex[m][n][p] = s.CExE(m, n, p) * Ex[m][n][p] + s.CExHz(m, n, p) * (Hz[m][n][p] - Hz[m][n - 1][p]) - s.CExHy(m, n, p)*(Hy[m][n][p] - Hy[m][n][p - 1]);

		*/


		parallel_for(0, M - 1, [this, &s](size_t m)
								{parallel_for(1, N - 1, [this, &s, m](size_t n)
														{parallel_for(1, P - 1, [this, &s, m, n](size_t p)
																				{Ex[m][n][p] = s.CExE(m, n, p) * Ex[m][n][p] + s.CExHz(m, n, p) * (Hz[m][n][p] - Hz[m][n - 1][p])
																								- s.CExHy(m, n, p)*(Hy[m][n][p] - Hy[m][n][p - 1]); }); }); });

		return 0;
	}//End function updateElectricField



	void EMField3D::updateElectricField_onBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
		vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ)
	{
		int m, n, p;

		if (s.Is_periodic_in_yDirection())
		{
			//Update Ez on high y-boundary of unit cell and apply y-periodicity

			/***********************************************************************************************************************************************************************
			for (m = 1; m < M - 1; m++)
				for (p = 0; p < P - 1; p++)
				{
					Ez[m][N - 1][p] = s.CEzE(m, N - 1, p) * Ez[m][N - 1][p] + s.CEzHy(m, N - 1, p) * (Hy[m][N - 1][p] - Hy[m - 1][N - 1][p]) - s.CEzHx(m, N - 1, p)*(HxOutsideY[m][p] - Hx[m][N - 1 - 1][p]);
					Ez[m][0][p] = Ez[m][N - 1][p];
				}//end for
				 //End for

			*/

			parallel_for(1, M - 1, [this, &s, &HxOutsideY](size_t m)
									{parallel_for(0, P - 1, [this, &s, &HxOutsideY, m](size_t p)
															{Ez[m][N - 1][p] = s.CEzE(m, N - 1, p) * Ez[m][N - 1][p] + s.CEzHy(m, N - 1, p) * (Hy[m][N - 1][p] - Hy[m - 1][N - 1][p])
																				- s.CEzHx(m, N - 1, p)*(HxOutsideY[m][p] - Hx[m][N - 1 - 1][p]);
															Ez[m][0][p] = Ez[m][N - 1][p]; }); });


			//Update Ex on high y-boundary of unit cell except for nodes also on high z-boundary and apply y-periodicity

			/***************************************************************************************************************************************************************************
			for (m = 0; m < M - 1; m++)
				for (p = 1; p < P - 1; p++)
				{
					Ex[m][N - 1][p] = s.CExE(m, N - 1, p) * Ex[m][N - 1][p] + s.CExHz(m, N - 1, p) * (HzOutsideY[m][p] - Hz[m][N - 1 - 1][p]) - s.CExHy(m, N - 1, p)*(Hy[m][N - 1][p] - Hy[m][N - 1][p - 1]);
					Ex[m][0][p] = Ex[m][N - 1][p];
				}//End for
				 //End for

			*/

			parallel_for(0, M - 1, [this, &s, &HzOutsideY](size_t m)
									{parallel_for(1, P - 1, [this, &s, &HzOutsideY, m](size_t p)
															{Ex[m][N - 1][p] = s.CExE(m, N - 1, p) * Ex[m][N - 1][p] + s.CExHz(m, N - 1, p) * (HzOutsideY[m][p] - Hz[m][N - 1 - 1][p])
																				- s.CExHy(m, N - 1, p)*(Hy[m][N - 1][p] - Hy[m][N - 1][p - 1]);
															Ex[m][0][p] = Ex[m][N - 1][p]; }); });

		}//End if

		if (s.Is_periodic_in_zDirection())
		{
			//Update Ex on high z-boundary of unit cell including nodes also on high y-boundary and apply z-periodicity and apply y-periodiciy if structure also y-periodic

			/***********************************************************************************************************************************************************************
			for (m = 0; m < M - 1; m++)
			{
				for (n = 1; n < N - 1; n++)
				{
					Ex[m][n][P - 1] = s.CExE(m, n, P - 1) * Ex[m][n][P - 1] + s.CExHz(m, n, P - 1) * (Hz[m][n][P - 1] - Hz[m][n - 1][P - 1]) - s.CExHy(m, n, P - 1)*(HyOutsideZ[m][n] - Hy[m][n][P - 1 - 1]);
					Ex[m][n][0] = Ex[m][n][P - 1];
				}//end for


				if (s.Is_periodic_in_yDirection())
				{
					Ex[m][N - 1][P - 1] = s.CExE(m, N - 1, P - 1) * Ex[m][N - 1][P - 1] + s.CExHz(m, N - 1, P - 1) * (HzOutsideY[m][P - 1] - Hz[m][n - 1 ][P - 1])
						- s.CExHy(m, n, P - 1)*(HyOutsideZ[m][n] - Hy[m][n][P - 1 - 1]);
					Ex[m][0][0] = Ex[m][N - 1][P - 1];
				}//end if

			}//End for
			*/


			parallel_for(0, M - 1, [this, &s, &HyOutsideZ](size_t m)
									{parallel_for(1, N - 1, [this, &s, &HyOutsideZ, m](size_t n)
															{Ex[m][n][P - 1] = s.CExE(m, n, P - 1) * Ex[m][n][P - 1] + s.CExHz(m, n, P - 1) * (Hz[m][n][P - 1] - Hz[m][n - 1][P - 1])
																				- s.CExHy(m, n, P - 1)*(HyOutsideZ[m][n] - Hy[m][n][P - 1 - 1]);
															Ex[m][n][0] = Ex[m][n][P - 1]; }); });

			if (s.Is_periodic_in_yDirection())
				parallel_for(0, M - 1, [this, &s, &HyOutsideZ, &HzOutsideY](size_t m)
										{Ex[m][N - 1][P - 1] = s.CExE(m, N - 1, P - 1) * Ex[m][N - 1][P - 1] + s.CExHz(m, N - 1, P - 1) * (HzOutsideY[m][P - 1] - Hz[m][N - 2][P - 1])
																- s.CExHy(m, N - 1, P - 1)*(HyOutsideZ[m][N - 1] - Hy[m][N - 1][P - 1 - 1]);
										Ex[m][0][0] = Ex[m][N - 1][P - 1]; });



			 //Update Ey on high z-boundary of unit cell and apply z-periodicity

			/*****************************************************************************************************************************************************************************
			for (m = 1; m < M - 1; m++)
				for (n = 0; n < N - 1; n++)
				{
					Ey[m][n][P - 1] = s.CEyE(m, n, P - 1) * Ey[m][n][P - 1] + s.CEyHx(m, n, P - 1) * (HxOutsideZ[m][n] - Hx[m][n][P - 1 - 1]) - s.CEyHz(m, n, P - 1)*(Hz[m][n][P - 1] - Hz[m - 1][n][P - 1]);
					Ey[m][n][0] = Ey[m][n][P - 1];
				}//End for
			*/

			parallel_for(1, M - 1, [this, &s, &HxOutsideZ](size_t m)
									{parallel_for(0, N - 1, [this, &s, &HxOutsideZ, m](size_t n)
															{Ey[m][n][P - 1] = s.CEyE(m, n, P - 1) * Ey[m][n][P - 1] + s.CEyHx(m, n, P - 1) * (HxOutsideZ[m][n] - Hx[m][n][P - 1 - 1])
																				- s.CEyHz(m, n, P - 1)*(Hz[m][n][P - 1] - Hz[m - 1][n][P - 1]);
															Ey[m][n][0] = Ey[m][n][P - 1]; }); });

		}//End if

	}//End function updateElectricField_onBoundaries_ofGridUnitCell



	void EMField3D::correctElectricField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy)
	{
		int n, p;

		//Correct Ez along TFSF boundary

		/**********************************************************************************************************************************************
		for (n = 1; n < N; n++)
			for (p = 0; p < P - 1; p++)
				Ez[mb][n][p] -= s.CEzHy(mb, n, p) * IFieldz.Hydir(mb - 1);
		*/


		parallel_for(1, N, [this, &s, &IFieldz](size_t n) 
							{parallel_for(0, P - 1, [this, &s, &IFieldz, n](size_t p)
													{Ez[mb][n][p] -= s.CEzHy(mb, n, p) * IFieldz.Hydir(mb - 1); }); });

		//Correct Ey along TFSF boundary

		/************************************************************************************************************************************
		for (n = 0; n < N - 1; n++)
			for (p = 1; p < P; p++)
				Ey[mb][n][p] += s.CEzHy(mb, n, p)*IFieldy.Hzdir(mb - 1);
		*/


		parallel_for(0, N - 1, [this, &s, &IFieldy](size_t n)
							{parallel_for(1, P, [this, &s, &IFieldy, n](size_t p)
													{Ey[mb][n][p] += s.CEzHy(mb, n, p)*IFieldy.Hzdir(mb - 1); }); });

	}//End function correctElectricField_alongTFSFxFrontBoundary



	void EMField3D::applySecondOrderABC_on_xBoundaries(vector<vector<vector<double>>> & EzOldLeft0, vector<vector<vector<double>>>& EzOldLeft1,
		vector<vector<vector<double>>>& EzOldRight0, vector<vector<vector<double>>>& EzOldRight1,
		vector<vector<vector<double>>> & EyOldLeft0, vector<vector<vector<double>>>& EyOldLeft1,
		vector<vector<vector<double>>>& EyOldRight0, vector<vector<vector<double>>>& EyOldRight1,
		double d, double C0, double C1, double C2)
	{
		int m;
		int n;
		int p;

		int mBack;

		/*//Simple ABC for Ez y-boundaries			Change to periodic boudary conditions
		for (m = 0; m < M; m++)
		{
		Ez[m][0] = Ez[m][1];
		Ez[m][N - 1] = Ez[m][N - 2];
		}//End for */

		//Second-order ABC for Ez left and right (i.e. above and below device)				
		for (n = 0; n < N; n++)
			for (p = 0; p < P - 1; p++)
			{
				Ez[0][n][p] = C0*(Ez[2][n][p] + EzOldLeft1[0][n][p]) + C1*(EzOldLeft0[0][n][p] + EzOldLeft0[2][n][p] - Ez[1][n][p] - EzOldLeft1[1][n][p]) + C2*EzOldLeft0[1][n][p] - EzOldLeft1[2][n][p];
				Ez[M - 1][n][p] = C0*(Ez[M - 3][n][p] + EzOldRight1[0][n][p]) + C1*(EzOldRight0[0][n][p] + EzOldRight0[2][n][p] - Ez[M - 2][n][p] - EzOldRight1[1][n][p])
					+ C2*EzOldRight0[1][n][p] - EzOldRight1[2][n][p];
			}//End for

			 //Update past field values
		for (n = 0; n < N; n++)
			for (p = 0; p < P - 1; p++)
			{
				for (m = 0; m < 3; m++)
				{
					EzOldLeft1[m][n][p] = EzOldLeft0[m][n][p];
					EzOldLeft0[m][n][p] = Ez[m][n][p];
				}//End for each m

				for (mBack = 0; mBack < 3; mBack++)
				{
					EzOldRight1[mBack][n][p] = EzOldRight0[mBack][n][p];
					EzOldRight0[mBack][n][p] = Ez[M - 1 - mBack][n][p];
				}//End for each mBack

			}//End For


			 //Second-order ABC for Ey left and right (i.e. above and below device)				
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P; p++)
			{
				Ey[0][n][p] = C0*(Ey[2][n][p] + EyOldLeft1[0][n][p]) + C1*(EyOldLeft0[0][n][p] + EyOldLeft0[2][n][p] - Ey[1][n][p] - EyOldLeft1[1][n][p]) + C2*EyOldLeft0[1][n][p] - EyOldLeft1[2][n][p];
				Ey[M - 1][n][p] = C0*(Ey[M - 3][n][p] + EyOldRight1[0][n][p]) + C1*(EyOldRight0[0][n][p] + EyOldRight0[2][n][p] - Ey[M - 2][n][p] - EyOldRight1[1][n][p])
					+ C2*EyOldRight0[1][n][p] - EyOldRight1[2][n][p];
			}//End for

			 //Update past field values
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P; p++)
			{
				for (m = 0; m < 3; m++)
				{
					EyOldLeft1[m][n][p] = EyOldLeft0[m][n][p];
					EyOldLeft0[m][n][p] = Ey[m][n][p];
				}//End for each m

				for (mBack = 0; mBack < 3; mBack++)
				{
					EyOldRight1[mBack][n][p] = EyOldRight0[mBack][n][p];
					EyOldRight0[mBack][n][p] = Ey[M - 1 - mBack][n][p];
				}//End for each mBack

			}//End For

			 //no need second-order ABCs for magnetic fields. 2/19/16

	}//End function applySecondOrderABC_on_xBoundaries

	

	void EMField3D::runningSum_for_Component_Phasors(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys,
		EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy,
		 const vector<int>& mp, int Np)
	{
		int xP, n, p, k, vInd;
		double wk;

		//Update phasor components on planes for reflected and transmitted power

		/**********************************************************************************************************************************************
		for (vInd = 0; vInd < numkVals; vInd++)
		{
			k = k0 + vInd*deltak;
			wk = (2 * PI*k / Npulse);
			for (xP = 0; xP < Number_of_power_monitors - 1; xP++)
			{
				for (n = 0; n < N - 1; n++)
					for (p = 0; p < P; p++)
					{
						Re_Eys(vInd)[xP][n][p] += (2.0 / Npulse)*Ey[mp[xP]][n][p] * cos(wk*q);
						Im_Eys(vInd)[xP][n][p] += (2.0 / Npulse)*Ey[mp[xP]][n][p] * (-sin(wk*q));

						Re_C_Hzs(vInd)[xP][n][p] += (2.0 / Npulse)*0.5*(Hz[mp[xP]][n][p] + Hz[mp[xP] - 1][n][p])*cos(wk*q);
						Im_C_Hzs(vInd)[xP][n][p] += (2.0 / Npulse)*0.5*(Hz[mp[xP]][n][p] + Hz[mp[xP] - 1][n][p])*sin(wk*q);

					}// End calculation for Eys and C_Hzs

				for (n = 0; n < N; n++)
					for (p = 0; p < P - 1; p++)
					{
						Re_Ezs(vInd)[xP][n][p] += (2.0 / Npulse)*Ez[mp[xP]][n][p] * cos(wk*q);
						Im_Ezs(vInd)[xP][n][p] += (2.0 / Npulse)*Ez[mp[xP]][n][p] * (-sin(wk*q));

						Re_C_Hys(vInd)[xP][n][p] += (2.0 / Npulse)*0.5*(Hy[mp[xP]][n][p] + Hy[mp[xP] - 1][n][p])*cos(wk*q);
						Im_C_Hys(vInd)[xP][n][p] += (2.0 / Npulse)*0.5*(Hy[mp[xP]][n][p] + Hy[mp[xP] - 1][n][p])*sin(wk*q);
					}// End calculation for Ezs and C_Hys

			}//End for each plane xP


			//Update phasor components on TFSF boundary
			for (n = 0; n < N - 1; n++)
				for (p = 0; p < P; p++)
				{
					Re_Eys(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.75*IFieldy.Eydir(mb) + 0.25*IFieldy.Eydir(mb - 1)) * cos(wk*q);
					Im_Eys(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.75*IFieldy.Eydir(mb) + 0.25*IFieldy.Eydir(mb - 1)) * (-sin(wk*q));

					Re_C_Hzs(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.25*IFieldy.Hzdir(mb) + 0.75*IFieldy.Hzdir(mb - 1))*cos(wk*q);
					Im_C_Hzs(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.25*IFieldy.Hzdir(mb) + 0.75*IFieldy.Hzdir(mb - 1))*sin(wk*q);

				}// End calculation for Eys and C_Hzs

			for (n = 0; n < N; n++)
				for (p = 0; p < P - 1; p++)
				{
					Re_Ezs(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.75*IFieldz.Ezdir(mb) + 0.25*IFieldz.Ezdir(mb - 1)) * cos(wk*q);
					Im_Ezs(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.75*IFieldz.Ezdir(mb) + 0.25*IFieldz.Ezdir(mb - 1)) * (-sin(wk*q));

					Re_C_Hys(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.25*IFieldz.Hydir(mb) + 0.75*IFieldz.Hydir(mb - 1))*cos(wk*q);
					Im_C_Hys(vInd)[Location_index_for_incident_power][n][p] += (2.0 / Npulse)*(0.25*IFieldz.Hydir(mb) + 0.75*IFieldz.Hydir(mb - 1))*sin(wk*q);
				}// End calculation for Ezs and C_Hys on TFSF boundary

		}//End for each frequency 
		*/


		parallel_for(0, numkVals, [&, this](size_t vInd) 
									{k = k0 + vInd*deltak; double wk = (2 * PI*k / Np); 
									parallel_for(0, Number_of_power_monitors - 1, [&, this, vInd](size_t xP) 
																					{parallel_for(0, N - 1, [&, this, vInd, xP](size_t n) 
																											{parallel_for(0, P, [&, this, vInd, xP, n](size_t p)
																																{Re_Eys(vInd)[xP][n][p] += (2.0 / Np)*Ey[mp[xP]][n][p] * cos(wk*q);
																																Im_Eys(vInd)[xP][n][p] += (2.0 / Np)*Ey[mp[xP]][n][p] * (-sin(wk*q));

																																Re_C_Hzs(vInd)[xP][n][p] 
																																	+= (2.0 / Np)*0.5*(Hz[mp[xP]][n][p] + Hz[mp[xP] - 1][n][p])*cos(wk*q);
																																Im_C_Hzs(vInd)[xP][n][p] 
																																	+= (2.0 / Np)*0.5*(Hz[mp[xP]][n][p] + Hz[mp[xP] - 1][n][p])*sin(wk*q); 
																																});// End calculation for Eys and C_Hzs
																											} );//End parallel_for over n

																					parallel_for(0, N, [&, this, vInd, xP](size_t n) 
																										{parallel_for(0, P - 1, [&, this, vInd, xP, n](size_t p) 
																																{Re_Ezs(vInd)[xP][n][p] += (2.0 / Np)*Ez[mp[xP]][n][p] * cos(wk*q);
																																Im_Ezs(vInd)[xP][n][p] += (2.0 / Np)*Ez[mp[xP]][n][p] * (-sin(wk*q));

																																Re_C_Hys(vInd)[xP][n][p] 
																																	+= (2.0 / Np)*0.5*(Hy[mp[xP]][n][p] + Hy[mp[xP] - 1][n][p])*cos(wk*q);
																																Im_C_Hys(vInd)[xP][n][p] 
																																	+= (2.0 / Np)*0.5*(Hy[mp[xP]][n][p] + Hy[mp[xP] - 1][n][p])*sin(wk*q); 
																																});// End calculation for Eys and C_Hzs 
																										});//End parallel_for over n 
																					});//End parallel_for over xP

									parallel_for(0, N - 1, [&, this, vInd](size_t n) 
															{parallel_for(0, P, [&, this, vInd, n](size_t p) 
																				{Re_Eys(vInd)[Location_index_for_incident_power][n][p]
																						+= (2.0 / Np)*(0.75*IFieldy.Eydir(mb) + 0.25*IFieldy.Eydir(mb - 1)) * cos(wk*q);
																				Im_Eys(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.75*IFieldy.Eydir(mb) + 0.25*IFieldy.Eydir(mb - 1)) * (-sin(wk*q));

																				Re_C_Hzs(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.25*IFieldy.Hzdir(mb) + 0.75*IFieldy.Hzdir(mb - 1))*cos(wk*q);
																				Im_C_Hzs(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.25*IFieldy.Hzdir(mb) + 0.75*IFieldy.Hzdir(mb - 1))*sin(wk*q); 
																				}// End calculation for Eys and C_Hzs
																		  ); 
															});//End parallel_for over n

									parallel_for(0, N, [&, this, vInd](size_t n)
															{parallel_for(0, P - 1, [&, this, vInd, n](size_t p)
																				{Re_Ezs(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.75*IFieldz.Ezdir(mb) + 0.25*IFieldz.Ezdir(mb - 1)) * cos(wk*q);
																				Im_Ezs(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.75*IFieldz.Ezdir(mb) + 0.25*IFieldz.Ezdir(mb - 1)) * (-sin(wk*q));

																				Re_C_Hys(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.25*IFieldz.Hydir(mb) + 0.75*IFieldz.Hydir(mb - 1))*cos(wk*q);
																				Im_C_Hys(vInd)[Location_index_for_incident_power][n][p] 
																						+= (2.0 / Np)*(0.25*IFieldz.Hydir(mb) + 0.75*IFieldz.Hydir(mb - 1))*sin(wk*q);
																				});// End calculation for Eys and C_Hzs
															});//End parallel_for over n

									});//End parallel_for over vInd

	}//End function runningSum_for_Component_Phasors



	void EMField3D::adjust_component_phasors_to_Yee_cell_face_centers(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys,
		 const vector<int>& mp)
	{
		int xP, n, p, k, vInd;

		for (vInd = 0; vInd < numkVals; vInd++)
		{
			for (xP = 0; xP < Number_of_power_monitors; xP++)
			{
				for (n = 0; n < N - 1; n++)
					for (p = 0; p < P - 1; p++)
					{
						Re_Ezs(vInd)[xP][n][p] = 0.5*(Re_Ezs(vInd)[xP][n][p] + Re_Ezs(vInd)[xP][n + 1][p]);
						Im_Ezs(vInd)[xP][n][p] = 0.5*(Im_Ezs(vInd)[xP][n][p] + Im_Ezs(vInd)[xP][n + 1][p]);

						Re_C_Hys(vInd)[xP][n][p] = 0.5*(Re_C_Hys(vInd)[xP][n][p] + Re_C_Hys(vInd)[xP][n + 1][p]);
						Im_C_Hys(vInd)[xP][n][p] = 0.5*(Im_C_Hys(vInd)[xP][n][p] + Im_C_Hys(vInd)[xP][n + 1][p]);


						Re_Eys(vInd)[xP][n][p] = 0.5*(Re_Eys(vInd)[xP][n][p] + Re_Eys(vInd)[xP][n][p + 1]);
						Im_Eys(vInd)[xP][n][p] = 0.5*(Im_Eys(vInd)[xP][n][p] + Im_Eys(vInd)[xP][n][p + 1]);

						Re_C_Hzs(vInd)[xP][n][p] = 0.5*(Re_C_Hzs(vInd)[xP][n][p] + Re_C_Hzs(vInd)[xP][n][p + 1]);
						Im_C_Hzs(vInd)[xP][n][p] = 0.5*(Im_C_Hzs(vInd)[xP][n][p] + Im_C_Hzs(vInd)[xP][n][p + 1]);

					}// End calculation for Ezs, C_Hys, Eys, C_Hzs

			}//End for each plane xP

		}//End for each frequency

	}//End function adjust_component_phasors_to_Yee_cell_face_centers



	void EMField3D::compute_spectral_power_densities(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys, Array1D_Vector3D& Pave_x,
		 const vector<int>& mp)
	{
		int xP, n, p, k, vInd;

		for (vInd = 0; vInd < numkVals; vInd++)
		{
			for (xP = 0; xP < Number_of_power_monitors; xP++)
			{
				for (n = 0; n < N - 1; n++)
					for (p = 0; p < P - 1; p++)
					{
						Pave_x(vInd)[xP][n][p] = 0.5*(Re_Eys(vInd)[xP][n][p] * Re_C_Hzs(vInd)[xP][n][p] - Im_Eys(vInd)[xP][n][p] * Im_C_Hzs(vInd)[xP][n][p])
							- 0.5*(Re_Ezs(vInd)[xP][n][p] * Re_C_Hys(vInd)[xP][n][p] - Im_Ezs(vInd)[xP][n][p] * Im_C_Hys(vInd)[xP][n][p]);

					}//End for each n,p

			}//End for each plane xP

		}//End for each frequency

	}//End function compute_spectral_power_densities



	void EMField3D::compute_spectral_power_flows(Array1D_Vector3D& Pave_x, Structure3DWithGrid& str, const vector<int>& mp)
	{
		int xP, n, p, k, vInd;

		for (vInd = 0; vInd < numkVals; vInd++)
		{
			for (xP = 0; xP < Number_of_power_monitors; xP++)
			{
				for (n = 0; n < N - 1; n++)
					for (p = 0; p < P - 1; p++)
					{
						Pflow(vInd)[xP] += Pave_x(vInd)[xP][n][p] * str.delAx(mp[xP], n, p);

					}//End for each n,p

				//Corect for units of grid increments being in micrometers while Pave_x elements being in W/m^2
				Pflow(vInd)[xP] = Pflow(vInd)[xP] * 1E-12;

				//Change reference direction for reflected power
				Pflow(vInd)[Location_index_for_reflected_power] = -Pflow(vInd)[Location_index_for_reflected_power];

			}//End for each plane xP

		}//End for each frequency

	}//End function compute_spectral_power_flows


	double EMField3D::time_average_incident_spectral_power_at_frequency_index(int k)
	{
		int vInd = (k - k0) / deltak; 
		return Pflow(vInd)[Location_index_for_incident_power];
	}

	double EMField3D::time_average_reflected_spectral_power_at_frequency_index(int k)
	{
		int vInd = (k - k0) / deltak;
		return Pflow(vInd)[Location_index_for_reflected_power];
	}

	double EMField3D::time_average_transmitted_spectral_power_at_frequency_index(int k)
	{
		int vInd = (k - k0) / deltak;
		return Pflow(vInd)[Location_index_for_transmitted_power];
	}

	double EMField3D::spectralPower(int vIndex, int locationIndex) { return Pflow(vIndex)[locationIndex]; }

	int EMField3D::TFSFboundaryEIndex() const { return mb; }

	int EMField3D::nEzNodePlanes_in_xDir() const { return M; }

	double EMField3D::Ez_comp(int mInd, int nInd, int pInd) const						//Returns Ez[mInd][nInd][pInd]
	{
		return Ez[mInd][nInd][pInd];
	}

	int EMField3D::numPowerMonitors() const { return Number_of_power_monitors; }




	//*************************************************************************************************************************************
	//*************************************************************************************************************************************
	//Functions for class EMField1D_zPolarized

	EMField1D_zPolarized::EMField1D_zPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt) :Ez(gridToUse->size()), Hy(gridToUse->size())
	{
		int m = 0;

		theGrid = gridToUse;
		N = theGrid->size();
		mb = 0;
		mD = theGrid->nNds(0) - 1; 

		Einc = IncidentField;

		q = 0;

		delt = aDelt;						//delt is in microseconds and delx in micrometers.
											//Initialize field

		for (m = 0; m < N; m++)
		{
			Ez[m] = 0;
			Hy[m] = 0;
		}//End for

		Ez[mb] = Einc(0, mb*theGrid->delx(0));

	}//End Constructor



	EMField1D_zPolarized::EMField1D_zPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse) :Ez(gridToUse->size()), Hy(gridToUse->size())
	{
		int m = 0;

		theGrid = gridToUse;
		N = theGrid->size();
		mb = int(theGrid->nNds(0) / 2 + 0.5);
		mD = theGrid->nNds(0) - 1;												// mD must be greater than mb. So improved version should check for this.

		Einc = IncidentField;

		delt = theGrid->minDelx() / c0;	//delt based on Courant number of unity. delt is in microseconds and delx in micrometers.

		q = 0;

		//Initialize field
		for (m = 0; m < N; m++)
		{
			Ez[m] = 0;
			Hy[m] = 0;
		}//End for

		Ez[mb] = Einc(0, mb*theGrid->delx(0));

	}//End Constructor



	int EMField1D_zPolarized::propagateIn(Structure1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1D_zPolarized*, int timeIndex), int displayTime)
	{
		double loss;												//Need code to protect against case of mb == 0
		double Sc;
		int i;
		int m;
		int mCount;
		int mCountL;
		int nL;

		//Old:	//double* CEzE = new double[N];
		//double* CEzH = new double[N];
		//double* CHyE = new double[N];
		//double* CHyH = new double[N];
		vector<double> CEzE(N);
		vector<double> CEzH(N);
		vector<double> CHyE(N);
		vector<double> CHyH(N);


		//*******************************************************************************************************************
		//Compute coeficients for grid layer above device

		mCount = 0;

		for (m = 0; m < theGrid->nNds(0); m++)
		{
			CEzE[m] = 1;
			CEzH[m] = Imp0;    //delt / eps0*theGrid->delx(0);
			CHyE[m] = 1 / Imp0;    //delt / (mu0*theGrid->delx(0));
			CHyH[m] = 1;

			mCount += 1;
		}//End for

		 //Correction for m =  nNodes(0) - 1 in case of multilayer grid because then
		 //(nNodes(0) - 1) + 1/2 is location of top surface of device.
		if (theGrid->numL() > 1)
			CHyE[m - 1] = CHyE[m - 1] * theGrid->delx(0) / (theGrid->delx(1) / 2 + theGrid->delx(0) / 2);		//Not needed if deltx is the same for first device layer as for layer above device


		//Compute coefficients for each device layer. Grid has one layer below device
		for (i = 1; i < theGrid->numL() - 1; i++)
		{
			mCountL = mCount;
			for (m = mCount; m < mCountL + theGrid->nNds(i); m++)
			{
				loss = dev.sigma(i - 1)*(delt*1E-6) / (2 * dev.ep(i - 1));		//Note that in calculating the quantity loss, time step is converted to seconds
				CEzE[m] = (1 - loss) / (1 + loss);								//CEzE[m] is dimensionless
				CEzH[m] = delt / (dev.ep(i - 1)*theGrid->delx(i)) / (1 + loss);	//CEzH[m] is in mks units
				CHyE[m] = delt / (dev.mu(i - 1)*theGrid->delx(i));				//CHyE[m] is in mks units
				CHyH[m] = 1;

				mCount += 1;
			}//End for each node in layer i

			 //correction for magnetic field boundary node except last one
			if (i < theGrid->numL() - 1)
			{
				CHyE[m - 1] = CHyE[m - 1] * theGrid->delx(i) / (theGrid->delx(i + 1) / 2 + theGrid->delx(i) / 2);
			}

			//For m = N-1, CHyE[m] not used
		}//End for each layer i


		 //Compute coefficients for vacuum layer below device.

		nL = theGrid->numL();
		Sc = c0*delt / theGrid->delx(nL - 1);
		mCountL = mCount;
		for (m = mCount; m < mCountL + theGrid->nNds(nL - 1); m++)
		{
			CEzE[m] = 1;
			CEzH[m] = Imp0*Sc;
			CHyH[m] = 1;
			CHyE[m] = Sc / Imp0;

			mCount += 1;
		}//End for each node in layer below the device


		 //*******************************************************************************************************************
		 //Time stepping

		for (q = 0; q < maxTime; q++)
		{
			//Simple ABC for Hy[N-1]
			Hy[N - 1] = Hy[N - 2];   //This should work provided the bottom layer is lossless and local Sc = 1. See Schneider

									 //Update magnetic field values
			for (m = 0; m < N - 1; m++)
				Hy[m] = CHyH[m] * Hy[m] + CHyE[m] * (Ez[m + 1] - Ez[m]);

			//Correction for Hy[mb-1]
			Hy[mb - 1] = Hy[mb - 1] - CHyE[mb - 1] * Einc(q*delt, mb*theGrid->delx(0));

			//Simple ABC for Ez[0]
			Ez[0] = Ez[1];

			//Update electric field values
			for (m = 1; m < N; m++)
				Ez[m] = CEzE[m] * Ez[m] + CEzH[m] * (Hy[m] - Hy[m - 1]);

			//Correction for Ez[mb]
			Ez[mb] = Ez[mb] + CEzH[mb] * (1 / Imp0)*Einc((q + 1.0 / 2.0)*delt, (mb - 1.0 / 2.0)*theGrid->delx(0));

			if (q % displayTime == 0)
			{
				EMField1D_zPolarized* PointerToMe = this;
				displayFunctionToCall(PointerToMe, q);
			}


		}//End Time stepping


		return 0;
	}//End propagateIn


	void EMField1D_zPolarized::compute_coefficients_for_medium_with_lossy_layer_underneath(Medium& am, vector<double>& CEzE, vector<double>& CEzH, vector<double>& CHyE, vector<double>& CHyH)
	{
		double loss;
		int m;
		int mCount;
		int mCountL;
		int nL;

			mCount = 0;
			for (m = 0; m < theGrid->nNds(0); m++)
			{
				loss = am.sigma()*(delt*1E-6) / (2 * am.ep());			//Note that here, time step is converted to seconds.

				CEzE[m] = (1 - loss) / (1 + loss);							//CEzE[m] is dimensionless
				CEzH[m] = delt / (am.ep()*theGrid->delx(0)) / (1 + loss);	//CEzH[m] is in mks units
				CHyE[m] = delt / (am.mu()*theGrid->delx(0));				//CHyE[m] is in mks units
				CHyH[m] = 1;

				mCount += 1;
			}//End for

			 //Compute coefficients for lossy layer below medium to terminate grid.
			 //Lossy layer with impedance matched to medium is intended to prevent reflection at interface and
			 //minimize reflection from bottom (right side) of grid.
			 //Should work if the medium is lossless (conductivity = 0). Otherwise need to look into case of medium is itself lossy 1/24/16.
			 //See Sec. 3.11 and 3.12 of Schneider.

			nL = theGrid->numL();
			mCountL = mCount;
			for (m = mCount; m < mCountL + theGrid->nNds(nL - 1); m++)
			{
				if (loss == 0)
					loss = 0.02*(m - mCount)/(double(mCountL + theGrid->nNds(nL - 1)) - 1.0);				// Max loss for the lossy bottom layer is the same as for the medium or the last
																//layer of the device unless it is zero. 
																//If the loss for the medium is zero (e.g. vaccum) loss for the lossy bottom layer set to value
																//used by Schneider for lossy layer terminating grid in program  3.8.
																//Note that the permitivity and permeability are the same as for the medium above.

				CEzE[m] = (1 - loss) / (1 + loss);									//CEzE[m] is dimensionless
				CEzH[m] = delt / (am.ep()* theGrid->delx(nL - 1)) / (1 + loss);		//CEzH[m] is in mks units
				CHyH[m] = (1 - loss) / (1 + loss);
				CHyE[m] = delt / (am.mu()*theGrid->delx(nL - 1)) / (1 + loss);						//CHyE[m] is in mks units

				mCount += 1;
			}//End for each node in lossy layer below the medium

	}//End function


	int EMField1D_zPolarized::Update(vector<double>& CEzE, vector<double>& CEzH, vector<double>& CHyE, vector<double>& CHyH)
	{

		double loss;
		int m;
		int mCount;
		int mCountL;
		int nL;


		if (theGrid->numL() == 2 && mb == 0)
		{
			 //*******************************************************************************************************************
			 //Update fields

			 //Update magnetic field values
			for (m = 0; m < N - 1; m++)
				Hy[m] = CHyH[m] * Hy[m] + CHyE[m] * (Ez[m + 1] - Ez[m]);

			Ez[0] = Einc(q*delt, 0.0*theGrid->delx(0));

			//Update electric field values
			for (m = 1; m < N; m++)
				Ez[m] = CEzE[m] * Ez[m] + CEzH[m] * (Hy[m] - Hy[m - 1]);

			q += 1;

			return 0;
		}
		else
		{
			//Error condition: this update function should not have been called
			return 1;
		}//End if

		
	}//End function Update of class EMField1D_zPolarized


	 //Accessors

	double EMField1D_zPolarized::Ezdir(int m) const { return Ez[m]; }		//Returns Ez[m]
	double  EMField1D_zPolarized::Hydir(int m) const { return Hy[m]; }		//Returns Hy[m]
	int EMField1D_zPolarized::nENodes() const { return N; }					//Returns the total number of electric field nodes in the grid
	int EMField1D_zPolarized::TFSFboundaryEIndex() const { return mb; }



	//***************************************************************************************************************************
	//Functions for EMField1D_yPolarized

	EMField1D_yPolarized::EMField1D_yPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt) :Ey(gridToUse->size()), Hz(gridToUse->size())
	{
		int m = 0;

		theGrid = gridToUse;
		N = theGrid->size();
		mb = 0;
		mD = theGrid->nNds(0) - 1;		//Change this.

		Einc = IncidentField;

		q = 0;


		delt = aDelt;								//delt is in microseconds and delx in micrometers.

													//Initialize field

		for (m = 0; m < N; m++)
		{
			Ey[m] = 0;
			Hz[m] = 0;
		}//End for

		Ey[mb] = Einc(0, mb*theGrid->delx(0));

	}//End Constructor EMField1D_yPolarized



	EMField1D_yPolarized::EMField1D_yPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse) :Ey(gridToUse->size()), Hz(gridToUse->size())
	{
		int m = 0;

		theGrid = gridToUse;
		N = theGrid->size();
		mb = int(theGrid->nNds(0) / 2 + 0.5);
		mD = theGrid->nNds(0) - 1;												// mD must be greater than mb. So improved version should check for this.

		Einc = IncidentField;

		q = 0;
		delt = theGrid->minDelx() / c0;	//delt based on Courant number of unity. delt is in microseconds and delx in micrometers. 

										//Initialize field
		for (m = 0; m < N; m++)
		{
			Ey[m] = 0;
			Hz[m] = 0;
		}//End for

		Ey[mb] = Einc(0, mb*theGrid->delx(0));

	}//End Constructor EMField1D_yPolarized



	int EMField1D_yPolarized::propagateIn(Structure1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1D_yPolarized*, int timeIndex), int displayTime)
	{
		double loss;												//Need code to protect against case of mb == 0
		double Sc = 1;
		int i;
		int m;
		int mCount;
		int mCountL;
		int nL;

		//Old:	//double* CEzE = new double[N];
		//double* CEzH = new double[N];
		//double* CHyE = new double[N];
		//double* CHyH = new double[N];
		vector<double> CEyE(N);
		vector<double> CEyH(N);
		vector<double> CHzE(N);
		vector<double> CHzH(N);

		//*******************************************************************************************************************
		//Compute coeficients for grid layer above device

		mCount = 0;

		for (m = 0; m < theGrid->nNds(0); m++)
		{
			CEyE[m] = 1.0;
			CEyH[m] = Imp0*Sc;    //delt / eps0*theGrid->delx(0);
			CHzE[m] = Sc / Imp0;    //delt / (mu0*theGrid->delx(0));
			CHzH[m] = 1.0;

			mCount += 1;
		}//End for

		 //Correction for m =  nNodes(0) - 1 in case of multilayer grid because then 
		 //(nNodes(0) - 1) + 1/2 is location of top surface of device.
		if (theGrid->numL() > 1)
			CHzE[m - 1] = CHzE[m - 1] * theGrid->delx(0) / (theGrid->delx(1) / 2 + theGrid->delx(0) / 2);		//Not needed if deltx is the same for first device layer as for layer above device


																												//Compute coefficients for each device layer. Grid has one layer below device

		for (i = 1; i < theGrid->numL() - 1; i++)
		{
			mCountL = mCount;
			for (m = mCount; m < mCountL + theGrid->nNds(i); m++)
			{
				loss = dev.sigma(i - 1)*(delt*1E-6) / (2 * dev.ep(i - 1));		//Note that in calculating the quantity loss, time step is converted to seconds
				CEyE[m] = (1 - loss) / (1 + loss);								//CEyE[m] is dimensionless
				CEyH[m] = delt / (dev.ep(i - 1)*theGrid->delx(i)) / (1 + loss);	//CEzH[m] is in mks units
				CHzE[m] = delt / (dev.mu(i - 1)*theGrid->delx(i));				//CHyE[m] is in mks units
				CHzH[m] = 1;

				mCount += 1;
			}//End for each node in layer i

			 //correction for magnetic field boundary node except last one
			if (i < theGrid->numL() - 1)
			{
				CHzE[m - 1] = CHzE[m - 1] * theGrid->delx(i) / (theGrid->delx(i + 1) / 2 + theGrid->delx(i) / 2);
			}

			//For m = N-1, CHyE[m] not used
		}//End for each layer i


		 //Compute coefficients for vacuum layer below device.

		nL = theGrid->numL();
		Sc = c0*delt / theGrid->delx(nL - 1);
		mCountL = mCount;
		for (m = mCount; m < mCountL + theGrid->nNds(nL - 1); m++)
		{
			CEyE[m] = 1;
			CEyH[m] = Imp0*Sc;
			CHzH[m] = 1;
			CHzE[m] = Sc / Imp0;

			mCount += 1;
		}//End for each node in layer below the device


		 //*******************************************************************************************************************
		 //Time stepping

		for (q = 0; q < maxTime; q++)
		{
			//Simple ABC for Hz[N-1]
			Hz[N - 1] = Hz[N - 2];   //This should work provided the bottom layer is lossless and local Sc = 1. See Schneider

									 //Update magnetic field values
			for (m = 0; m < N - 1; m++)
				Hz[m] = CHzH[m] * Hz[m] - CHzE[m] * (Ey[m + 1] - Ey[m]);

			//Correction for Hz[mb-1]
			Hz[mb - 1] = Hz[mb - 1] + CHzE[mb - 1] * Einc(q*delt, mb*theGrid->delx(0));

			//Simple ABC for Ey[0]
			Ey[0] = Ey[1];

			//Update electric field values
			for (m = 1; m < N; m++)
				Ey[m] = CEyE[m] * Ey[m] - CEyH[m] * (Hz[m] - Hz[m - 1]);

			//Correction for Ey[mb]
			Ey[mb] = Ey[mb] - CEyH[mb] * (1 / Imp0)*Einc((q + 1.0 / 2.0)*delt, (mb - 1.0 / 2.0)*theGrid->delx(0));	//Check that Hzinc is -(1/Imp0)*Einc

			if (q % displayTime == 0)
			{
				EMField1D_yPolarized* PointerToMe = this;
				displayFunctionToCall(PointerToMe, q);
			}


		}//End Time stepping


		return 0;
	}//End propagateIn of class EMField1D_yPolarized


	void EMField1D_yPolarized::compute_coefficients_for_medium_with_lossy_layer_underneath(Medium& am, vector<double>& CEyE, vector<double>& CEyH, vector<double>& CHzE, vector<double>& CHzH)
	{
		double loss;
		int m;
		int mCount;
		int mCountL;
		int nL;

		//Compute coeficients for top layer of two-layer grid

		if (q == 0)
		{
			mCount = 0;
			for (m = 0; m < theGrid->nNds(0); m++)
			{
				loss = am.sigma()*(delt*1E-6) / (2 * am.ep());			//Note that here, time step is converted to seconds.

				CEyE[m] = (1 - loss) / (1 + loss);							//CEzE[m] is dimensionless
				CEyH[m] = delt / (am.ep()*theGrid->delx(0)) / (1 + loss);	//CEzH[m] is in mks units
				CHzE[m] = delt / (am.mu()*theGrid->delx(0));				//CHyE[m] is in mks units
				CHzH[m] = 1;

				mCount += 1;
			}//End for

			 //Correction for m =  nNodes(0) - 1 
			 //(nNodes(0) - 1) + 1/2 is location of interface with lossy layer.
			if (theGrid->numL() > 1)
				CHzE[m - 1] = CHzE[m - 1] * theGrid->delx(0) / (theGrid->delx(1) / 2 + theGrid->delx(0) / 2);		//Not needed if deltx is the same for both layers


																													//Compute coefficients for lossy layer below medium to terminate grid.
																													//Lossy layer with impedance matched to medium is intended to prevent reflection at interface and
																													//minimize reflection form bottom (right side) of grid.
																													//Should work if the medium is lossless (conductivity = 0). Otherwise need to look into case of medium is itself lossy 1/24/16.
																													//See Sec. 3.11 and 3.12 of Schneider.
			nL = theGrid->numL();
			mCountL = mCount;
			for (m = mCount; m < mCountL + theGrid->nNds(nL - 1); m++)
			{
				if (loss == 0)
					loss = 0.02;								//loss for the lossy bottom layer is the same as for the last 
																//layer of the device unless it is zero, and if zero set to value 
																//used by Schneider for lossy layer terminating grid in program  3.8
																//Note that the pernitivity and permeability are the same as for the medium above

				CEyE[m] = (1 - loss) / (1 + loss);									//CEyE[m] is dimensionless
				CEyH[m] = delt / (am.ep()* theGrid->delx(nL - 1)) / (1 + loss);		//CEyH[m] is in mks units
				CHzH[m] = (1 - loss) / (1 + loss);									//Check this. Why not CHzH = 1?
				CHzE[m] = delt / (am.mu()*theGrid->delx(nL - 1)) / (1 + loss);						//CHzE[m] is in mks units

				mCount += 1;
			}//End for each node in lossy layer below the medium

		}//End if

	}//End function compute_coefficients_for_medium_with_lossy_layer_underneath of class EMField1D_yPolarized


	int EMField1D_yPolarized::Update(vector<double>& CEyE, vector<double>& CEyH, vector<double>& CHzE, vector<double>& CHzH)
	{

		double loss;
		int m;
		int mCount;
		int mCountL;
		int nL;

		if (theGrid->numL() == 2 && mb == 0)
		{
			 //Update fields

			 //Update magnetic field values
			for (m = 0; m < N - 1; m++)
				Hz[m] = CHzH[m] * Hz[m] - CHzE[m] * (Ey[m + 1] - Ey[m]);

			
			Ey[0] = Einc(q*delt, 0.0*theGrid->delx(0));

			//Update electric field values
			for (m = 1; m < N; m++)
				Ey[m] = CEyE[m] * Ey[m] - CEyH[m] * (Hz[m] - Hz[m - 1]);

			q += 1;

			return 0;
		}
		else
		{
			//Error condition: this update function should not havebeen called
			return 1;
		}
		
	}//End function Update of class EMField1D_yPolarized 


	 //*******************************************************************************************************************
	 //Accessors

	double EMField1D_yPolarized::Eydir(int m) const { return Ey[m]; }		//Returns Ez[m]
	double  EMField1D_yPolarized::Hzdir(int m) const { return Hz[m]; }		//Returns Hy[m]
	int EMField1D_yPolarized::nENodes() const { return N; }				//Returns the total number of electric field nodes in the grid
	int EMField1D_yPolarized::TFSFboundaryEIndex() const { return mb; }



}//End namespace FDTD
