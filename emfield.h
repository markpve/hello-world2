#ifndef _EMFIELD_H
#define _EMFIELD_H


#include "grid.h"

#include <vector>

//#include <windows.h>
//#include <ppl.h>


//using namespace concurrency;

using namespace std;

namespace FDTD
{

	//***************************************************************************************************************************************************
	//Class representing one-dimensional electromagmetic field with electric field in z-direction and propagation in x-direction

	class EMField1DzPolar
	{
	public:

		EMField1DzPolar() {}

		EMField1DzPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt);				//Field with incident field assigned at mb = 0 at time q = 0
																													//delt assigned aDelt
																													//delt is to be in microseconds and delx in micrometers.

		EMField1DzPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse);
		//Field with incident field assigned at a TFSF boundary between
		//Ez electric field node mb = int(theGrid->nNds(0)/2 + 0.5) and mb-1/2
		//in top layer of grid at q = 0.
		//The grid provided must be based on a device in order to subsequently use the prodageIn function.
		//delt based on Courant number of unity and minimum delx for grid.
		//delt is in microseconds and delx in micrometers.

		int propagateIn(Device1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1DzPolar*, int timeIndex), int displayTime);
		//Simulates propagation of field in 1D device described by first argument.
		//Implements TFSF boundary above device with assumption that mb is as assigned in second constructor (.0).
		//The dev argument must fit the grid used by *this.
		// Calls displayFunctionToCall after update at every displayTime time steps.



		int Update(Medium& aM, double CEzE[], double  CEzH[], double CHyE[], double CHyH[], int displayCode);
		//Update fields in homogeneous medium with lossy layer behind it for one time step.
		//Calculates coefficients at q = 0 only and stores them in array arguments.
		//Grid used by *this must be two-layer with bottom layer intended for
		//lossy layer to terminate grid.
		//medium is in the index zero layer. Implements lossy layer in second grid layer.


		double Ezdir(int m) const;					//Returns Ez[m]
		double Hydir(int m) const;					//Returns Hy[m]
		int nENodes() const;						//Returns the total number of electric field nodes in the grid
		int TFSFboundaryEIndex() const;			//Returns index of electric field node to the right of and nearest TFSF boundary

	private:

		vector<double> Ez;
		vector<double> Hy;

		int q;	//Time in units of delt
		double delt;	//Units are microseconds

		double(*Einc)(double t, double x);	//Pointer to function representing incident field

		int N;   //Number of Ez electric field nodes in grid. Equals number of Hy magnetic field nodes
				 //since this code assumes grid starts on left with Ez node and ends on right with Hy node.

		int mb;  //Index of Ez electric field node to right of SF/TF boundary or index of node with specified source
		int mD;  //Top surfaace of device (if any) is at magnetic field node mD + 1/2

		Grid1D* theGrid;	//Pointer to grid shared with device

	};
	

	//*****************************************************************************************************************************************************************
	//Class representing one-dimensional electromagmetic field with electric field in y-direction, magnetic field in z-direction and propagation in x-direction

	class EMField1DyPolar
	{
	public:
		EMField1DyPolar() {}

		EMField1DyPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt);			//Field with incident field assigned at mb = 0 at time q = 0
																												//delt assigned aDelt
																												//delt is to be in microseconds and delx in micrometers.

		EMField1DyPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse);			//Field with incident field assigned at a TFSF boundary between 
																								//Ey electric field node mb = int(theGrid->nNds(0)/2 + 0.5) and mb-1/2
																								//in top layer of grid at q = 0. 
																								//The grid provided must be based on a device in order to subsequently use the prodageIn function.  
																								//delt based on Courant number of unity and minimum delx for grid.
																								//delt is in microseconds and delx in micrometers.

		int propagateIn(Device1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1DyPolar*, int timeIndex), int displayTime);
		//Simulates propagation of field in 1D device described by first argument.
		//Implements TFSF boundary above device with assumption that mb is as assigned in second constructor (.0). 
		//The dev argument must fit the grid used by *this. 
		//Calls displayFunctionToCall after update at every displayTime time steps.



		int Update(Medium& aM, double CEzE[], double  CEzH[], double CHyE[], double CHyH[], int displayCode);			//Update fields in homogeneous medium with lossy layer behind it for one time step. 
																														//Calculates coefficients at q = 0 only and stores them in array arguments.
																														//Grid used by *this must be two-layer with bottom layer intended for 
																														//lossy layer to terminate grid.
																														//medium is in the index zero layer. Implements lossy layer in second grid layer.


		double Eydir(int m) const;					//Returns Ey[m]
		double Hzdir(int m) const;					//Returns Hz[m]
		int nENodes() const;					//Returns the total number of Ey  nodes in the grid
		int TFSFboundaryEIndex() const;			//Returns index of electric field node to the right of and nearest TFSF boundary

	private:

		vector<double> Ey;
		vector<double> Hz;

		int q;	//Time in units of delt
		double delt;	//Units are microseconds

		double(*Einc)(double t, double x);	//Pointer to function representing incident field

		int N;   //Number of Ey electric field nodes in grid. Equals number of Hz magnetic field nodes 
				 //since this code assumes grid starts on left with Ey node and ends on right with Hz node.

		int mb;  //Index of Ey electric field node to right of SF/TF boundary or index of node with specified source
		int mD;  //Top surfaace of device (if any) is at magnetic field node mD + 1/2

		Grid1D* theGrid;	//Pointer to grid shared with device

	};


	//***********************************************************************************************************************************
	//Class representing Three-dimensional electromagmetic field. 


class EMField3D
{
public:
	EMField3D() {}
	EMField3D(double(*IncidentEzField)(double t, double x), double(*IncidentEyField)(double t, double x), Grid3D* gridToUse);

	//int propagateIn(Device1D& dev, int maxTime, int displayCode);
	int propagateIn(Structure3DWithGrid& str, int maxTime, int k0, int deltak, int numkVals, int displayCode);
						//Simulate propagation of electomagnetic wave in unit cell of structure that is periodic in y and z-directions
						//maxTime is the number of time steps. k0 is lowest frequency index for which spectral power flow is to be calculated
						//deltak is change in succcessive frequencies for which spectral power flow is to be calculated
						//numkVals is the number of frequency values


private:
	vector <double> FieldInit;

	vector <vector <vector<double>>> Ex;		//Electric field
	vector <vector <vector<double>>> Ey;		
	vector <vector <vector<double>>> Ez;
	vector <vector <vector<double>>> Hx;		//Magnetic field 
	vector <vector <vector<double>>> Hy;
	vector <vector <vector<double>>> Hz;

	int q;	//Time in units of delt
	double delt;

	double(*Eincz)(double t, double x);	//Pointer to function representing incident Ez field propagating in x-direction.
	double(*Eincy)(double t, double x);	//Pointer to function representing incident Ey field propagating in x-direction.

	int M;		//M is the number of Ez nodes across 3D grid in x - direction. Also equals number of Hx maagnetic field nodes
				//and number of Ey nodes across 3D grid in x-direction since this code assumes planes of Ez, Hx, and Hy nodes on x-boundaries.
					
				//Number of Hz nodes across grid in x-direction equals M - 1. Number of Ex nodes across grid in x-direction equals M - 1.
				//Number of Hy magnetic firld nodes across grid in x-direction equals M - 1.

	int N;		//N is the number of Ez nodes across 3D grid in y-direction. Also equals the number of Ex nodes and number of Hy nodes 
				//across 3D grid in y-direction since this code assumes planes of Ez, Hy and ex nodes on y boundaries.
				
				//Number of Hx magnetic field nodes across grid in y-direction is N - 1. 
				//Number of Hz magnetic field nodes across grid in y-direction is N - 1.
				//Number of Ey nodes across grid in y-direction is N - 1.

	int P;		//Number of Ex nodes across grid in z-direction. Also equals number of Ey nodes and number of Hz nodes across grid in z-direction 
				//since this code assumes planes of Ex, Ey, and Hz nodes on z-boundaries.


	int mb;  //x-Index of Ez field nodes to right of TF/SF boundary
	int mD;  //Top surfaace of structure is at Hy magnetic field nodes with x-Index mD + 1/2

	Grid3D* theGrid;	//Pointer to grid shared with device

	const double Sc = 1 / sqrt(3.0);	//Choose delt = Sc*delMin/c0 so that local Courant number will not exceed 1/sqrt(3.0).

	double ScExterior;			// c0*delt/delx (= Sc for layer above and below device. Assumes delx the same above and below the device)

	int updateMagneticField(Structure3DWithGrid& str, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy);

	int updateMagneticField_outsideBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
															vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ);			
																	//Used only if the structure is the unit cell of a larger periodic structure 
																	//This instance for normal incidence only

	int updateElectricField(Structure3DWithGrid& str, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy);

	int updateElectricField_onBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
															vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ);
																	//Used only if the structure is the unit cell of a larger periodic structure 
																	//This instance for normal incidence only

	void correctMagneticField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy);

	void correctElectricField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy);

	void computeCoefficientsForSecondOrderABCs(double& d, double& C0, double& C1, double& C2);

	void applySecondOrderABC_on_xBoundaries(vector<vector<vector<double>>> & EzOldLeft0, vector<vector<vector<double>>>& EzOldLeft1,
							vector<vector<vector<double>>>& EzOldRight0, vector<vector<vector<double>>> EzOldRight1,
							vector<vector<vector<double>>> & EyOldLeft0, vector<vector<vector<double>>>& EyOldLeft1,
							vector<vector<vector<double>>>& EyOldRight0, vector<vector<vector<double>>> EyOldRight1,
							double d, double C0, double C1, double C2);

	void runningSum_for_Component_Phasors(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs, 
											Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys, 
											int k0, int deltak, int numkVals, const vector<int>& mp, const int NT);

	void adjust_component_phasors_to_Yee_cell_face_centers(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys,
		int k0, int deltak, int numkVals, const vector<int>& mp, const int NT);
										//After this function Re_Ezs(k)[xP][n][p] etc. will be value at the center of the x=const face of Yee cell with indexes mp[xp], n, p.
										//The actual mathematical index of the face center is mp[xp], n+1/2, p+1/2.

	void compute_spectral_power_densities(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys, Array1D_Vector3D& Pave_x,
		int k0, int deltak, int numkVals, const vector<int>& mp, const int NT);

	void compute_spectral_power_flows(Array1D_Vector3D Pave_x, vector<double> P, Structure3DWithGrid& str, int k0, int deltak, int numkVals);

};



EMField3D::EMField3D(double(*IncidentEzField)(double t, double x), double(*IncidentEyField)(double t, double x), Grid3D* gridToUse) :
	M(gridToUse->sizeX()), N(gridToUse->sizeY()), P(gridToUse->sizeZ()),
	Ex(M - 1, vector<vector<double>>(N, vector<double>(P, 0.0))),
	Ey(M, vector<vector<double>>(N - 1, vector<double>(P, 0.0))),
	Ez(M, vector<vector<double>>(N, vector<double>(P - 1, 0.0))),
	Hx(M, vector<vector<double>>(N - 1, vector<double>(P - 1, 0.0))),
	Hy(M - 1, vector<vector<double>>(N, vector<double>(P - 1, 0.0))),
	Hz(M - 1, vector<vector<double>>(N - 1, vector<double>(P, 0.0)))

{
	int i = 0; int n = 0; int p = 0;

	theGrid = gridToUse;

	mb = int(theGrid->nNdsXL(0) / 2 + 0.5);	//Note that nNds(0) must be at least two in order for a SF/TF boundary to be applied.
	mD = theGrid->nNdsXL(0) - 1;				//Should just assign value of mD from a Grid3D accessor. Must have mD >= mb + 2

	Eincz = IncidentEzField;
	Eincy = IncidentEyField;

	q = 0;
	delt = theGrid->minDelta()*Sc / c0;	//delt is in microseconds and delx in micrometers.

	ScExterior = c0*delt / theGrid->deltaxL(0);

	//Initialize plane  of Ez and Ey nodes with x-index mb to incident field values
	for (n = 0; n < N; n++)
	{
		for (p = 0; p < P - 1; p++)
		{
			Ez[mb][n][p] = Eincz(0, mb*theGrid->deltaxL(0));
			Ey[mb][n][p] = Eincy(0, mb*theGrid->deltaxL(0));
		}//End for

		Ey[mb][n][P - 1] = Eincy(0, mb*theGrid->deltaxL(0));
	}//End for


}//End Constructor



int EMField3D::propagateIn(Structure3DWithGrid& str, int maxTime, int k0, int deltak, int numkVals, int displayCode)
{
	int i;
	int j;
	int k;
	int m;
	int n;
	int p;
	int NT = maxTime;

	//**************************************************
	//objects needed for time-average power calculations 

	const int num_xPlanes = 2;		//Number of x = const planes for which time-average power calculated in addition to incident power

	vector <int> xIndexes_for_power_calculation(2, 0);		//Vector to store locations of x=const planes for which reflected power and transmitted power are to be calculated.

	xIndexes_for_power_calculation[0] = mb - 2;
	xIndexes_for_power_calculation[1] = str.grid()->numEzPlanesthrough(str.grid()->numLx() - 2);

	//Arrays for real and imaginary parts of component phasors needed for time-average power calculations

	Array1D_Vector3D Re_Eys(numkVals, 2, N - 1, P);
	Array1D_Vector3D Im_Eys(numkVals, 2, N - 1, P);
	Array1D_Vector3D Re_C_Hzs(numkVals, 2, N - 1, P);
	Array1D_Vector3D Im_C_Hzs(numkVals, 2, N - 1, P);

	Array1D_Vector3D Re_Ezs(numkVals, 2, N, P - 1);
	Array1D_Vector3D Im_Ezs(numkVals, 2, N, P - 1);
	Array1D_Vector3D Re_C_Hys(numkVals, 2, N, P - 1);
	Array1D_Vector3D Im_C_Hys(numkVals, 2, N, P - 1);

	Array1D_Vector3D Pave_x(numkVals, 2, N, P);

	vector<double> Pxtot(3, 0.0);


	//*************************************************************************
	//Objects needed for analytical absorbing boundary conditions at front and rear x = const boundaries

	double C0;		//Coefficients for analytical absorbing boundary conditions
	double C1;
	double C2;
	double d;


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


	//*******************************************************************************
	//objects needed for auxiliary one-dimensional incident field simulation to be used for corrections along the x=const TFSF boundary

	//Instantiate a two-layer Grid object. The first layer is for a vacuum and the second is for a lossy layer for grid termination.
	Grid1D auxGrid(mD, theGrid->deltaxL(0));

	//To simulate the incident field, instantiate an EMField1DzPolar using the auxiliary two-layer Grid1D .
	EMField1DzPolar IFieldz(Eincz, &auxGrid, delt);

	//To simulate the incident field, instantiate an EMField1DyPolar using the auxiliary two-layer Grid1D .
	EMField1DyPolar IFieldy(Eincy, &auxGrid, delt);

	//Instantiate a medium object reprersenting vaccuum for the auxiliary 1D fields
	Medium aM(0, 1, 1);

	//Allocate memory for 1D coefficient arrays
	double* CEzE1D = new double[mD];
	double* CEzH1D = new double[mD];
	double* CHyE1D = new double[mD];
	double* CHyH1D = new double[mD];

	double* CEyE1D = new double[mD];
	double* CEyH1D = new double[mD];
	double* CHzE1D = new double[mD];
	double* CHzH1D = new double[mD];


	//*******************************************************************************************
	//Objects needed for boundary conditions for periodicity in the y and z directions

	//Allocate memory for magnetic field components just outside the high y-boundary and just outside the high z-boundary
	//Used only if the structure is the unit cell of a larger periodic structure

	vector <vector<double>> HxOutsideY;			//Hx at nodes jsut outside high y-boundary. 
												//Used only if the structure is the unit cell of a larger structure periodic in the y-direction 
	vector <vector<double>> HzOutsideY;			//Hz at nodes jsut outside high y-boundary. 
												//Used only if the structure is the unit cell of a larger structure periodic in the y-direction 

	vector <vector<double>> HxOutsideZ;			//Hx at nodes jsut outside high z-boundary. 
												//Used only if the structure is the unit cell of a larger structure periodic in the z-direction 

	vector <vector<double>> HyOutsideZ;			//Hy at nodes jsut outside high z-boundary. 
												//Used only if the structure is the unit cell of a larger structure periodic in the z-direction

	


	//***************************************************************************************************************************
	//Time stepping

	for (q = 0; q < maxTime; q++)
	{

		updateMagneticField(str, IFieldz, IFieldy);

		correctMagneticField_alongTFSFxFrontBoundary(str, IFieldz, IFieldy);

		//for normal incidence
		if(str.Is_a_unitCell())
			updateMagneticField_outsideBoundaries_ofGridUnitCell(str, HxOutsideY,  HzOutsideY, HxOutsideZ,  HyOutsideZ);

		//Update auxilliary one-dimensional incident fields to use for correction to nodes adjacent to TF/SF boundary. 
		IFieldz.Update(aM, CEzE1D, CEzH1D, CHyE1D, CHyH1D, displayCode);
		IFieldy.Update(aM, CEyE1D, CEyH1D, CHzE1D, CHzH1D, displayCode);

		updateElectricField(str, IFieldz, IFieldy);

		//For normal incidence
		if (str.Is_a_unitCell())
			updateElectricField_onBoundaries_ofGridUnitCell(str, HxOutsideY, HzOutsideY, HxOutsideZ, HyOutsideZ);		//If the structure is not periodic in both y and z directions, need more code.

		correctElectricField_alongTFSFxFrontBoundary(str, IFieldz, IFieldy);

		applySecondOrderABC_on_xBoundaries(EzOldLeft0, EzOldLeft1, EzOldRight0, EzOldRight1, EyOldLeft0, EyOldLeft1, EyOldRight0, EyOldRight1, d, C0, C1, C2);

		runningSum_for_Component_Phasors(Re_Eys, Im_Eys, Re_C_Hzs, Im_C_Hzs,
			Re_Ezs, Im_Ezs, Re_C_Hys, Im_C_Hys,
			k0, deltak, numkVals, xIndexes_for_power_calculation, NT);

	}//End Time stepping

	adjust_component_phasors_to_Yee_cell_face_centers(Re_Eys, Im_Eys, Re_C_Hzs, Im_C_Hzs,
		Re_Ezs, Im_Ezs, Re_C_Hys, Im_C_Hys,
		k0, deltak, numkVals, xIndexes_for_power_calculation, NT);		//component phasors adjusted to Yee cell face centers on each xp plane

	compute_spectral_power_densities(Re_Eys, Im_Eys, Re_C_Hzs, Im_C_Hzs,
		Re_Ezs, Im_Ezs, Re_C_Hys, Im_C_Hys, Pave_x,
		k0,deltak, numkVals, xIndexes_for_power_calculation, NT);	//Time average x-component spectral power density calculated at Yee cell face centers on each xp plane

	



	return 0;
}//End propagateIn



 /*****************************************************************************************************************************/

void EMField3D::computeCoefficientsForSecondOrderABCs(double& d, double& C0, double& C1, double& C2)
{
	d = ScExterior + 2.0 + 1 / ScExterior;

	C0 = -(ScExterior - 2.0 + 1 / ScExterior) / d;
	C1 = -2.0*(ScExterior - 1 / ScExterior) / d;
	C2 = 4.0*(ScExterior + 1 / ScExterior) / d;
}//End function computeCoefficientsForSecondOrderABCs



int EMField3D::updateMagneticField(Structure3DWithGrid& s, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy)
{
	int m = 0;
	int n = 0;
	int p = 0;


	//Update Hx
	for (m = 0; m < M; m++)
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P - 1; p++)
			Hx[m][n][p] = s.CHxH(m, n, p)*Hx[m][n][p] - s.CHxEz(m, n, p)*(Ez[m][n+1][p] - Ez[m][n][p]) +s.CHxEy(m, n, p)*(Ey[m][n][p+1] - Ey[m][n][p]);

	//Update Hy
	for (m = 0; m < M - 1; m++)
		for (n = 0; n < N; n++)
			for (p = 0; p < P - 1; p++)
				Hy[m][n][p] = s.CHyH(m, n, p)*Hy[m][n][p] - s.CHyEx(m, n, p)*(Ex[m][n][p+1] - Ex[m][n][p]) + s.CHyEz(m, n, p)*(Ez[m+1][n][p] - Ez[m][n][p]);


	//Update Hz
	for (m = 0; m < M - 1; m++)
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P; p++)
				Hz[m][n][p] = s.CHzH(m, n, p) * Hz[m][n][p] - s.CHzEy(m, n, p) * (Ey[m+1][n][p] - Ey[m][n][p]) + s.CHzEx(m, n, p)*(Ex[m][n+1][p] - Ex[m][n][p]);

	return 0;
}//End function updateMagneticField


void EMField3D::correctMagneticField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy)
{
	int n, p;

	//Correct Hy along TFSF boundary
	for (n = 0; n < N; n++)
		for (p = 0; p < P - 1; p++)
			Hy[mb - 1][n][p] -= s.CHyEz(mb - 1, n, p)*IFieldz.Ezdir(mb);

	//Correct Hz along TFSF boundary
	for (n = 0; n < N - 1; n++)
		for (p = 0; p < P; p++)
			Hz[mb - 1][n][p] += s.CHzEy(mb - 1, n, p)*IFieldy.Eydir(mb);

}//End function correctMagneticField_alongTFSFxFrontBoundary


int EMField3D::updateMagneticField_outsideBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
	vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ)
{
	int m, n, p;

	if (s.Is_periodic_in_yDirection())
	{
		for (m = 0; m < M; m++)
			for (p = 0; p < P - 1; p++)
				HxOutsideY[m][p] = Hx[m][0][p];
			//End for
		//End for

		for (m = 0; m < M - 1; m++)
			for (p = 0; p < P; p++)
				HzOutsideY[m][p] = Hz[m][0][p];
			//End for
		//End for

	}//End if


	if (s.Is_periodic_in_zDirection())
	{
		for (m = 0; m < M; m++)
			for (n = 0; n < N - 1; n++)
				HxOutsideZ[m][n] = Hx[m][n][0];
			//End for
		//End for

		for (m = 0; m < M - 1; m++)
			for (n = 0; n < N; n++)
				HyOutsideZ[m][n] = Hy[m][n][0];
			//End for
		//End for

	}//End if

	return 0;

}//End functionupdateMagneticField_outsideBoundaries_ofGridUnitCell


int EMField3D::updateElectricField(Structure3DWithGrid& s, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy)
{
	int m;
	int n;
	int p;


	//Update Ez
	for (m = 1; m < M - 1; m++)
		for (n = 1; n < N - 1; n++)
			for (p = 0; p < P - 1; p++)
				Ez[m][n][p] = s.CEzE(m, n, p) * Ez[m][n][p] + s.CEzHy(m, n, p) * (Hy[m][n][p] - Hy[m - 1][n][p]) - s.CEzHx(m, n, p)*(Hx[m][n][p] - Hx[m][n-1][p]);

	//Update Ey
	for (m = 1; m < M - 1; m++)
		for (n = 0; n < N - 1; n++)
			for (p = 1; p < P - 1; p++)
				Ey[m][n][p] = s.CEyE(m, n, p) * Ey[m][n][p] + s.CEyHx(m, n, p) * (Hx[m][n][p] - Hx[m][n][p-1]) - s.CEyHz(m, n, p)*(Hz[m][n][p] - Hz[m-1][n][p]);

	//Update Ex
	for (m = 0; m < M - 1; m++)
		for (n = 1; n < N - 1; n++)
			for (p = 1; p < P - 1; p++)
				Ex[m][n][p] = s.CExE(m, n, p) * Ex[m][n][p] + s.CExHz(m, n, p) * (Hz[m][n][p] - Hz[m][n-1][p]) - s.CExHy(m, n, p)*(Hy[m][n][p] - Hy[m][n][p-1]);

	return 0;
}//End function updateElectricField


int EMField3D::updateElectricField_onBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
	vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ)
{
	int m, n, p;

	if (s.Is_periodic_in_yDirection())
	{	
		//Update Ez on high y-boundary of unit cell and apply y-periodicity
		for (m = 1; m < M - 1; m++)
			for (p = 0; p < P - 1; p++)
			{
				Ez[m][N - 1][p] = s.CEzE(m, N - 1, p) * Ez[m][N - 1][p] + s.CEzHy(m, N - 1, p) * (Hy[m][N - 1][p] - Hy[m - 1][N - 1][p]) - s.CEzHx(m, N - 1, p)*(HxOutsideY[m][p] - Hx[m][N - 1 - 1][p]);
				Ez[m][0][p] = Ez[m][N - 1][p];
			}//end for
		//End for
				

		//Update Ex on high y-boundary of unit cell except for nodes also on high z-boundary and apply y-periodicity
		for (m = 0; m < M - 1; m++)
			for (p = 1; p < P - 1; p++)
			{
				Ex[m][N - 1][p] = s.CExE(m, N - 1, p) * Ex[m][N - 1][p] + s.CExHz(m, N - 1, p) * (HzOutsideY[m][p] - Hz[m][N - 1 - 1][p]) - s.CExHy(m, n, p)*(Hy[m][n][p] - Hy[m][n][p - 1]);
				Ex[m][0][p] = Ex[m][N - 1][p];
			}//End for
		//End for																											

	}//End if
	
	if (s.Is_periodic_in_zDirection())
	{
		//Update Ex on high z-boundary of unit cell including nodes also on high y-boundary and apply z-periodicity and apply y-periodiciy if structure also y-periodic
		for (m = 0; m < M - 1; m++)
		{ 
			for (n = 1; n < N - 1; n++)
			{
				Ex[m][n][P - 1] = s.CExE(m, n, P - 1) * Ex[m][n][P - 1] + s.CExHz(m, n, P - 1) * (Hz[m][n][P - 1] - Hz[m][n - 1][P - 1]) - s.CExHy(m, n, P - 1)*(HyOutsideZ[m][n] - Hy[m][n][P - 1 - 1]);
				Ex[m][n][0] = Ex[m][n][P - 1];
			}//end for

			if(s.Is_periodic_in_yDirection())
			{
				Ex[m][N - 1][P - 1] = s.CExE(m, N - 1, P - 1) * Ex[m][N - 1][P - 1] + s.CExHz(m, N - 1, P - 1) * (HzOutsideY[m][P - 1] - Hz[m][n - 1][P - 1])
					- s.CExHy(m, n, P - 1)*(HyOutsideZ[m][n] - Hy[m][n][P - 1 - 1]);
				Ex[m][0][0] = Ex[m][N - 1][P - 1];
			}//end if

		}//End for
		
		 //Update Ey on high z-boundary of unit cell and apply z-periodicity
		for (m = 1; m < M - 1; m++)
			for (n = 0; n < N - 1; n++)
			{
				Ey[m][n][P - 1] = s.CEyE(m, n, P - 1) * Ey[m][n][P - 1] + s.CEyHx(m, n, P - 1) * (HxOutsideZ[m][n] - Hx[m][n][P - 1 - 1]) - s.CEyHz(m, n, P - 1)*(Hz[m][n][P - 1] - Hz[m - 1][n][P - 1]);
				Ey[m][n][0] = Ey[m][n][P - 1];
			}//End for

	}//End if


}//End function updateElectricField_onBoundaries_ofGridUnitCell


void EMField3D::correctElectricField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1DzPolar& IFieldz, EMField1DyPolar& IFieldy)
{
	int n, p;

	//Correct Ez along TFSF boundary
	for (n = 1; n < N; n++)
		for (p = 0; p < P - 1; p++)
			Ez[mb][n][p] -= s.CEzHy(mb, n, p) * IFieldz.Hydir(mb - 1);

	//Correct Ey along TFSF boundary
	for (n = 0; n < N - 1; n++)
		for (p = 1; p < P ; p++)
			Ey[mb][n][p] += s.CEzHy(mb, n, p)*IFieldy.Hzdir(mb - 1);


}//End function correctElectricField_alongTFSFxFrontBoundary


void EMField3D::applySecondOrderABC_on_xBoundaries(vector<vector<vector<double>>> & EzOldLeft0, vector<vector<vector<double>>>& EzOldLeft1,
									vector<vector<vector<double>>>& EzOldRight0, vector<vector<vector<double>>> EzOldRight1,
									vector<vector<vector<double>>> & EyOldLeft0, vector<vector<vector<double>>>& EyOldLeft1,
									vector<vector<vector<double>>>& EyOldRight0, vector<vector<vector<double>>> EyOldRight1,
									double d, double C0, double C1, double C2)
{
	int m;
	int n;
	int p;

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
		Ez[M - 1][n][p] = C0*(Ez[M - 3][n][p] + EzOldRight1[M - 1][n][p]) + C1*(EzOldRight0[M - 1][n][p] + EzOldRight0[M - 3][n][p] - Ez[M - 2][n][p] - EzOldRight1[M - 2][n][p])
						+ C2*EzOldRight0[M - 2][n][p] - EzOldRight1[M - 3][n][p];
	}//End for

	 //Update past field values
	for (n = 0; n < N; n++)
		for (p = 0; p < P - 1; p++)
	{
		EzOldLeft1[0][n][p] = EzOldLeft0[0][n][p];
		EzOldLeft1[1][n][p] = EzOldLeft0[1][n][p];
		EzOldLeft0[0][n][p] = Ez[0][n][p];
		EzOldLeft0[1][n][p] = Ez[1][n][p];

		EzOldRight1[0][n][p] = EzOldRight0[0][n][p];
		EzOldRight1[1][n][p] = EzOldRight0[1][n][p];
		EzOldRight0[0][n][p] = Ez[M - 1][n][p];
		EzOldRight0[1][n][p] = Ez[M - 2][n][p];
	}//End For


	 //Second-order ABC for Ey left and right (i.e. above and below device)				
	for (n = 0; n < N - 1; n++)
		for (p = 0; p < P; p++)
		{
			Ey[0][n][p] = C0*(Ey[2][n][p] + EyOldLeft1[0][n][p]) + C1*(EyOldLeft0[0][n][p] + EyOldLeft0[2][n][p] - Ey[1][n][p] - EyOldLeft1[1][n][p]) + C2*EyOldLeft0[1][n][p] - EyOldLeft1[2][n][p];
			Ey[M - 1][n][p] = C0*(Ey[M - 3][n][p] + EyOldRight1[M - 1][n][p]) + C1*(EyOldRight0[M - 1][n][p] + EyOldRight0[M - 3][n][p] - Ey[M - 2][n][p] - EyOldRight1[M - 2][n][p])
				+ C2*EyOldRight0[M - 2][n][p] - EyOldRight1[M - 3][n][p];
		}//End for

		 //Update past field values
	for (n = 0; n < N - 1; n++)
		for (p = 0; p < P; p++)
		{
			EyOldLeft1[0][n][p] = EyOldLeft0[0][n][p];
			EyOldLeft1[1][n][p] = EyOldLeft0[1][n][p];
			EyOldLeft0[0][n][p] = Ey[0][n][p];
			EyOldLeft0[1][n][p] = Ey[1][n][p];

			EyOldRight1[0][n][p] = EyOldRight0[0][n][p];
			EyOldRight1[1][n][p] = EyOldRight0[1][n][p];
			EyOldRight0[0][n][p] = Ey[M - 1][n][p];
			EyOldRight0[1][n][p] = Ey[M - 2][n][p];
		}//End For

	 //no need second-order ABCs for magnetic fields. 2/19/16

}//End function applySecondOrderABC



 //******************************************************************************************************************************************************



void EMField3D::runningSum_for_Component_Phasors(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs, 
													Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys, 
													int k0, int deltak, int numkVals, const vector<int>& mp, const int NT)
{
	int xP, n, p, k;
	double wk;

	for (k = k0; k < k0 + numkVals*deltak + 1; k += deltak)
	{
		wk = (2 * PI*k / NT);
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P; p++)
			{
				for (xP = 0; xP < 2; xP++)
				{
					Re_Eys(k)[xP][n][p] += (2.0/NT)*Ey[mp[xP]][n][p] * cos(wk*q);
					Im_Eys(k)[xP][n][p] += (2.0 / NT)*Ey[mp[xP]][n][p] * (-sin(wk*q));

					Re_C_Hzs(k)[xP][n][p] += (2.0 / NT)*0.5*(Hz[mp[xP]][n][p] + Hz[mp[xP] - 1][n][p])*cos(wk*q);
					Im_C_Hzs(k)[xP][n][p] += (2.0 / NT)*0.5*(Hz[mp[xP]][n][p] + Hz[mp[xP] - 1][n][p])*sin(wk*q);

				}//End calculation for each xP 
			}// End calculation for Eys and C_Hzs

		for (n = 0; n < N; n++)
			for (p = 0; p < P - 1; p++)
			{
				for (xP = 0; xP < 2; xP++)
				{
					Re_Ezs(k)[xP][n][p] += (2.0 / NT)*Ez[mp[xP]][n][p] * cos(wk*q);
					Im_Ezs(k)[xP][n][p] += (2.0 / NT)*Ez[mp[xP]][n][p] * (-sin(wk*q));

					Re_C_Hys(k)[xP][n][p] += (2.0 / NT)*0.5*(Hy[mp[xP]][n][p] + Hy[mp[xP] - 1][n][p])*cos(wk*q);
					Im_C_Hys(k)[xP][n][p] += (2.0 / NT)*0.5*(Hy[mp[xP]][n][p] + Hy[mp[xP] - 1][n][p])*sin(wk*q);

				}//End calculation for each xP 
			}// End calculation for Ezs and C_Hys

	}//End for each frequency 
		
	
}//End function runningSum_for_Component_Phasors


void EMField3D::adjust_component_phasors_to_Yee_cell_face_centers(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
	Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys,
	int k0, int deltak, int numkVals, const vector<int>& mp, const int NT)
{
	int xP, n, p, k;

	for (k = k0; k < k0 + numkVals*deltak + 1; k += deltak)
	{
		
		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P - 1; p++)
			{
				for (xP = 0; xP < 2; xP++)
				{
					Re_Ezs(k)[xP][n][p] = 0.5*(Re_Ezs(k)[xP][n][p] + Re_Ezs(k)[xP][n + 1][p]);
					Im_Ezs(k)[xP][n][p] = 0.5*(Im_Ezs(k)[xP][n][p] + Im_Ezs(k)[xP][n + 1][p]);
					
					Re_C_Hys(k)[xP][n][p] = 0.5*(Re_C_Hys(k)[xP][n][p] + Re_C_Hys(k)[xP][n + 1][p]);
					Im_C_Hys(k)[xP][n][p] = 0.5*(Im_C_Hys(k)[xP][n][p] + Im_C_Hys(k)[xP][n + 1][p]);

					
					Re_Eys(k)[xP][n][p] = 0.5*(Re_Eys(k)[xP][n][p] + Re_Eys(k)[xP][n][p + 1]);
					Im_Eys(k)[xP][n][p] = 0.5*(Im_Eys(k)[xP][n][p] + Im_Eys(k)[xP][n][p + 1]);

					Re_C_Hzs(k)[xP][n][p] = 0.5*(Re_C_Hzs(k)[xP][n][p] + Re_C_Hzs(k)[xP][n][p + 1]);
					Im_C_Hzs(k)[xP][n][p] = 0.5*(Im_C_Hzs(k)[xP][n][p] + Im_C_Hzs(k)[xP][n][p + 1]);

				}//End calculation for each xP 
			}// End calculation for Ezs, C_Hys, Eys, C_Hzs

	}//End for each frequency

}//End function adjust_component_phasors_to_Yee_cell_face_centers


void EMField3D::compute_spectral_power_densities(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
	Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys, Array1D_Vector3D& Pave_x,
	int k0, int deltak, int numkVals, const vector<int>& mp, const int NT)
{
	int xP, n, p, k;

	for (k = k0; k < k0 + numkVals*deltak + 1; k += deltak)
	{

		for (n = 0; n < N - 1; n++)
			for (p = 0; p < P - 1; p++)
			{
				for (xP = 0; xP < 2; xP++)
				{
					Pave_x(k)[xP][n][p] = 0.5*(Re_Eys(k)[xP][n][p] * Re_C_Hzs(k)[xP][n][p] - Im_Eys(k)[xP][n][p] * Im_C_Hzs(k)[xP][n][p])
						- 0.5*(Re_Ezs(k)[xP][n][p] * Re_C_Hys(k)[xP][n][p] - Im_Ezs(k)[xP][n][p] * Im_C_Hys(k)[xP][n][p]);

				}//End calculation for each xP 
			}//
	}//End for each frequency

}//End function compute_spectral_power_densities


void EMField3D::compute_spectral_power_flows(Array1D_Vector3D Pave_x, vector<double> PxTot, Structure3DWithGrid& str, int k0, int deltak, int numkVals)
{
	int j, k, n, p;
	int nNodes_so_far, pNodes_so_far = 0;

	double delA = 0;

	nNodes_so_far = 0;
	pNodes_so_far = 0;


	for (k = k0; k < k0 + numkVals*deltak + 1; k += deltak)
	{
		//Accumulate power total for first y-layer
		for (n = nNodes_so_far; n < nNodes_so_far + str.nNdsYL(0); n++)
		{
			//First z-layer
			for (p = pNodes_so_far; p < pNodes_so_far + str.nNdsZL(0); p++)
			{
				delA = str.deltayL(0)*str.deltazL(0);
				PxTot[1] = PxTot[1] + Pave_x(k)[1][n][p] * delA;

				pNodes_so_far += 1;
			}//

			 //Correct Yee cell face area for cell on z-interface
			if (str.numLz() > 1)
				PxTot[1] = PxTot[1] * str.deltayL(0)*0.5*(str.deltazL(0) + str.deltazL(1)) / delA;



		}// End for each n in first y-layer

	}
		
			

			

	for (j = 0; j < str.numLy() - 1; j++)
		for (k = 0; k < str.numLz() - 1; k++)
			for (n = nNodes_so_far; n < nNodes_so_far + str.nNdsYL(j); n++)
			{}
}



 EMField1DzPolar::EMField1DzPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt) :Ez(gridToUse->size()), Hy(gridToUse->size())
 {
 int m = 0;

 theGrid = gridToUse;
 N = theGrid->size();
 mb = 0;
 mD = N - 1;		//Change this.

 Einc = IncidentField;

 q = 0;


 delt = aDelt;								//delt is in microseconds and delx in micrometers.

 //Initialize field

 for (m = 0; m < N; m++)
 {
 Ez[m] = 0;
 Hy[m] = 0;
 }//End for

 Ez[mb] = Einc(0, mb*theGrid->delx(0));

 }//End Constructor






 EMField1DzPolar::EMField1DzPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse) :Ez(gridToUse->size()), Hy(gridToUse->size())
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
 Ez[m] = 0;
 Hy[m] = 0;
 }//End for

 Ez[mb] = Einc(0, mb*theGrid->delx(0));

 }//End Constructor




 int EMField1DzPolar::propagateIn(Device1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1DzPolar*, int timeIndex), int displayTime)
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
 EMField1DzPolar* PointerToMe = this;
 displayFunctionToCall(PointerToMe, q);
 }


 }//End Time stepping


 return 0;
 }//End propagateIn



 int EMField1DzPolar::Update(Medium& am, double CEzE[], double  CEzH[], double CHyE[], double CHyH[], int displayCode)
 {

 double loss;
 int m;
 int mCount;
 int mCountL;
 int nL;


 //*******************************************************************************************************************
 //Compute update coefficients on first time call

 //Compute coeficients for top layer of two-layer grid on first time call

 if (q == 0)
 {
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

 CEzE[m] = (1 - loss) / (1 + loss);									//CEzE[m] is dimensionless
 CEzH[m] = delt / (am.ep()* theGrid->delx(nL - 1)) / (1 + loss);		//CEzH[m] is in mks units
 CHyH[m] = (1 - loss) / (1 + loss);
 CHyE[m] = delt / (am.mu()*theGrid->delx(nL - 1));						//CHyE[m] is in mks units

 mCount += 1;
 }//End for each node in lossy layer below the medium

 }//End if


 //*******************************************************************************************************************
 //Update fields


 //Update magnetic field values
 for (m = 0; m < N - 1; m++)
 Hy[m] = CHyH[m] * Hy[m] + CHyE[m] * (Ez[m + 1] - Ez[m]);

 //Correction for Hy[mb-1] if needed
 if (mb > 0)
 Hy[mb - 1] = Hy[mb - 1] - CHyE[mb - 1] * Einc(q*delt, mb*theGrid->delx(0));

 if (mb > 0)
 Ez[0] = Ez[1];	//Simple ABC for Ez[0] in case mb > 0
 else
 Ez[0] = Einc(q*delt, 0.0*theGrid->delx(0));

 //Update electric field values
 for (m = 1; m < N; m++)
 Ez[m] = CEzE[m] * Ez[m] + CEzH[m] * (Hy[m] - Hy[m - 1]);

 //Correction for Ez[mb] in case of TFSF boundary
 if (mb > 0)
 Ez[mb] = Ez[mb] + CEzH[mb] * (1 / Imp0)*Einc((q + 1.0 / 2.0)*delt, (mb - 1.0 / 2.0)*theGrid->delx(0));


 q += 1;

 return 0;
 }//End Update


 //*******************************************************************************************************************
 //Accessors

 double EMField1DzPolar::Ezdir(int m) const { return Ez[m]; }		//Returns Ez[m]
 double  EMField1DzPolar::Hydir(int m) const { return Hy[m]; }		//Returns Hy[m]
 int EMField1DzPolar::nENodes() const { return N; }				//Returns the total number of electric field nodes in the grid
 int EMField1DzPolar::TFSFboundaryEIndex() const { return mb; }






EMField1DyPolar::EMField1DyPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt) :Ey(gridToUse->size()), Hz(gridToUse->size())
{
	int m = 0;

	theGrid = gridToUse;
	N = theGrid->size();
	mb = 0;
	mD = N - 1;		//Change this.

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

}//End Constructor



EMField1DyPolar::EMField1DyPolar(double(*IncidentField)(double t, double x), Grid1D* gridToUse) :Ey(gridToUse->size()), Hz(gridToUse->size())
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

}//End Constructor




int EMField1DyPolar::propagateIn(Device1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1DyPolar*, int timeIndex), int displayTime)
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
			EMField1DyPolar* PointerToMe = this;
			displayFunctionToCall(PointerToMe, q);
		}


	}//End Time stepping


	return 0;
}//End propagateIn



int EMField1DyPolar::Update(Medium& am, double CEyE[], double  CEyH[], double CHzE[], double CHzH[], int displayCode)
{

	double loss;
	int m;
	int mCount;
	int mCountL;
	int nL;


	//*******************************************************************************************************************
	//Compute update coefficients on first time call

	//Compute coeficients for top layer of two-layer grid on first time call

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

			CEyE[m] = (1 - loss) / (1 + loss);									//CEyE[m] is dimensionless
			CEyH[m] = delt / (am.ep()* theGrid->delx(nL - 1)) / (1 + loss);		//CEyH[m] is in mks units
			CHzH[m] = (1 - loss) / (1 + loss);									//Check this. Why not CHzH = 1?
			CHzE[m] = delt / (am.mu()*theGrid->delx(nL - 1));						//CHzE[m] is in mks units

			mCount += 1;
		}//End for each node in lossy layer below the medium

	}//End if


	 //*******************************************************************************************************************
	 //Update fields


	 //Update magnetic field values
	for (m = 0; m < N - 1; m++)
		Hz[m] = CHzH[m] * Hz[m] - CHzE[m] * (Ey[m + 1] - Ey[m]);

	//Correction for Hz[mb-1] if needed. TFSF boundary used if EMFieldTEz has been instantiated with mb > 0.
	if (mb > 0)
		Hz[mb - 1] = Hz[mb - 1] + CHzE[mb - 1] * Einc(q*delt, mb*theGrid->delx(0));

	if (mb > 0)
		Ey[0] = Ey[1];	//Simple ABC for Ey[0] in case mb > 0
	else
		Ey[0] = Einc(q*delt, 0.0*theGrid->delx(0));

	//Update electric field values
	for (m = 1; m < N; m++)
		Ey[m] = CEyE[m] * Ey[m] - CEyH[m] * (Hz[m] - Hz[m - 1]);

	//Correction for Ey[mb] if needed. TFSF boundary used if EMFieldTEz has been instantiated with mb > 0.
	if (mb > 0)
		Ey[mb] = Ey[mb] - CEyH[mb] * (1 / Imp0)*Einc((q + 1.0 / 2.0)*delt, (mb - 1.0 / 2.0)*theGrid->delx(0));	//Check that Hzinc is -(1/Imp0)*Einc


	q += 1;

	return 0;
}//End Update


 //*******************************************************************************************************************
 //Accessors

double EMField1DyPolar::Eydir(int m) const { return Ey[m]; }		//Returns Ez[m]
double  EMField1DyPolar::Hzdir(int m) const { return Hz[m]; }		//Returns Hy[m]
int EMField1DyPolar::nENodes() const { return N; }				//Returns the total number of electric field nodes in the grid
int EMField1DyPolar::TFSFboundaryEIndex() const { return mb; }





//*********************************************************************************************************************************
 


 

}

#endif

