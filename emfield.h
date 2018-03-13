#ifndef _EMFIELD_H
#define _EMFIELD_H


#include "grid.h"

#include <vector>



using namespace std;

namespace FDTD
{

	//***************************************************************************************************************************************************
	//Class representing one-dimensional electromagmetic field with electric field in z-direction and propagation in x-direction

	class EMField1D_zPolarized
	{
	public:

		EMField1D_zPolarized() {}

		EMField1D_zPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt);				//Field with incident field assigned at mb = 0 at time q = 0
																													//delt assigned aDelt
																													//delt is to be in microseconds and delx in micrometers.

		EMField1D_zPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse);
		//Field with incident field assigned at a TFSF boundary between
		//Ez electric field node mb = int(theGrid->nNds(0)/2 + 0.5) and mb-1/2
		//in top layer of grid at q = 0.
		//The grid provided must be based on a device in order to subsequently use the prodageIn function.
		//delt based on Courant number of unity and minimum delx for grid.
		//delt is in microseconds and delx in micrometers.

		int propagateIn(Structure1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1D_zPolarized*, int timeIndex), int displayTime);
		//Simulates propagation of field in 1D device described by first argument.
		//Implements TFSF boundary above device with assumption that mb is as assigned in second constructor (.0).
		//The dev argument must fit the grid used by *this.
		//Time step is calculated based on minimum deltax and slowest speed of propagation in any of the layers
		// Calls displayFunctionToCall after update at every displayTime time steps.

		void compute_coefficients_for_medium_with_lossy_layer_underneath(Medium& am, vector<double>& CEzE, vector<double>& CEzH, vector<double>& CHyE, vector<double>& CHyH);

		int Update(vector<double>& CEzE, vector<double>& CEzH, vector<double>& CHyE, vector<double>& CHyH);
		//Update fields in homogeneous medium with lossy layer behind it for one time step.
		//Grid used by *this must be two-layer with bottom layer intended for
		//lossy layer to terminate grid.
		//medium is in the index zero layer. Implements lossy layer in second grid layer.


		double Ezdir(int m) const;					//Returns Ez[m]
		double Hydir(int m) const;					//Returns Hy[m]
		int nENodes() const;						//Returns the total number of Ez x = const. node planes in the grid
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

	class EMField1D_yPolarized
	{
	public:
		EMField1D_yPolarized() {}

		EMField1D_yPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse, double aDelt);			//Field with incident field assigned at mb = 0 at time q = 0
																												//delt assigned aDelt
																												//delt is to be in microseconds and delx in micrometers.

		EMField1D_yPolarized(double(*IncidentField)(double t, double x), Grid1D* gridToUse);			//Field with incident field assigned at a TFSF boundary between 
																								//Ey electric field node mb = int(theGrid->nNds(0)/2 + 0.5) and mb-1/2
																								//in top layer of grid at q = 0. 
																								//The grid provided must be based on a device in order to subsequently use the prodageIn function.  
																								//delt based on Courant number of unity and minimum delx for grid.
																								//delt is in microseconds and delx in micrometers.

		int propagateIn(Structure1D& dev, int maxTime, int(*displayFunctionToCall)(const EMField1D_yPolarized*, int timeIndex), int displayTime);
		//Simulates propagation of field in 1D device described by first argument.
		//Implements TFSF boundary above device with assumption that mb is as assigned in second constructor (.0). 
		//The dev argument must fit the grid used by *this. 
		//Calls displayFunctionToCall after update at every displayTime time steps.

		void compute_coefficients_for_medium_with_lossy_layer_underneath(Medium& am, vector<double>& CEyE, vector<double>& CEyH, vector<double>& CHzE, vector<double>& CHzH);

		int Update(vector<double>& CEyE, vector<double>& CEyH, vector<double>& CHzE, vector<double>& CHzH);			//Update fields in homogeneous medium with lossy layer behind it for one time step. 
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
	//Default constructor

	EMField3D(double(*IncidentEzField)(double t, double x), double(*IncidentEyField)(double t, double x), Structure3DWithGrid& aStr, int anumkVals, int amb);
	
	//Constructs an EMField3D object representing an electromagnetic field that at t == 0 is described by IncidentEzField and IncidentEyField at the TFSF boundary between amb - 1 and amb 
	//and zero everywhere else in the computational domain.

	//The the size of the computational domain and the grid used to describe the field are obtained from aStr. 
	
	//numkVals must be the number of frequency values for which spectral power flows are to be computed.
	//amb must be the x-Index of Ez field nodes to right of TF/SF boundary

	
	int propagateIn(Structure3DWithGrid& str, double aTinc, int maxTime, int k0, int deltak, int(*displayFunctionToCall)(const vector<vector<vector<double>>>& F, int q),
					int displayTime, int(*displayFunctionToCall_for_1D)(const EMField1D_zPolarized*, int timeIndex));
	//Simulates propagation during t > 0 of the incident one-dimensional pulse or wave already defined at t == 0 at the TFSF boundary. 
	//The pulse or wave will be incident in the x-direction on, and propogates in, the unit cell str of a structure that is periodic in y and z-directions.
	//aTinc is the duration of the incident pulse or period of the harmonic source
	//maxTime is the exclusive upper bound of the time index. 
	//k0 is lowest frequency index for which spectral power flow is to be calculated.
	//deltak is difference in succcessive frequency indexs for which spectral power flow is to be calculated.

	//displayFunction to call is used to display the field values.
	//disPlayTime is the number of time steps untill next display.
	//displayFunctionToCall_for_1D is used to display the one-dimensional incident field above the structure

	int TFSFboundaryEIndex() const;				//Returns index of electric field node to the right of and nearest TFSF boundary

	int nEzNodePlanes_in_xDir() const;			//Returns the total number of Ez x = const. node planes in the grid	


	double time_average_incident_spectral_power_at_frequency_index(int k); 
	//k must be one of the frequency indexes for which spectral power computed. 
	//wk = PI*k/(Np/2). fk = k/Np, where Np is the pulse duration or period in number of time steps. 

	double time_average_reflected_spectral_power_at_frequency_index(int k);
	//Returns time_average spectral_power_at_frequency_index k reflected from the structure. 
	//k must be one of the frequency indexes for which spectral power computed. 
	//wk = PI*k/(Np/2). fk = k/Np, where Np is the pulse duration or period in number of time steps.
	//Note: The value returned will depend on the time duration of the simulation, which is not neccessaarily enough
	//to simulate all possible internal reflections in the structure.
																					

	double time_average_transmitted_spectral_power_at_frequency_index(int k);
	//Returns time_average spectral_power_at_frequency_index k transmitted through the structure. 
	//k must be one of the frequency indexes for which spectral power computed. 
	//wk = PI*k/(Np/2). fk = k/Np, where Np is the pulse duration or period in number of time steps.
	//Note: The value returned will depend on the time duration of the simulation, which is not neccessaarily enough
	//to simulate all possible internal reflections in the structure.

	double spectralPower(int vIndex, int locationIndex);

	double Ez_comp(int mInd, int nInd, int pInd) const;			//Returns Ez[mInd][nInd][pInd]

	int numPowerMonitors() const;				//returns //Number of x=const planes for which power flows calculated.

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

	int M;		//M is the number of Ez node planes across 3D grid in x - direction. Also equals number of Hx magnetic field node planes
				//and number of Ey node planes across 3D grid in x-direction since this code assumes planes of Ez, Hx, and Ey nodes on x-boundaries.
					
				//Number of Hz nodes across grid in x-direction equals M - 1. Number of Ex nodes across grid in x-direction equals M - 1.
				//Number of Hy magnetic firld nodes across grid in x-direction equals M - 1.

	int N;		//N is the number of Ez nodes across 3D grid in y-direction. Also equals the number of Ex nodes and number of Hy nodes 
				//across 3D grid in y-direction since this code assumes planes of Ez, Hy and Ex nodes on y boundaries.
				
				//Number of Hx magnetic field nodes across grid in y-direction is N - 1. 
				//Number of Hz magnetic field nodes across grid in y-direction is N - 1.
				//Number of Ey nodes across grid in y-direction is N - 1.

	int P;		//Number of Ex nodes across grid in z-direction. Also equals number of Ey nodes and number of Hz nodes across grid in z-direction 
				//since this code assumes planes of Ex, Ey, and Hz nodes on z-boundaries.


	int mb;		//x-Index of Ez field nodes to right of TF/SF boundary
	int mD;		//Top surfaace of structure is at Ez field node with x-Index mD

	Grid3D theGrid;	//Pointer to grid shared with device

	const double Sc = 1 / sqrt(3.0);	//Choose delt = Sc*delMin/c0 so that local Courant number will not exceed 1/sqrt(3.0).

	double ScExterior;					// c0*delt/delx (= Sc for layer above and below device. Assumes delx the same above and below the device)

	int numkVals;

	int k0;
	int deltak;

	int Number_of_power_monitors;		//Number of x=const planes for which power flows calculated

	int Location_index_for_incident_power;

	int Location_index_for_transmitted_power;

	int Location_index_for_reflected_power;

	Array1D_Vector1D Pflow;

	int updateMagneticField(Structure3DWithGrid& str, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy);

	int updateMagneticField_outsideBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
															vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ);			
																	//Used only if the structure is the unit cell of a larger periodic structure 
																	//This instance for normal incidence only

	int updateElectricField(Structure3DWithGrid& str, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy);

	void updateElectricField_onBoundaries_ofGridUnitCell(Structure3DWithGrid& s, vector <vector<double>>& HxOutsideY, vector <vector<double>>& HzOutsideY,
															vector <vector<double>>& HxOutsideZ, vector <vector<double>>& HyOutsideZ);
																	//Used only if the structure is the unit cell of a larger periodic structure 
																	//This instance for normal incidence only

	void correctMagneticField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy);

	void correctElectricField_alongTFSFxFrontBoundary(Structure3DWithGrid& s, EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy);

	void computeCoefficientsForSecondOrderABCs(double& d, double& C0, double& C1, double& C2);

	void applySecondOrderABC_on_xBoundaries(vector<vector<vector<double>>>& EzOldLeft0, vector<vector<vector<double>>>& EzOldLeft1,
							vector<vector<vector<double>>>& EzOldRight0, vector<vector<vector<double>>>& EzOldRight1,
							vector<vector<vector<double>>>& EyOldLeft0, vector<vector<vector<double>>>& EyOldLeft1,
							vector<vector<vector<double>>>& EyOldRight0, vector<vector<vector<double>>>& EyOldRight1,
							double d, double C0, double C1, double C2);

	void runningSum_for_Component_Phasors(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs, 
											Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys,
											EMField1D_zPolarized& IFieldz, EMField1D_yPolarized& IFieldy,
											const vector<int>& mp, int aNpulse);

	void adjust_component_phasors_to_Yee_cell_face_centers(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys,
		const vector<int>& mp);
										//After this function Re_Ezs(k)[xP][n][p] etc. will be value at the center of the x=const face of Yee cell with indexes mp[xp], n, p.
										//The actual mathematical index of the face center is mp[xp], n+1/2, p+1/2.

	void compute_spectral_power_densities(Array1D_Vector3D& Re_Eys, Array1D_Vector3D& Im_Eys, Array1D_Vector3D& Re_C_Hzs, Array1D_Vector3D& Im_C_Hzs,
		Array1D_Vector3D& Re_Ezs, Array1D_Vector3D& Im_Ezs, Array1D_Vector3D& Re_C_Hys, Array1D_Vector3D& Im_C_Hys, Array1D_Vector3D& Pave_x,
		const vector<int>& mp);

	void compute_spectral_power_flows(Array1D_Vector3D& Pave_x, Structure3DWithGrid& str, const vector<int>& mp);

};


}//End namespace FDTD


#endif

