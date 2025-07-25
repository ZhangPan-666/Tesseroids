// This cpp file is an example serial calculation script for 
// calculating magnetic potential field and its derivatives

// Make sure the parameter file "TFM.ForPar" is in the directory where the program is executed
// The resulting file name is Result_TFM.dat

// The core subfunction called by this script is named Tesseroid_MagneticEstimate. 
// For the definition of the input and output of this function, please refer to line 143.

#include "include/tesseroid.h"

int main(int argc, char** argv)
{
	int i, j, k;
	int NumPermeabilityModel, NumCoordinates;

	double AbsTol, RelTol;

	double* Fai1, * Fai2, * Lamda1, * Lamda2, * R1, * R2, * Mx, * My, * Mz;
	double* Longitude, * Latitude, * Radius;
	double* V,
		* Vx, * Vy, * Vz,
		* Vxx, * Vxy, * Vyy, * Vzx, * Vzy, * Vzz;

	char ParametersFile[1024] = "TFM.ForPar";
	char ResultsFile[1024] = "Result_TFM.dat";

	FILE* fp_Parameters;
	FILE* fp_Results;

	fp_Parameters = fopen(ParametersFile, "rb+");

	if (fp_Parameters == NULL) {
		printf("\nFailed to open parameter file ! \n");
		printf("Calculation program exit: calculation failed ! \n");
		return 0;
	}

	fread(&AbsTol, 8, 1, fp_Parameters);
	fread(&RelTol, 8, 1, fp_Parameters);

	fread(&NumPermeabilityModel, 4, 1, fp_Parameters);

	Fai1 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	Fai2 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	Lamda1 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	Lamda2 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	R1 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	R2 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	Mx = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	My = (double*)malloc(sizeof(double) * NumPermeabilityModel);
	Mz = (double*)malloc(sizeof(double) * NumPermeabilityModel);

	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&Fai1[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&Fai2[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&Lamda1[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&Lamda2[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&R1[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&R2[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&Mx[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&My[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumPermeabilityModel; i++)
	{
		fread(&Mz[i], 8, 1, fp_Parameters);
	}

	fread(&NumCoordinates, 4, 1, fp_Parameters);

	Longitude = (double*)malloc(sizeof(double) * NumCoordinates);
	Latitude = (double*)malloc(sizeof(double) * NumCoordinates);
	Radius = (double*)malloc(sizeof(double) * NumCoordinates);

	for (i = 0; i < NumCoordinates; i++)
	{
		fread(&Longitude[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumCoordinates; i++)
	{
		fread(&Latitude[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumCoordinates; i++)
	{
		fread(&Radius[i], 8, 1, fp_Parameters);
	}

	fclose(fp_Parameters);

	V = (double*)malloc(sizeof(double) * NumCoordinates);

	Vx = (double*)malloc(sizeof(double) * NumCoordinates);
	Vy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vz = (double*)malloc(sizeof(double) * NumCoordinates);

	Vxx = (double*)malloc(sizeof(double) * NumCoordinates);
	Vxy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vyy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vzx = (double*)malloc(sizeof(double) * NumCoordinates);
	Vzy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vzz = (double*)malloc(sizeof(double) * NumCoordinates);

	for (i = 0; i < NumCoordinates; i++)
	{

		V[i] = 0.0;

		Vx[i] = 0.0;
		Vy[i] = 0.0;
		Vz[i] = 0.0;

		Vxx[i] = 0.0;
		Vxy[i] = 0.0;
		Vyy[i] = 0.0;
		Vzx[i] = 0.0;
		Vzy[i] = 0.0;
		Vzz[i] = 0.0;

	}

	// "Tesseroid_MagneticEstimate" is a function that calculates the magnetic field and its derivatives of spherical prisms in 
	// serial mode.
	// 
	// Input Parameters:
	// 
	// NumPermeabilityModel: Total number of the permeability models, data type: integer, length: 1, units: none
	// Fai1: Starting latitudinal boundaries, data type: double, length: NumPermeabilityModel, units: degree
	// Fai2: Ending latitudinal boundaries, data type: double, length: NumPermeabilityModel, units: degree
	// Lamda1: Starting longitudinal boundaries, data type: double, length: NumPermeabilityModel, units: degree
	// Lamda2: Ending longitudinal boundaries, data type: double, length: NumPermeabilityModel, units: degree
	// R1: Starting radial boundaries, data type: double, length: NumPermeabilityModel, units: user-defined
	// R2: Ending radial boundaries, data type: double, length: NumPermeabilityModel, units: user-defined
	// Mx: Magnetization intensity (northern component) of the spherical prisms, data type: double, length: NumPermeabilityModel, units: user-defined
	// My: Magnetization intensity (eastern component) of the spherical prisms, data type: double, length: NumPermeabilityModel, units: user-defined
	// Mz: Magnetization intensity (radial component) of the spherical prisms, data type: double, length: NumPermeabilityModel, units: user-defined
	// NumCoordinates: Total number of the computational points, data type: integer, length: 1, units: none
	// Latitude: Latitude of the the computational points, data type: double, length: NumCoordinates, units: degree
	// Longitude: Longitude of the the computational points, data type: double, length: NumCoordinates, units: degree
	// Radius: Radius of the the computational points, data type: double, length: NumCoordinates, units: user-defined
	// AbsTol: Termination condition for integration: absolute error, data type: double, length 1, units: none
	// RelTol: Termination condition for integration: relative error, data type: double, length 1, units: none
	//
	// Output Parameters:
	// V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy, Vzz: 
	// These parameters are the result of the forward computation
	// data type: double, length: NumCoordinates, units: depending on the unit of the user input parameter

	mista_math::Tesseroid_MagneticEstimate(NumPermeabilityModel, Fai1, Fai2, Lamda1, Lamda2, R1, R2, Mx, My, Mz,
		NumCoordinates, Latitude, Longitude, Radius,
		AbsTol, RelTol,
		V,
		Vx, Vy, Vz,
		Vxx, Vxy, Vyy, Vzx, Vzy, Vzz);

	fp_Results = fopen(ResultsFile, "wb+");

	fwrite(&NumCoordinates, 4, 1, fp_Results);

	fwrite(Longitude, 8, NumCoordinates, fp_Results);
	fwrite(Latitude, 8, NumCoordinates, fp_Results);
	fwrite(Radius, 8, NumCoordinates, fp_Results);

	fwrite(V, 8, NumCoordinates, fp_Results);
	fwrite(Vx, 8, NumCoordinates, fp_Results);
	fwrite(Vy, 8, NumCoordinates, fp_Results);
	fwrite(Vz, 8, NumCoordinates, fp_Results);
	fwrite(Vxx, 8, NumCoordinates, fp_Results);
	fwrite(Vxy, 8, NumCoordinates, fp_Results);
	fwrite(Vyy, 8, NumCoordinates, fp_Results);
	fwrite(Vzx, 8, NumCoordinates, fp_Results);
	fwrite(Vzy, 8, NumCoordinates, fp_Results);
	fwrite(Vzz, 8, NumCoordinates, fp_Results);

	fclose(fp_Results);

	return 1;
}













