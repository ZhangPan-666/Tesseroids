// This cpp file is an example serial calculation script for 
// calculating gravitational potential field and its derivatives
// Make sure the parameter file "TFG.ForPar" is in the directory where the program is executed
// The resulting file name is Result_TFG.dat
// The core subfunction called by this script is named Tesseroid_GravityEstimate. 
// For the definition of the input and output of this function, please refer to line 153.

#include "include/tesseroid.h"

int main(int argc, char** argv)
{
	int i, j, k;
	int NumDensityModel, NumCoordinates;

	double AbsTol, RelTol;

	double* Fai1, * Fai2, * Lamda1, * Lamda2, * R1, * R2, * rho;
	double* Longitude, * Latitude, * Radius;
	double* V,
		* Vx, * Vy, * Vz,
		* Vxx, * Vxy, * Vyy, * Vzx, * Vzy, * Vzz,
		* Vxxx, * Vxxy, * Vxxz, * Vxyz, * Vyyx, * Vyyy, * Vyyz, * Vzzx, * Vzzy, * Vzzz;

	char ParametersFile[1024] = "TFG.ForPar";
	char ResultsFile[1024] = "Result_TFG.dat";

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

	fread(&NumDensityModel, 4, 1, fp_Parameters);

	Fai1 = (double*)malloc(sizeof(double) * NumDensityModel);
	Fai2 = (double*)malloc(sizeof(double) * NumDensityModel);
	Lamda1 = (double*)malloc(sizeof(double) * NumDensityModel);
	Lamda2 = (double*)malloc(sizeof(double) * NumDensityModel);
	R1 = (double*)malloc(sizeof(double) * NumDensityModel);
	R2 = (double*)malloc(sizeof(double) * NumDensityModel);
	rho = (double*)malloc(sizeof(double) * NumDensityModel);

	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&Fai1[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&Fai2[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&Lamda1[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&Lamda2[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&R1[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&R2[i], 8, 1, fp_Parameters);
	}
	for (i = 0; i < NumDensityModel; i++)
	{
		fread(&rho[i], 8, 1, fp_Parameters);
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

	Vxxx = (double*)malloc(sizeof(double) * NumCoordinates);
	Vxxy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vxxz = (double*)malloc(sizeof(double) * NumCoordinates);
	Vxyz = (double*)malloc(sizeof(double) * NumCoordinates);
	Vyyx = (double*)malloc(sizeof(double) * NumCoordinates);
	Vyyy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vyyz = (double*)malloc(sizeof(double) * NumCoordinates);
	Vzzx = (double*)malloc(sizeof(double) * NumCoordinates);
	Vzzy = (double*)malloc(sizeof(double) * NumCoordinates);
	Vzzz = (double*)malloc(sizeof(double) * NumCoordinates);

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

		Vxxx[i] = 0.0;
		Vxxy[i] = 0.0;
		Vxxz[i] = 0.0;
		Vxyz[i] = 0.0;
		Vyyx[i] = 0.0;
		Vyyy[i] = 0.0;
		Vyyz[i] = 0.0;
		Vzzx[i] = 0.0;
		Vzzy[i] = 0.0;
		Vzzz[i] = 0.0;
	}

	// "Tesseroid_GravityEstimate" is a function that calculates the gravitational field and its derivatives of spherical prisms in 
	// serial mode.
	// 
	// Input Parameters:
	// 
	// NumDensityModel: Total number of the density models, data type: integer, length: 1, units: none
	// Fai1: Starting latitudinal boundaries, data type: double, length: NumDensityModel, units: degree
	// Fai2: Ending latitudinal boundaries, data type: double, length: NumDensityModel, units: degree
	// Lamda1: Starting longitudinal boundaries, data type: double, length: NumDensityModel, units: degree
	// Lamda2: Ending longitudinal boundaries, data type: double, length: NumDensityModel, units: degree
	// R1: Starting radial boundaries, data type: double, length: NumDensityModel, units: user-defined
	// R2: Ending radial boundaries, data type: double, length: NumDensityModel, units: user-defined
	// rho: Density of the spherical prisms, data type: double, length: NumDensityModel, units: user-defined
	// NumCoordinates: Total number of the computational points, data type: integer, length: 1, units: none
	// Latitude: Latitude of the the computational points, data type: double, length: NumCoordinates, units: degree
	// Longitude: Longitude of the the computational points, data type: double, length: NumCoordinates, units: degree
	// Radius: Radius of the the computational points, data type: double, length: NumCoordinates, units: user-defined
	// AbsTol: Termination condition for integration: absolute error, data type: double, length 1, units: none
	// RelTol: Termination condition for integration: relative error, data type: double, length 1, units: none
	//
	// Output Parameters:
	// V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy, Vzz, Vxxx, Vxxy, Vxxz, Vxyz, Vyyx, Vyyy, Vyyz, Vzzx, Vzzy, Vzzz: 
	// These parameters are the result of the forward computation
	// data type: double, length: NumCoordinates, units: depending on the unit of the user input parameter

	mista_math::Tesseroid_GravityEstimate(NumDensityModel, Fai1, Fai2, Lamda1, Lamda2, R1, R2, rho,
		NumCoordinates, Latitude, Longitude, Radius,
		AbsTol, RelTol,
		V,
		Vx, Vy, Vz,
		Vxx, Vxy, Vyy, Vzx, Vzy, Vzz,
		Vxxx, Vxxy, Vxxz, Vxyz, Vyyx, Vyyy, Vyyz, Vzzx, Vzzy, Vzzz);

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
	fwrite(Vxxx, 8, NumCoordinates, fp_Results);
	fwrite(Vxxy, 8, NumCoordinates, fp_Results);
	fwrite(Vxxz, 8, NumCoordinates, fp_Results);
	fwrite(Vxyz, 8, NumCoordinates, fp_Results);
	fwrite(Vyyx, 8, NumCoordinates, fp_Results);
	fwrite(Vyyy, 8, NumCoordinates, fp_Results);
	fwrite(Vyyz, 8, NumCoordinates, fp_Results);
	fwrite(Vzzx, 8, NumCoordinates, fp_Results);
	fwrite(Vzzy, 8, NumCoordinates, fp_Results);
	fwrite(Vzzz, 8, NumCoordinates, fp_Results);

	fclose(fp_Results);

	return 1;
}













