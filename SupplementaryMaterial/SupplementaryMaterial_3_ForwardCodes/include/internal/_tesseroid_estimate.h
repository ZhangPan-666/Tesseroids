// This cpp header file contains some of the subfunctions needed for
// calculating gravitational and magnetic fields of a spherical prism
// For more details, see Zhang P.et al.,
// "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

#pragma once
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <array>
#include <algorithm>
#include <mpi.h>

#include "_integral2.h"
#include "_tesseroid_integralkernel.h"


namespace mista_math
{

inline void MPIDistribute(int Num, int* DataList, int* VectorShiftList)
{
	int i, j;
	int rank, nprocs;
	int DataPerBatch, NumRemain;

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	NumRemain = Num % nprocs;
	DataPerBatch = (Num - NumRemain) / nprocs;

	for (i = 1; i <= nprocs; i++)
	{
		if (i <= NumRemain)
		{
			DataList[i - 1] = DataPerBatch + 1;
		}
		else
		{
			DataList[i - 1] = DataPerBatch;
		}
	}

	for (i = 1; i <= nprocs; i++)
	{
		VectorShiftList[i - 1] = 0;

		if (i == 1)
		{
			continue;
		}

		for (j = 1; j <= (i - 1); j++)
		{
			VectorShiftList[i - 1] += DataList[j - 1];
		}

	}

}

inline std::vector<std::array<double, 4>> splitTesseroid(double Fai1, double Fai2,
	double Lamda1, double Lamda2,
	double P_Fai, double P_Lamda,
	double dFai, double dLamda) {
	if (P_Fai < Fai1 || P_Fai > Fai2 || P_Lamda < Lamda1 || P_Lamda > Lamda2) {
		return { {Fai1, Fai2, Lamda1, Lamda2} };
	}

	double half_dFai = dFai / 2.0;
	double half_dLamda = dLamda / 2.0;

	double small_Fai1 = std::max(P_Fai - half_dFai, Fai1);
	double small_Fai2 = std::min(P_Fai + half_dFai, Fai2);
	double small_Lamda1 = std::max(P_Lamda - half_dLamda, Lamda1);
	double small_Lamda2 = std::min(P_Lamda + half_dLamda, Lamda2);

	std::vector<std::array<double, 4>> split_rectangles;

	if (small_Fai2 < Fai2) {
		split_rectangles.push_back({ small_Fai2, Fai2, Lamda1, Lamda2 });
	}

	if (small_Fai1 > Fai1) {
		split_rectangles.push_back({ Fai1, small_Fai1, Lamda1, Lamda2 });
	}

	if (small_Lamda1 > Lamda1) {
		split_rectangles.push_back({ small_Fai1, small_Fai2, Lamda1, small_Lamda1 });
	}

	if (small_Lamda2 < Lamda2) {
		split_rectangles.push_back({ small_Fai1, small_Fai2, small_Lamda2, Lamda2 });
	}

	return split_rectangles;
}

inline void Tesseroid_GravityPointEstimate(double Fai1, double Fai2, double Lamda1, double Lamda2, double R1, double R2,
	double fai, double lamda, double r,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz,
	double* Vxxx, double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy, double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz)
{
	int i;

	double pi, RadianCorrection;
	double VSphericalShell,
		VxSphericalShell, VySphericalShell, VzSphericalShell,
		VxxSphericalShell, VxySphericalShell, VyySphericalShell, VzxSphericalShell, VzySphericalShell, VzzSphericalShell,
		VxxxSphericalShell, VxxySphericalShell, VxxzSphericalShell, VxyzSphericalShell, VyyxSphericalShell,
		VyyySphericalShell, VyyzSphericalShell, VzzxSphericalShell, VzzySphericalShell, VzzzSphericalShell;
	double VSplit,
		VxSplit, VySplit, VzSplit,
		VxxSplit, VxySplit, VyySplit, VzxSplit, VzySplit, VzzSplit,
		VxxxSplit, VxxySplit, VxxzSplit, VxyzSplit, VyyxSplit, VyyySplit, VyyzSplit, VzzxSplit, VzzySplit, VzzzSplit;

	std::function<Eigen::Matrix<double, 14, 14>(const Eigen::Matrix<double, 14, 14>&, const Eigen::Matrix<double, 14, 14>&)>
		func_V;
	std::function<Eigen::Matrix<double, 14, 14>(const Eigen::Matrix<double, 14, 14>&, const Eigen::Matrix<double, 14, 14>&)>
		func_Vx, func_Vy, func_Vz;
	std::function<Eigen::Matrix<double, 14, 14>(const Eigen::Matrix<double, 14, 14>&, const Eigen::Matrix<double, 14, 14>&)>
		func_Vxx, func_Vxy, func_Vyy, func_Vzx, func_Vzy, func_Vzz;
	std::function<Eigen::Matrix<double, 14, 14>(const Eigen::Matrix<double, 14, 14>&, const Eigen::Matrix<double, 14, 14>&)>
		func_Vxxx, func_Vxxy, func_Vxxz, func_Vxyz, func_Vyyx, func_Vyyy, func_Vyyz, func_Vzzx, func_Vzzy, func_Vzzz;

	pi = 4.0 * atan(1.0);
	RadianCorrection = (2.0 * pi * pi) / (180.0 * 360.0);

	if (lamda < -180.0) {
		lamda = 360.0 + lamda;
	}
	else if (lamda > 180.0) {
		lamda = -360.0 + lamda;
	}

	func_V = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelV(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vx = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVx(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vz = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVz(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vxx = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVxx(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vxy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVxy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vyy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVyy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vzx = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVzx(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vzy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVzy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vzz = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVzz(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vxxx = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVxxx(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vxxy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVxxy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vxxz = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVxxz(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vxyz = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVxyz(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vyyx = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVyyx(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vyyy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVyyy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vyyz = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVyyz(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vzzx = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVzzx(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vzzy = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVzzy(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};
	func_Vzzz = [&](const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI) {
		return Tesseroid_IntegralkernelVzzz(R2, R1, FaiI, LamdaI, r, fai, lamda);
		};

	if (abs(fai) == 90.0) {

		*V = integral2t(func_V, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vx = integral2t(func_Vx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vy = integral2t(func_Vy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vz = integral2t(func_Vz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vxx = integral2t(func_Vxx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vxy = integral2t(func_Vxy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vyy = integral2t(func_Vyy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vzx = integral2t(func_Vzx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vzy = integral2t(func_Vzy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vzz = integral2t(func_Vzz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vxxx = integral2t(func_Vxxx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vxxy = integral2t(func_Vxxy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vxxz = integral2t(func_Vxxz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vxyz = integral2t(func_Vxyz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vyyx = integral2t(func_Vyyx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vyyy = integral2t(func_Vyyy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vyyz = integral2t(func_Vyyz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vzzx = integral2t(func_Vzzx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vzzy = integral2t(func_Vzzy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
		*Vzzz = integral2t(func_Vzzz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);

		return;
	}

	if ((fai > Fai1) && (fai < Fai2) && (lamda > Lamda1) && (lamda < Lamda2)) {

		std::vector<std::array<double, 4>> subTesseroids;

		subTesseroids = splitTesseroid(-90.0, 90.0, -180.0, 180.0,
			(Fai1 + Fai2) / 2.0, (Lamda1 + Lamda2) / 2.0, Fai2 - Fai1, Lamda2 - Lamda1);

		int NumDensityModel = subTesseroids.size();

		if (r >= R2) {
			double M0 = 4.0 * pi / 3.0 * (pow(R2, 3) - pow(R1, 3));
			VSphericalShell = M0 / r;
			VzSphericalShell = -M0 / (r * r);
			VxxSphericalShell = -M0 / (r * r * r);
			VyySphericalShell = -M0 / (r * r * r);
			VzzSphericalShell = 2 * M0 / (r * r * r);
			VxxzSphericalShell = 3 * M0 / (r * r * r * r);
			VyyzSphericalShell = 3 * M0 / (r * r * r * r);
			VzzzSphericalShell = -6 * M0 / (r * r * r * r);
		}
		else if (r <= R1) {
			VSphericalShell = 2 * pi * (pow(R2, 2) - pow(R1, 2));
			VzSphericalShell = 0;
			VxxSphericalShell = 0;
			VyySphericalShell = 0;
			VzzSphericalShell = 0;
			VxxzSphericalShell = 0;
			VyyzSphericalShell = 0;
			VzzzSphericalShell = 0;
		}
		else {
			VSphericalShell = 2 * pi * (pow(R2, 2) - pow(r, 2) / 3.0 - 2 * pow(R1, 3) / (3.0 * r));
			VzSphericalShell = -4 * pi / 3.0 * (r - pow(R1, 3) / (r * r));
			VxxSphericalShell = -4 * pi / 3.0 * (1 - pow(R1, 3) / (r * r * r));
			VyySphericalShell = VxxSphericalShell;
			VzzSphericalShell = -4 * pi / 3.0 * (1 + 2 * pow(R1, 3) / (r * r * r));
			VxxzSphericalShell = -4 * pi * pow(R1, 3) / (r * r * r * r);
			VyyzSphericalShell = VxxzSphericalShell;
			VzzzSphericalShell = 8 * pi * pow(R1, 3) / (r * r * r * r);
		}

		VxSphericalShell = 0;
		VySphericalShell = 0;

		VxySphericalShell = 0;
		VzxSphericalShell = 0;
		VzySphericalShell = 0;

		VxxxSphericalShell = 0;
		VxxySphericalShell = 0;
		VxyzSphericalShell = 0;
		VyyxSphericalShell = 0;
		VyyySphericalShell = 0;
		VzzxSphericalShell = 0;
		VzzySphericalShell = 0;

		VSphericalShell /= (r * r * RadianCorrection);
		VzSphericalShell /= (r * RadianCorrection);
		VxxSphericalShell /= RadianCorrection;
		VyySphericalShell /= RadianCorrection;
		VzzSphericalShell /= RadianCorrection;
		VxxzSphericalShell *= (r / RadianCorrection);
		VyyzSphericalShell *= (r / RadianCorrection);
		VzzzSphericalShell *= (r / RadianCorrection);

		VSplit = 0;

		VxSplit = 0;
		VySplit = 0;
		VzSplit = 0;

		VxxSplit = 0;
		VxySplit = 0;
		VyySplit = 0;
		VzxSplit = 0;
		VzySplit = 0;
		VzzSplit = 0;

		VxxxSplit = 0;
		VxxySplit = 0;
		VxxzSplit = 0;
		VxyzSplit = 0;
		VyyxSplit = 0;
		VyyySplit = 0;
		VyyzSplit = 0;
		VzzxSplit = 0;
		VzzySplit = 0;
		VzzzSplit = 0;

		for (i = 0; i < NumDensityModel; i++) {

			VSplit = VSplit + integral2t(func_V, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);

			VxSplit = VxSplit + integral2t(func_Vx, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VySplit = VySplit + integral2t(func_Vy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzSplit = VzSplit + integral2t(func_Vz, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);

			VxxSplit = VxxSplit + integral2t(func_Vxx, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VxySplit = VxySplit + integral2t(func_Vxy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VyySplit = VyySplit + integral2t(func_Vyy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzxSplit = VzxSplit + integral2t(func_Vzx, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzySplit = VzySplit + integral2t(func_Vzy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzzSplit = VzzSplit + integral2t(func_Vzz, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);

			VxxxSplit = VxxxSplit + integral2t(func_Vxxx, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VxxySplit = VxxySplit + integral2t(func_Vxxy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VxxzSplit = VxxzSplit + integral2t(func_Vxxz, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VxyzSplit = VxyzSplit + integral2t(func_Vxyz, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VyyxSplit = VyyxSplit + integral2t(func_Vyyx, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VyyySplit = VyyySplit + integral2t(func_Vyyy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VyyzSplit = VyyzSplit + integral2t(func_Vyyz, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzzxSplit = VzzxSplit + integral2t(func_Vzzx, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzzySplit = VzzySplit + integral2t(func_Vzzy, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);
			VzzzSplit = VzzzSplit + integral2t(func_Vzzz, subTesseroids[i][0], subTesseroids[i][1], subTesseroids[i][2], subTesseroids[i][3], AbsTol, RelTol);

		}

		*V = VSphericalShell - VSplit;

		*Vx = VxSphericalShell - VxSplit;
		*Vy = VySphericalShell - VySplit;
		*Vz = VzSphericalShell - VzSplit;

		*Vxx = VxxSphericalShell - VxxSplit;
		*Vxy = VxySphericalShell - VxySplit;
		*Vyy = VyySphericalShell - VyySplit;
		*Vzx = VzxSphericalShell - VzxSplit;
		*Vzy = VzySphericalShell - VzySplit;
		*Vzz = VzzSphericalShell - VzzSplit;

		*Vxxx = VxxxSphericalShell - VxxxSplit;
		*Vxxy = VxxySphericalShell - VxxySplit;
		*Vxxz = VxxzSphericalShell - VxxzSplit;
		*Vxyz = VxyzSphericalShell - VxyzSplit;
		*Vyyx = VyyxSphericalShell - VyyxSplit;
		*Vyyy = VyyySphericalShell - VyyySplit;
		*Vyyz = VyyzSphericalShell - VyyzSplit;
		*Vzzx = VzzxSphericalShell - VzzxSplit;
		*Vzzy = VzzySphericalShell - VzzySplit;
		*Vzzz = VzzzSphericalShell - VzzzSplit;

		return;

	}

	*V = integral2t(func_V, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vx = integral2t(func_Vx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vy = integral2t(func_Vy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vz = integral2t(func_Vz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vxx = integral2t(func_Vxx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vxy = integral2t(func_Vxy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vyy = integral2t(func_Vyy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vzx = integral2t(func_Vzx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vzy = integral2t(func_Vzy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vzz = integral2t(func_Vzz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vxxx = integral2t(func_Vxxx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vxxy = integral2t(func_Vxxy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vxxz = integral2t(func_Vxxz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vxyz = integral2t(func_Vxyz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vyyx = integral2t(func_Vyyx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vyyy = integral2t(func_Vyyy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vyyz = integral2t(func_Vyyz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vzzx = integral2t(func_Vzzx, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vzzy = integral2t(func_Vzzy, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);
	*Vzzz = integral2t(func_Vzzz, Fai1, Fai2, Lamda1, Lamda2, AbsTol, RelTol);

	return;
}

inline void Tesseroid_GravityEstimate(int NumTesseroid, double* Fai1, double* Fai2, double* Lamda1, double* Lamda2, double* R1, double* R2, double* rho,
	int NumPoint, double* Latitude, double* Longitude, double* Radius,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz,
	double* Vxxx, double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy, double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz)
{
	int i, j;
	double pi, RadianCorrection;
	double Reportbar;
	double VTemp,
		VxTemp, VyTemp, VzTemp,
		VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp,
		VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp;

	time_t start, end;

	printf("\nEstimating 0-, 1-, 2- and 3- gradient tensor components of the tesseroid(s)  \n");
	printf("Serial computing mode \n\n");

	printf("\n%i tesseroid(s) included !\n", NumTesseroid);
	printf("%i coordinate(s) included !\n", NumPoint);
	printf("Absolute tolerance is set to: %e\n", AbsTol);
	printf("Relative tolerance is set to: %e\n\n", RelTol);

	pi = 4 * atan(1);
	RadianCorrection = (2 * pi * pi) / (180 * 360);

	Reportbar = 0.05;

	time(&start);

	printf("\nTotal progress :   |||||||||||||||||||| \nCurrent progress : ");
	for (i = 0; i < NumPoint; i++) {
		for (j = 0; j < NumTesseroid; j++) {

			Tesseroid_GravityPointEstimate(Fai1[j], Fai2[j], Lamda1[j], Lamda2[j], R1[j], R2[j],
				Latitude[i], Longitude[i], Radius[i],
				AbsTol, RelTol,
				&VTemp,
				&VxTemp, &VyTemp, &VzTemp,
				&VxxTemp, &VxyTemp, &VyyTemp, &VzxTemp, &VzyTemp, &VzzTemp,
				&VxxxTemp, &VxxyTemp, &VxxzTemp, &VxyzTemp, &VyyxTemp, &VyyyTemp, &VyyzTemp, &VzzxTemp, &VzzyTemp, &VzzzTemp);

			V[i] = V[i] + VTemp * rho[j];

			Vx[i] = Vx[i] + VxTemp * rho[j];
			Vy[i] = Vy[i] + VyTemp * rho[j];
			Vz[i] = Vz[i] + VzTemp * rho[j];

			Vxx[i] = Vxx[i] + VxxTemp * rho[j];
			Vxy[i] = Vxy[i] + VxyTemp * rho[j];
			Vyy[i] = Vyy[i] + VyyTemp * rho[j];
			Vzx[i] = Vzx[i] + VzxTemp * rho[j];
			Vzy[i] = Vzy[i] + VzyTemp * rho[j];
			Vzz[i] = Vzz[i] + VzzTemp * rho[j];

			Vxxx[i] = Vxxx[i] + VxxxTemp * rho[j];
			Vxxy[i] = Vxxy[i] + VxxyTemp * rho[j];
			Vxxz[i] = Vxxz[i] + VxxzTemp * rho[j];
			Vxyz[i] = Vxyz[i] + VxyzTemp * rho[j];
			Vyyx[i] = Vyyx[i] + VyyxTemp * rho[j];
			Vyyy[i] = Vyyy[i] + VyyyTemp * rho[j];
			Vyyz[i] = Vyyz[i] + VyyzTemp * rho[j];
			Vzzx[i] = Vzzx[i] + VzzxTemp * rho[j];
			Vzzy[i] = Vzzy[i] + VzzyTemp * rho[j];
			Vzzz[i] = Vzzz[i] + VzzzTemp * rho[j];

		}

		if ((double)i / (double)NumPoint >= Reportbar)
		{
			printf("|");
			Reportbar = Reportbar + 0.05;
		}
	}

	printf("|");

	time(&end);

	printf("\nTime for calculation: %.15lf seconds\n\n", difftime(end, start));

	for (i = 0; i < NumPoint; i++) {

		V[i] = V[i] * Radius[i] * Radius[i] * RadianCorrection;

		Vx[i] = Vx[i] * Radius[i] * RadianCorrection;
		Vy[i] = Vy[i] * Radius[i] * RadianCorrection;
		Vz[i] = Vz[i] * Radius[i] * RadianCorrection;

		Vxx[i] = Vxx[i] * RadianCorrection;
		Vxy[i] = Vxy[i] * RadianCorrection;
		Vyy[i] = Vyy[i] * RadianCorrection;
		Vzx[i] = Vzx[i] * RadianCorrection;
		Vzy[i] = Vzy[i] * RadianCorrection;
		Vzz[i] = Vzz[i] * RadianCorrection;

		Vxxx[i] = Vxxx[i] / Radius[i] * RadianCorrection;
		Vxxy[i] = Vxxy[i] / Radius[i] * RadianCorrection;
		Vxxz[i] = Vxxz[i] / Radius[i] * RadianCorrection;
		Vxyz[i] = Vxyz[i] / Radius[i] * RadianCorrection;
		Vyyx[i] = Vyyx[i] / Radius[i] * RadianCorrection;
		Vyyy[i] = Vyyy[i] / Radius[i] * RadianCorrection;
		Vyyz[i] = Vyyz[i] / Radius[i] * RadianCorrection;
		Vzzx[i] = Vzzx[i] / Radius[i] * RadianCorrection;
		Vzzy[i] = Vzzy[i] / Radius[i] * RadianCorrection;
		Vzzz[i] = Vzzz[i] / Radius[i] * RadianCorrection;

	}

	printf("\nEnd of calculation, function: Tesseroid_GravityEstimate, Exit!\n\n");

}

inline void Tesseroid_MagneticEstimate(int NumTesseroid, double* Fai1, double* Fai2, double* Lamda1, double* Lamda2, double* R1, double* R2, double* Mx, double* My, double* Mz,
	int NumPoint, double* Latitude, double* Longitude, double* Radius,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz)
{
	int i, j;
	double pi, RadianCorrection;
	double Reportbar;
	double Beita, cosBeita, sinBeita, cosDataLat, sinDataLat, cosLat, sinLat;
	double MxLocal, MyLocal, MzLocal;
	double VTemp,
		VxTemp, VyTemp, VzTemp,
		VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp,
		VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp;

	time_t start, end;

	printf("\nEstimating 0-, 1- and 2- gradient tensor components of the tesseroid(s)  \n");
	printf("Serial computing mode \n\n");

	printf("\n%i tesseroid(s) included !\n", NumTesseroid);
	printf("%i coordinate(s) included !\n", NumPoint);
	printf("Absolute tolerance is set to: %e\n", AbsTol);
	printf("Relative tolerance is set to: %e\n\n", RelTol);

	pi = 4 * atan(1);
	RadianCorrection = (2 * pi * pi) / (180 * 360);

	Reportbar = 0.05;

	time(&start);

	printf("\nTotal progress :   |||||||||||||||||||| \nCurrent progress : ");
	for (i = 0; i < NumPoint; i++) {
		for (j = 0; j < NumTesseroid; j++) {

			Beita = (Lamda2[j] + Lamda1[j]) / 2 - Longitude[i];
			cosBeita = cos(Beita / 180 * pi);
			sinBeita = sin(Beita / 180 * pi);
			cosDataLat = cos(Latitude[i] / 180 * pi);
			sinDataLat = sin(Latitude[i] / 180 * pi);
			cosLat = cos((Fai2[j] + Fai1[j]) / 2 / 180 * pi);
			sinLat = sin((Fai2[j] + Fai1[j]) / 2 / 180 * pi);

			MxLocal = (cosBeita * sinDataLat * sinLat + cosDataLat * cosLat) * Mx[j] +
				(sinBeita * sinDataLat) * My[j] +
				(cosBeita * sinDataLat * cosLat - cosDataLat * sinLat) * Mz[j];
			MyLocal = (-sinBeita * sinLat) * Mx[j] +
				cosBeita * My[j] -
				sinBeita * cosLat * Mz[j];
			MzLocal = (cosBeita * cosDataLat * sinLat - sinDataLat * cosLat) * Mx[j] +
				(sinBeita * cosDataLat) * My[j] +
				(cosBeita * cosDataLat * cosLat + sinDataLat * sinLat) * Mz[j];

			Tesseroid_GravityPointEstimate(Fai1[j], Fai2[j], Lamda1[j], Lamda2[j], R1[j], R2[j],
				Latitude[i], Longitude[i], Radius[i],
				AbsTol, RelTol,
				&VTemp,
				&VxTemp, &VyTemp, &VzTemp,
				&VxxTemp, &VxyTemp, &VyyTemp, &VzxTemp, &VzyTemp, &VzzTemp,
				&VxxxTemp, &VxxyTemp, &VxxzTemp, &VxyzTemp, &VyyxTemp, &VyyyTemp, &VyyzTemp, &VzzxTemp, &VzzyTemp, &VzzzTemp);

			V[i] = V[i] + MxLocal * VxTemp + MyLocal * VyTemp + MzLocal * VzTemp;

			Vx[i] = Vx[i] + MxLocal * VxxTemp + MyLocal * VxyTemp + MzLocal * VzxTemp;
			Vy[i] = Vy[i] + MxLocal * VxyTemp + MyLocal * VyyTemp + MzLocal * VzyTemp;
			Vz[i] = Vz[i] + MxLocal * VzxTemp + MyLocal * VzyTemp + MzLocal * VzzTemp;

			Vxx[i] = Vxx[i] + MxLocal * VxxxTemp + MyLocal * VxxyTemp + MzLocal * VxxzTemp;
			Vxy[i] = Vxy[i] + MxLocal * VxxyTemp + MyLocal * VyyxTemp + MzLocal * VxyzTemp;
			Vyy[i] = Vyy[i] + MxLocal * VyyxTemp + MyLocal * VyyyTemp + MzLocal * VyyzTemp;
			Vzx[i] = Vzx[i] + MxLocal * VxxzTemp + MyLocal * VxyzTemp + MzLocal * VzzxTemp;
			Vzy[i] = Vzy[i] + MxLocal * VxyzTemp + MyLocal * VyyzTemp + MzLocal * VzzyTemp;
			Vzz[i] = Vzz[i] + MxLocal * VzzxTemp + MyLocal * VzzyTemp + MzLocal * VzzzTemp;

		}

		if ((double)i / (double)NumPoint >= Reportbar)
		{
			printf("|");
			Reportbar = Reportbar + 0.05;
		}
	}

	printf("|");

	time(&end);

	printf("\nTime for calculation: %.15lf seconds\n\n", difftime(end, start));

	for (i = 0; i < NumPoint; i++) {

		V[i] = V[i] * Radius[i] * Radius[i] * RadianCorrection;

		Vx[i] = Vx[i] * Radius[i] * RadianCorrection;
		Vy[i] = Vy[i] * Radius[i] * RadianCorrection;
		Vz[i] = Vz[i] * Radius[i] * RadianCorrection;

		Vxx[i] = Vxx[i] * RadianCorrection;
		Vxy[i] = Vxy[i] * RadianCorrection;
		Vyy[i] = Vyy[i] * RadianCorrection;
		Vzx[i] = Vzx[i] * RadianCorrection;
		Vzy[i] = Vzy[i] * RadianCorrection;
		Vzz[i] = Vzz[i] * RadianCorrection;

	}

	printf("\nEnd of calculation, function: Tesseroid_MagneticEstimate, Exit!\n\n");

}


inline void Tesseroid_GravityEstimateOpenMP(int NumTesseroid, double* Fai1, double* Fai2, double* Lamda1, double* Lamda2, double* R1, double* R2, double* rho,
	int NumPoint, double* Latitude, double* Longitude, double* Radius,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz,
	double* Vxxx, double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy, double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz)
{
	int i, j;
	double pi, RadianCorrection;
	double Reportbar;
	double VTemp,
		VxTemp, VyTemp, VzTemp,
		VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp,
		VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp;

	time_t start, end;

	printf("\nEstimating 0-, 1-, 2- and 3- gradient tensor components of the tesseroid(s)  \n");
	printf("OpenMP parallel computing mode \n\n");

	printf("\n%i tesseroid(s) included !\n", NumTesseroid);
	printf("%i coordinate(s) included !\n", NumPoint);
	printf("Absolute tolerance is set to: %e\n", AbsTol);
	printf("Relative tolerance is set to: %e\n\n", RelTol);

	pi = 4 * atan(1);
	RadianCorrection = (2 * pi * pi) / (180 * 360);

	Reportbar = 0.1;

	time(&start);

#pragma omp parallel for default(shared) \
private(i, j, VTemp,VxTemp, VyTemp, VzTemp, \
			VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp, \
			VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp)
	for (i = 0; i < NumPoint; i++) {
		for (j = 0; j < NumTesseroid; j++) {

			Tesseroid_GravityPointEstimate(Fai1[j], Fai2[j], Lamda1[j], Lamda2[j], R1[j], R2[j],
				Latitude[i], Longitude[i], Radius[i],
				AbsTol, RelTol,
				&VTemp,
				&VxTemp, &VyTemp, &VzTemp,
				&VxxTemp, &VxyTemp, &VyyTemp, &VzxTemp, &VzyTemp, &VzzTemp,
				&VxxxTemp, &VxxyTemp, &VxxzTemp, &VxyzTemp, &VyyxTemp, &VyyyTemp, &VyyzTemp, &VzzxTemp, &VzzyTemp, &VzzzTemp);
#pragma omp atomic
			V[i] += VTemp * rho[j];
#pragma omp atomic
			Vx[i] += VxTemp * rho[j];
#pragma omp atomic
			Vy[i] += VyTemp * rho[j];
#pragma omp atomic
			Vz[i] += VzTemp * rho[j];

#pragma omp atomic
			Vxx[i] += VxxTemp * rho[j];
#pragma omp atomic
			Vxy[i] += VxyTemp * rho[j];
#pragma omp atomic
			Vyy[i] += VyyTemp * rho[j];
#pragma omp atomic
			Vzx[i] += VzxTemp * rho[j];
#pragma omp atomic
			Vzy[i] += VzyTemp * rho[j];
#pragma omp atomic
			Vzz[i] += VzzTemp * rho[j];

#pragma omp atomic
			Vxxx[i] += VxxxTemp * rho[j];
#pragma omp atomic
			Vxxy[i] += VxxyTemp * rho[j];
#pragma omp atomic
			Vxxz[i] += VxxzTemp * rho[j];
#pragma omp atomic
			Vxyz[i] += VxyzTemp * rho[j];
#pragma omp atomic
			Vyyx[i] += VyyxTemp * rho[j];
#pragma omp atomic
			Vyyy[i] += VyyyTemp * rho[j];
#pragma omp atomic
			Vyyz[i] += VyyzTemp * rho[j];
#pragma omp atomic
			Vzzx[i] += VzzxTemp * rho[j];
#pragma omp atomic
			Vzzy[i] += VzzyTemp * rho[j];
#pragma omp atomic
			Vzzz[i] += VzzzTemp * rho[j];

		}
	}

	time(&end);

	printf("\nTime for calculation: %.15lf seconds\n\n", difftime(end, start));
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < NumPoint; i++) {

		V[i] = V[i] * Radius[i] * Radius[i] * RadianCorrection;

		Vx[i] = Vx[i] * Radius[i] * RadianCorrection;
		Vy[i] = Vy[i] * Radius[i] * RadianCorrection;
		Vz[i] = Vz[i] * Radius[i] * RadianCorrection;

		Vxx[i] = Vxx[i] * RadianCorrection;
		Vxy[i] = Vxy[i] * RadianCorrection;
		Vyy[i] = Vyy[i] * RadianCorrection;
		Vzx[i] = Vzx[i] * RadianCorrection;
		Vzy[i] = Vzy[i] * RadianCorrection;
		Vzz[i] = Vzz[i] * RadianCorrection;

		Vxxx[i] = Vxxx[i] / Radius[i] * RadianCorrection;
		Vxxy[i] = Vxxy[i] / Radius[i] * RadianCorrection;
		Vxxz[i] = Vxxz[i] / Radius[i] * RadianCorrection;
		Vxyz[i] = Vxyz[i] / Radius[i] * RadianCorrection;
		Vyyx[i] = Vyyx[i] / Radius[i] * RadianCorrection;
		Vyyy[i] = Vyyy[i] / Radius[i] * RadianCorrection;
		Vyyz[i] = Vyyz[i] / Radius[i] * RadianCorrection;
		Vzzx[i] = Vzzx[i] / Radius[i] * RadianCorrection;
		Vzzy[i] = Vzzy[i] / Radius[i] * RadianCorrection;
		Vzzz[i] = Vzzz[i] / Radius[i] * RadianCorrection;

	}

	printf("\nEnd of calculation, function: Tesseroid_GravityEstimateOpenMP, Exit!\n\n");

}

inline void Tesseroid_MagneticEstimateOpenMP(int NumTesseroid, double* Fai1, double* Fai2, double* Lamda1, double* Lamda2, double* R1, double* R2, double* Mx, double* My, double* Mz,
	int NumPoint, double* Latitude, double* Longitude, double* Radius,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz)
{
	int i, j;
	double pi, RadianCorrection;
	double Reportbar;
	double Beita, cosBeita, sinBeita, cosDataLat, sinDataLat, cosLat, sinLat;
	double MxLocal, MyLocal, MzLocal;
	double VTemp,
		VxTemp, VyTemp, VzTemp,
		VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp,
		VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp;

	time_t start, end;

	printf("\nEstimating 0-, 1- and 2- gradient tensor components of the tesseroid(s)  \n");
	printf("OpenMP parallel computing mode \n\n");

	printf("\n%i tesseroid(s) included !\n", NumTesseroid);
	printf("%i coordinate(s) included !\n", NumPoint);
	printf("Absolute tolerance is set to: %e\n", AbsTol);
	printf("Relative tolerance is set to: %e\n\n", RelTol);

	pi = 4 * atan(1);
	RadianCorrection = (2 * pi * pi) / (180 * 360);

	Reportbar = 0.1;

	time(&start);

#pragma omp parallel for default(shared) \
private(i, j, \
			Beita, cosBeita, sinBeita, cosDataLat, sinDataLat, cosLat, sinLat, MxLocal, MyLocal, MzLocal, \
			VTemp,VxTemp, VyTemp, VzTemp, \
			VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp, \
			VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp)
	for (i = 0; i < NumPoint; i++) {
		for (j = 0; j < NumTesseroid; j++) {

			Beita = (Lamda2[j] + Lamda1[j]) / 2 - Longitude[i];
			cosBeita = cos(Beita / 180 * pi);
			sinBeita = sin(Beita / 180 * pi);
			cosDataLat = cos(Latitude[i] / 180 * pi);
			sinDataLat = sin(Latitude[i] / 180 * pi);
			cosLat = cos((Fai2[j] + Fai1[j]) / 2 / 180 * pi);
			sinLat = sin((Fai2[j] + Fai1[j]) / 2 / 180 * pi);

			MxLocal = (cosBeita * sinDataLat * sinLat + cosDataLat * cosLat) * Mx[j] +
				(sinBeita * sinDataLat) * My[j] +
				(cosBeita * sinDataLat * cosLat - cosDataLat * sinLat) * Mz[j];
			MyLocal = (-sinBeita * sinLat) * Mx[j] +
				cosBeita * My[j] -
				sinBeita * cosLat * Mz[j];
			MzLocal = (cosBeita * cosDataLat * sinLat - sinDataLat * cosLat) * Mx[j] +
				(sinBeita * cosDataLat) * My[j] +
				(cosBeita * cosDataLat * cosLat + sinDataLat * sinLat) * Mz[j];

			Tesseroid_GravityPointEstimate(Fai1[j], Fai2[j], Lamda1[j], Lamda2[j], R1[j], R2[j],
				Latitude[i], Longitude[i], Radius[i],
				AbsTol, RelTol,
				&VTemp,
				&VxTemp, &VyTemp, &VzTemp,
				&VxxTemp, &VxyTemp, &VyyTemp, &VzxTemp, &VzyTemp, &VzzTemp,
				&VxxxTemp, &VxxyTemp, &VxxzTemp, &VxyzTemp, &VyyxTemp, &VyyyTemp, &VyyzTemp, &VzzxTemp, &VzzyTemp, &VzzzTemp);

#pragma omp atomic
			V[i] += MxLocal * VxTemp + MyLocal * VyTemp + MzLocal * VzTemp;

#pragma omp atomic
			Vx[i] += MxLocal * VxxTemp + MyLocal * VxyTemp + MzLocal * VzxTemp;
#pragma omp atomic
			Vy[i] += MxLocal * VxyTemp + MyLocal * VyyTemp + MzLocal * VzyTemp;
#pragma omp atomic
			Vz[i] += MxLocal * VzxTemp + MyLocal * VzyTemp + MzLocal * VzzTemp;

#pragma omp atomic
			Vxx[i] += MxLocal * VxxxTemp + MyLocal * VxxyTemp + MzLocal * VxxzTemp;
#pragma omp atomic
			Vxy[i] += MxLocal * VxxyTemp + MyLocal * VyyxTemp + MzLocal * VxyzTemp;
#pragma omp atomic
			Vyy[i] += MxLocal * VyyxTemp + MyLocal * VyyyTemp + MzLocal * VyyzTemp;
#pragma omp atomic
			Vzx[i] += MxLocal * VxxzTemp + MyLocal * VxyzTemp + MzLocal * VzzxTemp;
#pragma omp atomic
			Vzy[i] += MxLocal * VxyzTemp + MyLocal * VyyzTemp + MzLocal * VzzyTemp;
#pragma omp atomic
			Vzz[i] += MxLocal * VzzxTemp + MyLocal * VzzyTemp + MzLocal * VzzzTemp;
		}
	}

	time(&end);

	printf("\nTime for calculation: %.15lf seconds\n\n", difftime(end, start));
#pragma omp parallel for default(shared) private(i)
	for (i = 0; i < NumPoint; i++) {

		V[i] = V[i] * Radius[i] * Radius[i] * RadianCorrection;

		Vx[i] = Vx[i] * Radius[i] * RadianCorrection;
		Vy[i] = Vy[i] * Radius[i] * RadianCorrection;
		Vz[i] = Vz[i] * Radius[i] * RadianCorrection;

		Vxx[i] = Vxx[i] * RadianCorrection;
		Vxy[i] = Vxy[i] * RadianCorrection;
		Vyy[i] = Vyy[i] * RadianCorrection;
		Vzx[i] = Vzx[i] * RadianCorrection;
		Vzy[i] = Vzy[i] * RadianCorrection;
		Vzz[i] = Vzz[i] * RadianCorrection;

	}

	printf("\nEnd of calculation, function: Tesseroid_MagneticEstimateOpenMP, Exit!\n\n");

}

inline void Tesseroid_GravityEstimateMPI(int NumTesseroid, double* Fai1, double* Fai2, double* Lamda1, double* Lamda2, double* R1, double* R2, double* rho,
	int NumPoint, double* Latitude, double* Longitude, double* Radius,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz,
	double* Vxxx, double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy, double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz)
{
	int i, j;
	int rank, nprocs;
	int Locals, Locale;
	int Position;

	int* DataList, * VectorShiftList;

	double pi, RadianCorrection;
	double t_Start, t_End;
	double VTemp,
		VxTemp, VyTemp, VzTemp,
		VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp,
		VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp;
	double* VLocal,
		* VxLocal, * VyLocal, * VzLocal,
		* VxxLocal, * VxyLocal, * VyyLocal, * VzxLocal, * VzyLocal, * VzzLocal,
		* VxxxLocal, * VxxyLocal, * VxxzLocal, * VxyzLocal, * VyyxLocal, * VyyyLocal, * VyyzLocal, * VzzxLocal, * VzzyLocal, * VzzzLocal;

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {

		printf("\nEstimating 0-, 1-, 2- and 3- gradient tensor components of the tesseroid(s)  \n");
		printf("MPI parallel computing mode \n");

		printf("\n%i tesseroid(s) included !\n", NumTesseroid);
		printf("%i coordinate(s) included !\n", NumPoint);
		printf("Absolute tolerance is set to: %e\n", AbsTol);
		printf("Relative tolerance is set to: %e\n\n", RelTol);

	}

	t_Start = MPI_Wtime();

	pi = 4 * atan(1);
	RadianCorrection = (2 * pi * pi) / (180 * 360);

	DataList = (int*)malloc(sizeof(int) * nprocs);
	VectorShiftList = (int*)malloc(sizeof(int) * nprocs);

	MPIDistribute(NumPoint, DataList, VectorShiftList);

	MPI_Barrier(MPI_COMM_WORLD);

	Locals = VectorShiftList[rank];
	Locale = VectorShiftList[rank] + DataList[rank] - 1;

	VLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	VxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	VxxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VxyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzzLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	VxxxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VxxyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VxxzLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VxyzLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyyxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyyyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyyzLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzzxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzzyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzzzLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	MPI_Barrier(MPI_COMM_WORLD);

	for (i = Locals; i <= Locale; i++) {

		Position = i - VectorShiftList[rank];

		VLocal[Position] = VxLocal[Position] = VyLocal[Position] = VzLocal[Position] =
			VxxLocal[Position] = VxyLocal[Position] = VyyLocal[Position] = VzxLocal[Position] = VzyLocal[Position] = VzzLocal[Position] =
			VxxxLocal[Position] = VxxyLocal[Position] = VxxzLocal[Position] = VxyzLocal[Position] = VyyxLocal[Position] =
			VyyyLocal[Position] = VyyzLocal[Position] = VzzxLocal[Position] = VzzyLocal[Position] = VzzzLocal[Position] = 0;

		for (j = 0; j < NumTesseroid; j++) {

			Tesseroid_GravityPointEstimate(Fai1[j], Fai2[j], Lamda1[j], Lamda2[j], R1[j], R2[j],
				Latitude[i], Longitude[i], Radius[i],
				AbsTol, RelTol,
				&VTemp,
				&VxTemp, &VyTemp, &VzTemp,
				&VxxTemp, &VxyTemp, &VyyTemp, &VzxTemp, &VzyTemp, &VzzTemp,
				&VxxxTemp, &VxxyTemp, &VxxzTemp, &VxyzTemp, &VyyxTemp, &VyyyTemp, &VyyzTemp, &VzzxTemp, &VzzyTemp, &VzzzTemp);

			VLocal[Position] += VTemp * rho[j];

			VxLocal[Position] += VxTemp * rho[j];
			VyLocal[Position] += VyTemp * rho[j];
			VzLocal[Position] += VzTemp * rho[j];

			VxxLocal[Position] += VxxTemp * rho[j];
			VxyLocal[Position] += VxyTemp * rho[j];
			VyyLocal[Position] += VyyTemp * rho[j];
			VzxLocal[Position] += VzxTemp * rho[j];
			VzyLocal[Position] += VzyTemp * rho[j];
			VzzLocal[Position] += VzzTemp * rho[j];

			VxxxLocal[Position] += VxxxTemp * rho[j];
			VxxyLocal[Position] += VxxyTemp * rho[j];
			VxxzLocal[Position] += VxxzTemp * rho[j];
			VxyzLocal[Position] += VxyzTemp * rho[j];
			VyyxLocal[Position] += VyyxTemp * rho[j];
			VyyyLocal[Position] += VyyyTemp * rho[j];
			VyyzLocal[Position] += VyyzTemp * rho[j];
			VzzxLocal[Position] += VzzxTemp * rho[j];
			VzzyLocal[Position] += VzzyTemp * rho[j];
			VzzzLocal[Position] += VzzzTemp * rho[j];

		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gatherv(VLocal, DataList[rank], MPI_DOUBLE, V, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gatherv(VxLocal, DataList[rank], MPI_DOUBLE, Vx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyLocal, DataList[rank], MPI_DOUBLE, Vy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzLocal, DataList[rank], MPI_DOUBLE, Vz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gatherv(VxxLocal, DataList[rank], MPI_DOUBLE, Vxx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VxyLocal, DataList[rank], MPI_DOUBLE, Vxy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyyLocal, DataList[rank], MPI_DOUBLE, Vyy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzxLocal, DataList[rank], MPI_DOUBLE, Vzx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzyLocal, DataList[rank], MPI_DOUBLE, Vzy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzzLocal, DataList[rank], MPI_DOUBLE, Vzz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gatherv(VxxxLocal, DataList[rank], MPI_DOUBLE, Vxxx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VxxyLocal, DataList[rank], MPI_DOUBLE, Vxxy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VxxzLocal, DataList[rank], MPI_DOUBLE, Vxxz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VxyzLocal, DataList[rank], MPI_DOUBLE, Vxyz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyyxLocal, DataList[rank], MPI_DOUBLE, Vyyx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyyyLocal, DataList[rank], MPI_DOUBLE, Vyyy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyyzLocal, DataList[rank], MPI_DOUBLE, Vyyz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzzxLocal, DataList[rank], MPI_DOUBLE, Vzzx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzzyLocal, DataList[rank], MPI_DOUBLE, Vzzy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzzzLocal, DataList[rank], MPI_DOUBLE, Vzzz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {

		for (i = 0; i < NumPoint; i++) {

			V[i] = V[i] * Radius[i] * Radius[i] * RadianCorrection;

			Vx[i] = Vx[i] * Radius[i] * RadianCorrection;
			Vy[i] = Vy[i] * Radius[i] * RadianCorrection;
			Vz[i] = Vz[i] * Radius[i] * RadianCorrection;

			Vxx[i] = Vxx[i] * RadianCorrection;
			Vxy[i] = Vxy[i] * RadianCorrection;
			Vyy[i] = Vyy[i] * RadianCorrection;
			Vzx[i] = Vzx[i] * RadianCorrection;
			Vzy[i] = Vzy[i] * RadianCorrection;
			Vzz[i] = Vzz[i] * RadianCorrection;

			Vxxx[i] = Vxxx[i] / Radius[i] * RadianCorrection;
			Vxxy[i] = Vxxy[i] / Radius[i] * RadianCorrection;
			Vxxz[i] = Vxxz[i] / Radius[i] * RadianCorrection;
			Vxyz[i] = Vxyz[i] / Radius[i] * RadianCorrection;
			Vyyx[i] = Vyyx[i] / Radius[i] * RadianCorrection;
			Vyyy[i] = Vyyy[i] / Radius[i] * RadianCorrection;
			Vyyz[i] = Vyyz[i] / Radius[i] * RadianCorrection;
			Vzzx[i] = Vzzx[i] / Radius[i] * RadianCorrection;
			Vzzy[i] = Vzzy[i] / Radius[i] * RadianCorrection;
			Vzzz[i] = Vzzz[i] / Radius[i] * RadianCorrection;

		}

	}

	t_End = MPI_Wtime();

	if (rank == 0) {
		printf("Computation Time (Units in Second) : %.2f", (t_End - t_Start));
		printf("\nEnd of calculation, function: Tesseroid_GravityEstimateMPI, Exit!\n\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	free(DataList);
	free(VectorShiftList);

	free(VLocal);
	free(VxLocal); free(VyLocal); free(VzLocal);
	free(VxxLocal); free(VxyLocal); free(VyyLocal); free(VzxLocal); free(VzyLocal); free(VzzLocal);
	free(VxxxLocal); free(VxxyLocal); free(VxxzLocal); free(VxyzLocal);
	free(VyyxLocal); free(VyyyLocal); free(VyyzLocal);
	free(VzzxLocal); free(VzzyLocal); free(VzzzLocal);

}

inline void Tesseroid_MagneticEstimateMPI(int NumTesseroid, double* Fai1, double* Fai2, double* Lamda1, double* Lamda2, double* R1, double* R2, double* Mx, double* My, double* Mz,
	int NumPoint, double* Latitude, double* Longitude, double* Radius,
	double AbsTol, double RelTol,
	double* V,
	double* Vx, double* Vy, double* Vz,
	double* Vxx, double* Vxy, double* Vyy, double* Vzx, double* Vzy, double* Vzz)
{
	int i, j;
	int rank, nprocs;
	int Locals, Locale;
	int Position;

	int* DataList, * VectorShiftList;

	double pi, RadianCorrection;
	double t_Start, t_End;
	double Beita, cosBeita, sinBeita, cosDataLat, sinDataLat, cosLat, sinLat;
	double MxLocal, MyLocal, MzLocal;
	double VTemp,
		VxTemp, VyTemp, VzTemp,
		VxxTemp, VxyTemp, VyyTemp, VzxTemp, VzyTemp, VzzTemp,
		VxxxTemp, VxxyTemp, VxxzTemp, VxyzTemp, VyyxTemp, VyyyTemp, VyyzTemp, VzzxTemp, VzzyTemp, VzzzTemp;
	double* VLocal,
		* VxLocal, * VyLocal, * VzLocal,
		* VxxLocal, * VxyLocal, * VyyLocal, * VzxLocal, * VzyLocal, * VzzLocal;

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {

		printf("\nEstimating 0-, 1- and 2- gradient tensor components of the tesseroid(s)  \n");
		printf("MPI parallel computing mode \n");

		printf("\n%i tesseroid(s) included !\n", NumTesseroid);
		printf("%i coordinate(s) included !\n", NumPoint);
		printf("Absolute tolerance is set to: %e\n", AbsTol);
		printf("Relative tolerance is set to: %e\n\n", RelTol);

	}

	t_Start = MPI_Wtime();

	pi = 4 * atan(1);
	RadianCorrection = (2 * pi * pi) / (180 * 360);

	DataList = (int*)malloc(sizeof(int) * nprocs);
	VectorShiftList = (int*)malloc(sizeof(int) * nprocs);

	MPIDistribute(NumPoint, DataList, VectorShiftList);

	MPI_Barrier(MPI_COMM_WORLD);

	Locals = VectorShiftList[rank];
	Locale = VectorShiftList[rank] + DataList[rank] - 1;

	VLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	VxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	VxxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VxyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VyyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzxLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzyLocal = (double*)malloc(sizeof(double) * DataList[rank]);
	VzzLocal = (double*)malloc(sizeof(double) * DataList[rank]);

	MPI_Barrier(MPI_COMM_WORLD);

	for (i = Locals; i <= Locale; i++) {

		Position = i - VectorShiftList[rank];

		VLocal[Position] = VxLocal[Position] = VyLocal[Position] = VzLocal[Position] =
			VxxLocal[Position] = VxyLocal[Position] = VyyLocal[Position] = VzxLocal[Position] = VzyLocal[Position] = VzzLocal[Position] = 0;

		for (j = 0; j < NumTesseroid; j++) {

			Beita = (Lamda2[j] + Lamda1[j]) / 2 - Longitude[i];
			cosBeita = cos(Beita / 180 * pi);
			sinBeita = sin(Beita / 180 * pi);
			cosDataLat = cos(Latitude[i] / 180 * pi);
			sinDataLat = sin(Latitude[i] / 180 * pi);
			cosLat = cos((Fai2[j] + Fai1[j]) / 2 / 180 * pi);
			sinLat = sin((Fai2[j] + Fai1[j]) / 2 / 180 * pi);

			MxLocal = (cosBeita * sinDataLat * sinLat + cosDataLat * cosLat) * Mx[j] +
				(sinBeita * sinDataLat) * My[j] +
				(cosBeita * sinDataLat * cosLat - cosDataLat * sinLat) * Mz[j];
			MyLocal = (-sinBeita * sinLat) * Mx[j] +
				cosBeita * My[j] -
				sinBeita * cosLat * Mz[j];
			MzLocal = (cosBeita * cosDataLat * sinLat - sinDataLat * cosLat) * Mx[j] +
				(sinBeita * cosDataLat) * My[j] +
				(cosBeita * cosDataLat * cosLat + sinDataLat * sinLat) * Mz[j];

			Tesseroid_GravityPointEstimate(Fai1[j], Fai2[j], Lamda1[j], Lamda2[j], R1[j], R2[j],
				Latitude[i], Longitude[i], Radius[i],
				AbsTol, RelTol,
				&VTemp,
				&VxTemp, &VyTemp, &VzTemp,
				&VxxTemp, &VxyTemp, &VyyTemp, &VzxTemp, &VzyTemp, &VzzTemp,
				&VxxxTemp, &VxxyTemp, &VxxzTemp, &VxyzTemp, &VyyxTemp, &VyyyTemp, &VyyzTemp, &VzzxTemp, &VzzyTemp, &VzzzTemp);

			VLocal[Position] += MxLocal * VxTemp + MyLocal * VyTemp + MzLocal * VzTemp;

			VxLocal[Position] += MxLocal * VxxTemp + MyLocal * VxyTemp + MzLocal * VzxTemp;
			VyLocal[Position] += MxLocal * VxyTemp + MyLocal * VyyTemp + MzLocal * VzyTemp;
			VzLocal[Position] += MxLocal * VzxTemp + MyLocal * VzyTemp + MzLocal * VzzTemp;

			VxxLocal[Position] += MxLocal * VxxxTemp + MyLocal * VxxyTemp + MzLocal * VxxzTemp;
			VxyLocal[Position] += MxLocal * VxxyTemp + MyLocal * VyyxTemp + MzLocal * VxyzTemp;
			VyyLocal[Position] += MxLocal * VyyxTemp + MyLocal * VyyyTemp + MzLocal * VyyzTemp;
			VzxLocal[Position] += MxLocal * VxxzTemp + MyLocal * VxyzTemp + MzLocal * VzzxTemp;
			VzyLocal[Position] += MxLocal * VxyzTemp + MyLocal * VyyzTemp + MzLocal * VzzyTemp;
			VzzLocal[Position] += MxLocal * VzzxTemp + MyLocal * VzzyTemp + MzLocal * VzzzTemp;

		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gatherv(VLocal, DataList[rank], MPI_DOUBLE, V, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gatherv(VxLocal, DataList[rank], MPI_DOUBLE, Vx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyLocal, DataList[rank], MPI_DOUBLE, Vy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzLocal, DataList[rank], MPI_DOUBLE, Vz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Gatherv(VxxLocal, DataList[rank], MPI_DOUBLE, Vxx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VxyLocal, DataList[rank], MPI_DOUBLE, Vxy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VyyLocal, DataList[rank], MPI_DOUBLE, Vyy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzxLocal, DataList[rank], MPI_DOUBLE, Vzx, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzyLocal, DataList[rank], MPI_DOUBLE, Vzy, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(VzzLocal, DataList[rank], MPI_DOUBLE, Vzz, DataList, VectorShiftList, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {

		for (i = 0; i < NumPoint; i++) {

			V[i] = V[i] * Radius[i] * Radius[i] * RadianCorrection;

			Vx[i] = Vx[i] * Radius[i] * RadianCorrection;
			Vy[i] = Vy[i] * Radius[i] * RadianCorrection;
			Vz[i] = Vz[i] * Radius[i] * RadianCorrection;

			Vxx[i] = Vxx[i] * RadianCorrection;
			Vxy[i] = Vxy[i] * RadianCorrection;
			Vyy[i] = Vyy[i] * RadianCorrection;
			Vzx[i] = Vzx[i] * RadianCorrection;
			Vzy[i] = Vzy[i] * RadianCorrection;
			Vzz[i] = Vzz[i] * RadianCorrection;

		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	t_End = MPI_Wtime();

	if (rank == 0) {
		printf("Computation Time (Units in Second) : %.2f", (t_End - t_Start));
		printf("\nEnd of calculation, function: Tesseroid_MagneticEstimateMPI, Exit!\n\n");
	}

	free(DataList);
	free(VectorShiftList);

	free(VLocal);
	free(VxLocal); free(VyLocal); free(VzLocal);
	free(VxxLocal); free(VxyLocal); free(VyyLocal); free(VzxLocal); free(VzyLocal); free(VzzLocal);

}


}