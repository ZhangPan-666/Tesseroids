// This cpp header file contains some of the subfunctions needed for
// calculating the surface integration kernels
// For more details, see Zhang P.et al.,
// "An alternative approach for accurately calculating gravitational and magnetic fields of a spherical prism."

#pragma once
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <Eigen/Dense>

namespace mista_math
{


inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelV(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi;
	double sin_half_Phi_sq, l2pow1, l1pow1;
	double cosPhi;
	
	Eigen::Matrix<double, 14, 14> KernelV;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);


	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);

			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelV(i, j) = -0.5 * cos_Fai * (2 * R2 + R2 * R2 - 2 * R1 - R1 * R1 + 2 * log((R2 - 1) / (R1 - 1)));
					continue;
				}
				else if (R1 < 1) {
					KernelV(i, j) = 0.5 * cos_Fai * (2 * R2 + R2 * R2 + 2 * R1 + R1 * R1 + 4 * log(r) + 2 * log((1 - R1) * (R2 - 1)));
					continue;
				}
				else {
					KernelV(i, j) = -0.5 * cos_Fai * (2 * R2 + R2 * R2 - 2 * R1 - R1 * R1 + 2 * log((R2 - 1) / (R1 - 1)));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelV(i, j) = 0.5 * cos_Fai * (-2 * R2 + R2 * R2 + 2 * R1 - R1 * R1 + 2 * log((R2 + 1) / (R1 + 1)));
				continue;
			}

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);

			cosPhi = cos(Phi);
			KernelV(i, j) = 0.5 * cos_Fai * (
				(3 * cosPhi + R2) * l2pow1 - (3 * cosPhi + R1) * l1pow1 +
				(1 - 3 * cosPhi * cosPhi) * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));
		}
	}

	return KernelV;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVx(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha;
	double sin_half_Phi_sq, l2pow1, l1pow1;
	double cscPhi, cosPhi, cos2Phi;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;
	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	Eigen::Matrix<double, 14, 14> KernelVx;

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVx(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVx(i, j) = cos_Fai * cos(Alpha) * (
				(0.5 * (cscPhi) * (1 - 3 * cos2Phi) * (l1pow1 - l2pow1) +
					0.5 * (-1.0 / tan(Phi) + 3 * cscPhi * cos(3 * Phi)) * (R2 * l1pow1 - R1 * l2pow1) +
					0.5 * (cscPhi) * (1 - cos2Phi) * (R2 * R2 * l1pow1 - R1 * R1 * l2pow1)) / (l2pow1 * l1pow1) -
				1.5 * sin(2 * Phi) * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));
		}
	}

	return KernelVx;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha;
	double sin_half_Phi_sq, l2pow1, l1pow1;
	double cscPhi, cosPhi, cos2Phi;

	h2 = r - R2;
	hRatio2 = h2 / r;
	h1 = r - R1;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	Eigen::Matrix<double, 14, 14> KernelVy;

	for (int i = 0; i < 14; ++i) {
		for (int j = 0; j < 14; ++j) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVy(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			sin_half_Phi_sq = pow(sin(Phi / 2), 2);

			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVy(i, j) = cos_Fai * sin(Alpha) * (
				(0.5 * (cscPhi) * (1 - 3 * cos2Phi) * (l1pow1 - l2pow1) +
					0.5 * (-1.0 / tan(Phi) + 3 * cscPhi * cos(3 * Phi)) * (R2 * l1pow1 - R1 * l2pow1) +
					0.5 * (cscPhi) * (1 - cos2Phi) * (R2 * R2 * l1pow1 - R1 * R1 * l2pow1)) / (l2pow1 * l1pow1) -
				1.5 * sin(2 * Phi) * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));
		}
	}

	return KernelVy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVz(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi;
	double sin_half_Phi_sq, l2pow1, l1pow1;
	double cosPhi;

	Eigen::Matrix<double, 14, 14> KernelVz;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVz(i, j) = -cos_Fai *
						((R2 - R1) / ((1 - R2) * (1 - R1)) + (R2 - R1) + 2 * log((1 - R2) / (1 - R1)));
					continue;
				}
				else if (R1 < 1) {
					KernelVz(i, j) = cos_Fai *
						((2 - R2 - R1) / ((1 - R2) * (1 - R1)) + (R2 + R1) + 4 * log(r) + 2 * log((1 - R1) * (R2 - 1)));
					continue;
				}
				else {
					KernelVz(i, j) = cos_Fai *
						((R2 - R1) / ((1 - R2) * (1 - R1)) + (R2 - R1) + 2 * log((1 - R2) / (1 - R1)));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVz(i, j) = cos_Fai *
					((R1 - R2) / ((1 + R2) * (1 + R1)) - (R2 - R1) + 2 * log((1 + R2) / (1 + R1)));
				continue;
			}

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);

			cosPhi = cos(Phi);

			KernelVz(i, j) = cos_Fai * (
				(3 * cosPhi * (l1pow1 - l2pow1) +
					(1 - 6 * cosPhi * cosPhi) * (R2 * l1pow1 - R1 * l2pow1) +
					cosPhi * (R2 * R2 * l1pow1 - R1 * R1 * l2pow1)) / (l2pow1 * l1pow1) +
				(1 - 3 * cosPhi * cosPhi) * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));

			continue;

		}
	}

	return KernelVz;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVxx(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3;
	double cscPhi, cosPhi, cos2Phi;

	Eigen::Matrix<double, 14, 14> KernelVxx;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVxx(i, j) = cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) -
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						log((1 - R2) / (1 - R1)));
					continue;
				}
				else if (R1 < 1) {
					KernelVxx(i, j) = -cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) +
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						2 * log(r) +
						log((R2 - 1) / (1 - R1)));
					continue;
				}
				else {
					KernelVxx(i, j) = -cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) -
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						log((R2 - 1) / (R1 - 1)));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVxx(i, j) = -cos_Fai * (
					((3 + 4 * R2) / (2 * pow(1 + R2, 2))) -
					((3 + 4 * R1) / (2 * pow(1 + R1, 2))) +
					log((1 + R2) / (1 + R1)));
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVxx(i, j) = cos_Fai * (
				pow(cos(Alpha), 2) * pow(cscPhi, 2) * (
					((-5 * cosPhi + 3 * pow(cosPhi, 3.0)) * (l1pow3 - l2pow3) +
						(-3 + 15 * pow(cosPhi, 2.0) - 6 * pow(cosPhi, 2.0) * cos2Phi) * (R2 * l1pow3 - R1 * l2pow3) +
						(-9 * pow(cosPhi, 3.0) + 3 * pow(cosPhi, 2.0) * cos(3 * Phi)) * (R2 * R2 * l1pow3 - R1 * R1 * l2pow3) +
						(-4 + 10 * pow(cosPhi, 2.0) - 4 * pow(cosPhi, 2.0) * cos2Phi) * (R2 * R2 * R2 * l1pow3 - R1 * R1 * R1 * l2pow3)) / (l2pow3 * l1pow3)) +
				((1.0 / tan(Phi)) * (cscPhi) * (l1pow1 - l2pow1) +
					(1 - pow((1.0 / tan(Phi)), 2.0)) * (R2 * l1pow1 - R1 * l2pow1)) / (l2pow1 * l1pow1) +
				(1 - 3 * Tx * Tx) * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));

			continue;
		}
	}

	return KernelVxx;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVxy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3;
	double cscPhi, cosPhi, cos2Phi;

	Eigen::Matrix<double, 14, 14> KernelVxy;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVxy(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVxy(i, j) = cos_Fai * 0.5 * sin(2 * Alpha) * (
				cscPhi * cscPhi * (
					(cosPhi * (-5 + 3 * pow(cosPhi, 2)) * (l1pow3 - l2pow3) +
						3 * (-1 + 7 * pow(cosPhi, 2) - 4 * pow(cosPhi, 4)) * (R2 * l1pow3 - R1 * l2pow3) +
						6 * pow(cosPhi, 3) * (-2 + cos2Phi) * (R2 * R2 * l1pow3 - R1 * R1 * l2pow3) +
						(3 * cos2Phi - cos(4 * Phi)) * (R2 * R2 * R2 * l1pow3 - R1 * R1 * R1 * l2pow3)) / (l2pow3 * l1pow3))) -
				cos_Fai * Tx * Ty * 3 * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1));

		}
	}

	return KernelVxy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVyy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3;
	double cscPhi, cosPhi, cos2Phi;

	Eigen::Matrix<double, 14, 14> KernelVyy;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVyy(i, j) = cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) -
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						log((1 - R2) / (1 - R1)));
					continue;
				}
				else if (R1 < 1) {
					KernelVyy(i, j) = -cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) +
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						2 * log(r) +
						log((R2 - 1) / (1 - R1)));
					continue;
				}
				else {
					KernelVyy(i, j) = -cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) -
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						log((R2 - 1) / (R1 - 1)));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVyy(i, j) = -cos_Fai * (
					((3 + 4 * R2) / (2 * pow(1 + R2, 2))) -
					((3 + 4 * R1) / (2 * pow(1 + R1, 2))) +
					log((1 + R2) / (1 + R1)));
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVyy(i, j) = cos_Fai * (
				pow(sin(Alpha), 2) * pow(cscPhi, 2) * (
					((-5 * cosPhi + 3 * pow(cosPhi, 3.0)) * (l1pow3 - l2pow3) +
						(-3 + 15 * pow(cosPhi, 2.0) - 6 * pow(cosPhi, 2.0) * cos2Phi) * (R2 * l1pow3 - R1 * l2pow3) +
						(-9 * pow(cosPhi, 3.0) + 3 * pow(cosPhi, 2.0) * cos(3 * Phi)) * (R2 * R2 * l1pow3 - R1 * R1 * l2pow3) +
						(-4 + 10 * pow(cosPhi, 2.0) - 4 * pow(cosPhi, 2.0) * cos2Phi) * (R2 * R2 * R2 * l1pow3 - R1 * R1 * R1 * l2pow3)) / (l2pow3 * l1pow3)) +
				((1.0 / tan(Phi)) * (cscPhi) * (l1pow1 - l2pow1) +
					(1 - pow((1.0 / tan(Phi)), 2.0)) * (R2 * l1pow1 - R1 * l2pow1)) / (l2pow1 * l1pow1) +
				(1 - 3 * Ty * Ty) * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));

			continue;
		}
	}

	return KernelVyy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVzx(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3;
	double cscPhi, cosPhi, cos2Phi;

	Eigen::Matrix<double, 14, 14> KernelVzx;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVzx(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;
			cscPhi = 1.0 / sin(Phi);

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVzx(i, j) = cos_Fai * (
				0.5 * cos(Alpha) * (
					(cscPhi * (1 - 3 * cos2Phi) * (l1pow3 - l2pow3) +
						cscPhi * 6 * cos(3 * Phi) * (R2 * l1pow3 - R1 * l2pow3) +
						cscPhi * 3 * (1 - 2 * cos2Phi - cos(4 * Phi)) * (R2 * R2 * l1pow3 - R1 * R1 * l2pow3) +
						(2 * (2 * cos(3 * Phi) * cscPhi - (1.0 / tan(Phi)))) * (R2 * R2 * R2 * l1pow3 - R1 * R1 * R1 * l2pow3)) / (l2pow3 * l1pow3))
				- 3 * cosPhi * Tx * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));

			continue;
		}
	}

	return KernelVzx;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVzy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3;
	double cscPhi, cosPhi, cos2Phi;

	Eigen::Matrix<double, 14, 14> KernelVzy;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVzy(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVzy(i, j) = cos_Fai * (
				0.5 * sin(Alpha) * (
					(cscPhi * (1 - 3 * cos2Phi) * (l1pow3 - l2pow3) +
						cscPhi * 6 * cos(3 * Phi) * (R2 * l1pow3 - R1 * l2pow3) +
						cscPhi * 3 * (1 - 2 * cos2Phi - cos(4 * Phi)) * (R2 * R2 * l1pow3 - R1 * R1 * l2pow3) +
						(2 * (2 * cos(3 * Phi) * cscPhi - (1.0 / tan(Phi)))) * (R2 * R2 * R2 * l1pow3 - R1 * R1 * R1 * l2pow3)) / (l2pow3 * l1pow3))
				- 3 * cosPhi * Ty * log((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));

			continue;
		}
	}

	return KernelVzy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVzz(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3;
	double cscPhi, cosPhi, cos2Phi;

	Eigen::Matrix<double, 14, 14> KernelVzz;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVzz(i, j) = -2 * cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) -
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						log((1 - R2) / (1 - R1)));
					continue;
				}
				else if (R1 < 1) {
					KernelVzz(i, j) = 2 * cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) +
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						2 * log(r) +
						log((R2 - 1) * (1 - R1)));
					continue;
				}
				else {
					KernelVzz(i, j) = 2 * cos_Fai * (
						((3 - 4 * R2) / (2 * pow(1 - R2, 2))) -
						((3 - 4 * R1) / (2 * pow(1 - R1, 2))) +
						log((R2 - 1) / (R1 - 1)));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVzz(i, j) = 2 * cos_Fai * (
					((3 + 4 * R2) / (2 * pow(1 + R2, 2))) -
					((3 + 4 * R1) / (2 * pow(1 + R1, 2))) +
					log((1 + R2) / (1 + R1)));
				continue;
			}

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVzz(i, j) = cos_Fai * (
				(3 * cosPhi * (l1pow3 - l2pow3) +
					(-5 - 6 * cos2Phi) * (R2 * l1pow3 - R1 * l2pow3) +
					2 * cosPhi * (4 + 3 * cos2Phi) * (R2 * R2 * l1pow3 - R1 * R1 * l2pow3) +
					2 * (-1 - 2 * cos2Phi) * (R2 * R2 * R2 * l1pow3 - R1 * R1 * R1 * l2pow3)) / (l2pow3 * l1pow3) +
				(1 - 3 * cosPhi * cosPhi) * logl((cosPhi - R2 + l2pow1) / (cosPhi - R1 + l1pow1)));

			continue;
		}
	}

	return KernelVzz;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVxxx(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cscPhi, cosPhi, cos2Phi, cosAlphapow2;

	Eigen::Matrix<double, 14, 14> KernelVxxx;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVxxx(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);
			cosAlphapow2 = pow(cos(Alpha), 2);

			KernelVxxx(i, j) = -cos_Fai * cos(Alpha) * pow(cscPhi, 3) * (
				(cosAlphapow2 * 8 * (
					(1 - R2 * cosPhi) * l1pow5 - (1 - R1 * cosPhi) * l2pow5) -
					cosAlphapow2 * 32 * cosPhi * (
						R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * (1 - R1 * cosPhi) * l2pow5) +
					cosAlphapow2 * 4 * (5 + 7 * pow(cosPhi, 2)) * (
						pow(R2, 2) * (1 - R2 * cosPhi) * l1pow5 - pow(R1, 2) * (1 - R1 * cosPhi) * l2pow5) +
					cosAlphapow2 * 4 * cosPhi * (-9 + cos2Phi) * (
						pow(R2, 3) * (1 - R2 * cosPhi) * l1pow5 - pow(R1, 3) * (1 - R1 * cosPhi) * l2pow5) +
					cosAlphapow2 * (15 - 10 * pow(cosPhi, 2) + 3 * pow(cosPhi, 4)) * (
						pow(R2, 4) * (1 - R2 * cosPhi) * l1pow5 - pow(R1, 4) * (1 - R1 * cosPhi) * l2pow5)) / (l2pow5 * l1pow5) +
				(-6 * (
					(1 - R2 * cosPhi) * l1pow3 - (1 - R1 * cosPhi) * l2pow3) +
					12 * cosPhi * (
						R2 * (1 - R2 * cosPhi) * l1pow3 - R1 * (1 - R1 * cosPhi) * l2pow3) +
					1.5 * (-5 + cos2Phi) * (
						pow(R2, 2) * (1 - R2 * cosPhi) * l1pow3 - pow(R1, 2) * (1 - R1 * cosPhi) * l2pow3)) / (l2pow3 * l1pow3));

		}
	}

	return KernelVxxx;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVxxy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cscPhi, cosPhi, cos2Phi, cosAlphapow2;

	Eigen::Matrix<double, 14, 14> KernelVxxy;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVxxy(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);
			cosAlphapow2 = pow(cos(Alpha), 2);

			KernelVxxy(i, j) = -cos_Fai * sin(Alpha) * pow(cscPhi, 3) * (
				cosAlphapow2 * (8 *
					((1 - R2 * cosPhi) * l1pow5 - (1 - R1 * cosPhi) * l2pow5) -
					32 * cosPhi * (
						R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * (1 - R1 * cosPhi) * l2pow5) +
					4 * (5 + 7 * pow(cosPhi, 2)) * (
						pow(R2, 2) * (1 - R2 * cosPhi) * l1pow5 - pow(R1, 2) * (1 - R1 * cosPhi) * l2pow5) +
					4 * cosPhi * (-9 + cos2Phi) * (
						pow(R2, 3) * (1 - R2 * cosPhi) * l1pow5 - pow(R1, 3) * (1 - R1 * cosPhi) * l2pow5) +
					(15 - 10 * pow(cosPhi, 2) + 3 * pow(cosPhi, 4)) * (
						pow(R2, 4) * (1 - R2 * cosPhi) * l1pow5 - pow(R1, 4) * (1 - R1 * cosPhi) * l2pow5)) / (l2pow5 * l1pow5) +
				(-2 * ((1 - R2 * cosPhi) * l1pow3 - (1 - R1 * cosPhi) * l2pow3) +
					4 * cosPhi * (R2 * (1 - R2 * cosPhi) * l1pow3 - R1 * (1 - R1 * cosPhi) * l2pow3) +
					(-3 + pow(cosPhi, 2)) * (pow(R2, 2) * (1 - R2 * cosPhi) * l1pow3 - pow(R1, 2) * (1 - R1 * cosPhi) * l2pow3)) / (l2pow3 * l1pow3));

		}
	}

	return KernelVxxy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVxxz(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cscPhi, cosPhi, cos2Phi, cosAlphapow2;

	Eigen::Matrix<double, 14, 14> KernelVxxz;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVxxz(i, j) = cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) -
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
				else if (R1 < 1) {
					KernelVxxz(i, j) = -cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) +
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
				else {
					KernelVxxz(i, j) = -cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) -
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVxxz(i, j) = -cos_Fai * (
					((1 + 3 * R2 + 3 * R2 * R2) / (pow(1 + R2, 3))) -
					((1 + 3 * R1 + 3 * R1 * R1) / (pow(1 + R1, 3))));
				continue;
			}

			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);

			KernelVxxz(i, j) = cos_Fai * (
				(1 * (pow(R2, 3) * l1pow5 - pow(R1, 3) * l2pow5) -
					2 * cosPhi * (pow(R2, 4) * l1pow5 - pow(R1, 4) * l2pow5) +
					(1 - 3 * pow(Tx, 2)) * (pow(R2, 5) * l1pow5 - pow(R1, 5) * l2pow5)) / (l2pow5 * l1pow5));

			continue;
		}
	}

	return KernelVxxz;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVxyz(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {

	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;

	Eigen::Matrix<double, 14, 14> KernelVxyz;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVxyz(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			KernelVxyz(i, j) = -3 * cos_Fai * Tx * Ty * (pow(R2, 5) * l1pow5 - pow(R1, 5) * l2pow5) / (l2pow5 * l1pow5);

			continue;
		}
	}

	return KernelVxyz;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVyyx(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cscPhi, cosPhi, cos2Phi, sinAlphapow2;

	Eigen::Matrix<double, 14, 14> KernelVyyx;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVyyx(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);
			sinAlphapow2 = pow(sin(Alpha), 2);

			KernelVyyx(i, j) = -cos_Fai * cos(Alpha) * pow(cscPhi, 3) * (
				sinAlphapow2 * (8 *
					((1 - R2 * cosPhi) * l1pow5 - (1 - R1 * cosPhi) * l2pow5) -
					32 * cosPhi * (
						R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * (1 - R1 * cosPhi) * l2pow5) +
					4 * (5 + 7 * cosPhi * cosPhi) * (
						R2 * R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * R1 * (1 - R1 * cosPhi) * l2pow5) +
					4 * cosPhi * (-9 + cos2Phi) * (
						R2 * R2 * R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * R1 * R1 * (1 - R1 * cosPhi) * l2pow5) +
					(15 - 10 * cosPhi * cosPhi + 3 * cosPhi * cosPhi * cosPhi * cosPhi) * (
						R2 * R2 * R2 * R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * R1 * R1 * R1 * (1 - R1 * cosPhi) * l2pow5)) / (l2pow5 * l1pow5) +
				(-2 * ((1 - R2 * cosPhi) * l1pow3 - (1 - R1 * cosPhi) * l2pow3) +
					4 * cosPhi * (R2 * (1 - R2 * cosPhi) * l1pow3 - R1 * (1 - R1 * cosPhi) * l2pow3) +
					(-3 + cosPhi * cosPhi) * (R2 * R2 * (1 - R2 * cosPhi) * l1pow3 - R1 * R1 * (1 - R1 * cosPhi) * l2pow3)) / (l2pow3 * l1pow3));

		}
	}

	return KernelVyyx;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVyyy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cscPhi, cosPhi, cos2Phi, sinAlphapow2;

	Eigen::Matrix<double, 14, 14> KernelVyyy;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVyyy(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow3 = l2pow1 * l2pow1 * l2pow1;
			l1pow3 = l1pow1 * l1pow1 * l1pow1;
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cscPhi = 1 / sin(Phi);
			cosPhi = cos(Phi);
			cos2Phi = cos(2 * Phi);
			sinAlphapow2 = pow(sin(Alpha), 2);

			KernelVyyy(i, j) = -cos_Fai * sin(Alpha) * pow(cscPhi, 3) * (
				sinAlphapow2 * (8 *
					((1 - R2 * cosPhi) * l1pow5 - (1 - R1 * cosPhi) * l2pow5) -
					32 * cosPhi * (
						R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * (1 - R1 * cosPhi) * l2pow5) +
					4 * (5 + 7 * cosPhi * cosPhi) * (
						R2 * R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * R1 * (1 - R1 * cosPhi) * l2pow5) +
					4 * cosPhi * (-9 + cos2Phi) * (
						R2 * R2 * R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * R1 * R1 * (1 - R1 * cosPhi) * l2pow5) +
					(15 - 10 * cosPhi * cosPhi + 3 * pow(cosPhi, 4)) * (
						R2 * R2 * R2 * R2 * (1 - R2 * cosPhi) * l1pow5 - R1 * R1 * R1 * R1 * (1 - R1 * cosPhi) * l2pow5)) / (l2pow5 * l1pow5) +
				(-6 * ((1 - R2 * cosPhi) * l1pow3 - (1 - R1 * cosPhi) * l2pow3) +
					12 * cosPhi * (R2 * (1 - R2 * cosPhi) * l1pow3 - R1 * (1 - R1 * cosPhi) * l2pow3) +
					1.5 * (-5 + cos2Phi) * (R2 * R2 * (1 - R2 * cosPhi) * l1pow3 - R1 * R1 * (1 - R1 * cosPhi) * l2pow3)) / (l2pow3 * l1pow3));
		}
	}

	return KernelVyyy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVyyz(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cosPhi;

	Eigen::Matrix<double, 14, 14> KernelVyyz;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVyyz(i, j) = cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) -
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
				else if (R1 < 1) {
					KernelVyyz(i, j) = -cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) +
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
				else {
					KernelVyyz(i, j) = -cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) -
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVyyz(i, j) = -cos_Fai * (
					((1 + 3 * R2 + 3 * R2 * R2) / (pow(1 + R2, 3))) -
					((1 + 3 * R1 + 3 * R1 * R1) / (pow(1 + R1, 3))));
				continue;
			}

			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cosPhi = cos(Phi);

			KernelVyyz(i, j) = cos_Fai * (
				(1.0 * (R2 * R2 * R2 * l1pow5 - R1 * R1 * R1 * l2pow5) -
					2.0 * cosPhi * (R2 * R2 * R2 * R2 * l1pow5 - R1 * R1 * R1 * R1 * l2pow5) +
					(1.0 - 3.0 * Ty * Ty) * (R2 * R2 * R2 * R2 * R2 * l1pow5 - R1 * R1 * R1 * R1 * R1 * l2pow5)) / (l2pow5 * l1pow5));

			continue;
		}
	}

	return KernelVyyz;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVzzx(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cosPhi;

	Eigen::Matrix<double, 14, 14> KernelVzzx;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVzzx(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cosPhi = cos(Phi);

			KernelVzzx(i, j) = 3.0 * cos_Fai * Tx * (
				(R2 * R2 * R2 * R2 * (1.0 - R2 * cosPhi) * l1pow5 - R1 * R1 * R1 * R1 * (1.0 - R1 * cosPhi) * l2pow5) /
				(l2pow5 * l1pow5));

			continue;
		}
	}

	return KernelVzzx;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVzzy(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cosPhi;

	Eigen::Matrix<double, 14, 14> KernelVzzy;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps) || Phi >(M_PI / 2 + acos_eps)) {
				KernelVzzy(i, j) = 0;
				continue;
			}

			Alpha = atan2(sin_Lamda_diff * cos_Fai, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;
			Ty = cos_Fai * sin_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cosPhi = cos(Phi);

			KernelVzzy(i, j) = 3.0 * cos_Fai * Ty * (
				(R2 * R2 * R2 * R2 * (1.0 - R2 * cosPhi) * l1pow5 - R1 * R1 * R1 * R1 * (1.0 - R1 * cosPhi) * l2pow5) /
				(l2pow5 * l1pow5));

			continue;
		}
	}

	return KernelVzzy;
}

inline Eigen::Matrix<double, 14, 14> Tesseroid_IntegralkernelVzzz(double R2, double R1,
	const Eigen::Matrix<double, 14, 14>& FaiI, const Eigen::Matrix<double, 14, 14>& LamdaI,
	double r, double fai, double lamda) {


	const double rad = M_PI / 180.0;
	const double acos_eps = acos(1e-5);

	double h2, h1, hRatio2, hRatio1;
	double Fai_rad, Lamda_diff_rad, cos_Lamda_diff, sin_Lamda_diff;
	double num, den, Phi, Alpha, Tx, Ty;
	double sin_half_Phi_sq, l2pow1, l1pow1, l2pow3, l1pow3, l2pow5, l1pow5;
	double cosPhi;

	Eigen::Matrix<double, 14, 14> KernelVzzz;

	h2 = r - R2;
	h1 = r - R1;
	hRatio2 = h2 / r;
	hRatio1 = h1 / r;

	R2 = 1 - hRatio2;
	R1 = 1 - hRatio1;
    double cos_fai = cos(fai * rad);
    double sin_fai = sin(fai * rad);

	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 14; j++) {

			Fai_rad = FaiI(i, j) * rad;
			Lamda_diff_rad = (LamdaI(i, j) - lamda) * rad;
			double cos_Fai = cos(Fai_rad);
			double sin_Fai = sin(Fai_rad);
			cos_Lamda_diff = cos(Lamda_diff_rad);
			sin_Lamda_diff = sin(Lamda_diff_rad);


			num = hypot(cos_Fai * sin_Lamda_diff, cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff);
			den = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_Lamda_diff;
			Phi = atan2(num, den);

			if (Phi < (M_PI / 2 - acos_eps)) {
				if (R2 < 1) {
					KernelVzzz(i, j) = -2 * cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) -
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
				else if (R1 < 1) {
					KernelVzzz(i, j) = 2 * cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) +
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
				else {
					KernelVzzz(i, j) = 2 * cos_Fai * (
						((1 - 3 * R2 + 3 * R2 * R2) / (pow(1 - R2, 3))) -
						((1 - 3 * R1 + 3 * R1 * R1) / (pow(1 - R1, 3))));
					continue;
				}
			}

			if (Phi > (M_PI / 2 + acos_eps)) {
				KernelVzzz(i, j) = 2 * cos_Fai * (
					((1 + 3 * R2 + 3 * R2 * R2) / (pow(1 + R2, 3))) -
					((1 + 3 * R1 + 3 * R1 * R1) / (pow(1 + R1, 3))));
				continue;
			}

			Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_Lamda_diff;

			sin_half_Phi_sq = pow(sin(Phi / 2), 2);
			l2pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio2) + hRatio2 * hRatio2);
			l1pow1 = sqrt(2 * (2 * sin_half_Phi_sq) * (1 - hRatio1) + hRatio1 * hRatio1);
			l2pow5 = l2pow1 * l2pow1 * l2pow1 * l2pow1 * l2pow1;
			l1pow5 = l1pow1 * l1pow1 * l1pow1 * l1pow1 * l1pow1;

			cosPhi = cos(Phi);

			KernelVzzz(i, j) = cos_Fai * (
				(-2.0 * (R2 * R2 * R2 * l1pow5 - R1 * R1 * R1 * l2pow5) +
					4.0 * cosPhi * (R2 * R2 * R2 * R2 * l1pow5 - R1 * R1 * R1 * R1 * l2pow5) +
					(1.0 - 3.0 * cosPhi * cosPhi) * (R2 * R2 * R2 * R2 * R2 * l1pow5 - R1 * R1 * R1 * R1 * R1 * l2pow5)) / (l2pow5 * l1pow5));

			continue;
		}
	}

	return KernelVzzz;
}

}