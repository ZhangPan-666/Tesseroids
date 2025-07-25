#pragma once

#include <iostream>
#include <array>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <queue>
#include <Eigen/Dense>

namespace mista_math
{


// Definition of Gauss-Kronrod (3,7)
struct Gauss3Kronrod7Const
{
    const std::array<double, 7> Nodes = {
        -0.9604912687080202, -0.7745966692414834, -0.4342437493468026,
        0.0,
        0.4342437493468026, 0.7745966692414834, 0.9604912687080202
    };
    const Eigen::Matrix<double, 1, 7> LowWeights = {0, 5.0 / 9.0, 0, 8.0 / 9.0, 0, 5.0 / 9.0, 0};
    const Eigen::Matrix<double, 1, 7> HighWeights = {
        0.1046562260264672, 0.2684880898683334, 0.4013974147759622,
        0.4509165386584744,
        0.4013974147759622, 0.2684880898683334, 0.1046562260264672
    };

    const Eigen::Array<double, 14, 1> NArray = {
        0.009877182822994962, 0.05635083268962915, 0.14143906266329936, 0.25, 0.35856093733670064, 0.44364916731037085, 0.490122817177005, 
        0.509877182822995, 0.5563508326896291, 0.6414390626632993, 0.75, 0.8585609373367007, 0.9436491673103709, 0.990122817177005
    };
    const Eigen::Array<double, 6, 1> VtstIdx ={25,  49,  73,  97,  121,  145};

};



class SubRectangle {
public:
    double q_;           // Estimated integral value
    double e_;           // Estimation error
    double L_;           // Left boundary
    double R_;           // Right boundary
    double B_;           // Bottom boundary
    double T_;           // Top boundary
    double adjer_;       // Adjusted error

    // Constructor
    SubRectangle(double q, double e, double L, double R, double B, double T, double adjer)
        : q_(q), e_(e), L_(L), R_(R), B_(B), T_(T), adjer_(adjer) {}

    bool operator<(const SubRectangle& other) const {
        return adjer_ < other.adjer_;
    }
};


class RectangleInfoManager {
private:
    std::priority_queue<SubRectangle, std::vector<SubRectangle>> pq;
    double err_ok = 0.0;
    double errbnd = 0.0;
    double sum_adjer_in_pq = 0.0; 

public:
    void saveRectInfo(const std::array<double, 4>& Qsub, const std::array<double, 4>& esub,
                      double thetaL, double thetaR, double phiB, double phiT,
                      double TOL, double AREA, double EPS100, double ADJUST) {
        double dthetad2 = (thetaR - thetaL) / 2;
        double thetaM = thetaL + dthetad2;
        double dphid2 = (phiT - phiB) / 2;
        double phiM = phiB + dphid2;
        double localtol = TOL * dthetad2 * dphid2 / AREA;
        localtol = std::max(std::abs(localtol), EPS100 * std::abs(Qsub[0]+Qsub[1]+Qsub[2]+Qsub[3]));

        // Define the boundary information of the subrectangle
        struct Boundary {
            double L, R, B, T;
        };
        std::array<Boundary, 4> boundaries = {
            thetaL, thetaM, phiB, phiM,
            thetaM, thetaR, phiB, phiM,
            thetaL, thetaM, phiM, phiT,
            thetaM, thetaR, phiM, phiT
        };

        // Handle each sub-rectangle
        for (size_t i = 0; i < 4; ++i) {
            double adjer_i = ADJUST * esub[i];
            if (adjer_i > localtol) {
                pq.push(SubRectangle(Qsub[i], esub[i], boundaries[i].L, boundaries[i].R,
                                     boundaries[i].B, boundaries[i].T, adjer_i));
                sum_adjer_in_pq += adjer_i;
            } else {
                err_ok += adjer_i;
            }
        }

        errbnd = err_ok + sum_adjer_in_pq;
    }

    // Get the sub-rectangle with the largest error
    SubRectangle getNextRectangle() {
        if (pq.empty()) {
            return SubRectangle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        SubRectangle rect = pq.top();
        pq.pop();
        sum_adjer_in_pq -= rect.adjer_;
        errbnd = err_ok + sum_adjer_in_pq;
        return rect;
    }

    bool isEmpty() const {
        return pq.empty();
    }

    double getErrOk() const {
        return err_ok;
    }

    double getErrBnd() const {
        return errbnd;
    }
};


/**
 * @brief Integral partitions are computed for a given rectangular region and 
 * the Gauss - Kronrod method is used to estimate the integral value and error.
 * 
 * This function divides the specified rectangular region [thetaL, thetaR] x [phiB, phiT] into subregions,
 * and calculate the integral approximation and error estimates for each subregion, 
 * it will also be processed accordingly to the first evaluation flags and boundary conditions.
 * 
 * @tparam N The number of integral nodes, the default is 14.
 * @tparam The type of the constant constructor for the KernalT integration rule, which defaults to Gauss3Kronrod7Const.
 * @param thetaL The rectangular region is at the left boundary of the theta dimension.
 * @param thetaR The rectangular region is at the right boundary of the theta dimension.
 * @param phiB The rectangular region is at the lower boundary of the phi dimension.
 * @param phiT The rectangular region is at the upper boundary of the phi dimension.
 * @param Qsub Output parameter storing the integral approximation for the four subregions.
 * @param esub Output parameter that stores the estimation errors for the four subregions.
 * @param FIRSTFUNEVAL Reference parameter indicating whether this is the first time the function is evaluated.
 * @param NFE Reference parameter to record the total number of times the function was evaluated.
 * @param XMAX The maximum value of the integration area in the x-direction.
 * @param XMIN The minimum value of the integration area in the x-direction.
 * @param YMAX The maximum value of the integration area in the y-direction.
 * @param YMIN The minimum value of the integration area in the y-direction.
 * @param FUN A function handle to be integrated that accepts as input two NxN matrices and returns an NxN matrix.
 * @param ATOL Absolute error tolerance for controlling the accuracy of the integral.
 * @param RTOL Relative error tolerance, used to control the precision of the integral.
 */
template <size_t N = 14, class KernalT = Gauss3Kronrod7Const>
inline void PartitionIntegralQuad(
    double thetaL, double thetaR, double phiB, double phiT,
    Eigen::Matrix<double, 1, 4>& Qsub, Eigen::Matrix<double, 1, 4>& esub,
    bool& FIRSTFUNEVAL, int& NFE,
    double XMAX, double XMIN, double YMAX, double YMIN,
    const std::function<Eigen::Matrix<double, N, N>(const Eigen::Matrix<double, N, N>&, const Eigen::Matrix<double, N, N>&)>& FUN
) {

    const KernalT g3k7;
    const Eigen::Matrix<double, 1, N/2>& WT3 = g3k7.LowWeights;
    const Eigen::Matrix<double, 1, N/2>& WT7 = g3k7.HighWeights;
    const Eigen::Array<double, N, 1>& NARRAY = g3k7.NArray;

    double dtheta = thetaR - thetaL;
    Eigen::Array<double, N, 1> theta = thetaL + NARRAY * dtheta;
    Eigen::Array<double, N, 1> x = 0.5 * (XMAX + XMIN) + 0.5 * (XMAX - XMIN) * theta.cos();

    if (!FIRSTFUNEVAL && (x[0] == XMAX || x[N - 1] == XMIN)) {
        return;
    }

    Eigen::Matrix<double, N, N> X = x.replicate(1, N).transpose();

    double dphi = phiT - phiB;
    Eigen::Array<double, N, 1> phi = phiB + NARRAY * dphi;

    Eigen::Array<double, N, N> Y;
    double dydt = (YMAX - YMIN);

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			Y(i, j) = YMIN + (0.5 + 0.5 * cos(phi(i))) * dydt;
		}
	}

    bool returnFlag = false;
    for (size_t i = 0; i < N; ++i) {
        if ((Y(0, i) == YMAX) || (Y(N - 1, i) == YMIN)) {
            returnFlag = true;
            break;
        }
    }

    if (!FIRSTFUNEVAL && returnFlag) {
        return;
    }

    Eigen::Matrix<double, N, N> Z = FUN(X, Y);
    ++NFE;
   
    Eigen::Matrix<double, N, N> temp;
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			Z(i, j) *= 0.25 * (XMAX - XMIN) * sin(phi[i]) * (dydt * sin(theta[j]));
		}
	}

    Eigen::Matrix<double, N/2, 2*N> Z_part;
    for (size_t i = 0; i < N/2; i++) {
        for (size_t j = 0; j < N; j++)
        {
            Z_part(i, j) = Z(i, j);
        }
    }
    for (size_t i = N/2; i < N; i++) {
        for (size_t j = 0; j < N; j++)
        {
            Z_part(i-N/2, j+N) = Z(i, j);
        }
    }

    double r = (dtheta / 4) * (dphi / 4);
    Eigen::Matrix<double, 1, 2*N> tmp = WT7 * Z_part;
    Eigen::Matrix<double, N/2, 4> tmp_reshaped = tmp.reshaped(N/2, 4);
    Qsub = WT7 * tmp_reshaped * r;

    tmp = WT3 * Z_part;
    tmp_reshaped = tmp.reshaped(N/2, 4);
    esub = (WT3  * (tmp_reshaped) * r - Qsub);
    esub = ((esub.array()).abs()).matrix();

}


// function integral2t
inline double integral2t(
    const std::function<Eigen::Matrix<double, 14, 14>(
        const Eigen::Matrix<double, 14, 14>&,
        const Eigen::Matrix<double, 14, 14>&)> FUN,
    double XMIN, double XMAX, double YMIN, double YMAX,
    double ATOL = 1e-10, double RTOL = 1e-6, int maxFunEvals = 1000
) {
    bool FIRSTFUNEVAL = true;
    int NFE = 0;

    double thetaL = 0;
    double thetaR = M_PI;
    double phiB = 0;
    double phiT = M_PI;

    double AREA = (thetaR - thetaL) * (phiT - phiB);
    Eigen::Matrix<double, 1, 4> Qsub;
    Eigen::Matrix<double, 1, 4> esub;

    // Calculate the initial approximation
    PartitionIntegralQuad<14, Gauss3Kronrod7Const>(thetaL, thetaR, phiB, phiT, Qsub, esub, FIRSTFUNEVAL, NFE, XMAX, XMIN, YMAX, YMIN, FUN);
    double Q = Qsub.sum();
    double EPS100 = 100 * std::numeric_limits<double>::epsilon();
    if (RTOL < EPS100) {
        RTOL = EPS100;
    }
    double rtold8 = std::max(RTOL / 8, EPS100);
    double atold8 = ATOL / 8;
    double TOL = EPS100 * std::abs(Q);

    RectangleInfoManager manager;
    double ADJUST = 1.0;
    manager.saveRectInfo({Qsub(0), Qsub(1), Qsub(2), Qsub(3)}, {esub(0), esub(1), esub(2), esub(3)}, thetaL, thetaR, phiB, phiT, TOL, AREA, EPS100, ADJUST);

    while (!manager.isEmpty() && manager.getErrBnd() > TOL) {
        if (NFE >= maxFunEvals) {
            break;
        }
        SubRectangle rect = manager.getNextRectangle();
        PartitionIntegralQuad<14, Gauss3Kronrod7Const>(rect.L_, rect.R_, rect.B_, rect.T_, Qsub, esub, FIRSTFUNEVAL, NFE, XMAX, XMIN, YMAX, YMIN, FUN);
        double Newq = Qsub.sum();
        ADJUST = std::min(1.0, std::abs(rect.q_ - Newq) / rect.e_);
        Q = Q + (Newq - rect.q_);
        TOL = std::max(atold8, rtold8 * std::abs(Q));
        manager.saveRectInfo({Qsub(0), Qsub(1), Qsub(2), Qsub(3)}, {esub(0), esub(1), esub(2), esub(3)}, rect.L_, rect.R_, rect.B_, rect.T_, TOL, AREA, EPS100, ADJUST);
    }

    return Q;
}


} // namespace mista_math
