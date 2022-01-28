// ------------------------------------------------------------------
// Copyright (C) 2018 Lawrence Buckingham.
//
// This file is part of the supplementary material which accompanies 
// the paper:
//	Lawrence Buckingham, Shlomo Geva, and James M. Hogan. 2018.
//	Protein database search using compressed k-mer vocabularies.
//	In 23rd Australasian Document Computing Symposium (ADCS '18), 
//	December 11--12, 2018, Dunedin, New Zealand
//	https://doi.org/10.1145/3291992.3291997
//
// This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3, or (at your 
// option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License 
// along with this program; see the file "license.txt".  If not, see
// <http://www.gnu.org/licenses/>.
// ------------------------------------------------------------------

#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201103L
#endif

#include <cmath>

#include "Distribution.hpp"
#include "Exception.hpp"
#include "Constants.h"
#include "JohnCook.h"

namespace QutBio{
    class NormalDistribution : Distribution {
    private:
        double mu, sigma;
    public:
        NormalDistribution(double mu = 0, double sigma = 1) : mu(mu), sigma(sigma) {}

        double Cdf(double t) {
            return Cdf( t, mu, sigma );
        };

        static double Cdf(double t, double mu, double sigma) {
            return ( 1 + erf((t-mu)/(sigma * SQRT_2)) ) / 2;
        };

        double Pdf(double t) {
            // https://en.wikipedia.org/wiki/Normal_distribution
            double twoSigmaSquared = 2 * sigma * sigma;
            double x = t - mu;
            return exp(-x * x / twoSigmaSquared) / sqrt( M_PI * twoSigmaSquared );
        };

        static double Pdf(double t, double mu, double sigma) {
            // https://en.wikipedia.org/wiki/Normal_distribution
            double twoSigmaSquared = 2 * sigma * sigma;
            double x = t - mu;
            return exp(-x * x / twoSigmaSquared) / sqrt( M_PI * twoSigmaSquared );
        };

        double InverseCdf(double p) {
            double z = JohnCook::JC::NormalCDFInverse( p );
            return z * sigma + mu;
        }

        double Mean(void) {
            return mu;
        }

        double StdDev(void) {
            return sigma;
        }
    };
}
