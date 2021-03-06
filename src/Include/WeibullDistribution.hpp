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
#include "Util.hpp"

namespace QutBio {
	class WeibullDistribution : Distribution {
	private:
		double scale, shape;
	public:
		WeibullDistribution( double scale = 1, double shape = 1 ) : scale( fabs( scale ) ), shape( shape ) {}

		double Cdf( double t ) {
			return Cdf( t, scale, shape );
		};

		static double Cdf( double t, double scale, double shape ) {
			// https://en.wikipedia.org/wiki/Weibull_distribution
			return t < 0 ? 0 : Util::OneMinusExpX( -pow( t / scale, shape ) );
		};

		double Pdf( double t ) {
			return Pdf( t, scale, shape );
		};

		static double Pdf( double t, double scale, double shape ) {
			// https://en.wikipedia.org/wiki/Weibull_distribution
			double a = t / scale;
			double b = pow( a, shape - 1 );
			double c = pow( b, shape );
			return t < 0 ? 0 : b * exp( -c ) * shape / scale;
		};

		double InverseCdf( double p ) {
			return scale * pow( -log( 1 - p ), -shape );
		}

		double Mean( void ) {
			return scale * gamma( 1 + 1 / shape );
		}

		double StdDev( void ) {
			double g1 = gamma( 1 + 2 / shape );
			double g2 = gamma( 1 + 1 / shape );
			return scale * sqrt( g1 + g2 * g2 );
		}

		/**
		 *	Uses linear regression to fit this distribution to an empirical
		 *	(or exact) cumulative distribution function F, evaluated at points
		 *	x.
		 *
		 *	If this works properly, then this.Cdf(x[i]) = F[i] when the operation.
		 *	Is complete.
		 *
		 *	I found that fitting the tails was problematic, so we ignore F[i] < 0.01
		 *	and F[i] > 0.99.
		 *
		 *	Arguments:
		 *	x	List of observation times / positions / distances.
		 *	F	Parallel list of CDF values.
		 */

		void FitToCdf( vector<double> & x, vector<double> & F ) {
			vector<double> logX;
			vector<double> logMinusR;

			for ( size_t i = 0; i < x.size() && i < F.size(); i++ ) {
				if ( x[i] <= 0 || F[i] < 0.01 || F[i] > 0.99 ) continue;

				logX.push_back( log( x[i] ) );
				// logMinusR.push_back( log(-Util::LogOnePlusX(-F[i])));
				logMinusR.push_back( log( -log( 1 - F[i] ) ) );

				//cerr << "\t\t" "x[" << i << "] = " << x[i] << "\n\n";
				//cerr << "\t\t" "F[" << i << "] = " << F[i] << "\n\n";
				//cerr << "\t\t" "logX[" << i << "] = " << logX[i] << "\n";
				//cerr << "\t\t" "logMinusR[" << i << "] = " << logMinusR[i] << "\n\n";
			}

			double a, b;
			Util::LinFit( logX, logMinusR, logX.size(), a, b );

			shape = a;
			scale = exp( -b / a );

			//cerr << "shape = " << shape << "\n";
			//cerr << "scale = " << scale << "\n";
		}

		double Scale() { return scale; }

		double Shape() { return shape; }
	};
}
