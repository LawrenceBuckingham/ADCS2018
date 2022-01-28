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

namespace QutBio {
	class Distribution {
	public:
		virtual ~Distribution() {}

		virtual double Cdf( double t ) = 0;
		virtual double Pdf( double t ) = 0;
		virtual double InverseCdf( double t ) = 0;
		virtual double Mean( void ) = 0;
		virtual double StdDev( void ) = 0;

		// Gets (min, max) such that min = inf{t:Cdf(t) > 0} and
		// max = sup{t:Cdf(t) < 1}.
		void GetSupport( double & min, double & max, double epsilon = 1e-10 ) {
			double hi, lo;

			// Get min:
			hi = Mean();
			lo = hi - 100 * StdDev();

			while ( fabs(hi - lo) > epsilon ) {
				min = ( lo + hi ) / 2;
				//( cerr << "lo = " << lo << "\n" ).flush();
				//( cerr << "hi = " << hi << "\n" ).flush();
				//( cerr << "fabs(hi - lo) = " << fabs( hi - lo ) << "\n" ).flush();
				//( cerr << "min = " << min << "\n" ).flush();

				if ( Cdf( min ) <= 0 ) {
					lo = min;
				}
				else {
					hi = min;
				}
			}

			min = ( lo + hi ) / 2;

			// Get max:
			lo = Mean();
			hi = lo + 100 * StdDev();

			while ( fabs( hi - lo ) > epsilon ) {
				max = ( lo + hi ) / 2;
				//(cerr << "max = " << max << "\n").flush();

				if ( Cdf( max ) >= 1 ) {
					hi = max;
				}
				else {
					lo = max;
				}
			}

			max = ( lo + hi ) / 2;
		}
	};

	class ScaledDistribution : public Distribution {
		double scale;
		Distribution &baseDistribution;
	public:
		ScaledDistribution(
			double scale,
			Distribution &baseDistribution
			//
		) : scale( scale ), baseDistribution( baseDistribution ) {}

		virtual ~ScaledDistribution() {}

		virtual double Cdf( double t ) {
			return baseDistribution.Cdf( t / scale );
		}

		virtual double Pdf( double t ) {
			return baseDistribution.Pdf( t / scale );
		}

		virtual double InverseCdf( double t ) {
			return baseDistribution.InverseCdf(t) * scale;
		}

		virtual double Mean( void ) {
			return scale * baseDistribution.Mean();
		}

		virtual double StdDev( void ) {
			return fabs(scale) * baseDistribution.StdDev();
		}

	};

}
