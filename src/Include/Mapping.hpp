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

#include "DoubleArrayExtensions.hpp"

namespace Lmvq {

	using DAX = DoubleArrayExtensions;

	class Mapping {
	public:
		virtual void Map(double *x, double *y) = 0;

		virtual void Map(double **x, double **y, size_t n ) {
			for ( size_t i = 0; i < n; i++ ) {
				Map(x[i], y[i]);
			}
		}

		virtual uint InputDimension(void) = 0;

		virtual uint OutputDimension(void) = 0;

		virtual bool IsDefinedAt(double *x) {
			return true;
		}

		virtual bool IsDifferentiableAt(double *x) {
			return false;
		}

		double GetAverageClassificationError(double ** x, double ** y, size_t n) {
			return GetClassificationError(x, y, 0) / n;
		}

		virtual double GetClassificationError(
			double ** x,
			double ** y,
			size_t n,
			double * individual_error = 0
			) {
			int C = OutputDimension();
			double e = 0;

			double * y_predicted = new double[C];

			for (size_t t = 0; t < n; t++ ) {
				Classify(x[t], y_predicted);

				double this_error = DAX::CompareTo( y[t], C, y_predicted ) == 0 ? 0 : 1;

				e += this_error;

				if ( individual_error ) {
					individual_error[t] = this_error;
				}
			}

			return e;
		}

		double GetAverageEdge(
			double ** x,
			double ** y,
			size_t n,
			double * individual_edge = 0
			) {
			return GetTotalEdge(x, y, n, individual_edge) / n;
		}

		virtual double GetTotalEdge(
			double ** x,
			double ** y,
			size_t n,
			double * individual_edge = 0
			) {
			int C = OutputDimension();
			double    e = 0;

			double* y_predicted = new double[C];

			for (size_t t = 0; t < n; t++ ) {
				Map(x[t], y_predicted);

				double this_edge = DAX::Dot(y[t], n, y_predicted);

				e += this_edge;

				if ( individual_edge ) {
					individual_edge[t] = this_edge;
				}
			}

			delete[] y_predicted;
			return e;
		}

		virtual void Classify(
			double * x,
			double * y
			) {
			DAX::HardLimit(y, OutputDimension());
		}

		virtual void Classify(
			double ** x,
			double ** y,
			size_t n
			) {
			for (size_t i = 0; i < n; i++ ) {
				Classify(x[i], y[i]);
			}
		}
	};
}
