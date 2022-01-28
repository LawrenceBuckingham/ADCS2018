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
	#include "Histogram.hpp"
	#include "SimilarityMatrix.hpp"

	class KmerDistributions {
	public:
		static Histogram<double> GetOneMerDistanceDistribution(
			const SimilarityMatrix & blosum62,
			const Histogram<char> & symbolDist,
			Histogram<double> & oneMerDistances
			) {
			string alphabet = blosum62.Symbols();
			auto idx = alphabet.find('*');

			if ( idx != string::npos ) {
				alphabet.erase(idx, 1);
			}

			oneMerDistances.GetOneMerHistogram<char, char>(
				symbolDist,
				[blosum62](char x, char y) { return blosum62.MaxValue() - blosum62.Similarity(x, y); }
			);

			return oneMerDistances;
		}

		static Histogram<double> GetOneMerSimilarityDistribution(
			const SimilarityMatrix & blosum62,
			const Histogram<char> & symbolDist,
			Histogram<double> & oneMerDistances
			) {
			string alphabet = blosum62.Symbols();
			auto idx = alphabet.find('*');

			if ( idx != string::npos ) {
				alphabet.erase(idx, 1);
			}

			oneMerDistances.GetOneMerHistogram<char, char>(
				symbolDist,
				[blosum62](char x, char y) { return blosum62.Similarity(x, y); }
			);

			return oneMerDistances;
		}

		static void GetHousdorffAverageFragmentDistributions(
			int maxK,
			int fragLength,
			const Histogram<double> & oneMerDistances,
			map<int, DiscreteDistribution> & hausdorffFragmentDistributions
			) {
			Histogram<double> kmerDistances = oneMerDistances;

			for ( int k = 2; k <= maxK; k++ ) {
				Histogram<double> newHistogram;
				kmerDistances.DoConvolution(oneMerDistances, newHistogram);
				kmerDistances = newHistogram;
				DiscreteDistribution discrete;
				DiscreteDistribution minDist;
				discrete.SetPmf(kmerDistances);
				discrete.GetMinimumDistribution(fragLength, minDist);

				Histogram<double> currentSum = minDist.Pmf();

				for ( int f = 2; f <= fragLength; f++ ) {
					Histogram<double> newSum;
					currentSum.DoConvolution(minDist.Pmf(), newSum);
					newSum.Cleanup([](double key, double value) { return value <= 0; });
					currentSum = newSum;
				}

				Histogram<double> averagePmf;

				for ( auto p : currentSum.data ) {
					averagePmf.data[p.first / fragLength] = p.second;
				}

				DiscreteDistribution averageDistribution;
				averageDistribution.SetPmf(averagePmf);

				DiscreteDistribution hausdorffAverageDistribution;
				averageDistribution.GetMaximumDistribution(2, hausdorffAverageDistribution);
				hausdorffAverageDistribution.Cleanup();
				hausdorffFragmentDistributions[k] = hausdorffAverageDistribution;
			}
		}

		static void GetHausdorffAverageFragmentDistributions(
			int maxK,
			int fragLength,
			const SimilarityMatrix & similarityMatrix,
			const Histogram<char> & symbolDist,
			map<int, DiscreteDistribution> & hausdorffFragmentDistributions
			) {
			Histogram<double> oneMerDistances;
			GetOneMerDistanceDistribution(similarityMatrix, symbolDist, oneMerDistances);
			GetHousdorffAverageFragmentDistributions(maxK, fragLength, oneMerDistances, hausdorffFragmentDistributions);
		}
	};
}