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

#include <random>

namespace QutBio {
	struct UniformRealRandom {
		std::mt19937 rng;
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution;

		UniformRealRandom(unsigned long seed) : rng(seed), generator(rng()), distribution(0.0, 1.0) {
		}

		double operator()() {
			return distribution(generator);
		}
	};

	template<typename T = int>
	struct UniformIntRandom {
		std::mt19937 rng;
		std::default_random_engine generator;
		std::uniform_int_distribution<T> distribution;

		UniformIntRandom(unsigned long seed, T min, T max) : rng(seed), generator(rng()), distribution( min, max ) {
		}

		/**
		 *	Call without arguments to generate a random value in the stored interval, [min,max].
		 */
		T operator()() {
			auto x = distribution(generator);
			return x;
		}

		/**
		 *	Call with explicit values of min and max to ad-hoc sample a different interval than 
		 *	that which is stored.
		 */
		T operator()( T min, T max ) {
			std::uniform_int_distribution<T> d( min, max);
			return d(generator);
		}
	};
}
