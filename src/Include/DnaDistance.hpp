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

#include "DistanceType.hpp"

namespace QutBio {

	// Making this an object instead of a static function will let me
	// use it in the KmerClusterAL template class.
	class DnaDistance {
#if INSTRUMENT_DNA_DIST
		size_t callCounter = 0;
#endif
	public:

#define MASK1 (0x5555555555555555ULL)
#define MASK2 (MASK1<<1)

		// DNA kmers are always encoded in a single 64 bit word,
		// so we can just calculate the Hamming distance by bit 
		// operations.
		Distance GetDistance1(KmerWord x, KmerWord y) const {
			uint64_t a = x ^ y;
			uint64_t b = a & MASK1;
			uint64_t c = (a & MASK2) >> 1;
			uint64_t d = b | c;
			return POPCOUNT(d);
		}

		// Hot spot: lines 38, 39
		Distance operator()(KmerWord *x, KmerWord *y, uint kmerLength) const {
#if INSTRUMENT_DNA_DIST
			callCounter++;
#endif
			Distance dist = 0;

			for (size_t i = 0; i * (CHAR_BIT * sizeof(KmerWord) / 2) < kmerLength; i++) {
				uint64_t a = x[i] ^ y[i];
				uint64_t b = a & MASK1;
				uint64_t c = (a & MASK2) >> 1;
				uint64_t d = b | c;
				dist += POPCOUNT(d);
			}

			return dist;
		}

#undef MASK1
#undef MASK2

#if INSTRUMENT_DNA_DIST
		size_t CallCounter() {
			return callCounter;
		}
#endif

	};
}
