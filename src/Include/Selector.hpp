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

#include <cstdlib>
#include <sstream>

#include "HBRandom.hpp"


namespace QutBio {
	/// <summary>
	/// Selects "numberWanted" choices from a candidate pool containing "outOf" items 
	/// according to a uniform distribution.
	/// </summary>
	class Selector {
		size_t    numberWanted = 0;
		size_t    outOf = 0;
		size_t	stillWanted = 0;
		size_t	remaining = 0;
		UniformRealRandom & rand;

	public:
		Selector(UniformRealRandom & rand, size_t numberWanted, size_t outOf) :
			numberWanted(numberWanted), outOf(outOf), rand(rand) {
			Reset();
			Validate();
		}

		/// <summary>
		/// This property should be called in sequence until the number wanted 
		/// has been obtained, which is guaranteed to happen by the time the 
		/// the property has been evaluated "out of" times.
		/// <para>Side effects: outOf and numberWanted are updated as part of the decision process.</para>
		/// </summary>
		bool SelectThis(void) {
			Validate();

			if ( stillWanted == 0 || remaining == 0 ) return false;

			double probability = (double) stillWanted / remaining;
			double random = rand();

			remaining--;

			if ( random <= probability ) {
				stillWanted--;
				return true;
			}
			else {
				return false;
			}
		}

		/// <summary>
		/// Gets the number of items that have yet to be chosen.
		/// </summary>
		size_t StillWanted(void) {
			return stillWanted;
		}

		/// <summary>
		/// Gets the number of candidate items available for selection.
		/// </summary>
		size_t Remaining(void) {
			return remaining;
		}

		/// <summary>
		/// Gets the original number wanted.
		/// </summary>
		size_t NumberWanted(void) {
			return numberWanted;
		}

		/// <summary>
		/// Gets the original size of the pool of candidates.
		/// </summary>
		size_t OutOf(void) {
			return outOf;
		}

		/// <summary>
		/// Restore the initial values of Remaining and StillWanted in 
		/// preparation for a new run of selections.
		/// </summary>
		void Reset(void) {
			remaining = outOf;
			stillWanted = numberWanted;
		}

	private:
		void Validate() {
			if ( numberWanted > outOf ) {
				stringstream str;
				str << "selector: numberWanted > outOf\n\tnumberWanted = " << numberWanted << "\n\toutOf = " << outOf;
				throw Exception(str.str(), __FILE__, __LINE__);
			}

			if ( stillWanted > remaining ) {
				stringstream str;
				str << "selector: stillWanted > remaining\n\tstillWanted = " << stillWanted << "\n\tremaining = " << remaining;
				throw Exception(str.str(), __FILE__, __LINE__);
			}

			if ( remaining > outOf ) {
				stringstream str;
				str << "selector: remaining > outOf\n\tremaining = " << remaining << "\n\toutOf = " << outOf;
				throw Exception(str.str(), __FILE__, __LINE__);
			}
		}
	};
}
