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
	typedef unsigned char byte;
	typedef unsigned int uint;
	typedef unsigned long ulong;

	template <typename T1, typename T2, typename T3>
	struct Triple {
		T1 item1;
		T2 item2;
		T3 item3;

		Triple() {}

		Triple(const T1 & item1, const T2 & item2, const T3 & item3) :
			item1(item1), item2(item2), item3(item3) {}

		friend bool operator<(const Triple & lhs, const Triple & rhs) {
			if (lhs.item1 < rhs.item1) return true;
			if (rhs.item1 < lhs.item1) return false;
			if (lhs.item2 < rhs.item2) return true;
			if (rhs.item2 < lhs.item2) return false;
			if (lhs.item3 < rhs.item3) return true;
			if (rhs.item3 < lhs.item3) return false;
			return false;
		}
	};
}
