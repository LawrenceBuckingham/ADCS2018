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

#include "EnumBase.hpp"

#include <string>
#include <mutex>
#include <vector>

using namespace std;

namespace QutBio {
	typedef class DistanceType * pDistanceType;

	class DistanceType : public EnumBase {
	private:
		DistanceType(string literal, int value) : EnumBase(literal, value) {}
		static std::mutex m;

	public:
		static DistanceType * HalperinEtAl() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("HalperinEtAl", 0);
			return &value;
		}

		static DistanceType * UngappedEdit() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("UngappedEdit", 1);
			return &value;
		}

		static DistanceType * BlosumDistance() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("BlosumDistance", 2);
			return &value;
		}

		/**
		 *	<summary>
		 *		A custom distance which is defined by a similarity matrix loaded at run time.
		 *	</summary>
		 */
		static DistanceType * Custom() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("Custom", 3);
			return &value;
		}

		static vector<EnumBase *> Values() {
			std::vector<EnumBase *> result(4);
			result[0] = HalperinEtAl();
			result[1] = UngappedEdit();
			result[2] = BlosumDistance();
			result[3] = Custom();
			return result;
		}
	};

	typedef DistanceType *pDistanceType;
}
