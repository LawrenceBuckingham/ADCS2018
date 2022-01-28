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

#include <unordered_map>

#include "Assert.hpp"
#include "Delegates.hpp"

using namespace QutBio;

namespace QutBio {
	/**
	<summary>
		A lookup table mapping some key of type K to _pointer_to_ values of type V, which owns the stored pointers.
		The values are dynamically allocated when a value is added by a factory function which is responsible for
		correct initialisation.
		When the collection is deleted the pointers are destroyed.
	<summary>
		*/

	template<typename K, typename V>
	class LookupTable {
	private:
		unordered_map<K, V *> map;
	public:
		LookupTable() {}

		~LookupTable() {
			for ( auto kv : map ) {
				delete kv.second;
			}
		}

		LookupTable( const LookupTable & other ) = delete;

		LookupTable & operator = ( const LookupTable & other ) = delete;

		void Add( const K & key, Func<V *> factory ) {
			assert_false( map.find( key ) != map.end() );
			V * value = factory();
			map[key] = value;
		}

		V * operator[] ( const K & key ) {
			return map[key];
		}

		void ForEach( Action2<K&, V *> action ) {
			for ( auto kv : map ) {
				action( kv.first, kv.second );
			}
		}
	};
}
