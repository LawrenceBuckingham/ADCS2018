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

#include <cstring>
#include <cstdio>
#include <ios>
#include <functional>

namespace QutBio {
	// Class representing a sub-string of a raw char array. 
	class Substring {
	private:
		// The address of the first character in the instance. 
		// NB: strlen(string) > length with high probability
		const char * chars;

		// The number of characters in the kmer.
		size_t length;

	public:
		Substring(const char * str, size_t start, size_t length) : chars(str + start), length(length) {}

		virtual ~Substring() {}

		const char * Chars() const { return chars; }

		size_t Length() const { return length; }

		size_t size() const { return length; }

		const char & operator[](size_t index) const {
			return chars[index];
		}

		friend bool operator==(const Substring & lhs, const Substring & rhs) {
			return strncmp(lhs.chars, rhs.chars, lhs.length) == 0;
		}

		friend bool operator!=(const Substring & lhs, const Substring & rhs) {
			return strncmp(lhs.chars, rhs.chars, lhs.length) != 0;
		}

		friend bool operator<(const Substring & lhs, const Substring & rhs) {
			return strncmp(lhs.chars, rhs.chars, lhs.length) < 0;
		}

		friend std::ostream & operator<<(std::ostream & str, const Substring & s) {
			for (size_t i = 0; i < s.length; i++) {
				str << s.chars[i];
			}

			return str;
		}

		struct Hash {
			size_t operator()(const Substring & __val) const noexcept {
				//fprintf( stderr, "Hashing ");
				//__val.fprint( stderr ); 
				//fputc( '\n', stderr );
				return _Hash_impl::hash((void *) __val.chars, __val.length);
			}
		};

		void fprint( FILE * stream ) const {
			for (size_t i = 0; i < length; i++) {
				fputc( chars[i], stream );
			}
		}
	};

	typedef Substring * pSubstring;
}
