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

#include <cstdint>

namespace QutBio {

	typedef struct BitRep {
		uint64_t lo, hi;
	} BitRep;

	#define BITS_PER_WORD (sizeof(uint64_t)<<3)

	class CharMap {
	public:
		BitRep bits[128];

		CharMap() {
			memset(bits, 0, sizeof (bits));
		}

		~CharMap() {
		}

		static const CharMap & Blosum62QueryEncoding() {
			static CharMap * encoding = 0;

			if (!encoding) {
				encoding = new CharMap();
				encoding->bits['a'].lo = 4196281838917878893ull;
				encoding->bits['r'].lo = 7650073181085339229ull;
				encoding->bits['n'].lo = 16820669735176575068ull;
				encoding->bits['d'].lo = 13974388523644329108ull;
				encoding->bits['c'].lo = 3754921625820924652ull;
				encoding->bits['q'].lo = 2966207013620391484ull;
				encoding->bits['e'].lo = 3615058547148921981ull;
				encoding->bits['g'].lo = 7081679552086086861ull;
				encoding->bits['h'].lo = 11541817753105046620ull;
				encoding->bits['i'].lo = 1989466549711871335ull;
				encoding->bits['l'].lo = 2034503094600777063ull;
				encoding->bits['k'].lo = 8693214589846654589ull;
				encoding->bits['m'].lo = 8078335720694865167ull;
				encoding->bits['f'].lo = 1334197304103321434ull;
				encoding->bits['p'].lo = 2323611251589552409ull;
				encoding->bits['s'].lo = 3613470385805040733ull;
				encoding->bits['t'].lo = 3560253876602510204ull;
				encoding->bits['w'].lo = 134217727ull;
				encoding->bits['y'].lo = 2454539073132404596ull;
				encoding->bits['v'].lo = 10528300240591231349ull;
				encoding->bits['b'].lo = 16242056208945323541ull;
				encoding->bits['z'].lo = 2417730536003701791ull;
				encoding->bits['x'].lo = 4840157387973009236ull;
				encoding->bits['A'].lo = 4196281838917878893ull;
				encoding->bits['R'].lo = 7650073181085339229ull;
				encoding->bits['N'].lo = 16820669735176575068ull;
				encoding->bits['D'].lo = 13974388523644329108ull;
				encoding->bits['C'].lo = 3754921625820924652ull;
				encoding->bits['Q'].lo = 2966207013620391484ull;
				encoding->bits['E'].lo = 3615058547148921981ull;
				encoding->bits['G'].lo = 7081679552086086861ull;
				encoding->bits['H'].lo = 11541817753105046620ull;
				encoding->bits['I'].lo = 1989466549711871335ull;
				encoding->bits['L'].lo = 2034503094600777063ull;
				encoding->bits['K'].lo = 8693214589846654589ull;
				encoding->bits['M'].lo = 8078335720694865167ull;
				encoding->bits['F'].lo = 1334197304103321434ull;
				encoding->bits['P'].lo = 2323611251589552409ull;
				encoding->bits['S'].lo = 3613470385805040733ull;
				encoding->bits['T'].lo = 3560253876602510204ull;
				encoding->bits['W'].lo = 134217727ull;
				encoding->bits['Y'].lo = 2454539073132404596ull;
				encoding->bits['V'].lo = 10528300240591231349ull;
				encoding->bits['B'].lo = 16242056208945323541ull;
				encoding->bits['Z'].lo = 2417730536003701791ull;
				encoding->bits['X'].lo = 4840157387973009236ull;
			}

			return (*encoding);
		}

		static const CharMap & Blosum62SubjectEncoding() {
			static CharMap * encoding = 0;

			if (!encoding) {
				encoding = new CharMap();
				encoding->bits['a'].lo = 2863761771407970925ull;
				encoding->bits['r'].lo = 7651199062198035261ull;
				encoding->bits['n'].lo = 14505852547472661084ull;
				encoding->bits['d'].lo = 3595913551146720277ull;
				encoding->bits['c'].lo = 3755053567216261860ull;
				encoding->bits['q'].lo = 3006730097971289629ull;
				encoding->bits['e'].lo = 12874384598663773244ull;
				encoding->bits['g'].lo = 7658265648044020940ull;
				encoding->bits['h'].lo = 11541819024448920664ull;
				encoding->bits['i'].lo = 269102453885837161ull;
				encoding->bits['l'].lo = 584915626282040166ull;
				encoding->bits['k'].lo = 6558930587529087837ull;
				encoding->bits['m'].lo = 8073852185476959501ull;
				encoding->bits['f'].lo = 1334828286049501018ull;
				encoding->bits['p'].lo = 7007073522020817209ull;
				encoding->bits['s'].lo = 4262410801802746462ull;
				encoding->bits['t'].lo = 8316072681063168622ull;
				encoding->bits['w'].lo = 134217727ull;
				encoding->bits['y'].lo = 2455735375069421426ull;
				encoding->bits['v'].lo = 17516751889262022129ull;
				encoding->bits['b'].lo = 7054334882014501973ull;
				encoding->bits['z'].lo = 2390568716419798137ull;
				encoding->bits['x'].lo = 5930836213530205298ull;
				encoding->bits['A'].lo = 2863761771407970925ull;
				encoding->bits['R'].lo = 7651199062198035261ull;
				encoding->bits['N'].lo = 14505852547472661084ull;
				encoding->bits['D'].lo = 3595913551146720277ull;
				encoding->bits['C'].lo = 3755053567216261860ull;
				encoding->bits['Q'].lo = 3006730097971289629ull;
				encoding->bits['E'].lo = 12874384598663773244ull;
				encoding->bits['G'].lo = 7658265648044020940ull;
				encoding->bits['H'].lo = 11541819024448920664ull;
				encoding->bits['I'].lo = 269102453885837161ull;
				encoding->bits['L'].lo = 584915626282040166ull;
				encoding->bits['K'].lo = 6558930587529087837ull;
				encoding->bits['M'].lo = 8073852185476959501ull;
				encoding->bits['F'].lo = 1334828286049501018ull;
				encoding->bits['P'].lo = 7007073522020817209ull;
				encoding->bits['S'].lo = 4262410801802746462ull;
				encoding->bits['T'].lo = 8316072681063168622ull;
				encoding->bits['W'].lo = 134217727ull;
				encoding->bits['Y'].lo = 2455735375069421426ull;
				encoding->bits['V'].lo = 17516751889262022129ull;
				encoding->bits['B'].lo = 7054334882014501973ull;
				encoding->bits['Z'].lo = 2390568716419798137ull;
				encoding->bits['X'].lo = 5930836213530205298ull;
			}

			return (*encoding);
		}
	};
}
