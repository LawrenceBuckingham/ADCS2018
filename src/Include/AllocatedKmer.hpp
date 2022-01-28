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
#include <atomic>

#include "SimilarityMatrix.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"

namespace QutBio {
#if USE_OMP
#include <omp.h>
	class AllocatedKmer : public Kmer {
		bool isAllocated{ 0 };
		omp_lock_t lock;
	public:
		AllocatedKmer(const Substring & substring) : Kmer(substring) {
			omp_init_lock(&lock);
		}

		AllocatedKmer(const AllocatedKmer & other) : Kmer(other) {
			omp_init_lock(&lock);
			isAllocated = other.isAllocated;
		}

		void operator=(const AllocatedKmer & other) {
			Kmer::operator=(other);
			omp_init_lock(&lock);
			isAllocated = other.isAllocated;
		}

		virtual ~AllocatedKmer() {
			omp_destroy_lock(&lock);
		}

		bool IsAllocated() {
			return isAllocated;
		}

		void Allocate() {
			isAllocated = true;
		}

		void Lock() {
			omp_set_lock(&lock);
		}

		void Unlock() {
			omp_unset_lock(&lock);
		}
	};
#else
	struct AllocatedKmer : Kmer {
		void * connection;

		AllocatedKmer(const Substring & substring) : Kmer(substring) {}

		virtual ~AllocatedKmer() {}

		void Lock() {}

		void Unlock() {}
	};
#endif

}
