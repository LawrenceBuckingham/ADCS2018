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

#include <unordered_map>
#include <vector>

#include "Kmer.hpp"
#include <cstdio>
#include <chrono>
#include <thread>

using namespace std;

namespace QutBio
{

class KmerIndex : public unordered_map<Substring, Kmer *, Substring::Hash>
{
	// TODO: this is not a valid use of inheritance. Convert the map to a member variable.

  protected:
	vector<Kmer *> allKmers;

  public:
	using BaseType = unordered_map<Substring, Kmer *, Substring::Hash>;

	/*
		**	Summary:
		**		Constructs a KmerIndex.
		**
		**	Parameters:
		**		dataset:    a list of sequences;
		**		kmerLength: the word length for tiling;
		*/
	KmerIndex(
		EncodedFastaSequence **dataset,
		size_t length,
		size_t kmerLength //
	)
	{
		for (size_t i = 0; i < length; i++)
		{
			auto seq = dataset[i];
			const char *residues{seq->Sequence().c_str()};
			const size_t seqLen{seq->Sequence().length()};

			if ( seqLen < kmerLength ) continue;

			uint kmerCount = seq->KmerCount( kmerLength );

			for (size_t kmerPos = 0; kmerPos < kmerCount; kmerPos++)
			{
				Substring s(residues, kmerPos, kmerLength);
				auto item = this->find(s);

				if (item == this->end())
				{
					Kmer * kmer = new Kmer(s);
					kmer->Add(seq, kmerPos);
					std::pair<Substring, Kmer *> p(s, kmer);
					BaseType::insert(p);

					// fprintf(stderr, "Inserted %s:%zu\n", seq->Id().c_str(), kmerPos );
				}
				else
				{
					item->second->Add(seq, kmerPos);

					// fprintf(stderr, "Updated %s:%zu\n", seq->Id().c_str(), kmerPos );
				}

				// std::this_thread::sleep_for( std::chrono::seconds( 1 ) );
			}
		}
	}

	KmerIndex(
		const vector<EncodedFastaSequence *> dataset,
		size_t kmerLength //
	) : KmerIndex( (EncodedFastaSequence **)( dataset.data() ), dataset.size(), kmerLength ) // 
	{}

	KmerIndex( 
		const vector<Subsequence> & substrings,
		size_t kmerLength
		//
	) {
		for ( auto & substring: substrings ){
			auto seq = substring.source;
			const char *residues{ seq->Sequence().c_str() };
			uint kmerCount = seq->KmerCount( kmerLength );

			for ( size_t kmerPos = substring.start; kmerPos < kmerCount && kmerPos + kmerLength <= substring.start + substring.length; kmerPos++ ){
				Substring s( residues, kmerPos, kmerLength );
				auto item = this->find( s );

				if ( item == this->end() ){
					Kmer * kmer = new Kmer( s );
					kmer->Add( seq, kmerPos );
					std::pair<Substring, Kmer *> p( s, kmer );
					BaseType::insert( p );
				}
				else{
					item->second->Add( seq, kmerPos );
				}
			}
		}
	}

	KmerIndex( 
		const Subsequence * substrings,
		size_t length,
		size_t kmerLength
		//
	) {
		for ( size_t i = 0; i < length; i++ ){
			auto & substring = substrings[i];
			auto seq = substring.source;
			const char *residues{ seq->Sequence().c_str() };
			uint kmerCount = seq->KmerCount(kmerLength);

			for ( size_t kmerPos = substring.start; kmerPos < kmerCount && kmerPos + kmerLength <= substring.start + substring.length; kmerPos++ ) {
				Substring s( residues, kmerPos, kmerLength );
				auto item = this->find( s );

				if ( item == this->end() ){
					Kmer * kmer = new Kmer( s );
					kmer->Add( seq, kmerPos );
					std::pair<Substring, Kmer *> p( s, kmer );
					BaseType::insert( p );
				}
				else{
					item->second->Add( seq, kmerPos );
				}
			}
		}
	}

	/**
		**	Summary:
		**		Destructor.
		*/
	virtual ~KmerIndex() {
		for (auto &p : *this) {
			delete p.second;
			p.second = 0;
		}
	}

	/**
		**	Summary:
		**		Updates an internal cache of kmer pointers and then returns it by reference.
		*/
	vector<Kmer *> &GetKmers()
	{
		if (allKmers.size() != this->size())
		{
			allKmers.clear();

			for (auto &pair : *this)
			{
				allKmers.push_back(pair.second);
			}
		}

		return allKmers;
	}

	size_t size()
	{
		return BaseType::size();
	}
};
} // namespace QutBio
