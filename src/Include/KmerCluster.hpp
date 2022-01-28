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

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <functional>
#include <numeric>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <csignal>

#include "FastaSequence.hpp"
#include "Kmer.hpp"
#include "KmerClusterPrototype.hpp"
#include "KmerDistanceCache.hpp"
#include "SimilarityMatrix.hpp"
#include "Selector.hpp"
#include "SimilarityMatrix.hpp"
#include "KmerIndex.hpp"

#if defined(SHOW_PROGRESS)
#define PROGRESS(x) x
#else
#define PROGRESS(x)
#endif

#if USE_OMP
#include <omp.h>
#endif

namespace QutBio {
	/**
		**	<summary>
		**		Template class representing cluster of Kmers with a "central" prototype.
		**	</summary>
		**	<typeparam name=DistanceFunction>
		**		A class which provides a function
		**			Distance GetDistance(EncodedKmer, EncodedKmer, size_t);
		**		which returns the distance between two kmers.
		**	</typeparam>
		*/
	template <typename DistanceFunction, typename KmerType>
	struct KmerCluster {
		// using DistanceFunction = function<Distance (KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength )>;
		using Cluster = KmerCluster<DistanceFunction, KmerType>;

		Kmer prototype;
		vector<KmerType> kmers;
		size_t expectedSize;
		size_t index;
		const DistanceFunction &distanceFunction;
		vector<vector<KmerType *>> kmersPerThread;

#if USE_OMP
		omp_lock_t lock;
#endif

		static bool *interrupted() {
			static bool interrupted = false;
			return &interrupted;
		}

		static void sigintHandler( int arg ) {
			auto i = interrupted();
			*i = true;
		}

		KmerCluster(
			Kmer prototype,
			size_t expectedSize,
			const DistanceFunction &distanceFunction ) : prototype( prototype ),
			expectedSize( expectedSize ),
			distanceFunction( distanceFunction ) {

#if USE_OMP
			omp_init_lock( &lock );
#endif
		}

		KmerCluster( const KmerCluster & other ) = delete;
		KmerCluster & operator=( const KmerCluster & other ) = delete;

		//KmerCluster(PrototypeSequenceType * protoSeq, size_t expectedSize, const DistanceFunction & distanceFunction) :
		//	prototype(Substring(protoSeq->Sequence().c_str(), 0, protoSeq->Sequence().length())),
		//	expectedSize(expectedSize),
		//	distanceFunction(distanceFunction)
		//	//
		//{
		//	prototype->Add(protoSeq, 0);
		//}

		virtual ~KmerCluster() {
#if USE_OMP
			omp_destroy_lock( &lock );
#endif
		}

		friend ostream & operator << ( ostream &str, const Cluster & cluster ) {
			str << "Cluster"
				<< "," << cluster.kmers.size() << "," << cluster.prototype << endl;

			for ( auto & kmer : cluster.kmers ) {
				str << kmer << '\n';
			}

			return str;
		}

		/**
		 *	Gets the number of kmer instances assigned to the cluster.
		 *	This will generally be greater than the number of kmers, because each
		 *	kmer may appear (identically) in multiple sequences or positions within
		 *	a sequence.
		 */

		size_t InstanceCount() {
			size_t count = std::accumulate( kmers.begin(), kmers.end(), 0,
				[]( size_t accumulated, KmerType &kmer ) { return accumulated + kmer.Instances().size(); } );
			return count;
		}

		/*
		**	Appends a kmer to the list of kmers attached to this cluster.
		**	Not OMP thread-safe.
		*/

		void Add( const KmerType &kmer ) {
			kmers.push_back( kmer );
		}

		/*
		**	Appends a kmer to the list of kmers attached to this cluster.
		**	OMP thread-safe.
		*/

		void AddParallel( const KmerType &kmer ) {
#if USE_OMP
			omp_set_lock( &lock );
#endif
			kmers.push_back( kmer );
#if USE_OMP
			omp_unset_lock( &lock );
#endif
		}

		void Add( const vector<Kmer> &kmers ) {
			for ( auto &kmer : kmers ) {
				Add( kmer );
			}
		}

		void Add( const vector<Kmer *> &kmers ) {
			for ( auto kmer : kmers ) {
				Add( *kmer );
			}
		}

		void AddParallel( const vector<Kmer> &kmers ) {
#if USE_OMP
			omp_set_lock( &lock );
#endif
			for ( auto &kmer : kmers ) {
				Add( kmer );
			}
#if USE_OMP
			omp_unset_lock( &lock );
#endif
		}

		void AddParallel( const vector<Kmer *> &kmers ) {
#if USE_OMP
			omp_set_lock( &lock );
#endif
			for ( auto &kmer : kmers ) {
				Add( *kmer );
			}
#if USE_OMP
			omp_unset_lock( &lock );
#endif
		}

		size_t Size() const {
			size_t total = 0;

			for ( auto k : kmers ) {
				total += k.Instances().size();
			}

			return total;
		}

		vector<Kmer>::iterator begin() {
			return kmers.begin();
		}

		vector<Kmer>::iterator end() {
			return kmers.end();
		}

		void AllocateKmersToThreads( size_t numThreads ) {
			kmersPerThread.resize( numThreads );

			for ( size_t i = 0; i < kmers.size(); i++ ) {
				kmersPerThread[i % numThreads].push_back( &kmers[i] );
			}
		}

		virtual double DistanceTo( const Kmer &kmer ) {
			EncodedKmer thisCode = prototype.PackedEncoding();
			EncodedKmer kmerCode = kmer.PackedEncoding();
			Distance dist = distanceFunction( thisCode, kmerCode, prototype.Substr().Length() );
			return dist;
		}

		virtual double DistanceTo( EncodedFastaSequence *seq, size_t kmerPosition ) {
			EncodedKmer thisCode = prototype.PackedEncoding();
			EncodedKmer kmerCode = seq->GetEncodedKmer( kmerPosition );
			Distance dist = distanceFunction( thisCode, kmerCode, prototype.Substr().Length() );
			return dist;
		}

		virtual double DistanceTo( EncodedKmer encodedKmer ) {
			EncodedKmer thisCode = prototype.PackedEncoding();
			Distance dist = distanceFunction( thisCode, encodedKmer, prototype.Substr().Length() );
			return dist;
		}

		static void DoIncrementalClusteringParallel(
			uint K,
			double threshold,
			size_t alphaSize,
			const DistanceFunction &symbolCodeDist,
			size_t requestedClusterCount,
			vector<KmerCluster *> &clusters,
			vector<KmerType *> &kmers,
			size_t firstClusterIndex,
			size_t &firstUnallocIndex,
			UniformRealRandom &rand,
			function<KmerClusterPrototype::Pointer( Kmer *kmer )> createPrototype
			//
		) {
			if ( clusters.size() == 0 && requestedClusterCount == 0 ) {
				cerr << "clusters.size() == 0 && requestedClusterCount == 0\n\n";
				throw Exception( "clusters.size() == 0 && requestedClusterCount == 0", FileAndLine );
			}

			// TODO: Re-factor so that the use of requestedClusterCount and firstClusterIndex is not needed.

			if ( requestedClusterCount > 0 ) {
				size_t available = kmers.size() - firstUnallocIndex;
				size_t wanted = std::min( requestedClusterCount, available );

				for ( size_t i = firstUnallocIndex; i < firstUnallocIndex + wanted; i++ ) {
					KmerClusterPrototype *proto = createPrototype( kmers[i] );
					Kmer kmer = proto->SingletonKmer();
					clusters.push_back( new KmerCluster( kmer, 0, symbolCodeDist ) );
				}
			}

#pragma omp parallel for
			for ( size_t i = firstUnallocIndex; i < kmers.size(); i++ ) {
				KmerType *kmer = kmers[i];
				auto kmerEncoding = kmers[i]->PackedEncoding();
				KmerCluster *attractor = 0;
				auto attractorDistance = numeric_limits<Distance>::max();

				for ( size_t j = firstClusterIndex; j < clusters.size(); j++ ) {
					// Try to map kmer to a pre-existing cluster if possible, otherwise search the
					// newly added clusters.
					auto prototypeEncoding = clusters[j]->prototype.PackedEncoding();

					// Always go with the first cluster that falls within the threshold.
					// This way, any well-conserved kmers will be quickly mapped into clusters.
					// Hopefully, these are in the majority.
					auto dist = symbolCodeDist( kmerEncoding, prototypeEncoding, K );

					if ( dist <= threshold ) {
						attractor = clusters[j];
						attractorDistance = dist;
						break;
					}
				}

				if ( attractor ) {
					kmer->DistanceFromPrototype( attractorDistance );
					attractor->AddParallel( *kmer );
#pragma omp critical
					{
						if ( i > firstUnallocIndex ) {
							std::swap( kmers[firstUnallocIndex], kmers[i] );
						}
						firstUnallocIndex++;
					}
				}
			}
		}

		static void DoExhaustiveIncrementalClustering(
			KmerIndex &kmerIndex,
			uint wordLength,
			double threshold,
			size_t alphaSize,
			const DistanceFunction &distanceFunction,
			UniformRealRandom &rand,
			size_t clusterIncrement,
			function<KmerClusterPrototype *( Kmer *kmer )> createPrototype,
			vector<KmerCluster *> &clusters
			//
		) {
			size_t initialClusterCount = clusters.size();
			size_t N = kmerIndex.size();
			size_t firstUnallocIndex = 0;
			vector<KmerType *> allKmers = kmerIndex.GetKmers();

			// Move any kmers that cannot possibly land in a cluster to a place where they will not
			// mess things up. We don't want them in any clusters anyway.
			// This applies to data which uses a distance derived from a similarity matrix rather
			// than a metric.
			for ( size_t i = 0; i < N; i++ ) {
				auto kmer = allKmers[i];
				auto encoding = kmer->PackedEncoding();
				Distance selfMatchDistance = distanceFunction( encoding, encoding, wordLength );

				if ( selfMatchDistance > threshold ) {
					if ( i > firstUnallocIndex ) {
						std::swap( allKmers[i], allKmers[firstUnallocIndex] );
						firstUnallocIndex++;
					}
				}

				i++;
			}

			// Shuffle the kmers to avoid pathological selection.
			for ( size_t i = firstUnallocIndex; i < N; i++ ) {
				auto newLoc = firstUnallocIndex + (size_t) ( rand() * ( N - firstUnallocIndex ) );
				std::swap( allKmers[i], allKmers[newLoc] );
			}

			size_t firstClusterIndex = 0;
			size_t increment = clusters.size() > 0 ? 0 : clusterIncrement;

			while ( firstUnallocIndex < N ) {
				( cerr << "\r" << ( N - firstUnallocIndex ) << " unassigned kmers.                               " ).flush();
				auto previous = firstUnallocIndex;

				// On the first pass, see if we can extend the contents of previously discovered clusters.

				DoIncrementalClusteringParallel(
					wordLength,
					threshold,
					alphaSize,
					distanceFunction,
					increment,
					clusters,
					allKmers,
					firstClusterIndex,
					firstUnallocIndex,
					rand,
					createPrototype );

				if ( firstUnallocIndex == previous ) {
					// cerr << "No further clustering possible. Finishing.\n";
					break;
				}

				increment = clusterIncrement;
				firstClusterIndex = clusters.size();
			}

			( cerr << "Adding " << ( clusters.size() - initialClusterCount ) << " new clusters...\n" ).flush();
		}

		//-----------------------------------------------------------------------------------------
		//	Version that assigns a fixed range of kmers to each thread.
		//-----------------------------------------------------------------------------------------

		static void DoIncrementalClusteringParallel(
			uint K,
			double threshold,
			size_t alphaSize,
			const DistanceFunction &symbolCodeDist,
			size_t requestedClusterCount,
			vector<KmerCluster *> &clusters,
			vector<KmerType *> &kmers,
			size_t firstClusterIndex,
			vector<size_t> &firstUnallocIndex,
			UniformRealRandom &rand,
			function<KmerClusterPrototype::Pointer( Kmer *kmer )> createPrototype,
			size_t numThreads
			//
		) {
			if ( clusters.size() == 0 && requestedClusterCount == 0 ) {
				cerr << "clusters.size() == 0 && requestedClusterCount == 0\n\n";
				throw Exception( "clusters.size() == 0 && requestedClusterCount == 0", FileAndLine );
			}

			size_t N = kmers.size();

			// TODO: Re-factor so that the use of requestedClusterCount and firstClusterIndex is not needed.

			if ( requestedClusterCount > 0 ) {
				for ( size_t threadId = 0; threadId < numThreads; threadId++ ) {
					// Don't parallel this... createPrototype is not thread safe.

					size_t endIdx = ( threadId + 1 ) * N / numThreads;
					size_t available = endIdx - firstUnallocIndex[threadId];

					size_t desired = ( threadId + 1 ) * requestedClusterCount / numThreads - threadId * requestedClusterCount / numThreads;
					size_t wanted = std::min( desired, available );

					if ( false ) {
						( cerr << "Thread " << threadId << ": wanted = " << wanted << "\n" ).flush();
						( cerr << "Thread " << threadId << ": firstUnallocIndex[threadId] = " << firstUnallocIndex[threadId] << "\n" ).flush();
					}

					for ( size_t i = firstUnallocIndex[threadId]; i < firstUnallocIndex[threadId] + wanted; i++ ) {
						KmerClusterPrototype *proto = createPrototype( kmers[i] );
						Kmer kmer = proto->SingletonKmer();
						KmerCluster *cluster = new KmerCluster( kmer, 0, symbolCodeDist );

						{
							if ( false )
								( cerr << "Thread " << threadId << ": adding a cluster\n" ).flush();
							clusters.push_back( cluster );
						}
					}
				}
			}

#if USE_OMP
#pragma omp parallel
#endif
			{
#if USE_OMP
				size_t threadId = omp_get_thread_num();
#else
				size_t threadId = 0;
#endif

				// Rely on firstUnallocIndex for the beginning of the unprocessed part of the range.
				size_t endIdx = ( threadId + 1 ) * N / numThreads;

				for ( size_t i = firstUnallocIndex[threadId]; i < endIdx; i++ ) {
					KmerType *kmer = kmers[i];
					KmerWord *kmerEncoding = kmers[i]->PackedEncoding();
					KmerCluster *attractor = 0;
					auto attractorDistance = numeric_limits<Distance>::max();

					for ( size_t j = firstClusterIndex; j < clusters.size(); j++ ) {
						// Try to map kmer to a pre-existing cluster if possible, otherwise search the
						// newly added clusters.
						auto prototypeEncoding = clusters[j]->prototype.PackedEncoding();

						// Always go with the first cluster that falls within the threshold.
						// This way, any well-conserved kmers will be quickly mapped into clusters.
						// Hopefully, these are in the majority.
						auto dist = symbolCodeDist( kmerEncoding, prototypeEncoding, K );

						if ( dist <= threshold ) {
							attractor = clusters[j];
							attractorDistance = dist;
							break;
						}
					}

					if ( attractor ) {
						kmer->DistanceFromPrototype( attractorDistance );
						attractor->AddParallel( *kmer );

						if ( i > firstUnallocIndex[threadId] ) {
							std::swap( kmers[firstUnallocIndex[threadId]], kmers[i] );
							firstUnallocIndex[threadId]++;
						}
					}
				}
			}
		}

		static size_t GetUnallocated( size_t numThreads, size_t N, vector<size_t> &firstUnallocIndex ) {
			size_t numAllocated = 0;

			for ( size_t i = 0; i < numThreads; i++ ) {
				size_t beginIdx = i * N / numThreads;
				numAllocated += firstUnallocIndex[i] - beginIdx;
			}

			return numAllocated;
		}

		static void DoExhaustiveIncrementalClustering(
			KmerIndex &kmerIndex,
			uint wordLength,
			double threshold,
			size_t alphaSize,
			const DistanceFunction &distanceFunction,
			UniformRealRandom &rand,
			size_t clusterIncrement,
			function<KmerClusterPrototype *( Kmer *kmer )> createPrototype,
			vector<KmerCluster *> &clusters,
			size_t numThreads
			//
		) {
			size_t initialClusterCount = clusters.size();
			size_t N = kmerIndex.size();
			vector<KmerType *> allKmers = kmerIndex.GetKmers();
			vector<size_t> firstUnallocIndex( numThreads );

#if USE_OMP
#pragma omp parallel
#endif
			{
#if USE_OMP
				size_t threadId = omp_get_thread_num();
#else
				size_t threadId = 0;
#endif
				// Set up a range of kmers which will be processed solely by the current thread.
				size_t beginIdx = threadId * N / numThreads;
				size_t endIdx = ( threadId + 1 ) * N / numThreads;

				// High-water mark for this thread is initially set to the beginning of the
				// range.
				firstUnallocIndex[threadId] = beginIdx;

				// Move any kmers that cannot possibly land in a cluster to a place where they will not
				// mess things up. We don't want them in any clusters anyway.
				// This applies to data which uses a distance derived from a similarity matrix rather
				// than a metric.
				for ( size_t i = beginIdx; i < endIdx; i++ ) {
					auto kmer = allKmers[i];
					auto encoding = kmer->PackedEncoding();
					Distance selfMatchDistance = distanceFunction( encoding, encoding, wordLength );

					if ( selfMatchDistance > threshold ) {
						if ( i > firstUnallocIndex[threadId] ) {
							std::swap( allKmers[i], allKmers[firstUnallocIndex[threadId]] );
							firstUnallocIndex[threadId]++;
						}
					}
				}

				// Shuffle the kmers to avoid pathological selection.
				for ( size_t i = firstUnallocIndex[threadId]; i < endIdx; i++ ) {
					auto newLoc = firstUnallocIndex[threadId] + (size_t) ( rand() * ( endIdx - firstUnallocIndex[threadId] ) );
					std::swap( allKmers[i], allKmers[newLoc] );
				}
			}

			size_t firstClusterIndex = 0;
			size_t increment = clusters.size() > 0 ? 0 : clusterIncrement;

			size_t numAllocated = GetUnallocated( numThreads, N, firstUnallocIndex );

			( cerr << "N = " << N
				<< "\nnumAllocated = " << numAllocated
				<< "\n" )
				.flush();

			while ( numAllocated < N ) {
				( cerr << "\r" << ( N - numAllocated ) << " unassigned kmers.                               " ).flush();
				auto previous = numAllocated;

				// On the first pass, see if we can extend the contents of previously discovered clusters.

				DoIncrementalClusteringParallel(
					wordLength,
					threshold,
					alphaSize,
					distanceFunction,
					increment,
					clusters,
					allKmers,
					firstClusterIndex,
					firstUnallocIndex,
					rand,
					createPrototype,
					numThreads );

				numAllocated = GetUnallocated( numThreads, N, firstUnallocIndex );

				if ( numAllocated == previous ) {
					// cerr << "No further clustering possible. Finishing.\n";
					break;
				}

				increment = clusterIncrement;
				firstClusterIndex = clusters.size();
			}

			( cerr << "Adding " << ( clusters.size() - initialClusterCount ) << " new clusters...\n" ).flush();
		}

		static void InitialiseClusters(
			PointerList<EncodedFastaSequence> &protos,
			size_t wordLength,
			const DistanceFunction &dist,
			vector<KmerCluster *> &clusters ) {
			for ( auto proto_ : protos ) {
				pKmerClusterPrototype proto = (pKmerClusterPrototype) proto_;
				Kmer protoKmer = proto->SingletonKmer();
				clusters.push_back( new KmerCluster( protoKmer, 0, dist ) );
			}
		}
	};
} // namespace QutBio
