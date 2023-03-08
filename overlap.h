#pragma once
#include <string>
#include <cstdio>
#include <vector>
#include <unordered_map>
#include <memory>
#include <set>
#include <list>
#include <thread>
#include "seqan/seq_io.h"
#include "seqan/align.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graphviz.hpp"

namespace cwd {
	using namespace std;
	template<typename T1, typename T2, typename T3 = std::hash<T1>>
	using hash = std::unordered_multimap<T1, T2, T3>;
	using dnaPos_t = ushort;
	using kmer_t = std::string;

	typedef struct {
		uint readID;
		dnaPos_t begin;
	} hashValue_t;

	struct kmer_hash {
		uint operator()(const kmer_t& kmer) const
		{
			uint hash = 0;
			for (auto dntp : kmer)
			{
				switch (char(dntp))
				{
				case 'A': hash <<= 2;break;
				case 'T': hash = (hash << 2) + 1;break;
				case 'C': hash = (hash << 2) + 2;break;
				case 'G': hash = (hash << 2) + 3;break;
				default:
					break;
				}
			}
			return hash;
		}
	};

	typedef seqan::StringSet<seqan::Dna5String> seqData_t;
	// typedef std::tuple<uint, uint, uint> hashValue_t; // (readId, begin, end)
	typedef hash<kmer_t, hashValue_t> kmerHashTable_t;
	typedef std::pair<uint, uint> numPair_t;
	typedef numPair_t posPair_t;

	typedef struct {
		dnaPos_t SP1;
		dnaPos_t SP2;
		dnaPos_t EP1;
		dnaPos_t EP2;
		bool orient = true;
	} overlapInfo_t;

	typedef struct alignInfo_t {
		bool orient = true;
		dnaPos_t SP1;
		dnaPos_t SP2;

		alignInfo_t(bool orient_, dnaPos_t SP1_, dnaPos_t SP2_)
		{
			orient = orient_;
			SP1 = SP1_;
			SP2 = SP2_;
		}
	} alignInfo_t;

	typedef struct assemblyInfo_t {
		uint r1, r2;
		dnaPos_t SP1;
		dnaPos_t EP1;
		dnaPos_t SP2;
		dnaPos_t EP2;
		bool orient = true;
		//float precision = 0;
	} assemblyInfo_t;

	//using vertex_descriptor = boost::graph_traits<AGraph>::vertex_descriptor;

	uint maxKmerFrequency(std::ifstream& kmerFrequency);
	void filterKmer(kmerHashTable_t& kmerHashTable, const std::string& kfFileName);
	unique_ptr<hash<uint, alignInfo_t>> findSameKmer(kmerHashTable_t& kmerHashTable, seqData_t & seq, uint r);
	bool findSmallerSameKmer(seqData_t& seq, uint r, uint t, uint SKMER_LEN, int s, int s2, int d, bool orient);
	kmerHashTable_t* createKmerHashTable(const seqData_t& seq, bool isFull = false);
	vector<shared_ptr<vector<alignInfo_t>>> chainFromStart(seqData_t& seq, vector<alignInfo_t>& cks, int k, int ks, int alpha, int beta, double gamma, int r, int t);
	vector<assemblyInfo_t> finalOverlap(vector<shared_ptr<vector<alignInfo_t>>>& chain_v, uint len1, uint len2, uint r, uint i, int chainLen, int ovLen);
	void loadSeqData(const std::string& seqFileName, seqan::StringSet<seqan::CharString>& ID, seqData_t& seq);
	void outputOverlapInfo(uint r, uint i, vector<shared_ptr<vector<alignInfo_t>>>& chain_v, seqData_t& seq, seqan::StringSet<seqan::CharString> & ID, ofstream& outFile, int minSize, int chainLen, int ovLen);
	void mainProcess(kmerHashTable_t& kmerHashTable, seqData_t& seq, seqan::StringSet<seqan::CharString> & ID, uint block1, uint block2, ofstream& outFile, int chainLen, int ovLen);
	kmer_t revComp(const kmer_t& kmer);

	template<typename T, typename R>
	vector<alignInfo_t> getCommonKmerSet(T range, R read, const uint KMER_LEN)
	{
		vector<alignInfo_t> cks;
		//unique_ptr<hash<kmer_t, alignInfo_t>> commonKmerSet( new hash<kmer_t, alignInfo_t>() );
		std::set<kmer_t> keySet;
		while (range.first != range.second)
		{
			uint startPos1 = (range.first->second).SP1;
			uint startPos2 = (range.first->second).SP2;
			bool orient = range.first->second.orient;
			kmer_t kmer = { seqan::begin(read) + startPos1, seqan::begin(read) + startPos1 + KMER_LEN };
			if (keySet.insert(kmer).second)
			{
				//commonKmerSet->insert({ kmer, {orient, startPos1, startPos2} });
				cks.emplace_back( orient, startPos1, startPos2 );
			}
			++range.first;
		}
		//for (auto kmer : keySet)
		//{
		//	if (commonKmerSet->count(kmer) > 1)
		//	{
		//		commonKmerSet->erase(kmer);
		//	}
		//}
		//auto iend = chain->begin();
		// cout << CKS.size();
		//for (auto& kmer : commonKmerSet)
		//{
		//	cks.push_back(kmer.second);
		//}

		return cks;
	}

	//bool isConnected(AGraph& g, vertex_descriptor a, vertex_descriptor b);
	//void DFS(cwd::AGraph& g, vertex_descriptor i, vector<bool>& visited);
	std::set<size_t> finalOverlap2(vector<shared_ptr<vector<alignInfo_t>>>& chain_v, uint len1, uint len2, uint r, uint i, int chainLen, int ovLen);
	void mainProcess2(cwd::kmerHashTable_t& kmerHashTable, seqData_t& seq, seqan::StringSet<seqan::CharString>& ID, int block1, int block2, ofstream& outFile, int chainLen, int ovLen, std::set<size_t> & dump);
	void readPAF(const string& file, int minOverlapLen);
}