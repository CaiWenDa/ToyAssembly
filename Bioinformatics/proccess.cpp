#include "proccess.h"

using namespace std;
using namespace seqan;
using namespace cwd;

kmerHashTable_t& cwd::createKmerHashTable(const seqData_t& seq)
{
	kmerHashTable_t* kmerHashTable = new kmerHashTable_t();
	uint readID = 0;
	for (auto& read : seq)
	{
		uint pos = 0;
		for (auto i = begin(read); i < end(read) - KMER_LEN; i += KMER_STEP, pos += KMER_STEP)
		{
			kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
			hashValue_t kmerInfo{ readID, pos }; //, pos + KMER_LEN };
			kmerHashTable->insert({ kmer, kmerInfo });
		}
		readID++;
	}
	cout << "Created HashTable!\n";
	return *kmerHashTable;
}

unique_ptr<list<alignInfo_t>> cwd::chainFromStart(hash<kmer_t, alignInfo_t>& CKS, int k, int ks, int alpha, double beta, double gamma)
{
	unique_ptr<list<alignInfo_t>> chain( new list<alignInfo_t>() );
	chain->push_back(CKS.begin()->second);
	vector<alignInfo_t> cks;
	// cout << CKS.size();
	for (auto& kmer : CKS)
	{
		cks.push_back(kmer.second);
	}
	sort(cks.begin(), cks.end(), [](alignInfo_t a, alignInfo_t b) { return a.SP1 < b.SP1; });
	for (auto ix = cks.begin(), nextx = next(ix); ix != cks.end() && nextx != cks.end();)
	{
		uint d1 = 0;
		uint d2 = 0;
		if (ix->SP2 < nextx->SP2)
		{
			d1 = nextx->SP1 - ix->SP1;
			d2 = nextx->SP2 - ix->SP2;
			if ((d1 < alpha && d2 < alpha) and (double(max(d1, d2) - min(d1, d2)) / max(d1, d2) < gamma))
			{
				chain->push_back(*ix);
				chain->push_back(*nextx);
				ix = next(nextx);
				nextx = next(ix);
			}
			else
			{
				ix = nextx;
				nextx++;
			}
		}
		else
		{
			nextx++;
			continue;
		}
		
	}
	return chain;
}

uint cwd::maxKmerFrequency(ifstream& kmerFrequency)
{
	vector<uint> vFrequency;
	if (kmerFrequency.is_open())
	{
		kmer_t kmer;
		uint frequency;
		uint maxFrequency = 0;
		while (kmerFrequency >> kmer >> frequency)
		{
			vFrequency.push_back(frequency);
			if (maxFrequency < frequency)
			{
				maxFrequency = frequency;
			}
		}
		auto count = new uint[maxFrequency];
		memset(count, 0, maxFrequency);
		for (auto& x : vFrequency)
		{
			count[x]++;
		}
		double s = 0;
		for (uint i = 0; i < maxFrequency; i++)
		{
			// cout << count[i] << " ";
			s += count[i];
			if (s > 0.9 * maxFrequency)
			{
				if (i > 2) return i;
				else return 3;
			}
		}
	}
	return 0;
}

overlapInfo_t cwd::finalOverlap(list<alignInfo_t>& chain, uint len1, uint len2)
{
	uint P1 = chain.begin()->SP1; // P1
	uint Q1 = min_element(chain.begin(), chain.end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//chain.begin()->SP2; // Q1
	uint Pnk = prev(chain.end())->SP1 + KMER_LEN; // Pn + k
	uint Qnk = max_element(chain.begin(), chain.end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//prev(chain.end())->SP2 + KMER_LEN; // Qn + k;

	uint ovl_str1, ovl_str2, ovl_end1, ovl_end2;
	ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk + KMER_LEN;
/*	if (P1 > Q1 && len1 - Pnk <= len2 - Qnk)
	{
		ovl_str1 = P1 - Q1;// -1;
		ovl_end1 = len1 - 1;
		ovl_str2 = 0;
		ovl_end2 = Qnk + len1 - Pnk;// -1;
	}
	else if (P1 <= Q1 && len1 - Pnk <= len2 - Qnk)
	{
		ovl_str1 = 0;
		ovl_end1 = len1 - 1;
		ovl_str2 = Q1 - P1;// -1;
		ovl_end2 = Qnk + len1 - Pnk;// -1;
	}
	else if (P1 > Q1 && len1 - Pnk > len2 - Qnk)
	{
		ovl_str1 = P1 - Q1;
		ovl_end1 = Pnk + len2 - Qnk;// -1;
		ovl_str2 = 0;
		ovl_end2 = len2 - 1;
	}
	else if (P1 <= Q1 && len1 - Pnk >= len2 - Qnk)
	{
		ovl_str1 = 0;
		ovl_end1 = Pnk + len2 - Qnk;// - 1;
		ovl_str2 = Q1 - P1;// -1;
		ovl_end2 = len2 - 1;
	}
*/
	return { ovl_str1, ovl_str2, ovl_end1, ovl_end2 , chain.begin()->orient};
}

unique_ptr<cwd::hash<uint, alignInfo_t>> cwd::findSameKmer(kmerHashTable_t& kmerHashTable, Dna5String& seq)
{
	//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
	unique_ptr<cwd::hash<uint, alignInfo_t>> kmerSet( new cwd::hash<uint, alignInfo_t> ); //[length(seq)];
	auto& read1 = seq;
	kmer_t kmer1;
	for (auto i = begin(read1); i < end(read1) - KMER_LEN; i += KMER_STEP)
	{
		kmer1 = { i, i + KMER_LEN };
		auto rangeP = kmerHashTable.equal_range(kmer1);
		auto rangeN = kmerHashTable.equal_range(revComp(kmer1));
		auto range = rangeP;
		bool orient = true;
		if (distance(rangeP.first, rangeP.second) < distance(rangeN.first, rangeN.second))
		{
			orient = false;
			range = rangeN;
		}
		while (range.first != range.second)
		{
			uint readID = range.first->second.readID;
			uint startPos1 = distance(begin(read1), i);
			uint startPos2 = range.first->second.begin;
			if (readID != 0)
				kmerSet->insert({ readID, {orient, startPos1, startPos2 } });
			++range.first; // Increment begin iterator
		}
	}
	return kmerSet;
}

seqData_t* cwd::loadSeqData(const string& seqFileName, StringSet<CharString>& ID, seqData_t& seq)
{
	StringSet<CharString> id;
	//seqData_t* seq = new StringSet<Dna5String>();
	SeqFileIn seqFileIn(seqFileName.c_str());
	readRecords(id, seq, seqFileIn);
	cout << "Read seqFile\n";
	for (auto& str : id)
	{
		StringSet<CharString	> split;
		strSplit(split, str);
		erase(split[1], 0, 3);
		appendValue(ID, split[1]);
	}

	return &seq;
}

void cwd::filterKmer(kmerHashTable_t& kmerHashTable, const string& kfFileName)
{
	ifstream kmerFrequency(kfFileName);
	if (kmerFrequency.is_open())
	{
		kmer_t kmer;
		uint frequency;
		uint maxFrequency = maxKmerFrequency(kmerFrequency);
		cout << "Kmer Frequncy Range : [ " << 2 << " , " << maxFrequency << " ]\n";
		kmerFrequency.clear();
		kmerFrequency.seekg(0);
		while (kmerFrequency >> kmer >> frequency)
		{
			if (frequency < 2 || frequency > maxFrequency)
			{
				kmerHashTable.erase(kmer);
			}
		}
	}
}

void cwd::outputOverlapInfo(uint r, uint & i, std::__cxx11::list<alignInfo_t>& chain, seqData_t& seq, StringSet<CharString> & ID, ofstream& outFile)
{
	//cout << "seq: " << r << " and: " << i << " ";
	//cout << " has OVERLAP! " << endl;
	outFile << r << "," << i << ",";
	auto ovl = finalOverlap(chain, length(seq[r]), length(seq[i]));
	//cout << ovl.SP1 << " , ";
	//cout << ovl.EP1 << " , ";
	//cout << ovl.SP2 << " , ";
	//cout << ovl.EP2 << endl;
	outFile << ovl.orient << ",";
	outFile << ovl.SP1 << ",";
	outFile << ovl.EP1 << ",";
	outFile << ovl.SP2 << ",";
	outFile << ovl.EP2 << "," << length(seq[r]) << "," << length(seq[i]) << endl;
}

void cwd::mainProcess(cwd::kmerHashTable_t& kmerHashTable, seqData_t& seq, StringSet<CharString> & ID)
{
	// 取出表中的一行 ，放到新的表 commonKmerSet 中，然后再去除重复的 kmer
	ofstream outFile("result_sampled_.txt", ios_base::out);
	for (uint r = 0; r < length(seq); r++)
	{
		//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
		auto kmerSet = findSameKmer(kmerHashTable, seq[r]);
		for (uint i = r + 1; i < length(seq); i++)
		{
			if (ID[r] == "41869" && ID[i] == "225759")
				cout << r << ", " << i << endl;
			auto range = kmerSet->equal_range(i);
			if (range.first == range.second) continue;
			auto commonKmerSet = getCommonKmerSet(range, seq[r]);
			//cout << commonKmerSet.size() << endl;
			if (commonKmerSet->size() > 0)
			{
				auto chain = chainFromStart(*commonKmerSet, 15, 13, 200, 0.5, 0.2);
				if (chain->size() > 2)
				{
					outputOverlapInfo(r, i, *chain, seq, ID, outFile);
				}
			}
		}
		seqan::clear(seq[r]);
	}
}

kmer_t cwd::revComp(const kmer_t& kmer)
{
	kmer_t rvcKmer = kmer;
	for (auto& dNtp : rvcKmer)
	{
		switch (dNtp)
		{
		case 'A': dNtp = 'T'; break;
		case 'T': dNtp = 'A'; break;
		case 'C': dNtp = 'G'; break;
		case 'G': dNtp = 'C'; break;
		default:
			break;
		}
	}
	std::reverse(rvcKmer.begin(), rvcKmer.end())	;
	// cout << rvcKmer;
	return rvcKmer;
}
