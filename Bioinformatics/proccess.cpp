#include "proccess.h"
#include "assembler.hpp"
#include <boost/graph/copy.hpp>
#include "boost/graph/depth_first_search.hpp"
using namespace std;
using namespace seqan;
using namespace cwd;
void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient = true);


mutex fileMutex;
mutex overlapMutex;
uint KMER_STEP = 1;
const int KMER_LIMIT = 51;
const int OVL_TIP_LEN = 100;
const int CHAIN_LEN = 2;
const double DETECT_RATIO = 0.4;
vector<shared_ptr<list<assemblyInfo_t>>> assemblyChain;
shared_ptr<AGraph> assemblyGraph;
vector<assemblyInfo_t> overlap;
seqData_t assemblySeq;
shared_ptr<SubGraph> tmpGraph;
vector<ComponentGraph> comps;

kmerHashTable_t* cwd::createKmerHashTable(const seqData_t& seq, bool isFull)
{
	kmerHashTable_t* kmerHashTable = new kmerHashTable_t();
	uint readID = 0;
	std::default_random_engine dre;
	std::uniform_int_distribution<int> di(1, KMER_LIMIT);
	cout << "Creating HashTable...\n";
	for (auto& read : seq)
	{
		int dlen = length(read) * DETECT_RATIO;
		uint pos = 0;
		if (!isFull)
		{
			auto headEnd = begin(read) + dlen - KMER_LEN;
			auto tailBegin = end(read) - dlen - KMER_LEN;
			for (auto i = begin(read); i < headEnd; i += KMER_STEP, pos += KMER_STEP)
			{
				KMER_STEP = di(dre);
				kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
				hashValue_t kmerInfo{ readID, pos }; //, pos + KMER_LEN };
				kmerHashTable->insert({ kmer, kmerInfo });
			}

			uint pos2 = distance(begin(read), tailBegin);
			for (auto i = tailBegin; i < end(read) - KMER_LEN; i += KMER_STEP, pos2 += KMER_STEP)
			{
				KMER_STEP = di(dre);
				kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
				hashValue_t kmerInfo{ readID, pos2 }; //, pos + KMER_LEN };
				kmerHashTable->insert({ kmer, kmerInfo });
			}
		}
		else
		{
			for (auto i = begin(read); i < end(read); i += KMER_STEP, pos += KMER_STEP)
			{
				KMER_STEP = di(dre);
				kmer_t kmer = { i, i + KMER_LEN }; // *kmer = string(left:b, right:s)
				hashValue_t kmerInfo{ readID, pos }; //, pos + KMER_LEN };
				kmerHashTable->insert({ kmer, kmerInfo });
			}
		}
		readID++;
	}
	cout << "HashTable has been created!\n";
	return kmerHashTable;
}

vector<shared_ptr<list<alignInfo_t>>> cwd::chainFromStart(seqData_t& seq, vector<alignInfo_t>& cks, int k, int ks, int alpha, int beta, double gamma, int r, int t)
{
	shared_ptr<list<alignInfo_t>> chain = make_shared<list<alignInfo_t>>();
	vector<decltype(chain)> chain_v;
	// chain->push_back(CKS.begin()->second);
	sort(cks.begin(), cks.end(), [](alignInfo_t& a, alignInfo_t& b) { return a.SP1 < b.SP1; });
	chain->push_back(*cks.begin());
	for (auto ix = cks.begin(), nextx = next(ix); ix != cks.end() && nextx != cks.end();)
	{
		uint d1 = 0;
		uint d2 = 0;
		if (ix->SP2 < nextx->SP2 && ix->orient && nextx->orient || ix->SP2 > nextx->SP2 && !ix->orient && !nextx->orient)
		{
			if (ix->SP2 < nextx->SP2)
			{
				d1 = nextx->SP1 - ix->SP1;
				d2 = nextx->SP2 - ix->SP2;
			}
			else
			{
				d1 = nextx->SP1 - ix->SP1;
				d2 = ix->SP2 - nextx->SP2;
			}
			if ((d1 < alpha && d2 < alpha))// and (double(max(d1, d2) - min(d1, d2)) / max(d1, d2) < gamma))
			{
				//chain->push_back(*ix);
				chain->push_back(*nextx);
				//ix = next(nextx);
				//nextx = next(ix);
				ix = nextx;
				nextx++;
			}
			else if (d1 < beta && d2 < beta)
			{
				if ((double)d1 / d2 > 0.9 || (double)d2 / d1 > 0.9)
				{
					int s = min(ix->SP1, nextx->SP1), s2 = min(ix->SP2, nextx->SP2);
					int min_d = min(d1, d2);
					if (findSmallerSameKmer(seq, r, t, ks, s, s2, min_d,ix->orient && nextx->orient))
					{
						chain->push_back(*nextx);
						ix = nextx;
						nextx++;
					}
					else
					{
						if (chain->size() > CHAIN_LEN)
						{
							chain_v.push_back(chain);
							chain = make_shared<list<alignInfo_t>>();
							//chain->push_back(F_END);
							//iend = prev(chain->end());
							chain->push_back(*nextx);
						}
						else
							chain->erase(chain->begin(), chain->end());
						//chain->erase(next(iend), chain->end());
						ix = nextx;
						nextx++;
					}
				}
			}
			else
			{
				if (chain->size() > CHAIN_LEN)
				{
					chain_v.push_back(chain);
					chain = make_shared<list<alignInfo_t>>();
					//chain->push_back(F_END);
					//iend = prev(chain->end());
					chain->push_back(*nextx);
				}
				else
					//chain->erase(next(iend), chain->end());
					chain->erase(chain->begin(), chain->end());
				ix = nextx;
				nextx++;
				continue;
			}
		}
		else
		{
			if (chain->size() > CHAIN_LEN)
				chain_v.push_back(chain);
			chain = make_shared<list<alignInfo_t>>();
			//chain = decltype(chain)(new list<alignInfo_t>());
			//chain->push_back(F_END);//
			//iend = prev(chain->end());
			ix = nextx;
			nextx++;
			continue;
		}
	}
	//chain->push_back(F_END);
	if (chain->size() > CHAIN_LEN)
		chain_v.push_back(chain);
	return chain_v;
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

vector<assemblyInfo_t> cwd::finalOverlap(vector<shared_ptr<list<alignInfo_t>>>& chain_v, uint len1, uint len2, uint r, uint i, int chainLen, int ovLen)
{
	//auto r1 = chain.begin();
	//auto r2 = find_if(r1, chain.end(), [](alignInfo_t& a) { return a.orient == false && a.SP1 == numeric_limits<uint>::max() && a.SP2 == numeric_limits<uint>::max(); });
	vector<assemblyInfo_t> res;
	//auto ch = *max_element(chain_v.begin(), chain_v.end(), [](shared_ptr<list<alignInfo_t>>& a, shared_ptr<list<alignInfo_t>>& b) { return a->size() < b->size();});
	//auto ch = chain_v[0];
	for (auto& ch : chain_v)
	{
		if (ch->size() > chainLen)
		{
			uint P1 = ch->begin()->SP1; // P1
			uint Q1 = min_element(ch->begin(), ch->end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//chain.begin()->SP2; // Q1
			uint Pnk = ch->rbegin()->SP1 + KMER_LEN; // Pn + k
			uint Qnk = max_element(ch->begin(), ch->end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//prev(chain.end())->SP2 + KMER_LEN; // Qn + k;

			uint ovl_str1, ovl_str2, ovl_end1, ovl_end2;
			ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk + KMER_LEN;
			assemblyInfo_t a;
			a.r1 = r;
			a.r2 = i;
			a.SP1 = ovl_str1;
			a.SP2 = ovl_str2;
			a.EP1 = ovl_end1;
			a.EP2 = ovl_end2;
			a.orient = ch->begin()->orient;
			auto len = max(a.EP1 - a.SP1, a.EP2 - a.SP2);
			if (len > ovLen /*&& len < min(len1, len2) * 1*/)
			{
				res.push_back(a);
			}
			else
			{
				//res.push_back(a);
			}
		}
		//r1 = next(r2);
		//r2 = find_if(r1, chain_v.end(), [](alignInfo_t a) { return a.orient == false && a.SP1 == F_END.SP1 && a.SP2 == F_END.SP2; });
	}
	
	//uint P1 = chain.begin()->SP1; // P1
	//uint Q1 = min_element(chain.begin(), chain.end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//chain.begin()->SP2; // Q1
	//uint Pnk = prev(chain.end())->SP1 + KMER_LEN; // Pn + k
	//uint Qnk = max_element(chain.begin(), chain.end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//prev(chain.end())->SP2 + KMER_LEN; // Qn + k;

	//uint ovl_str1, ovl_str2, ovl_end1, ovl_end2;
	//ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk + KMER_LEN;
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
	return res;
}

unique_ptr<cwd::hash<uint, alignInfo_t>> cwd::findSameKmer(kmerHashTable_t& kmerHashTable, seqData_t & seq, uint r)
{
	//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
	auto kmerSet = make_unique<cwd::hash<uint, alignInfo_t>>(); //[length(seq)];
	auto& read1 = seq[r];
	kmer_t kmer1;
	for (auto i = begin(read1); i < end(read1) - KMER_LEN; i++)
	{
		kmer1 = { i, i + KMER_LEN };
		auto rangeP = kmerHashTable.equal_range(kmer1);
		auto rangeN = kmerHashTable.equal_range(revComp(kmer1));
		kmer_t rkmer1 = revComp(kmer1);
		auto range = rangeP;
		bool orient = true;
		int d = distance(range.first, range.second);
		//{
//			std::default_random_engine dre;
//			std::uniform_int_distribution<int> di(0, 30);
//			int x = di(dre);
//			kmer1[x] = revComp(kmer1.c_str() + x)[0];
//			rangeP = kmerHashTable.equal_range(kmer1);
//			rangeN = kmerHashTable.equal_range(revComp(kmer1));
//			range = rangeP;
		//}
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
			//kmer_t search = { begin(seq[readID]) + startPos2, begin(seq[readID]) + startPos2 + KMER_LEN };
			//if (search == kmer1 || search == rkmer1)
			kmerSet->insert({ readID, {orient, startPos1, startPos2 } });
			//else
			//{
			//	kmerSet->insert({ readID, {orient, startPos1, startPos2 } });
			//}
			++range.first; // Increment begin iterator
		}
	}
	return kmerSet;
}

void cwd::loadSeqData(const string& seqFileName, StringSet<CharString>& ID, seqData_t& seq)
{
	StringSet<CharString> id;
	SeqFileIn seqFileIn(seqFileName.c_str());
	cout << "Reading seqFile...\n";
	readRecords(id, seq, seqFileIn);
	cout << "seqFile has been read.\n";
	//ofstream dict("dict3.txt", ios_base::out);
	//int i = 0;
	//for (auto& str : id)
	//{
	//	StringSet<CharString> split;
	//	strSplit(split, str);
	//	//erase(split[1], 0, 3);
	//	//appendValue(ID, split[1]);
	//	dict << split[0] << " " << i++ << endl;
	//}
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

void cwd::outputOverlapInfo(uint r, uint i, vector<shared_ptr<list<alignInfo_t>>>& chain_v, seqData_t& seq, StringSet<CharString> & ID, ofstream& outFile, int minSize, int chainLen, int ovLen)
{
	auto v_ovl = finalOverlap(chain_v, length(seq[r]), length(seq[i]), r, i, chainLen, ovLen);

	//for (auto& ovl : v_ovl)
	//{
	//	if (ovl.EP1 - ovl.SP1 > 600 && ovl.EP2 - ovl.SP2 > 600)
	//	{
	//		outFile << boost::format("%u, %u, %u, %u, %u, %u, %u, %u, %u\n")
	//			% r % i % ovl.orient % ovl.SP1 % ovl.EP1 % ovl.SP2 % ovl.EP2 % length(seq[r]) % length(seq[i]);
	//	}
	//}
	
	//vector<assemblyInfo_t> v_ass;
	//transform(v_ovl.begin(), v_ovl.end(), back_inserter(v_ass),
	//		[=](overlapInfo_t& a) {
	//			return assemblyInfo_t{ r, i, a.SP1, a.EP1, a.SP2, a.EP2, a.orient };
	//		});
	
	overlap.insert(overlap.end(), v_ovl.begin(), v_ovl.end());
}

void cwd::mainProcess(cwd::kmerHashTable_t& kmerHashTable, seqData_t& seq, StringSet<CharString> & ID, int block1, int block2, ofstream& outFile, int chainLen, int ovLen)
{
	// 取出表中的一行 ，放到新的表 commonKmerSet 中，然后再去除重复的 kmer
	for (uint r = block1; r < block2; r++)
	{
		//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
		auto kmerSet = findSameKmer(kmerHashTable, seq, r);
		for (uint i = r + 1; i < length(seq); i++)
		{
			//if (r == 6423 && i == 6859)
			//	cout << string{ begin(seq[r]), begin(seq[r]) + 5347 } << "\n" << string{ begin(seq[r]), begin(seq[r]) + 5347 } << endl;
			//else continue;
			auto range = kmerSet->equal_range(i);
			if (range.first == range.second)
			{
				continue;
			}
			else
			{
				auto commonKmerSet = getCommonKmerSet(range, seq[r]);
				//cout << commonKmerSet.size() << endl;
				if (commonKmerSet.size() > 0)
				{
					auto chain_v = chainFromStart(seq, commonKmerSet, KMER_LEN, 15, 300, 500, 0.2, r, i);
					if (chain_v.size() > 0)
					{
						lock_guard<mutex> lock(fileMutex);
						outputOverlapInfo(r, i, chain_v, seq, ID, outFile, 600, chainLen, ovLen);
					}
				}
			}
		}
		//seqan::clear(seq[r]);
	}
	cout << boost::format("mapped %d - %d reads.\n") % block1 % block2;
	//std::sort(info.begin(), info.end(), 
	//	[](assemblyInfo_t& a, assemblyInfo_t& b)
	//	{
	//		if (a.r1 < b.r1)
	//		{
	//			return true;
	//		}
	//		else if (a.r1 == b.r1)
	//		{
	//			return (a.r2 < b.r2);
	//		}
	//		else
	//		{
	//			return false;
	//		}
	//	});
	//assembler();
}

void cwd::assembler(const seqData_t& seq, ofstream & seqOut)
{
	//ofstream outAssembly("toyAssembly.csv");
	//for (auto& chain : assemblyChain)
	//{
	//	if (chain->size() < 3)
	//	{
	//		continue;
	//	}
	//	for (auto& ovl : *chain)
	//	{
	//		outAssembly << boost::format("%u, %u, %u, %u, %u, %u\n") 
	//			% ovl.r1 % ovl.r2 % ovl.SP1 % ovl.EP1 % ovl.SP2 % ovl.EP2;
	//		
	//	}
	//	//cout << "--------------------------------------------------------------\n";
	//	outAssembly << endl;
	//}
	//assemblyChain.clear();
	//outAssembly.close();

	ofstream outAssembly, outPath;
	outPath.open("toyAssembly_path.txt");
	outAssembly.open("toyAssembly_graph.txt");
	int i = 0;
	for (auto& aGraph : comps)
	{
	//	boost::dynamic_properties dp;
	//	dp.property("node_id", boost::get(&AVertex::r, *aGraph));
	//	dp.property("node_id", boost::get(&AVertex::r, *aGraph));
	//	dp.property("label", boost::get(vertex_property_t(), *aGraph));
	//	dp.property("label", boost::get(edge_property_t(), *aGraph));
	//	boost::write_graphviz(outAssembly, *aGraph, dp);
		AGraph g;
		copy_graph(aGraph, g);
		boost::write_graphviz(outAssembly, g, boost::make_label_writer(boost::get(vertex_property_t(), aGraph)),
		make_edge_writer(boost::get(&AEdge::adj, g)));
		outAssembly << endl;
		outPath << "Graph: " << i++;
		auto ps = findPath(g, seq, outPath);
		//std::set<int> sp(p.begin(), p.end());
		for (auto& p : ps)
		{
			generateContig(g, p, seq, seqOut, i++);
		}
		//set_union(un.begin(), un.end(), sp.begin(), sp.end(), inserter(un, un.begin()));
	}
	comps.clear();
	overlap.clear();
	outAssembly.close();
	outPath.close();
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
		default: dNtp = 'n';break;
		}
	}
	std::reverse(rvcKmer.begin(), rvcKmer.end());
	// cout << rvcKmer;
	return rvcKmer;
}

bool cwd::findSmallerSameKmer(seqData_t& seq, uint r, uint t, uint SKMER_LEN, int s, int s2, int d, bool orient)
{
	//unique_ptr<cwd::hash<uint, alignInfo_t>> kmerSet(new cwd::hash<uint, alignInfo_t>);
	auto& read1 = seq[r];
	auto& read2 = seq[t];
	
	string gap1 = { begin(read1) + s, begin(read1) + s + d };
	string gap2 = { begin(read2) + s2, begin(read2) + s2 + d };
	/*kmer_t kmer1, kmer2;
	int count = 0;
	for (auto i = begin(read1) + s; i < begin(read1) + e - SKMER_LEN; i++)
	{
		for (auto j = begin(read2) + s; j < begin(read2) + e - SKMER_LEN;)
		{
			int a = distance(begin(read2), j), b = distance(begin(read1), i);
			if (a > b)
			{
				t = a - b;
				if (t >= SKMER_LEN)
				{
					break;
				}
			}
			else
			{
				t = b - a;
				if (t >= SKMER_LEN)
				{
					j += t;
				}
			}
			kmer1 = { i, i + SKMER_LEN };
			kmer2 = { j, j + SKMER_LEN };
			if (kmer1 == kmer2)
			{
				count++;
				j += SKMER_LEN;
				i += SKMER_LEN;
			}
			else
			{
				j++;
			}
		}
	}
	if (count > 1)
	{
		float score = float(count) / (float(e - s ) / SKMER_LEN);
		if (score > 0.5)
		{
			return true;
		}
	}*/
	return hamming(gap1, gap2) < 0.5f;
}

std::string cwd::getCurrentDate()
{
	time_t nowtime;
	nowtime = time(NULL); //获取日历时间   
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d-%H-%M", localtime(&nowtime));
	return tmp;
}

float cwd::jaccard(string& a, string& b)
{
	if (a == "" && b == "")
	{
		return 1.0f;
	}
	// 都为空相似度为 1
	if (a == "" || b == "")
	{
		return 0.0f;
	}
	std::set<int> aChar(a.begin(), a.end());
	std::set<int> bChar(b.begin(), b.end());
	std::set<int> in;
	std::set<int> un;
	set_intersection(aChar.begin(), aChar.end(), bChar.begin(), bChar.end(), inserter(in, in.begin()));
	set_union(aChar.begin(), aChar.end(), bChar.begin(), bChar.end(), inserter(un, un.begin()));
	int inter = in.size(); // 交集数量
	int uni = un.size(); // 并集数量
	if (inter == 0) return 0;
	return ((float)inter) / (float)uni;
}

float cwd::hamming(string& a, string& b)
{
	if (a == "" || b == "")
	{
		return 0.0f;
	}
	if (a.length() != b.length())
	{
		return 0.0f;
	}

	int disCount = 0;
	for (int i = 0; i < a.length(); i++)
	{
		if (a.at(i) != b.at(i))
		{
			disCount++;
		}
	}
	return (float)disCount / (float)a.length();
}

void cwd::toyAssembly(seqData_t& seq, int block1, int block2)
{
	using vertex_descriptor = boost::graph_traits<AGraph>::vertex_descriptor;
	vertex_descriptor src, dst;
	uint idx = 0;
	assemblyGraph = make_shared<AGraph>();
	tmpGraph = make_shared<SubGraph>();
	auto e_prop = boost::get(&AEdge::weight, *assemblyGraph);
	for (int i = block1; i < block2; i++)
	{
		auto& ovl = overlap[i];
		if (true)
		{
			//int j = length(seq[i]);
			bool flag = false;
			uint r = ovl.r1;
			uint i = ovl.r2;
			uint len_r = length(seq[r]);
			uint len_i = length(seq[i]);
			double weight = -1.0 * (len_r + len_i - ovl.EP1 + ovl.SP1);
			bool isHeadTail = false, isTailHead = false, isHeadHead = false, isTailTail = false;
			if (
				(isTailHead = len_r - ovl.EP1 < OVL_TIP_LEN && ovl.SP2 < OVL_TIP_LEN && ovl.orient) ||
				(isHeadTail = len_i - ovl.EP2 < OVL_TIP_LEN && ovl.SP1 < OVL_TIP_LEN && ovl.orient) ||
				//(isTailHead = len_r - ovl.EP1 < OVL_TIP_LEN && len_i - ovl.EP2 < OVL_TIP_LEN && !ovl.orient) ||
				//(isHeadTail = ovl.SP2 < OVL_TIP_LEN && ovl.SP1 < OVL_TIP_LEN && !ovl.orient) ||
				(isHeadHead = ovl.SP1 < OVL_TIP_LEN && ovl.SP2 < OVL_TIP_LEN && !ovl.orient) ||
				(isTailTail = len_r - ovl.EP1 < OVL_TIP_LEN && len_i - ovl.EP2 < OVL_TIP_LEN && !ovl.orient)
				) //TODO direction
			{
				if (0)
				{
					for (auto& chain : assemblyChain)
					{
						auto ix = find_if(chain->begin(), chain->end(),
							[=](assemblyInfo_t& a) {
							return a.r1 == r || a.r1 == i || a.r2 == r || a.r2 == i;
						});
						if (ix != chain->end())
						{
							//chain->push_back(ovl);
							flag = true;
							chain->insert(next(ix), ovl);
							//break;
						}
						//if (chain->rbegin()->r1 == r || chain->rbegin()->r1 == i)
						//{
						//	chain->push_back({ r, i, ovl.SP1, ovl.EP1, ovl.SP2, ovl.EP2, true });
						//	flag = true;
						//	break;
						//}
						//else if (chain->rbegin()->r2 == r || chain->rbegin()->r2 == i)
						//{
						//	chain->push_back({ i, r, ovl.SP1, ovl.EP1, ovl.SP2, ovl.EP2, true });
						//	flag = true;
						//	break;
						//}
					}
					if (!flag)
					{
						auto n_chain = make_shared<list<assemblyInfo_t>>();
						n_chain->push_back(ovl);
						assemblyChain.push_back(n_chain);
					}
				}
				
				auto& aGraph = assemblyGraph;
				boost::graph_traits<AGraph>::vertex_iterator vi, vi_end;
				boost::tie(vi, vi_end) = boost::vertices(*aGraph);
				auto prop = boost::get(vertex_property_t(), *aGraph);
				auto vx = find_if(vi, vi_end,
					[=](vertex_descriptor ix) {
						AVertex y = prop[ix];
						return y.r == r || y.r == i;
				});
				if (vx != vi_end)
				{
					flag = true;
					AVertex v = { ovl.r1, ovl.SP1, ovl.EP1 };
					AVertex v2 = { ovl.r2, ovl.SP2, ovl.EP2 };
					if (prop[*vx].r == i)
					{
						auto vx2 = find_if(vi, vi_end,
							[=](vertex_descriptor ix) {
							AVertex y = prop[ix];
							return y.r == r;
						});
						if (vx2 != vi_end)
						{
							auto e = boost::edge(*vx, *vx2, *aGraph);
							if (e.first.m_eproperty)
							{
								auto edg = e.first;
								if (e_prop[edg] > weight)
								{
									e_prop[edg] = weight;
								}
							}

							else
							{
								if (isTailHead)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								}
								else if (isHeadTail)
								{
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								}
								else if (isHeadHead)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								}
								else if (isTailTail)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								}
								else
								{
									cout << "none\n";
								}

							}
							//add to an overlap vector
						}
						else
						{
							src = boost::add_vertex(v, *aGraph);
							if (isHeadTail)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								// *vx -> src
							}
							else if (isHeadHead)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								boost::add_edge(src, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
							}
							else if (isTailHead)
							{
								boost::add_edge(src, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
							}
							else if (isTailTail)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								boost::add_edge(src, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
							}
							else
							{
								cout << "none\n";
							}
							//TODO: condition of direction
						}
					}
					else
					{
						auto vx2 = find_if(vi, vi_end,
							[=](vertex_descriptor ix) {
							AVertex y = prop[ix];
							return y.r == i;
						});
						if (vx2 != vi_end)
						{
							auto e = boost::edge(*vx, *vx2, *aGraph);
							if (e.first.m_eproperty)
							{
								auto edg = e.first;
								if (e_prop[edg] > weight)
								{
									e_prop[edg] = weight;
								}
							}

							else
							{
								if (isTailHead)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								}
								else if (isHeadTail)
								{
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								}
								else if (isHeadHead)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								}
								else if (isTailTail)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								}
								else
								{
									cout << "none\n";
								}
							}
							//add to an overlap vector
						}
						else
						{
							dst = boost::add_vertex(v2, *aGraph);
							if (isTailHead)
							{
								boost::add_edge(*vx, dst, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								// *vx -> dst
							}
							else if (isHeadTail)
							{
								boost::add_edge(dst, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
							}
							else if (isHeadHead)
							{
								boost::add_edge(*vx, dst, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								boost::add_edge(dst, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
							}
							else if (isTailTail)
							{
								boost::add_edge(*vx, dst, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								boost::add_edge(dst, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
							}
							else
							{
								cout << "none\n";
							}
							//TODO: condition of direction
						}
					}
				}

				if (!flag)
				{
					auto& aGraph = assemblyGraph;
					AVertex v = { ovl.r1, ovl.SP1, ovl.EP1 };
					AVertex v2 = { ovl.r2, ovl.SP2, ovl.EP2 };
					src = boost::add_vertex(v, *aGraph);
					dst = boost::add_vertex(v2, *aGraph);
					if (isTailHead)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
					}
					else if (isHeadTail)
					{
						boost::add_edge(dst, src, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
					}
					else if (isHeadHead)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
						boost::add_edge(dst, src, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
					}
					else if (isTailTail)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
						boost::add_edge(dst, src, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
					}
					else
					{
						cout << "none\n";
					}
				}
			}
		}
	}
	connected_components_subgraphs(*assemblyGraph);
}

void cwd::connected_components_subgraphs(AGraph const& g)
{
	using namespace boost;
	boost::copy_graph(g, *tmpGraph);
	std::shared_ptr<vector<unsigned long>> mapping = std::make_shared<vector<unsigned long>>(num_vertices(g));
	size_t num = connected_components(*tmpGraph, mapping->data());

	auto& component_graphs = comps;

	for (size_t i = 0; i < num; i++)
		component_graphs.emplace_back(g,
			[mapping, i, &g](AGraph::edge_descriptor e) {
		return mapping->at(source(e, g)) == i
			|| mapping->at(target(e, g)) == i;
	},
			[mapping, i](AGraph::vertex_descriptor v) {
		return mapping->at(v) == i;
	});

	//return component_graphs;
}

bool cwd::isConnected(AGraph& g, vertex_descriptor a, vertex_descriptor b)
{
	//vector<bool> visited(boost::num_vertices(g));
	SubGraph gs;
	boost::copy_graph(g, gs);
	std::vector<boost::default_color_type> color_map(boost::num_vertices(gs));
	auto color = boost::make_iterator_property_map(color_map.begin(), boost::get(boost::vertex_index_t(), gs));
	boost::depth_first_search(gs, boost::default_dfs_visitor(), color, a);
	//auto xx = color[b];

	//DFS(g, a, visited);
	return color[b] != boost::default_color_type::white_color;
}

void cwd::DFS(cwd::AGraph& g, vertex_descriptor i, vector<bool>& visited)
{
	AGraph::out_edge_iterator ei, eend;
	auto ep = boost::out_edges(i, g);
	ei = ep.first;
	eend = ep.second;
	for (auto e = ei; e != eend; e++)
	{
		auto v = boost::target(*e, g);
		visited[v] = true;
		DFS(g, v, visited);
	}
}

std::set<size_t> cwd::finalOverlap2(vector<shared_ptr<list<alignInfo_t>>>& chain_v, uint len1, uint len2, uint r, uint i, int chainLen, int ovLen)
{
	std::set<size_t> res;
	for (auto& ch : chain_v)
	{
		if (ch->size() > chainLen)
		{
			uint P1 = ch->begin()->SP1; // P1
			uint Q1 = min_element(ch->begin(), ch->end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//chain.begin()->SP2; // Q1
			uint Pnk = ch->rbegin()->SP1 + KMER_LEN; // Pn + k
			uint Qnk = max_element(ch->begin(), ch->end(), [](alignInfo_t& a, alignInfo_t& b) {return a.SP2 < b.SP2;})->SP2;//prev(chain.end())->SP2 + KMER_LEN; // Qn + k;

			uint ovl_str1, ovl_str2, ovl_end1, ovl_end2;
			ovl_str1 = P1, ovl_str2 = Q1, ovl_end1 = Pnk, ovl_end2 = Qnk + KMER_LEN;
			assemblyInfo_t a;
			a.r1 = r;
			a.r2 = i;
			a.SP1 = ovl_str1;
			a.SP2 = ovl_str2;
			a.EP1 = ovl_end1;
			a.EP2 = ovl_end2;
			a.orient = ch->begin()->orient;
			auto len = max(a.EP1 - a.SP1, a.EP2 - a.SP2);
			if (len >= min(len1, len2) * 0.95)
			{
				res.insert(r);
				res.insert(i);
			}
		}
	}

	return res;
}

void cwd::mainProcess2(cwd::kmerHashTable_t& kmerHashTable, seqData_t& seq, StringSet<CharString>& ID, int block1, int block2, ofstream& outFile, int chainLen, int ovLen, std::set<size_t> & dump)
{
	// 取出表中的一行 ，放到新的表 commonKmerSet 中，然后再去除重复的 kmer
	for (uint r = block1; r < block2; r++)
	{
		//每一个读数一个表，用 ReadID 作为索引，记录 readx 与 readID 之间的相同的 kmer
		auto kmerSet = findSameKmer(kmerHashTable, seq, r);
		for (uint i = r + 1; i < length(seq); i++)
		{
			auto range = kmerSet->equal_range(i);
			if (range.first == range.second)
			{
				continue;
			}
			else
			{
				auto commonKmerSet = getCommonKmerSet(range, seq[r]);
				//cout << commonKmerSet.size() << endl;
				if (commonKmerSet.size() > 0)
				{
					auto chain_v = chainFromStart(seq, commonKmerSet, KMER_LEN, 15, 300, 500, 0.2, r, i);
					if (chain_v.size() > 0)
					{
						lock_guard<mutex> lock(fileMutex);
						auto v_ovl = finalOverlap2(chain_v, length(seq[r]), length(seq[i]), r, i, chainLen, ovLen);
						set_union(dump.begin(), dump.end(), v_ovl.begin(), v_ovl.end(), inserter(dump, dump.begin()));
					}
				}
			}
		}
	}
	cout << boost::format("mapped %d - %d reads.\n") % block1 % block2;
}
