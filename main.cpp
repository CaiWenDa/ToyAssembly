﻿#include "assembly.h"
#include "utility.h"
#include "memory_info.h"
#include <cstdlib>
#include <malloc.h>
#include <getopt.h>
#include <iostream>
#include <sys/stat.h>

extern int optind, opterr, optopt;
extern char* optargi;

void usage()
{
	std::cerr <<
		"usage: alphaHiASM \n\
		--file - file_1[file_2 ...]\n\
		--output - dir PATH\n\
		--genomeSize SIZE\n\
		--threads int\n\
		[--minOverlap SIZE]\n\
		[--help]\n\
		[--paf inOverlapFile]\n\
		[--overlap outOverlapFile]\n\
		[--kmerLen kmer len]\n\
		[--step detect step]\n";

}

bool parseOption(int argc, char* argv[],
	std::string& seqFileName,
	std::string& asmFileName,
	int& genomeSize,
	int& thread_i,
	int& minOverlapLen,
	std::string& paf,
	std::string& overlap,
	uint& KMER_LEN,
	uint& KMER_STEP)
{
	int index = 0;
	int c = 0; //用于接收选项
	//定义长选项
	static struct option long_options[] =
	{
		{"help", no_argument, NULL, 'h'},
		{"file", required_argument, NULL, 'f'},
		{"output", required_argument, NULL, 'o'},
		{"genomeSize", required_argument, NULL, 'g'},
		{"minOverlap", required_argument, NULL, 0},
		{"theads", required_argument, NULL, 't'},
		{"paf", required_argument, NULL, 0},
		{"overlap", required_argument, NULL, 0},
		{"step", required_argument, NULL, 0},
		{"kmerLen", required_argument, NULL, 'k'}
	};

	/*循环处理参数*/
	while (EOF != (c = getopt_long(argc, argv, "hf:o:g:t:k:", long_options, &index)))
	{
		using std::cerr;
		using std::endl;
		switch (c)
		{
		case 'h':
			usage();
			break;
		case 'f':
			seqFileName = optarg;
			break;
			//-n选项必须要参数
		case 'o':
			asmFileName = optarg;
			break;
		case 'g':
			genomeSize = cwd::parseGenomeSize(optarg);
			break;
		case 't':
			thread_i = atoi(optarg);
			break;
		case 'k':
			KMER_LEN = atoi(optarg);
			break;
			//表示选项不支持
		case '?':
			cerr << "unknown option: " << static_cast<char>(optopt) << endl;
			break;
		case 0:
			if (!strcmp(long_options[index].name, "minOverlap"))
				minOverlapLen = atoi(optarg);
			if (!strcmp(long_options[index].name, "paf"))
				paf = optarg;
			if (!strcmp(long_options[index].name, "overlap"))
				overlap = optarg;
			if (!strcmp(long_options[index].name, "step"))
				KMER_STEP = atoi(optarg);
		default:
			break;
		}
	}

	if (argc <= 4)
	{
		usage();
		return false;
	}
	if (seqFileName.empty() || asmFileName.empty() || genomeSize == 0)
	{
		usage();
		return false;
	}

	return true;
}

void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient = true);

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace seqan;
	using namespace cwd;
	int minOverlapLen = 2000;
	extern int thread_i;
	extern int genomeSize;
	extern uint KMER_STEP;
	extern uint KMER_LEN;
	//string seqFileName = "/home/caiwenda/dmel.trimmedReads_20x.fasta";
	string seqFileName;//argv[2];//"/home/caiwenda/SRR11292120_sample.fastq";
	string asmFileName;//argv[3];//"result_chm13_debug713.fasta";
	//string asmFileName = "tmp.fasta";
	string outFilePre = getCurrentDate();
	//string seqFileName = "result_ecoli_debug.fasta";
	string kfFileName = "/home/caiwenda/ecoli_kmerFrequency.txt";
	//string seqFileName = "dmel.trimmedReads_10x.fasta";
	// string kfFileName = "/publicdata/Reads/HiFi/D.mel/kmer31.txt";
	string paf;
	string overlapFile;
	string graphFileName;// = "~/toyAssembly_graph.txt";
	if (!parseOption(argc, argv, seqFileName, asmFileName, genomeSize,
		thread_i, minOverlapLen, paf, overlapFile, KMER_LEN, KMER_STEP))
	{
		return 1;
	}
	cerr << "SYSTEM INFORMATION: \n";
	cerr << "Total RAM: "
		<< getMemorySize() / 1024 / 1024 / 1024 << " Gb\n";
	cerr << "Available RAM: "
		<< getFreeMemorySize() / 1024 / 1024 / 1024 << " Gb\n";
	cerr << "Total CPUs: " << std::thread::hardware_concurrency() << endl;
	cerr << endl;
	cerr << "PARAMETERS: \n";
	cerr << "seqFile: " << seqFileName << endl;
	cerr << "genomeSize: " << genomeSize << endl;
	cerr << "minOverlap: " << minOverlapLen << endl;
	cerr << "kmerSize: " << KMER_LEN << endl;
	cerr << "kmerStep: " << KMER_STEP << endl;
	cerr << "thread(s): " << thread_i << endl;
	cerr << endl;
	cerr << getCurrentTime() << "\t>>>STAGE: Reading file(s)\n";
	// cout << "frequencyFile : " << kfFileName << endl;
	seqData_t seq;
	StringSet<CharString> ID;
	clock_t start = clock();
	//ofstream outFile("result-" + getCurrentDate() + ".csv", ios_base::out);
	ofstream outFile(overlapFile, ios_base::out);
	if (!outFile.is_open())
	{
		cerr << "Cannot open file: " << overlapFile << endl;
		return 1;
	}
	ofstream seqOut(asmFileName, ios_base::out);
	if (!seqOut.is_open())
	{
		cerr << "Cannot open file: " << asmFileName << endl;
		return 1;
	}
	try
	{
		struct stat info;
		stat(seqFileName.c_str(), &info);
		auto size = info.st_size;
		loadSeqData(seqFileName, ID, seq);
		for (int i = optind; i < argc; i++)
		{
			loadSeqData(argv[i], ID, seq);
			stat(argv[i], &info);
			size += info.st_size;
		}
		cerr << "Reads: " << seqan::length(seq) << endl;
		cerr << "Estimate reads coverage: " << float(size) / genomeSize << endl;
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		getchar();
		return 1;
	}
	seqan::clear(ID);
	malloc_trim(0);

	if (!paf.empty())
	{
		cerr << "Reading PAF file...\n";
		try
		{
			readPAF(paf, minOverlapLen);
		}
		catch (exception& e)
		{
			cerr << e.what() << endl;
			getchar();
			return 1;
		}
	}
	else
	{
		cerr << getCurrentTime() << "\t>>>STAGE: Hashing\n";
		auto kmerHashTable = createKmerHashTable(seq, true);
		//filterKmer(*kmerHashTable, kfFileName);
		uint block1 = 0;
		uint block2 = 0;
		uint b_size = length(seq) / thread_i;
		//ofstream seqOut;
		vector<thread> threadPool;
		cerr << getCurrentTime() << "\t>>>STAGE: Detecting Overlap...\n";
#if 1
		for (size_t i = 0; i < thread_i; i++)
		{
			block2 += b_size;
			threadPool.emplace_back(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, minOverlapLen);
			block1 = block2;
		}
		if (length(seq) % thread_i != 0)
		{
			threadPool.emplace_back(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, minOverlapLen);
		}
		for (auto& th : threadPool)
		{
			th.join();
		}
		kmerHashTable->clear();
		delete kmerHashTable;
		threadPool.clear();
		malloc_trim(0);

	}
	//boost::thread_group tg;
	//tg.create_thread(bind(createOverlapGraph, ref(seq), block1, block2));
	cerr << getCurrentTime() << "\t>>>STAGE: Creating Overlap Graph...\n";
	extern vector<assemblyInfo_t> overlap;
	if (0 and !graphFileName.empty())
	{
		cerr << "Reading Graph file...\n";
		vector<AGraph> v_g;
		readOverlapGraph(graphFileName, v_g);
	}
	else
	{
		createOverlapGraph(seq, 0, overlap.size());
	}
	cerr << getCurrentTime() << "\t>>>STAGE: Assembling reads...\n";
	//clear(seq);
	//loadSeqData(seqFileName, ID, seq);
	assembler(seq, seqOut);

	if (false)
	{
		std::set<size_t> all;
		for (size_t i = 0; i < length(seq); i++)
		{
			all.insert(i);
		}
		extern std::set<size_t> delReads;
		std::set<size_t> una;
		set_difference(all.begin(), all.end(), delReads.begin(), delReads.end(), inserter(una, una.end()));
		if (!una.empty())
		{
			for (auto& i : una)
			{
				auto assembly = seq[i];
				seqOut << boost::format(">contig_%d; length=%d; reads=%d; type=dna\n") % (length(seq) + i) % length(assembly) % 1;
				auto a = begin(assembly);
				for (; a < prev(end(assembly), 1000); std::advance(a, 1000))
				{
					seqOut << string{ a, next(a, 1000) } << endl;
				}
				seqOut << string{ a, end(assembly) } << endl;
				seqOut << endl;
			}
		}
	}

	//ifstream overlapFile;
	//assemblyInfo_t ovl;
	//uint len1, len2;
	//while (overlapFile >> ovl.r1 >> ovl.r2 >> ovl.orient >> ovl.SP1 >> ovl.EP1 >> ovl.SP2 >> ovl.EP2 >> len1 >> len2)
	//{
	//	overlap.emplace_back(ovl);
	//}
#endif

#if 0
	seqOut.close();
	seqOut.open("result_ecoli_debug3.fasta");
	extern seqData_t assemblySeq;
	extern int OVL_TIP_LEN;
	OVL_TIP_LEN = 200000;
	ovLen = 1000;
	while (length(assemblySeq) > 1)
	{
		seq = assemblySeq;
		clear(assemblySeq);
		kmerHashTable = createKmerHashTable(seq);
		block1 = 0;
		block2 = 0;
		b_size = length(seq) / thread_i;
		//outFile.open("result-" + getCurrentDate() + ".csv", ios_base::out);
		//ofstream outFile;
		cout << "Detecting Overlap...\n";
		for (size_t i = 0; i < thread_i; i++)
		{
			block2 += b_size;
			threadPool.push_back(thread(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, ovLen));
			block1 = block2;
		}
		if (length(seq) % thread_i != 0)
		{
			threadPool.push_back(thread(mainProcess,
				ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, ovLen));
		}
		for (auto& th : threadPool)
		{
			th.join();
		}
		kmerHashTable->clear();
		threadPool.clear();
		malloc_trim(0);
		//boost::thread_group tg;
		//tg.create_thread(bind(createOverlapGraph, ref(seq), block1, block2));
		cout << "Creating Overlap Graph...\n";
		extern vector<assemblyInfo_t> overlap;
		createOverlapGraph(seq, 0, overlap.size());
		cout << "Assembling reads...\n";
		assembler(seq, seqOut);
	}
#endif

#if 0
	seqOut.close();
	seqan::clear(seq);
	loadSeqData(asmFileName, ID, seq);
	seqan::clear(ID);
	kmerHashTable = createKmerHashTable(seq, false);
	block1 = 0;
	block2 = 0;
	b_size = length(seq) / thread_i;
	std::set<size_t> dump;
	threadPool.clear();
	cout << "Detecting Overlap...\n";
	for (size_t i = 0; i < thread_i; i++)
	{
		block2 += b_size;
		threadPool.push_back(thread(mainProcess2,
			ref(*kmerHashTable), ref(seq), ref(ID), block1, block2, ref(outFile), 2, 600, ref(dump)));
		block1 = block2;
	}
	if (length(seq) % thread_i != 0)
	{
		threadPool.push_back(thread(mainProcess2,
			ref(*kmerHashTable), ref(seq), ref(ID), block1, length(seq), ref(outFile), 2, 600, ref(dump)));
	}
	for (auto& th : threadPool)
	{
		th.join();
	}
	kmerHashTable->clear();
	delete kmerHashTable;
	threadPool.clear();
	malloc_trim(0);

	seqOut.open("remove_dup.fasta");
	std::set<size_t> all;
	for (size_t i = 0; i < length(seq); i++)
	{
		all.insert(i);
	}
	std::set<size_t> una;
	set_difference(all.begin(), all.end(), dump.begin(), dump.end(), inserter(una, una.end()));
	for (auto& i : una)
	{
		seqOut << boost::format(">contig_%d length=%d; reads=%d; type=dna\n") % i % length(seq[i]) % 1;
		auto a = begin(seq[i]);
		for (; a < prev(end(seq[i]), 1000); std::advance(a, 1000))
		{
			seqOut << string{ a, next(a, 1000) } << endl;
		}
		seqOut << string{ a, end(seq[i]) } << endl;
		seqOut << endl;
	}
#endif

	cerr << getCurrentTime() << "\tDone!\n";
	cerr << "CPU time: " << (clock() - start) / CLOCKS_PER_SEC << " sec(s)\n";
	auto peakMemByte = getPeakRSS();
	if (peakMemByte >= 1024 * 1024 * 1024)
	{
		cerr << "Peak RAM usage: " << peakMemByte / 1024 / 1024 / 1024 << " Gb\n";
	}
	else if (peakMemByte >= 1024 * 1024)
	{
		cerr << "Peak RAM usage: " << peakMemByte / 1024 / 1024 << " Mb\n";
	}
	else
	{
		cerr << "Peak RAM usage: " << peakMemByte / 1024 << " Kb\n";
	}
	outFile.close();
	seqOut.close();
	getchar();
	return 0;
}

void compSeqInRange(cwd::seqData_t& seq, uint r1, uint r2, uint start1, uint start2, uint end1, uint end2, uint len, bool orient)
{
	using namespace std;
	if (!orient)
	{
		cerr << string{ begin(seq[r1]) + start1, begin(seq[r1]) + start1 + len } << endl;
		cerr << cwd::revComp(string{ begin(seq[r2]) + end2 - len, begin(seq[r2]) + end2 }) << endl;

	}
	else
	{
		cerr << string{ begin(seq[r1]) + start1, begin(seq[r1]) + start1 + len } << endl;
		cerr << string{ begin(seq[r2]) + start2, begin(seq[r2]) + start2 + len } << endl;
	}
}