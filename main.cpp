#include "assembly.h"
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
	string seqFileName;
	string asmFileName;
	string outFilePre = getCurrentDate();
	string paf;
	string overlapFile = "overlap.csv";
	string graphFileName;
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
	seqData_t seq;
	StringSet<CharString> ID;
	clock_t start = clock();
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
		uint block1 = 0;
		uint block2 = 0;
		uint b_size = length(seq) / thread_i;
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