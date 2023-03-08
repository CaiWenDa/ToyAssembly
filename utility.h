#pragma once
#include "overlap.h"

namespace cwd
{
	void compSeqInRange(const seqan::Dna5String& r1, const seqan::Dna5String& r2, int r, uint start1, uint start2, uint end1, uint end2, uint len, bool orient);
	float hamming(string& a, string& b);
	uint parseGenomeSize(const char* str);
	float jaccard(string& a, string& b);
	std::string getCurrentDate();
	std::string getCurrentTime();

}


