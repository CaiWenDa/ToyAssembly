#include "utility.h"
#include <ctime>
#include <iostream>
#include <algorithm>

using namespace cwd;

void cwd::compSeqInRange(const seqan::Dna5String& r1, const seqan::Dna5String& r2, int r, uint start1, uint start2, uint end1, uint end2, uint len, bool orient)
{
	using namespace std;
	if (len > 30)
	{
		if (!orient)
		{
			string s1 = { begin(r1) + start1, begin(r1) + start1 + len };
			string s2 = cwd::revComp(string{ begin(r2) + end2 - len, begin(r2) + end2 });
			cerr << string{ begin(r1),  begin(r1) + 30 } + "...";
			cerr << string{ s1.begin(), s1.begin() + 30 } + "..." + string{ s1.end() - 30, s1.end() };
			cerr << "..." + string{ end(r1) - 30,  end(r1) } << endl;
			cerr << string{ begin(r2),  begin(r2) + 30 } + "...";
			cerr << string{ s2.begin(), s2.begin() + 30 } + "..." + string{ s2.end() - 30, s2.end() };
			cerr << "..." + string{ end(r2) - 30,  end(r2) } << endl;
		}

		else
		{
			string s1 = { begin(r1) + start1, begin(r1) + start1 + len };
			string s2 = { begin(r2) + start2, begin(r2) + start2 + len };
			cerr << string{ begin(r1),  begin(r1) + 30 } + "...";
			cerr << string{ s1.begin(), s1.begin() + 30 } + "..." + string{ s1.end() - 30, s1.end() };
			cerr << "..." + string{ end(r1) - 30,  end(r1) } << endl;
			cerr << string{ begin(r2),  begin(r2) + 30 } + "...";
			cerr << string{ s2.begin(), s2.begin() + 30 } + "..." + string{ s2.end() - 30, s2.end() };
			cerr << "..." + string{ end(r2) - 30,  end(r2) } << endl;
		}

	}
}

std::string cwd::getCurrentDate()
{
	time_t nowtime;
	nowtime = time(NULL); //获取日历时间   
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d-%H-%M", localtime(&nowtime));
	return tmp;
}

std::string cwd::getCurrentTime()
{
	time_t nowtime;
	nowtime = time(NULL); //获取日历时间   
	char tmp[64];
	strftime(tmp, sizeof(tmp), "[%Y-%m-%d %H:%M:%S]", localtime(&nowtime));
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

int convert(const string& paf)
{
	using namespace std;
	FILE* fp = fopen(paf.c_str(), "r");
	ifstream dict("");
	string Lname, id;
	map<string, string> dic;
	while (dict >> Lname >> id)
	{
		dic.insert({ Lname, id });
	}
	ofstream output("");
	string tName;
	int a, b, c, d, f, g, h, i, j;
	char e;
	char LName[100] = { 0 };
	char tname[100] = { 0 };
	//while (fscanf(fp, "%d,%[^,],%d,%d,%d,%c,%[^,],%d,%d,%d", &a, LName, &b, &c, &d, &e, tname, &f, &g, &h) != EOF)
	//    output << dic[string(LName)] << "," << dic[string(tname)] << "," << b << "," << c << "," << d << "," << e << ","
		//<< f << "," << g << "," << h << endl;
	while (fscanf(fp, "%[^,],%d,%d,%d,%c,%[^,],%d,%d,%d,%d,%d,%d\n", LName, &a, &b, &c, &e, tname, &d, &f, &g, &h, &i, &j) != EOF)
	{
		cout << dic[string(LName)] << endl;
		output << dic[string(LName)] << "," << a << "," << b << "," << c << "," << e << ","
			<< dic[string(tname)] << "," << d << "," << f << "," << g << "," << h << "," << i << "," << j << endl;
	}

	return 0;
}

uint cwd::parseGenomeSize(const char* str)
{
	string arg = str;
	uint result = 0;
	char exp = *arg.rbegin();
	string num_str = { arg.begin(), arg.end() - 1 };
	float num = atof(num_str.c_str());
	if (exp == 'k' or exp == 'K')
	{
		result = num * 1000;
	}
	else if (exp == 'm' or exp == 'M')
	{
		result = num * 1000 * 1000;
	}
	else if (exp == 'g' or exp == 'G')
	{
		result = num * 1000 * 1000 * 1000;
	}
	else
	{
		return 0;
	}
	return result;
}