#define MAIN_FILE
#include <commonMacro.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

typedef std::pair<std::string, char> SymUnit;
typedef std::pair<std::string, int> SymMap;

bool CompSym(const SymUnit &sym1, const SymUnit &sym2)
{
	return sym1.first < sym2.first;
}

void SortTheSymbol(const std::string seshatSymFName = "../ND-seshat/Grammar/symbol.types")
{
	std::fstream fs(seshatSymFName, std::ios::in);

	if (!fs.is_open())
		HL_CERR("Failed to open the file " + seshatSymFName);

	std::vector<SymUnit> vseshatSym;
	int seshatN;
	fs >> seshatN;
	for (size_t i = 0; i < seshatN; i++)
	{
		std::string symName;
		char type;
		fs >> symName >> type;
		vseshatSym.push_back(SymUnit(symName, type));
	}
	fs.close();

	std::sort(vseshatSym.begin(), vseshatSym.end(), CompSym);
	fs.open(seshatSymFName, std::ios::out);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " + seshatSymFName);

	fs << vseshatSym.size() << std::endl;
	for (size_t i = 0; i < seshatN; i++)
	{
		fs << vseshatSym[i].first << " " << vseshatSym[i].second << std::endl;
	}
	fs.close();
}

int main(int argc, char *argv[])
{
	std::string seshatSymFName = "../ND-seshat/Grammar/symbol_nd.types";
	std::fstream fs(seshatSymFName, std::ios::in);

	if (!fs.is_open())
		HL_CERR("Failed to open the file " + seshatSymFName);

	std::vector<SymUnit> vseshatSym;
	int seshatN;
	fs >> seshatN;
	for (size_t i = 0; i < seshatN; i++)
	{
		std::string symName;
		char type;
		fs >> symName >> type;
		vseshatSym.push_back(SymUnit(symName, type));
	}
	fs.close();
	
	std::vector<SymMap> vSymMap;

	std::string NDCharMapFName = "../ND-seshat/VOC2007/charmap.txt";
	fs.open(NDCharMapFName, std::ios::in);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " + NDCharMapFName);

	while (!fs.eof())
	{
		std::string symName;
		int index = -1;
		fs >> symName;

		//exclude the empty line
		if (symName == "\n" || symName.empty()) continue;

		for (size_t i = 0; i < seshatN; i++)
		{
			if (symName == vseshatSym[i].first)
			{
				index = i;
				break;
			}
		}
		vSymMap.push_back(SymUnit(symName, index));
	}
	fs.close();

	std::string NDCharMapFName_ = "../ND-seshat/VOC2007/charmap_.txt";
	fs.open(NDCharMapFName_, std::ios::out);
	if (!fs.is_open())
		HL_CERR("Failed to open the file " + NDCharMapFName_);

	fs << vSymMap.size() << std::endl;
	for (size_t i = 0; i < vSymMap.size(); i++)
	{
		fs << vSymMap[i].first << " " << vSymMap[i].second << std::endl;
	}
	fs.close();

	return 0;
}