#define MAIN_FILE
#include <commonMacro.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "grammar.h"
#include "symSet.h"
#include "sample.h"
#include "meParser.h"

std::string testFile = "VOC2007/Annotations/hwf_0000007.xml";
std::string VOC2007ImgPath = "VOC2007/JPEGImages/";
std::string VOC2007CharMapFName = "VOC2007/charmap_.txt";
bool IsUniformFile = false;
bool IsShowGraph = !false;
bool IsShowSample = !false;
bool withGT = true;

inline bool getSymbolMap(const std::string &filename, std::map<std::string, int> &symbolMap)
{
	if (!symbolMap.empty()) symbolMap.clear();
	std::fstream fs(filename, std::ios::in);
	if (!fs.is_open())
		HL_CERR_RETURN_FALSE("Failed to open the file " + filename);
	int charMapSize;
	fs >> charMapSize;
	for (size_t i = 0; i < charMapSize; i++)
	{
		std::string symbolName;
		int mapID;
		fs >> symbolName >> mapID;
		if (symbolMap.find(symbolName) == symbolMap.end())
		{
			symbolMap[symbolName] = mapID;
		}
		else
			HL_CERR_RETURN_FALSE("Symbol is duplicated in charmap file " + filename);
	}
	fs.close();
}

bool checkLatex(std::string &latex1, std::string &latex2)
{
	std::string s1, s2;
	for (size_t i = 0; i < latex1.size(); i++)
		if (latex1[i] != ' ')s1.push_back(latex1[i]);

	for (size_t i = 0; i < latex2.size(); i++)
		if (latex2[i] != ' ')s2.push_back(latex2[i]);

	std::vector<std::vector<std::string>> replaceMap = {
		{ "{\\int}", "int" },
		{ "{\\log}", "log" },
		{ "\\int", "int" },
		{ "\\log", "log" },
		{ "{log}", "log" },
		{ "\\\\", "\\" },
		{ "\\{", "{" },
		{ "\\}", "}" },
		//{ "\\lg", "lg" }
	};

	size_t pos = std::string::npos;
	for (size_t i = 0; i < replaceMap.size(); i++)
	{
		std::string &replaceStr = replaceMap[i][0], &aimStr = replaceMap[i][1];
		while ((pos = s1.find(replaceStr)) != std::string::npos)
			s1.replace(pos, replaceStr.size(), aimStr);

		while ((pos = s2.find(replaceStr)) != std::string::npos)
			s2.replace(pos, replaceStr.size(), aimStr);
	}

	std::cout << "s1     = " << s1 << std::endl;
	std::cout << "s2(GT) = " << s2 << std::endl;
	return s1 == s2;
}

int parseCmdArgs(int argc, char** argv)
{
	for (int i = 1; i < argc; i++)
	{
		if (std::string(argv[i]) == "-f")
		{
			testFile = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-voc2007imgdir")
		{
			VOC2007ImgPath = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-uniform")
		{
			if (std::string(argv[i + 1]) == "yes")
				IsUniformFile = true;
			else if (std::string(argv[i + 1]) == "no")
				IsUniformFile = false;
			i++;
		}
		else if (std::string(argv[i]) == "-withGT")
		{
			if (std::string(argv[i + 1]) == "yes")
				withGT = true;
			else if (std::string(argv[i + 1]) == "no")
				withGT = false;
			i++;
		}
		else if (std::string(argv[i]) == "-showSample")
		{
			if (std::string(argv[i + 1]) == "yes")
				IsShowSample = true;
			else if (std::string(argv[i + 1]) == "no")
				IsShowSample = false;
			i++;
		}
		else if (std::string(argv[i]) == "-showgraph")
		{
			if (std::string(argv[i + 1]) == "yes")
				IsShowGraph = true;
			else if (std::string(argv[i + 1]) == "no")
				IsShowGraph = false;
			i++;
		}
		else
		{
			return 0;
		}
	}

	return 1;
}

int main(int argc, char *argv[])
{
	parseCmdArgs(argc, argv);

	std::shared_ptr<SymSet> pSymSet = std::make_shared<SymSet>();
	pSymSet->load("Grammar/symbol_nd.types");
	
	std::shared_ptr<Grammar> pG = std::make_shared<Grammar>();
	pG->reSetup("Grammar/mathexp.gram", pSymSet);

	char gmmfile[1024] = "Grammar/sparels.gmm";
	std::shared_ptr<GMM> pGMM = std::make_shared<GMM>(gmmfile);

	MeParser meparser(pSymSet, pG, pGMM);

	std::shared_ptr<Sample> sample = std::make_shared<Sample>(pSymSet);

	std::map<std::string, int> symbolMap;
	getSymbolMap(VOC2007CharMapFName, symbolMap);
	
	std::string gtLatex;
	if (IsUniformFile)
	{
		sample->LoadFromUnifromFile(testFile, gtLatex, symbolMap, withGT);
	}
	else
	{
		sample->LoadFromVOC2007XML(testFile, VOC2007ImgPath, gtLatex);
	}
	
	//sample->ShowSample();
	std::string latexResult;
	IntevalTime(latexResult = meparser.parse(sample););
	
	if (withGT)
	{
		if (checkLatex(latexResult, gtLatex))
		{
			std::cout << "True Sample" << std::endl;
		}
		else
		{
			std::cout << "False Sample" << std::endl;
		}
		std::cout << "GT Latex : " << gtLatex << std::endl;
	}

	if (IsShowSample)
	{
		sample->ShowSample("SampleDebug");
	}

	if (IsShowGraph)
	{
		cv::Mat graphImg;
		std::shared_ptr<RelationSet> pRelSet = meparser.getRelationSet();
		pRelSet->DrawInImage(sample->getRGBImg(), graphImg);
		size_t pos = testFile.rfind("\\");
		pos = pos == std::string::npos ? 0 : pos;
		std::string pureFileName;
		pureFileName.assign(testFile.begin() + pos, testFile.end() - 4);
		cv::imshow(pureFileName, graphImg);
		int keyValue = cv::waitKey(0);
	}
		
	return 0;
}