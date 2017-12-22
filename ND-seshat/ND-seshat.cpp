#define MAIN_FILE
#include <commonMacro.h>

#include "grammar.h"
#include "symSet.h"
#include "sample.h"
#include "meParser.h"

std::string testXML = "VOC2007/Annotations/hwf_0000007.xml";
std::string imgPath = "VOC2007/JPEGImages/";

int parseCmdArgs(int argc, char** argv)
{
	for (int i = 1; i < argc; i++)
	{
		if (std::string(argv[i]) == "-xml")
		{
			testXML = std::string(argv[i + 1]);
			i++;
		}
		else if (std::string(argv[i]) == "-imgdir")
		{
			imgPath = std::string(argv[i + 1]);
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
	sample->LoadFromVOC2007XML(testXML, imgPath);
	//sample->ShowSample();
	
	meparser.parse(sample);

	sample->ShowSample();
	return 0;
}