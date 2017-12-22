/*Copyright 2014 Francisco Alvaro

 This file is part of SESHAT.

    SESHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SESHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SESHAT.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _PRODUCTION_
#define _PRODUCTION_

#include <cstdio>
#include <string>
#include <list>
#include <vector>
#include <commonMacro.h>
//#include "cellcyk.h"

//Binary productions of the grammar (2D-PCFG)
class ProductionB {
protected:
	std::string outStr;

public:
	int S;
	int A, B;
	std::string sS, sA, sB;

	float prior;
	char merge_cen;

	ProductionB(int s, int a, int b) : S(s), A(a), B(b) {}
	ProductionB(int s, int a, int b, float pr, const std::string &out)
		: S(s), A(a), B(b), outStr(out)
	{
		prior = pr > 0.0 ? log(pr) : -FLT_MAX;
		setMerges('C');
	}

	ProductionB(int s, int a, int b, std::string &ss, std::string &sa, std::string &sb) : S(s), A(a), B(b), sS(ss), sA(sa), sB(sb) {}
	ProductionB(int s, int a, int b, std::string &ss, std::string &sa, std::string &sb, float pr, const std::string &out)
		: S(s), A(a), B(b), sS(ss), sA(sa), sB(sb), outStr(out)
	{
		prior = pr > 0.0 ? log(pr) : -FLT_MAX;
		setMerges('C');
	}

	~ProductionB() {}

	/*float solape(Hypothesis *a, Hypothesis *b);
	void printOut(Grammar *G, Hypothesis *H);*/
	void setMerges(char c)
	{
		merge_cen = c;
	}

	

	bool check_out()
	{
		if (outStr.find("$1") == std::string::npos &&
			outStr.find("$2") == std::string::npos)
			return false;
		return true;
	}

	std::string &get_outstr()
	{
		return outStr;
	}

	//Pure virtual functions
	/*virtual char tipo() = 0;
	virtual void print() = 0;
	virtual void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid) = 0;*/
};

//Production S -> A : B
class ProductionH : public ProductionB {

public:
	ProductionH(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionH(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}

	~ProductionH() {}

	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production: S -> A / B
class ProductionV : public ProductionB {

public:
	ProductionV(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionV(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionV() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production: S -> A /u B
class ProductionU : public ProductionB {

public:
	ProductionU(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionU(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionU() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production: S -> A /e B
class ProductionVe : public ProductionB {

public:
	ProductionVe(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionVe(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionVe() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};



//Production: S -> A sse B
class ProductionSSE : public ProductionB {

public:
	ProductionSSE(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionSSE(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionSSE() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};



//Production: S -> A ^ B
class ProductionSup : public ProductionB {

public:
	ProductionSup(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionSup(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionSup() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production: S -> A _ B
class ProductionSub : public ProductionB {

public:
	ProductionSub(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionSub(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionSub() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production: S -> A ins B
class ProductionIns : public ProductionB {

public:
	ProductionIns(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionIns(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionIns() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production: S -> A mroot B
class ProductionMrt : public ProductionB {

public:
	ProductionMrt(int s, int a, int b)
		: ProductionB(s, a, b) {}
	ProductionMrt(int s, int a, int b, float pr, const std::string &out)
		: ProductionB(s, a, b, pr, out) {}
	~ProductionMrt() {}
	/*void print();
	char tipo();
	void mergeRegions(Hypothesis *a, Hypothesis *b, Hypothesis *s);
	void print_mathml(Grammar *G, Hypothesis *H, FILE *fout, int *nid);*/
};


//Production S -> term ( N clases )
class ProductionT{

public:
	ProductionT(int s, int nclases) :
		S(s), N(nclases), vTexStr(nclases), vMlType(nclases, 'z'),
		vProbs(nclases, 0), vClases(nclases, false) {}

	~ProductionT() {}

	void setClase(int k, float pr, const std::string& tex, char mlt)
	{
		vClases[k] = true;
		if (vTexStr[k].empty())
		{
			vTexStr[k] = tex;
			vProbs[k] = pr > 0.0 ? log(pr) : -FLT_MAX;
			vMlType[k] = mlt;
		}
		else
			HL_WARNING("WARNING: Terminal " << k << " redefined with label " << tex);
	}

	bool getClase(int k)
	{
		return vClases[k];
	}

	float getPrior(int k)
	{
		return vProbs[k];
	}

	std::string getTeX(int k)
	{
		return vTexStr[k];
	}

	char getMLtype(int k)
	{
		return vMlType[k];
	}

	int  getNoTerm()
	{
		return S;
	}

	void print()
	{
		int nc = 0;

		for (int i = 0; i<N; i++)
			if (vClases[i])
				nc++;

		std::cout << S << " -> [" << nc << "clases]" << std::endl;
	}

private:
	int S;
	std::vector<bool> vClases;
	std::vector<std::string> vTexStr;
	std::vector<char> vMlType;
	std::vector<float> vProbs;
	int N;
};

#endif
