#pragma once
#include <algorithm>
#include "symSet.h"
#include "grammar.h"
#include "sample.h"
#include "tablecyk.h"
#include "logspace.h"
#include "sparel.h"

#define NB 1

//#define SHOW_PIPELINE

void PrintSymSeg(std::shared_ptr<Hypothesis> &H)
{
	if (!H.use_count())
		HL_CERR("The Null Pointer");

	if (!H->pt.use_count()) {
		PrintSymSeg(H->hleft);
		PrintSymSeg(H->hright);
	}
	else {
		std::string clatex = H->pt->getTeX(H->clase);
		std::cout << clatex << std::endl;
	}
}

void PrintLatex(std::shared_ptr<Hypothesis> &H, std::shared_ptr<Grammar> &pG) {
	if (!H.use_count())
		HL_CERR("The Null Pointer");

	if (!H->pt.use_count())
	{
		//H->prod->printOut(pG, H);
		std::string &outStr = H->prod->get_outstr();
		if (!outStr.empty()) {

			int pd1 = outStr.find("$1");
			int pd2 = outStr.find("$2");

			int i = 0;
			if (pd2 >= 0 && pd1 >= 0 && pd2 < pd1) 
			{
				while (outStr[i] != '$' || outStr[i + 1] != '2') 
				{
					putchar(outStr[i]);
					i++;
				}
				i += 2;

				PrintLatex(H->hright, pG);
				if (H->hright->clase < 0)

				while (outStr[i] != '$' || outStr[i + 1] != '1')
				{
					putchar(outStr[i]);
					i++;
				}
				i += 2;

				PrintLatex(H->hleft, pG);
			}
			else
			{
				if (pd1 >= 0)
				{
					while (outStr[i] != '$' || outStr[i + 1] != '1')
					{
						putchar(outStr[i]);
						i++;
					}
					i += 2;

					PrintLatex(H->hleft, pG);
				}
				if (pd2 >= 0)
				{
					while (outStr[i] != '$' || outStr[i + 1] != '2')
					{
						putchar(outStr[i]);
						i++;
					}
					i += 2;

					PrintLatex(H->hright, pG);
				}
			}

			while (outStr[i]) 
			{
				putchar(outStr[i]);
				i++;
			}
		}
		
	}
	else {
		std::string clatex = H->pt->getTeX(H->clase);
		std::cout << clatex;
	}
	//std::cout << std::endl;
}

void drawCell(cv::Mat &src, std::shared_ptr<CellCYK> &pCell)
{
	coo &box = pCell->box;
	cv::Rect roi(box.x, box.y, box.s - box.x, box.t - box.y);
	cv::rectangle(src, roi, RandomColor(), 1);
}

void drawCellWithColor(cv::Mat &src, std::shared_ptr<CellCYK> &pCell, cv::Scalar &color)
{
	coo &box = pCell->box;
	cv::Rect roi(box.x, box.y, box.s - box.x, box.t - box.y);
	cv::rectangle(src, roi, color, 3);
}

class MeParser
{
public:
	MeParser(const std::shared_ptr<SymSet> &pSymSet_,
			 const std::shared_ptr<Grammar> &pG_,
			 const std::shared_ptr<GMM> &pGMM_)
		: pSymSet(pSymSet_), pG(pG_), pGMM(pGMM_) 
	{
		clusterF = 0.15680604;
		ptfactor = 0.37846575;
		pbfactor = 0.14864657;
		rfactor = 0.63225349;
		qfactor = 1.88593577;
		InsPen = 2.11917745;
	}
	MeParser() {}
	~MeParser()	{}

	void reSetup(const std::shared_ptr<SymSet> &pSymSet_,
				 const std::shared_ptr<Grammar> &pG_,
				 const std::shared_ptr<GMM> &pGMM_)
	{
		pSymSet = pSymSet_;
		pG = pG_;
		pGMM = pGMM_;
	}

	void parse(std::shared_ptr<Sample> &M)
	{
		//Compute the normalized size of a symbol for sample M
		M->detRefSymbol();

		M->computeSegDistance(M->RX, M->RY);

		int N = M->getSegUnitSize();
		int K = pG->noTerminales.size();

		//Cocke-Younger-Kasami (CYK) algorithm for 2D-SCFG
		TableCYK tcyk(N, K);

		std::cout << "CYK table initialization:" << std::endl;
		initCYKterms(M, tcyk, N, K);

		std::vector<std::shared_ptr<LogSpace>> logspace(N + 1);
		std::list<std::shared_ptr<CellCYK>> c1setH, c1setV, c1setU, c1setI, c1setM, c1setS;
		//SpaRel SPR();
		pSPR = std::make_shared<SpaRel>(pGMM, M);

		//Init spatial space for size 1
		logspace[1] = std::make_shared<LogSpace>(tcyk.get(1), tcyk.size(1), M->RX, M->RY);

		std::cout << "\nCYK parsing algorithm" << std::endl;
		std::cout << "Size 1: Generated " << tcyk.size(1) << std::endl;

		//CYK algorithm main loop
		for (int talla = 2; talla <= N; talla++)
		{
			for (int a = 1; a < talla; a++)
			{
				int b = talla - a;
				for (std::shared_ptr<CellCYK> c1 = tcyk.get(a); c1.use_count(); c1 = c1->sig)
				{
					//Clear lists
					c1setH.clear();
					c1setV.clear();
					c1setU.clear();
					c1setI.clear();
					c1setM.clear();
					c1setS.clear();


					//Get the subset of regions close to c1 according to different spatial relations
					logspace[b]->getH(c1, c1setH); //Horizontal (right)
					logspace[b]->getV(c1, c1setV); //Vertical (down)
					logspace[b]->getU(c1, c1setU); //Vertical (up)
					logspace[b]->getI(c1, c1setI); //Inside (sqrt)
					logspace[b]->getM(c1, c1setM); //mroot (sqrt[i])

#ifdef SHOW_PIPELINE
					int time = 1;
					cv::Mat tmpImg1 = M->getRGBImg(), tmpImg2;
					drawCellWithColor(tmpImg1, c1, cv::Scalar(0, 0, 255));
					int winx = 10, winy = 10;
					tmpImg2 = tmpImg1.clone();
					cv::imshow("Current cell", tmpImg2);
					cv::moveWindow("Current cell", winx, winy);
					cv::waitKey(time);
#endif // SHOW_PIPELINE
					
					for (auto c2 = c1setH.begin(); c2 != c1setH.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE

						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsH, Grammar::H);
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsSup, Grammar::SUP);
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsSub, Grammar::SUB);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setH traverse", tmpImg2);
					cv::moveWindow("c1setH traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setV.begin(); c2 != c1setV.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsV, Grammar::V);
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsVe, Grammar::VE);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setV traverse", tmpImg2);
					cv::moveWindow("c1setV traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setU.begin(); c2 != c1setU.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsV, Grammar::V);
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsVe, Grammar::VE);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setU traverse", tmpImg2);
					cv::moveWindow("c1setU traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setI.begin(); c2 != c1setI.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsIns, Grammar::INS);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setI traverse", tmpImg2);
					cv::moveWindow("c1setI traverse", winx, winy);
					cv::waitKey(time);

					tmpImg2 = tmpImg1.clone();
#endif // SHOW_PIPELINE
					for (auto c2 = c1setM.begin(); c2 != c1setM.end(); c2++)
					{
#ifdef SHOW_PIPELINE
						drawCell(tmpImg2, *c2);
#endif // SHOW_PIPELINE
						traverseProductionB(c1, *c2, M, tcyk, talla, pG->prodsMrt, Grammar::MRT);
					}
#ifdef SHOW_PIPELINE
					cv::imshow("c1setM traverse", tmpImg2);
					cv::moveWindow("c1setM traverse", winx, winy);
					cv::waitKey(time);
#endif // SHOW_PIPELINE

					//Look for combining {x_subs} y {x^sups} in {x_subs^sups}
					//TODO :


					
				}//for (std::shared_ptr<CellCYK> c1 = tcyk.get(a); c1.use_count(); c1 = c1->sig)
			}//for (int a = 1; a < talla; a++)

			std::cout << "Size " << talla << ": Generated " << tcyk.size(talla) << std::endl;

			if (talla <= N)
			{
				//Create new logspace structure of size "talla"
				logspace[talla] = std::make_shared<LogSpace>(tcyk.get(talla), tcyk.size(talla), M->RX, M->RY);

#ifdef SHOW_PIPELINE
				int time = 0;
				cv::Mat tmpImg1 = M->getRGBImg(), tmpImg2;
				int winx = 10, winy = 10;
				for (std::shared_ptr<CellCYK> c = tcyk.get(talla); c.use_count(); c = c->sig)
				{
					drawCellWithColor(tmpImg1, c, RandomColor());
				}
				std::stringstream ioStr;
				ioStr << talla;
				cv::imshow("talla" + ioStr.str(), tmpImg1);
				cv::moveWindow("talla" + ioStr.str(), winx, winy);
				cv::waitKey(time);
#endif // SHOW_PIPELINE
			}

			
		}//for (int talla = 2; talla <= N; talla++)


		std::shared_ptr<Hypothesis> mlh;
		std::shared_ptr<CellCYK> mlc;
		int mlhIdx;
		tcyk.getMLInfo(mlh, mlc, mlhIdx);

		assert(mlh == mlc->vNoTerm[mlhIdx]);

		if (!mlh.use_count() || !mlc.use_count())
			HL_CERR("\nNo hypothesis found!!");

		std::cout << "\nMost Likely Hypothesis " << mlc->talla << " Segmentations" << std::endl;

		std::cout << "Math Symbols : " << std::endl;
		PrintSymSeg(mlh);
		std::cout << std::endl;

		std::cout << "Latex : " << std::endl;
		PrintLatex(mlh, pG);
		std::cout << std::endl;
	}

	std::shared_ptr<SymSet> pSymSet;
	std::shared_ptr<Grammar> pG;
	std::shared_ptr<GMM> pGMM;
	std::shared_ptr<SpaRel> pSPR;

	//ClusterF
	float clusterF;

	//ProductionTSF, ProductionBSF, RelationSF
	float ptfactor, pbfactor, rfactor;

	//SymbolSF, InsPenalty
	float qfactor, InsPen;

private:

	void initCYKterms(std::shared_ptr<Sample> &M, TableCYK &tcyk, int N, int K)
	{
		for (int i = 0; i<M->getSegUnitSize(); i++) {

			int cmy, asc, des;

			std::cout << "Segment " << i << std::endl;

			std::vector<int> clase(NB);
			std::vector<float> pr(NB);

			M->getSegUnitInfo(i, NB, clase, pr, cmy, asc, des);

			std::shared_ptr<CellCYK> pCell = std::make_shared<CellCYK>(K, N);

			pCell->setRegion(*M, i);

			bool insertar = false;
			for (size_t i = 0; i < pG->prodTerms.size(); i++)
			{
				std::shared_ptr<ProductionT> &prod = pG->prodTerms[i];
				int ntID = prod->getNoTerm();

				for (int k = 0; k < NB; k++)
				{
					if (pr[k] > 0.0 && prod->getClase(clase[k]) && prod->getPrior(clase[k]) > -FLT_MAX)
					{
						//the origin probability
						/*float prob = log(InsPen)
							+ ptfactor * prod->getPrior(clase[k])
							+ qfactor  * log(pr[k])
							+ dfactor  * log(duration->prob(clase[k], 1));*/
						float prob = log(InsPen) 
							+ ptfactor * prod->getPrior(clase[k]) 
							+ qfactor * log(pr[k]);

						if (pCell->vNoTerm[ntID].use_count())
						{
							if (pCell->vNoTerm[ntID]->pr > prob + prod->getPrior(clase[k]))
								continue;
							else
								pCell->vNoTerm[ntID].reset();
						}

						insertar = true;

						//Create new symbol
						pCell->vNoTerm[ntID] = std::make_shared<Hypothesis>(clase[k], prob, pCell->box);
						pCell->vNoTerm[ntID]->pt = prod;

						//Compute the vertical centroid according to the type of symbol
						int cen, type = pSymSet->symType(clase[k]);
						if (type == 0)       cen = cmy; //Normal
						else if (type == 1) cen = asc; //Ascendant
						else if (type == 2) cen = des; //Descendant
						else                cen = (pCell->box.t + pCell->box.y)*0.5; //Middle point

																	   //Vertical center
						pCell->vNoTerm[ntID]->lcen = cen;
						pCell->vNoTerm[ntID]->rcen = cen;
					}
				}
			}

			if (insertar)
			{
				for (int j = 0; j < K; j++) {
					if (pCell->vNoTerm[j].use_count()) {
						std::cout << pSymSet->strClase(pCell->vNoTerm[j]->clase) << " "
							<< pG->key2str(j) << " " << exp(pCell->vNoTerm[j]->pr) << std::endl;
						/*printf("%12s [%s] %g\n", pSymSet->strClase(pCell->vNoTerm[j]->clase),
						pG->key2str(j), exp(pCell->vNoTerm[j]->pr));*/
					}
				}

				//Add to parsing table (size=1)
				tcyk.add(1, pCell, -1, pG->esInit);
			}
		}
	}

	 std::shared_ptr<CellCYK> fusion(std::shared_ptr<Sample> &M, std::shared_ptr<ProductionB>& pd, std::shared_ptr<CellCYK>& pCA,
									 int ntIDA, std::shared_ptr<CellCYK>& pCB, int ntIDB, double prob)
	 {
		 std::shared_ptr<CellCYK> pCS;

		 if (!pCA->compatible(pCB) || pd->prior == -FLT_MAX)
			 return pCS;

		 //Penalty according to distance between Cell
		 float grpen = 1.0;
		 //TODO : induce the penalty function
		 if (clusterF > 0.0) {

			 grpen = M->group_penalty(pCA->ccc, pCB->ccc);
			 //If distance is infinity -> not visible
			 if (grpen >= M->INF_DIST)
				 return pCS;

			 //Compute penalty
			 grpen = std::max(-0.999999f, grpen);
			 grpen = 1.0 / (1.0 + grpen);
			 grpen = pow(grpen, clusterF);
			 //grpen = 1 * grpen;
			 //std::cout << "grpen = " << grpen << std::endl;
		 }

		 //Get nonterminal
		 int ps = pd->S;
		 int N = M->getSegUnitSize();
		 int K = pG->noTerminales.size();

		 pCS = std::make_shared<CellCYK>(K, N);

		 std::shared_ptr<Hypothesis> &A = pCA->vNoTerm[ntIDA], &B = pCB->vNoTerm[ntIDB];

		 //Compute the (log)probability
		 //prob = pbfactor * pd->prior + rfactor * log(prob * grpen) + A->pr + B->pr;
		 prob = pbfactor *pd->prior + rfactor * log(prob * grpen) + A->pr + B->pr;

		 //Copute resulting region
		 pCS->box.x = std::min(pCA->box.x, pCB->box.x);
		 pCS->box.y = std::min(pCA->box.y, pCB->box.y);
		 pCS->box.s = std::max(pCA->box.s, pCB->box.s);
		 pCS->box.t = std::max(pCA->box.t, pCB->box.t);

		 //Set the Cell covered
		 pCS->ccUnion(*pCA, *pCB);

		 int clase = -1;
		 if (!pd->check_out() && pSymSet->checkClase(pd->get_outstr()))
			 clase = pSymSet->keyClase(pd->get_outstr());

		 //Create hypothesis
		 pCS->vNoTerm[ps] = std::make_shared<Hypothesis>(clase, prob, pCS->box);

		 MergeRegionsCenter(pCA, ntIDA, pCB, ntIDB, pCS, ps, pd->merge_cen);

		 pCS->vNoTerm[ps]->hleft = A;
		 pCS->vNoTerm[ps]->hright = B;
		 pCS->vNoTerm[ps]->prod = pd;

		 if (clase >= 0)
		 {
			 for (size_t i = 0; i < pG->prodTerms.size(); i++)
			 {
				 std::shared_ptr<ProductionT> &prod = pG->prodTerms[i];

				 if (prod->getClase(clase) && prod->getPrior(clase) > -FLT_MAX)
				 {
					 pCS->vNoTerm[ps]->pt = prod;
					 break;
				 }
			 }
		 }

		 return pCS;
	 }

	 
	 void traverseProductionB(std::shared_ptr<CellCYK> &c1, 
							  std::shared_ptr<CellCYK> &c2,
							  std::shared_ptr<Sample> &M, TableCYK &tcyk, int talla,
							  std::vector<std::shared_ptr<ProductionB>> &vProds, Grammar::PBTYPE pType)
	 {
		 for (size_t i = 0; i < vProds.size(); i++)
		 {
			 std::shared_ptr<ProductionB> &pd = vProds[i];

			 if (pd->prior == -FLT_MAX) continue;
			 //Production S -> A B
			 int ps = pd->S;
			 int pa = pd->A;
			 int pb = pd->B;

			 if (c1->vNoTerm[pa].use_count() && c2->vNoTerm[pb].use_count()) {
				 double cdpr = 0.0;
				 switch (pType)
				 {
				 case Grammar::H:
					 cdpr = pSPR->getHorProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 case Grammar::SUP:
					 cdpr = pSPR->getSupProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 case Grammar::SUB:
					 cdpr = pSPR->getSubProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 case Grammar::V:
				 case Grammar::VE:
					 cdpr = pSPR->getVerProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 case Grammar::INS:
					 cdpr = pSPR->getInsProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 case Grammar::MRT:
					 cdpr = pSPR->getMrtProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 case Grammar::SSE:
					 //cdpr = pSPR->getSupProb(c1->vNoTerm[pa], c2->vNoTerm[pb]);
					 break;
				 default:
					 break;
				 }
				 //double cdpr = SPR.getHorProb(c1->noterm[pa], (*c2)->noterm[pb]);
				 if (cdpr <= 0.0) continue;

				 std::shared_ptr<CellCYK> pCell = fusion(M, pd, c1, pa, c2, pb, cdpr);
				 //CellCYK *cd = fusion(M, *it, c1->noterm[pa], (*c2)->noterm[pb], M->nStrokes(), cdpr);

				 if (!pCell.use_count()) continue;

				 if (pCell->vNoTerm[ps].use_count()) {
					 tcyk.add(talla, pCell, ps, pG->esInit); //Add to parsing table (size=talla)
				 }
				 else {
					 tcyk.add(talla, pCell, -1, pG->esInit); //Add to parsing table
				 }
			 }
		 }
	 }

};

