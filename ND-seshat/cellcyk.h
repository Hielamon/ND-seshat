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
#ifndef _CELLCYK_
#define _CELLCYK_
#include <cstdio>
#include <vector>
#include "hypothesis.h"
#include "sample.h"



class CellCYK{
public:
	//Bounding box spatial region coordinates
	coo box;
	//int x, y; //top-left
	//int s, t; //bottom-right

			  //Hypotheses for every non-terminals
	int nnt;
	std::vector<std::shared_ptr<Hypothesis>> vNoTerm;

	//Strokes covered in this cell
	int nc;
	std::vector<bool> ccc;
	int talla; //total number of strokes

			   //Next cell in linked list (CYK table of same size)
	std::shared_ptr<CellCYK> sig;

	//Methods
	CellCYK(int n, int ncc)
	{
		sig = NULL;
		nnt = n;
		nc = ncc;
		talla = 0;

		//Create (empty) hypotheses
		vNoTerm.resize(nnt);

		//Create (empty) strokes covered
		ccc.resize(nc, false);
	}
	~CellCYK() {}

	void setRegion(Sample &m, const int &segIdx)
	{
		SegUnit &seg = m.getSegUnit(segIdx);
		ccc[segIdx] = true;

		box.x = seg.ROI.x;
		box.y = seg.ROI.y;
		box.s = seg.ROI.br().x;
		box.t = seg.ROI.br().y;
	}

	//Comparison operator for logspace ordering
	bool operator<(const CellCYK &C)
	{
		if (box.x < C.box.x)
			return true;
		if (box.x == C.box.x) {
			if (box.y < C.box.y)
				return true;
			if (box.y == C.box.y) {
				if (box.s < C.box.s)
					return true;
				if (box.s == C.box.s)
					if (box.t < C.box.t)
						return true;
			}
		}
		return false;
	}

	//Set the covered strokes to the union of cells A and B
	void ccUnion(CellCYK &A, CellCYK &B)
	{
		for (int i = 0; i<nc; i++)
			ccc[i] = (A.ccc[i] || B.ccc[i]) ? true : false;
	}

	//Check if cell H covers the same strokes that this
	bool ccEqual(std::shared_ptr<CellCYK> &pCell)
	{
		if (talla != pCell->talla)
			return false;

		for (int i = 0; i<nc; i++)
			if (ccc[i] != pCell->ccc[i])
				return false;

		return true;
	}

	//Check if the intersection between the strokes of this cell and H is empty
	bool compatible(std::shared_ptr<CellCYK> &pCell)
	{
		for (int i = 0; i<nc; i++)
			if (ccc[i] && pCell->ccc[i])
				return false;

		return true;
	}

  
};

inline void MergeRegionsCenter(std::shared_ptr<CellCYK>& pCA, int ntIDA,
						 std::shared_ptr<CellCYK>& pCB, int ntIDB,
						 std::shared_ptr<CellCYK>& pCS, int ntIDS,
						 char merge_cen)
{
	std::shared_ptr<Hypothesis> &a = pCA->vNoTerm[ntIDA];
	std::shared_ptr<Hypothesis> &b = pCB->vNoTerm[ntIDB];
	std::shared_ptr<Hypothesis> &s = pCS->vNoTerm[ntIDS];
	switch (merge_cen) {
	case 'A': //Data Hypothesis a
		s->lcen = a->lcen;
		s->rcen = a->rcen;
		break;
	case 'B': //Data Hypothesis b
		s->lcen = b->lcen;
		s->rcen = b->rcen;
		break;
	case 'C': //Center point
		s->lcen = (pCA->box.y + pCA->box.t) / 2;
		s->rcen = (pCB->box.y + pCB->box.t) / 2;
		break;
	case 'M': //Mean of both centers
		s->lcen = (a->lcen + b->lcen) / 2; //a->lcen;
		s->rcen = (a->rcen + b->rcen) / 2; //b->rcen;
		break;
	default:
		HL_CERR("Error: Unrecognized option " << merge_cen << " in merge regions");
	}
}


#endif
