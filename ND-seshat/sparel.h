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
#ifndef _SPAREL_
#define _SPAREL_

class CellCYK;

#include <cstdio>
#include "hypothesis.h"
#include "cellcyk.h"
#include "gmm.h"
#include "sample.h"

//Aux functions
std::shared_ptr<Hypothesis> leftmost(std::shared_ptr<Hypothesis> &h) {
	if (!h.use_count())
		HL_CERR("Invalid Pointer");

	if (h->pt.use_count())
		return h;

	std::shared_ptr<Hypothesis> izq = leftmost(h->hleft);
	std::shared_ptr<Hypothesis> der = leftmost(h->hright);

	return izq->box.x < der->box.x ? izq : der;
}

std::shared_ptr<Hypothesis> rightmost(std::shared_ptr<Hypothesis> &h) {
	if (!h.use_count())
		HL_CERR("Invalid Pointer");

	if (h->pt.use_count())
		return h;

	std::shared_ptr<Hypothesis> izq = rightmost(h->hleft);
	std::shared_ptr<Hypothesis> der = rightmost(h->hright);

	return izq->box.s > der->box.s ? izq : der;
}

//Percentage of the area of region A that overlaps with region B
float solape(std::shared_ptr<Hypothesis> &a, std::shared_ptr<Hypothesis> &b) {
	int x = std::max(a->box.x, b->box.x);
	int y = std::max(a->box.y, b->box.y);
	int s = std::min(a->box.s, b->box.s);
	int t = std::min(a->box.t, b->box.t);

	if (s >= x && t >= y) {
		float aSolap = (s - x + 1.0)*(t - y + 1.0);
		float aTotal = (a->box.s - a->box.x + 1.0)*(a->box.t - a->box.y + 1.0);

		return aSolap / aTotal;
	}

	return 0.0;
}

class SpaRel {
public:
	const int NRELS = 6;
	const int NFEAT = 9;

private:
	std::shared_ptr<GMM> model;
	std::shared_ptr<Sample>mue;
	std::vector<float> probs;

	double compute_prob(std::shared_ptr<Hypothesis> &h1, std::shared_ptr<Hypothesis> &h2, int k)
	{
		//Set probabilities according to spatial constraints  

		if (k <= 2) {
			//Check left-to-right order constraint in Hor/Sub/Sup relationships
			std::shared_ptr<Hypothesis> rma = rightmost(h1);
			std::shared_ptr<Hypothesis> lmb = leftmost(h2);

			if (lmb->box.x < rma->box.x || lmb->box.s <= rma->box.s)
				return 0.0;
		}

		//Compute probabilities
		std::vector<float> sample(NFEAT);

		getFeas(h1, h2, sample, mue->RY);

		//Get spatial relationships probability from the model
		model->posterior(&sample[0], &probs[0]);

		//Slightly smooth probabilities because GMM classifier can provide
		//to biased probabilities. Thsi way we give some room to the
		//language model (the 2D-SCFG grammar)
		smooth(probs);

		return probs[k];
	}
	void smooth(std::vector<float> &post)
	{
		for (int i = 0; i<NRELS; i++)
			post[i] = (post[i] + 0.02) / (1.00 + NRELS*0.02);
	}

public:
	SpaRel(std::shared_ptr<GMM> &gmm, std::shared_ptr<Sample> &m)
	{
		model = gmm;
		mue = m;
		probs.resize(NRELS);
	}

	~SpaRel() {}

	void getFeas(std::shared_ptr<Hypothesis> &a, std::shared_ptr<Hypothesis> &b, std::vector<float> &sample, int ry)
	{
		//Normalization factor: combined height
		float F = std::max(a->box.t, b->box.t) - std::min(a->box.y, b->box.y) + 1;

		sample[0] = (b->box.t - b->box.y + 1) / F;
		sample[1] = (a->rcen - b->lcen) / F;
		sample[2] = ((a->box.s + a->box.x) / 2.0 - (b->box.s + b->box.x) / 2.0) / F;
		sample[3] = (b->box.x - a->box.s) / F;
		sample[4] = (b->box.x - a->box.x) / F;
		sample[5] = (b->box.s - a->box.s) / F;
		sample[6] = (b->box.y - a->box.t) / F;
		sample[7] = (b->box.y - a->box.y) / F;
		sample[8] = (b->box.t - a->box.t) / F;
	}

	double getHorProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		return compute_prob(ha, hb, 0);
	}
	double getSubProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		return compute_prob(ha, hb, 1);
	}
	double getSupProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		return compute_prob(ha, hb, 2);
	}
	double getVerProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb, bool strict = false)
	{
		//Pruning
		if (hb->box.y < (ha->box.y + ha->box.t) / 2
			|| abs((ha->box.x + ha->box.s) / 2 - (hb->box.x + hb->box.s) / 2) > 2.5*mue->RX
			|| (hb->box.x > ha->box.s || hb->box.s < ha->box.x))
			return 0.0;

		if (!strict)
			return compute_prob(ha, hb, 3);

		//Penalty for strict relationships
		float penalty = abs(ha->box.x - hb->box.x) / (3.0*mue->RX)
			+ abs(ha->box.s - hb->box.s) / (3.0*mue->RX);

		if (penalty > 0.95) penalty = 0.95;

		return (1.0 - penalty) * compute_prob(ha, hb, 3);
	}
	double getInsProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		if (solape(hb, ha) < 0.5 ||
			hb->box.x < ha->box.x || hb->box.y < ha->box.y)
			return 0.0;

		return compute_prob(ha, hb, 4);
	}
	double getMrtProb(std::shared_ptr<Hypothesis> &ha, std::shared_ptr<Hypothesis> &hb)
	{
		return compute_prob(ha, hb, 5);
	}
};

#endif
