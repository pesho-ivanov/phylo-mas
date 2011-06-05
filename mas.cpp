// Steel&Warnow_1993.pdf
// g++ mas.cpp -l newick -L newick_parser

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <assert.h>
using namespace std;

extern "C" {
#include "newick_parser/seqUtil.h"
#include "newick_parser/Newickform.h"
}

#include "myTree.h"

namespace MAST {
	namespace MASTprivate {
		// [p1][p2][q1][q1] means e-tree (p1,p2) in P and e-tree (q1,q2) in Q
		int dp[MAX_N][MAX_N][MAX_N][MAX_N]; // MST size for for each pair of subtrees TODO: change to hash_map
		vector<four> path[MAX_N][MAX_N][MAX_N][MAX_N]; // vector of subtrees of (p1,p2) used in the best solution for each pair of subtrees TODO: change to hash_map

		void matchSubtrees(Tree *P, Tree *Q, int &bestdp, vector<four> &bestpath, four pos, bool inv);
		void matchTreeWithSubtree(Tree *P, Tree *Q, int &bestdp, vector<four> &bestpath, four pos, bool inv);
		string subtree2string(Tree *P, Tree *Q, four pos, int lvl);
		string tree2string(Tree *P, Tree *Q, four half);
		int solveDP(Tree *P, Tree *Q, four pos);
	}
	
//public:
	string intersect(Tree *P, Tree *Q);
}

// index of v1 in V[v2]
int parentIndex(Tree *T, int v1, int v2)
{
	for(int i=0; i<T->V[v2].size(); i++)
		if(T->V[v2][i] == v1)
			return i;
	
	return -1;
}

// O(k!), k=max_degree; TODO: optimize next_permutation with equal numbers
void MAST::MASTprivate::matchSubtrees(Tree *P, Tree *Q, int &bestdp, vector<four> &bestpath, four pos, bool inv)
{
	int pchilds, qchilds;

	pchilds = P->V[pos.p2].size()-1; // minus one because of the parent
	qchilds = Q->V[pos.q2].size()-1;

	//cout << pchilds << " " << qchilds << endl;
	
	if (pchilds > qchilds) { // permute on P
		return matchSubtrees(Q, P, bestdp, bestpath, pos.inv(), true);
	} else { // permute on Q
		int c, i, j, p3, q3, qParentIndex, pParentIndex;
		four next;
		vector<int> perm(qchilds);
		vector<four> currpath;

		pParentIndex = parentIndex(P, pos.p1, pos.p2);
		qParentIndex = parentIndex(Q, pos.q1, pos.q2);
		assert(qParentIndex!=-1);

		//fprintf(stderr, "### pPar=%d qPar=%d\n", pParentIndex, qParentIndex);
		// miss the parent q1 in the permutation
		for(i=0; i<qchilds; i++)
			perm[i] = i + (i>=qParentIndex);
			
		if (pos==four(21,22,35,36)) {
			fprintf(stderr, "### pchilds=%d qchilds=%d | %d %d %d %d\n", pchilds, qchilds, P->V[pos.p2][0], P->V[pos.p2][1], Q->V[pos.q2][ perm[0] ], Q->V[pos.q2][ perm[1] ]);
		}

		do {
			int currsum = 0;
			currpath.clear();

			for(i=j=0; i<pchilds; j++) if(j!=pParentIndex) {
				p3 = P->V[pos.p2][ j ];
				q3 = Q->V[pos.q2][ perm[i] ];
				//for (j=0; j<qchilds; j++) if ((q3=Q->V[q2][ perm[j] ]) != q1) 
				
				if(!inv) {
					next = four(pos.p2, p3, pos.q2, q3);
					if (pos==four(21,22,35,36)) next.print();
					if(c=solveDP(P, Q, next)) {
						currsum += c;
						currpath.push_back(next);
					}
				}
				else {
					next = four(pos.q2, q3, pos.p2, p3);
					if(c=solveDP(Q, P, next)) { // return the right order <P, Q>
						currsum += c;
						currpath.push_back(next);
					}
				}
				
				i++;
			}

			if (currsum > bestdp) {
				bestdp = currsum;
				bestpath = currpath;
			}
		} while(next_permutation(perm.begin(), perm.end()));
	}
}

void MAST::MASTprivate::matchTreeWithSubtree(Tree *P, Tree *Q, int &bestdp, vector<four> &bestpath, four pos, bool inv)
{
	int c;
	four next;
	vector<int>::iterator it;

	for (it=Q->V[pos.q2].begin(); it!=Q->V[pos.q2].end(); it++)
		if (*it != pos.q1) {
			if (!inv) {
				next = four(pos.p1, pos.p2, pos.q2, *it);
				if (bestdp < (c=solveDP(P, Q, next))) {
					bestdp = c;
					bestpath.clear();
					bestpath.push_back(next);
				}
			} else {
				next = four(pos.q2, *it, pos.p1, pos.p2);
				if (bestdp < (c=solveDP(Q, P, next))) {
					bestdp = c;
					bestpath.clear();
					bestpath.push_back(next);
				}
			}
		}
}

int MAST::MASTprivate::solveDP(Tree *P, Tree *Q, four pos)
{
	int &bestdp = dp[pos.p1][pos.p2][pos.q1][pos.q2];
	vector<four> &bestpath = path[pos.p1][pos.p2][pos.q1][pos.q2];
	
	if (bestdp != -1) return bestdp;

	set<string> labelsP, labelsQ;
	
	labelsP = P->getLabels(pos.p1, pos.p2);
	labelsQ = Q->getLabels(pos.q1, pos.q2);

	if (labelsP.size() == 1 || labelsQ.size() == 1) {
		if (labelsP.size() == 1) {
			bestdp = labelsQ.count(*labelsP.begin());
			bestpath.clear(); // no children
		} else {
			bestdp = labelsP.count(*labelsQ.begin());

			// path through the vertex whose subtree contains the same label
			bestpath.clear();
			if (bestdp && labelsP.size()>1) {
				for (vector<int>::iterator it=P->V[pos.p2].begin(); it!=P->V[pos.p2].end(); it++)
					if (*it != pos.p1) {
						if (P->getLabels(pos.p2,*it).count(*labelsQ.begin())) {
							bestpath.push_back( four(pos.p2, *it, pos.q1, pos.q2) );
							break;
						}
					}
			}
		}
		
		return bestdp;
	}

	//cerr << q2 << endl;
	matchSubtrees(P, Q, bestdp, bestpath, pos, false);
	// TODO: remove comment
	matchTreeWithSubtree(P, Q, bestdp, bestpath, pos, false); // tree (p1,p2) with the subtrees (q2,*)
	matchTreeWithSubtree(Q, P, bestdp, bestpath, pos.inv(), true); // tree (q1,q2) with the subtrees (p2,*)

	return bestdp;
}

string MAST::MASTprivate::subtree2string(Tree *P, Tree *Q, four pos, int lvl=0)
{
	string str;
	newick_node *root = P->nodes[pos.p2];

	ostringstream out;

	fprintf(stderr, "%3d %*s dp[%d][%d][%d][%d]=%d -> %s %s\n", lvl, lvl, "",
			pos.p1, pos.p2, pos.q1, pos.q2, dp[pos.p1][pos.p2][pos.q1][pos.q2], root->taxon, Q->nodes[pos.q2]->taxon);

	if (P->V[pos.p2].size() == 1) // leaf
		out << root->taxon << ":" << root->dist;
	else {
		vector<four>::iterator it;
		vector<four> *childs;

		//childs = &path[pos.p1][pos.p2][pos.q1][pos.q2];

		while(1) {
			childs = &path[pos.p1][pos.p2][pos.q1][pos.q2];
			if (pos.p2 == (*childs)[0].p1)
				break;
			pos = (*childs)[0];
		}

		out << "(";
		for (it=childs->begin(); it!=childs->end(); it++) {
			if (it!=childs->begin()) out << ",";
			out << subtree2string(P, Q, *it, lvl+1);
		}

		if (root->taxon != NULL) {
			//out << "):bar" << root->dist;
			out << ")" << root->taxon << ":" << root->dist;
		} else
			out << "):" << root->dist;
	}

	return out.str();
}

string MAST::MASTprivate::tree2string(Tree *P, Tree *Q, four half)
{
	four otherHalf = four(half.p2, half.p1, half.q2, half.q1);
	fprintf(stderr, "half subtree = (%d,%d,%d,%d)\n", half.p1, half.p2, half.q1, half.q2);
	return "(" + subtree2string(P, Q, half) + "," + subtree2string(P, Q, otherHalf) + ");";
}

// Maximum Agreement Subtree
// returns string in newick's form
string MAST::intersect(Tree *P, Tree *Q)
{
	using namespace MASTprivate;
	
/*
	hash_map< pair<edge, edge>, int > solveDP;

	edge a(1,2), b(3,4);
	solveDP[ make_pair(a,b) ] = 5;
	
	cout << solveDP[ make_pair(a,b) ] << endl;
*/

	int p1, p2, q1, q2;
	vector<int>::iterator itP, itQ;

	int c, MAX = 0;
	four pos, halfBestSubtree;

	memset(dp, -1, sizeof(dp)); // TODO: change to hash_map

	for (pos.p1=0; pos.p1<P->n; pos.p1++)
		for (itP=P->V[pos.p1].begin(); itP!=P->V[pos.p1].end(); itP++) {
			pos.p2 = *itP;
			for (pos.q1=0; pos.q1<Q->n; pos.q1++)
				for (itQ=Q->V[pos.q1].begin(); itQ!=Q->V[pos.q1].end(); itQ++) {
					pos.q2 = *itQ;
					//printf ("%d %d %d %d\n", p1, *itP, q1, *itQ);
					
					if (MAX < (c=(solveDP(P, Q, pos) + solveDP(P, Q, pos.inv_edges())))) {
						MAX = c;
						halfBestSubtree = pos;
						//bestSubtree2 = four(*itP, p1, *itQ, q1);
					}
					
					//printf ("1. dp[%d][%d][%d][%d] = %d\n", p1, *itP, q1, *itQ, solveDP(P, Q, p1, *itP, q1, *itQ));
					//printf ("1. dp[%d][%d][%d][%d]+dp[%d][%d][%d][%d] = %d\n", p1, *itP, q1, *itQ, *itP, p1, *itQ, q1, c);
					//if(c==3) cerr << tree2string(P, Q, four(p1,*itP,q1,*itQ));*/
				}
		}

	cout << MAX << endl;
	fprintf(stderr, "Intersection leaves: %d\n", MAX);

	return tree2string(P, Q, halfBestSubtree);
}

void output(char *filename, string &TreeString)
{
	FILE *fout = fopen(filename, "w");

	fprintf(fout, "%s\n", TreeString.c_str());
}

int main(int argc, char **args)
{
	if (argc != 4) {
		printf ("Use arguments: <FirstInputTreeFile> <SecondInputTreeFile> <OutputAgreementTreeFile>\n");
		return 1;
	}

	// Initialize memory management procedure
	seqMemInit();

	Tree P(args[1], args[1]), Q(args[2], args[2]);
	string MASstring = MAST::intersect(&P, &Q); // TODO: return a root

	output(args[3], MASstring);
	
	//cerr << MASstring << endl;
	//cout << MASstring << endl;

	//char *MASs = MASstr.c_str();
	//MAS = parseTree(MASs);

	// End memory management procedutre and free all allocated space
	seqFreeAll();

	return 0;
}
