
/*
//#include <hash_map>
//using namespace __gnu_cxx;

const size_t	BASE	= ((1<<19)-1);

struct edge {
	size_t from, to;		// addresses of newick_node's
	edge () {}
	edge (size_t _from, size_t _to) : from(_from), to(_to) {}
	
	bool operator==(const edge &b) const { return from==b.from && to==b.to; }
};

namespace __gnu_cxx {
	template<> class hash < pair<edge, edge> > {
		public:
		size_t operator () (const pair<edge, edge> &x) const {
			size_t res;
			
			res =            x.first.from;
			res = res*BASE + x.second.from;
			res = res*BASE + x.first.to;
			res = res*BASE + x.second.to;
			
			return res;
		}
	};
}
*/

int 		VERBOSE = 2;
const int	MAX_N	= 64; // TODO: dynamic memory

struct four {
	int p1, p2, q1, q2;

	four() {}
	four(int _p1, int _p2, int _q1, int _q2) : p1(_p1), p2(_p2), q1(_q1), q2(_q2) {}

	bool operator==(four a) { return p1==a.p1 && p2==a.p2 && q1==a.q1 && q2==a.q2; }

	four inv() { return four(q1, q2, p1, p2); }	// inverses the order of the trees
	four inv_edges() { return four(p2, p1, q2, q1); } // inverses the direction of edges
	void print() { fprintf(stderr, "%d %d %d %d\n", p1, p2, q1, q2); }	
};

class Tree {
	//int subtreeSizeMem[MAX_N][MAX_N];			// Memoized!!! subtreeSizeMem[v1][v2] is about the size of the subtree rooted in v2
	set<string> labelsMem[MAX_N][MAX_N];		// Memoized!!! labelsMem[v1][v2] includes all the labels(taxons) in the subtree (v1,v2)
	int num;									// used in suppressDegreeTwoVerticesRec(); new number of v1

	newick_node *file2tree(char *pcInputFile);
	
	void suppressDegreeTwoVerticesRec(int v1, int v2, int prevnum, vector<int> Vnew[], newick_node *nodesnew[]);
	void suppressDegreeTwoVertices();			// contracts all edges neighbouring a degree two vertex

	int findleaf();
	void printTreeRec(int prev, int i, int lvl);
		
public:
	char treeName[128];							// for identification and debug
	
	int n;										// nodes in tree
	newick_node *nodes[MAX_N];					// back reference my Tree -> newick
	vector<int> V[MAX_N];						// useful representation V[v1][v2] describes the tree with root v2 and doesn't include v1

	Tree() {}
	Tree(char *_filename, char *_treeName, bool _suppressDegreeTwoVertices);

	set<string> getLabels(int v1, int v2);		// labels in subtree (v1,v2)
	void Newickform2myTree(newick_node *curr);	// convert to my representation
	// int subtreeSize(int v1, int v2);

	void printTree();									// for debug
	void printSubtreeLabels(int v1, int v2);	// for debug
};

Tree::Tree(char *_filename, char *_treeName, bool _suppressDegreeTwoVertices=1)
{
	n = 0;
	strcpy(treeName, _treeName);

	newick_node *root = file2tree(_filename);
	Newickform2myTree(root);

	if (_suppressDegreeTwoVertices)
		suppressDegreeTwoVertices();

	if (VERBOSE >= 2) {
		printTree();
	}
}

void eraseSpaces(char *str)
{
	char *l, *c;

	for(l=c=str; *c!='\0'; c++)
		if(!isspace(*c))
			*l++ = *c;
	*l = '\0';		
}

newick_node *Tree::file2tree(char *pcInputFile)
{
	int iLen, iMaxLen;
	char *pcTreeStr;
	newick_node *root;
	char acStrArray[256];

	FILE *f;
	int i;

	// Open tree file
	f = fopen(pcInputFile, "r");
	if (f == NULL)
	{
		printf("Cannot open input file %s. Please check file name again.\n", pcInputFile);
		seqFreeAll();
		exit(-1);
	}

	// Read in tree string
	pcTreeStr = NULL;
	iLen = 0;
	iMaxLen = 0;
	while (1)
	{
		memset(acStrArray, '\0', 256);
		fgets(acStrArray, 255, f);
		eraseSpaces(acStrArray);
		if (acStrArray[0] == '\0' && feof(f))
			break;
		inputString(acStrArray, &pcTreeStr, &iLen, &iMaxLen);
	}
	fclose(f);

	//cout << '\n' << pcTreeStr << '\n' << endl;
	root = parseTree(pcTreeStr);

	// Free occupied memory
	seqFree(pcTreeStr);

	return root;
}

/*
// not used
int Tree::subtreeSize(int v1, int v2)
{
	int &size = subtreeSizeMem[v1][v2];

	if (size != -1) return size;
	else size = 1;

	vector<int>::iterator it;
	for (it=V[v2].begin(); it!=V[v2].end(); it++)
		if (*it != v1)
			size += subtreeSize(v2, *it);

	return size;
}
*/

void Tree::printTreeRec(int prev, int i, int lvl)
{
	fprintf(stderr, "%*s%d(%s)\n", lvl, "", i, nodes[i]->taxon);

	for (int j=0; j<V[i].size(); j++)
		if (prev != V[i][j])
			printTreeRec(i, V[i][j], lvl+1);
}

int Tree::findleaf()
{
	for(int i=0; i<n; i++)
		if(V[i].size() == 1)
			return i;
	
	return 0;
}

void Tree::printTree()
{
	int leaf;

	fprintf(stderr, "TreeName: %s\n", treeName);

	leaf = findleaf();
	printTreeRec(0, V[leaf][0], 0);
	
	/*
	for(int i=0; i<n; i++) {
		fprintf(stderr, "%d(%Trs) --> ", i, nodes[i]->taxon);
		for (int j=0; j<V[i].size(); j++) {
			if (j) fprintf(stderr, ", ");
			fprintf(stderr, "%d", V[i][j]);
		}
		fprintf(stderr, "\n");
	}
	*/
	
	fprintf(stderr, "\n");
}

void Tree::printSubtreeLabels(int v1, int v2)
{
	fprintf(stderr, "%s: Labels in subtree (%d,%d):\n", treeName, v1, v2);

	for (set<string>::iterator it=labelsMem[v1][v2].begin(); it!=labelsMem[v1][v2].end(); it++)
		fprintf(stderr, "   %s\n", it->c_str());
}

void Tree::Newickform2myTree(newick_node *curr)
{
	newick_child *child;
	int i, currn=n;
	
	nodes[currn] = curr;
	n++;
	for(child=curr->child, i=0; i<curr->childNum; child=child->next, i++) {
		V[currn].push_back (n);
		V[n].push_back (currn);

		Newickform2myTree(child->node);
	}
}

// changes the number of the vertex 'vold' to 'vnew' in the neighbours of 'v1'
/*
void Tree::renumerateVertex(int v1, int vnum, int vnew)
{
	vector<int>::iterator it;
	for(it=V[v1].begin(); it!=V[v1].end(); it++)
		if (*it == vold) {
			*it = vnew;
			break;
		}
}
*/

// renumerates the numbering to the interval [0,nprim)
// TODO: changes the node array
// returns the new number of vertices
// current position: V[v1] and Vnew[currnum]
// invariant: v1 doesn't need to be deleted
void Tree::suppressDegreeTwoVerticesRec(int v1, int v2, int prevnum, vector<int> Vnew[], newick_node *nodesnew[])
{
	int v3, vprev, currnum;

	for(vprev=v1; V[v2].size() == 2; vprev=v2, v2=v3)
		v3 = (V[v2][0] != vprev) ? V[v2][0] : V[v2][1];

	if (vprev!=v1) { // some verteces skipped
		suppressDegreeTwoVerticesRec(vprev, v2, num, Vnew, nodesnew);
	} else {
		num++;
		Vnew[prevnum].push_back(num);
		Vnew[num].push_back(prevnum);

		currnum = num;
		nodesnew[currnum] = nodes[v2];

		//fprintf(stderr, "(%d,%d) prevnum=%d currnum=%d num=%d\n", v1, v2, prevnum, currnum, num);

		if (V[v2].size() == 1) { // leaf
			return;
		} else {
			for(int i=0; i<V[v2].size(); i++)
				if ((v3=V[v2][i]) != v1)
					suppressDegreeTwoVerticesRec(v2, v3, currnum, Vnew, nodesnew);
		}
	}
}

// assumes 0-th node as root
void Tree::suppressDegreeTwoVertices()
{
	int i, nprim, leaf;
	vector<int> Vnew[MAX_N];
	newick_node *nodesnew[MAX_N];

	leaf = findleaf();
	assert(V[leaf].size() == 1);
	//fprintf(stderr, "leaf=%d\n", leaf);

	num=0;
	nodesnew[0] = nodes[leaf];
	suppressDegreeTwoVerticesRec(leaf, V[leaf][0], 0, Vnew, nodesnew);

	for (i=0, nprim=n; i<n; i++) {
		if (Vnew[i].empty()) {
			nprim = i;
			break;
		}

		V[i] = Vnew[i];
		nodes[i] = nodesnew[i];
	}
	
	for (; i<n; i++) {
		V[i].clear();
		nodes[i] = NULL;
	}

	fprintf(stderr, "n=%d nprim=%d\n", n, nprim);
	n = nprim;
}

set<string> Tree::getLabels(int v1, int v2)
{
	set<string> &curr = labelsMem[v1][v2];
	if (!curr.empty()) return curr;

	if (V[v2].size() == 1) { // includes only the parent v1
		curr.insert(nodes[v2]->taxon);
	} else {
		vector<int>::iterator it;
		for (it=V[v2].begin(); it!=V[v2].end(); it++)
			if (*it != v1) {
				getLabels(v2, *it); // calculates labelsMem
				curr.insert(labelsMem[v2][*it].begin(), labelsMem[v2][*it].end());
			}
	}
	
	return curr;
}

