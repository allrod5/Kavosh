#include <fstream>
#include <string>
//#include <vector>
#include <exception>
#include <stdexcept>
#include <GraphData.hpp>
#include <nauty.h>
#include <stdio.h>

template <typename T>
struct GD::Cell {
	T label;
	Cell* next;
};

template <typename T>
GD::GraphList<T>::GraphList() {
	ini = new Cell<T>;
	ini->next = ini;
	cursor = ini;
}

template <typename T>
GD::GraphList<T>::~GraphList() {

}

template <typename T>
void GD::GraphList<T>::AddGraph() {
	Cell<T>* aux = new Cell<T>;
	cursor->next = aux;
}

void GD::Enumerate(TNGraph G, int k, /*GD::GraphList<*/int/*>* */kSubgraphs) {
	/*for(int i=0; G.IsNode(i); i++) {
		std::vector<bool> Visited(G.GetNodes());
		Visited[i] = true;
		std::vector<int> S(1);
		S[0] = i;
		GD::Explore(G, i, Visited, S, k-1, 1);
		Visited[i] = false;
	}
	*/
}
/*
void GD::Explore(TNGraph G, int root, std::vector<bool> &Visited, std::vector<int> &S, int Remainder, int depth) {
	if(Remainder==0) return;
	int* ValList[G.GetNodes()] = GD::Validate(G, Visited, S, root);

}

int* GD::Validate(TNGraph G, std::vector<bool> &Visited, std::vector<int> &S, int root) {
	int ValList[G.GetNodes()] = {0};
	for(int i=0; i<sizeof(S)/sizeof(*S); i++) {
		for(int j=S[i]+1; j<G.GetNodes(); j++) {
			if(G.IsEdge(S[i],j) && !Visited[j]) {
				Visited[j] = true;
				ValList[j] = 1;
			}
		}
	}
	return &ValList;
}*/

/*
graph* GD::ReadInput(char* file_path) {
	printf("setting up data...\n");
	DYNALLSTAT(int,lab1,lab1_sz);
    DYNALLSTAT(int,lab2,lab2_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
    DYNALLSTAT(graph,g1,g1_sz);
    DYNALLSTAT(graph,g2,g2_sz);
    DYNALLSTAT(graph,cg1,cg1_sz);
    DYNALLSTAT(graph,cg2,cg2_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    int n,m,i;
    n = 6301;
    size_t k;
    options.getcanon = TRUE;
    m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
    DYNALLOC1(int,lab2,lab2_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    DYNALLOC1(int,map,map_sz,n,"malloc");
    DYNALLOC2(graph,g1,g1_sz,n,m,"malloc");
    EMPTYGRAPH(g1,m,n);

    printf("reading file...\n");
    std::ifstream infile(file_path);
	std::string line;
	int a, b;
	while (std::getline(infile, line)){
		std::istringstream iss(line);
	    if (!(iss >> a >> b)) continue; // line is a comment
		ADDONEARC(g1,a,b,m);
	}
	printf("input read!\n");

	return g1;
}

void GD::ExtractData(PNGraph* G, int k, SubgraphList* sglist) {
	//enumerate
	GD::Enumerate(G, k, sglist);
	//classify
	GD::Classify(sglist);
}

void GD::Enumerate(PNGraph* G, int k, SubgraphList* sglist) {
	graph *g;
	//std::cout << "hey";
}

void GD::ExtractData(graph* G, int k, SubgraphList* sglist) {
	//enumerate
	printf("starting enumeration...\n");
	GD::Enumerate(G, k, sglist);
	//classify
	printf("starting classification...\n");
	GD::Classify(sglist);
}

GD::SubgraphList* GD::Classify(SubgraphList* sglist) {

}

void GD::DiscoverMotifs(SubgraphList* mlist) {

}
*/