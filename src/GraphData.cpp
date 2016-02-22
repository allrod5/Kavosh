#include <GraphData.hpp>
#include <GetTimeMs64.hpp>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

struct GD::Cell {
	std::vector<long unsigned int> vertices;
	std::vector<long unsigned int> label;
	Cell* next = NULL;
};

GD::GraphList::GraphList() {
	ini = new Cell;
	cursor = ini;
}

GD::GraphList::~GraphList() {
	cursor = ini;
	while(cursor) {
		ini = cursor->next;
		delete cursor;
		cursor = ini;
	}
}

void GD::GraphList::AddGraph() {
	Cell* aux = new Cell;
	cursor->next = aux;
	cursor = aux;
}

#ifdef DEBUG
	int LEVEL = 0;
	int COUNTER = 0;
#endif

void GD::PrintThat(GD::GraphList *kSubgraphs) {
	printf("PRINTING THAT\n");
	kSubgraphs->cursor = kSubgraphs->ini->next;
	while(kSubgraphs->cursor) {
		printf("%lu %lu %lu\n", kSubgraphs->cursor->vertices.at(0), kSubgraphs->cursor->vertices.at(1), kSubgraphs->cursor->vertices.at(2));
		kSubgraphs->cursor = kSubgraphs->cursor->next;
	}
}

void GD::Enumerate(TNGraph &G, long unsigned int k, GD::GraphList *kSubgraphs) {
	uint64 t0 = GetTimeMs64();
	std::vector<std::vector<long unsigned int>> Neighbors(G.GetMxNId());
	uint64 t1 = GetTimeMs64();
	GD::mapNeighbors(G, Neighbors);												// map neighbors for all vertices to avoid repeatly discovering the neighborhood of a vertex
	uint64 t2 = GetTimeMs64();
	#ifdef DEBUG
		if(DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1) printf("GetNodes: %lu ms\nMap Neighbors: %lu ms\n", t1-t0, t2-t1);
	#endif
	//for(long unsigned int i=0; G.IsNode(i); i++) { 								// iterate through all vertices of G
	for (TNGraph::TNodeI NI = G.BegNI(); NI < G.EndNI(); NI++) {
		int i = NI.GetId();
		#ifdef DEBUG
			if(DEBUG_LEVEL>=2) printf("Root: %d\n", i);
		#endif
		std::vector<bool> Visited(G.GetMxNId()); 								// list to keep track of visited vertices
		Visited[i] = true; 														// mark the current root as visited
		std::vector<long unsigned int> S(1);									// vertices selected in current level
		S[0] = i;																// only the root is selected in the first level
		std::vector<long unsigned int> subgraph(S);								// vector to store the vertices of the subgraph
		GD::Explore(G, Neighbors, i, Visited, S, subgraph, kSubgraphs, k-1);	// Explore all possible combinations that generete different trees
		Visited[i] = false;														// unmark current root as visited at the end
	}
	uint64 t5 = GetTimeMs64();
	#ifdef DEBUG
		if(DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1) {
			printf("Enumeration > Exploration loop: %lu ms\n", t5-t2);
			printf("\nSubgraphs of size %lu found: %d\n", k, COUNTER);
		}
		COUNTER = 0;
	#endif
}

void GD::Explore(TNGraph &G, std::vector<std::vector<long unsigned int>> &Neighbors, long unsigned int root,
	std::vector<bool> &Visited, std::vector<long unsigned int> &S, std::vector<long unsigned int> &subgraph,
	GD::GraphList *kSubgraphs, long unsigned int Remainder) {
	#ifdef DEBUG
		LEVEL++;
	#endif
	uint64 t0 = GetTimeMs64();
	if(Remainder==0) {								// if there are no more vertices to select the k-subgraph is formed (trivial case of the recursion)
		kSubgraphs->AddGraph();						// then add it to my
		kSubgraphs->cursor->vertices = subgraph;	// my data structure
		#ifdef DEBUG
			COUNTER++;									// count one more k-subgraph discovered
			if(DEBUG_LEVEL>=2) {
				for(int p=0; p<LEVEL; p++) printf("\t");
				printf("Subgraph: ");
				for(int p=0; p<subgraph.size(); p++) printf("%lu ", subgraph[p]);
				printf("\n");
			}
			LEVEL--;
		#endif
		return;
	}
	std::vector<long unsigned int> ValList;							// create a vector to store valid children at this level of the tree
	GD::Validate(Neighbors, Visited, S, root, ValList);				// fill the vector with the actual children
	uint64 t1 = GetTimeMs64();
	long unsigned int n = std::min((long unsigned int) ValList.size(), Remainder); // get the number of vertices to select at this level
	#ifdef DEBUG
		if(DEBUG_LEVEL>=4 || DEBUG_LEVEL==-1) printf("\tExplore > Validade: %lu ms\n", t1-t0);
		if(DEBUG_LEVEL>=2) {
			for(long unsigned int p=0; p<LEVEL; p++) printf("\t");
			printf("Valid Children: ");
			for(long unsigned int p=0; p<ValList.size(); p++) printf("%lu ", ValList[p]);
			printf("\n");
		}
	#endif
	for(long unsigned int i=1; i<=n; i++) {											// iterate from selecting a single vertex to selecting n vertices
		std::vector<long unsigned int> Index(ValList.size()-i, 0);					// create a vector to keep track of which vertices to delete from ValList in order to generate a combination
		std::vector<long unsigned int> C = GD::genComb(i, ValList, n, Index);		// generate the initial combination
		do {																									// repeat
			#ifdef DEBUG
				if(DEBUG_LEVEL>=2) {
					for(long unsigned int p=0; p<LEVEL; p++) printf("\t");
					printf("Selection: ");
					for(long unsigned int p=0; p<C.size(); p++) printf("%lu ", C[p]);
					printf("\n");
				}
			#endif
			subgraph.insert(subgraph.end(), C.begin(), C.end());												// compose subgraph with the current combination
			GD::Explore(G, Neighbors, root, Visited, C, subgraph, kSubgraphs, Remainder-i);						// recursive call to advance to next level of the implicit tree
			uint64 t2 = GetTimeMs64();
			GD::updateIndex(Index, i);																			// update Index to generate the next combination
			C = GD::genComb(i, ValList, n, Index);																// generate next combination
			uint64 t3 = GetTimeMs64();
			#ifdef DEBUG
				if(DEBUG_LEVEL>=4 || DEBUG_LEVEL==-1) printf("Get next combination time: %lu ms\n", t3-t2);
			#endif
			subgraph.erase(subgraph.end() - C.size(), subgraph.end());											// clear old combination from the subgraph "constructor"
		} while(!Index.empty() && Index.back()!=0);
	}
	for(long unsigned int i=0; i<ValList.size(); i++) Visited.at(ValList[i]) = false;	// unmark selected vertices at current level
	#ifdef DEBUG
		LEVEL--;
	#endif
}

void GD::mapNeighbors(TNGraph &G, std::vector<std::vector<long unsigned int>> &Neighbors) {
	#ifdef DEBUG
		if(DEBUG_LEVEL>=2) printf("Mapping neighborhood...\n");
	#endif
	for (TNGraph::TNodeI NI1 = G.BegNI(); NI1 < G.EndNI(); NI1++)
		for (TNGraph::TNodeI NI2 = G.BegNI(); NI2 < G.EndNI(); NI2++)
			if(NI1.IsNbrNId(NI2.GetId()))
				Neighbors[NI1.GetId()].push_back(NI2.GetId());
	/*for(long unsigned int i=0; i<G.GetNodes(); i++)
		for(long unsigned int j=0; j<G.GetNodes(); j++)
			if((G.GetNI(i)).IsNbrNId(j))
				Neighbors[i].push_back(j);*/
	#ifdef DEBUG
		if(DEBUG_LEVEL>=2) printf("Neighborhood mapped.\n");
	#endif
}

void GD::Validate(std::vector<std::vector<long unsigned int>>& Neighbors, std::vector<bool> &Visited, std::vector<long unsigned int> &S, long unsigned int root, std::vector<long unsigned int> &ValList) {
	for(long unsigned int i=0; i<S.size(); i++) {							// S is the selected children of the last level (or root)
		for(long unsigned int j=0; j<Neighbors[S[i]].size(); j++) {			// j is the label of possible children, which can't be greater than root
			if(Neighbors[S[i]][j]>root && !Visited[Neighbors[S[i]][j]]) {	// if not visited and it's in the neighborhood of S[i]
				Visited[Neighbors[S[i]][j]] = true;							// mark as visited
				ValList.push_back(Neighbors[S[i]][j]);						// include as valid children
			}
		}
	}
}

std::vector<long unsigned int> GD::genComb(long unsigned int i, std::vector<long unsigned int> ValList, long unsigned int n, std::vector<long unsigned int> &Index, long unsigned int level) {
	if(!Index.empty())
		while(ValList.size()!=i)
			ValList.erase(ValList.begin() + Index[level++]);
	return ValList;
}

void GD::updateIndex(std::vector<long unsigned int> &Index, long unsigned int n) {
	if(Index.empty()) return;					// nothing to update if Index is empty
	long unsigned int i = Index.size() - 1;		// go to last position of Index
	Index[i]++;									// increase it
	for(i; i>0; i--) {							// ascend Index
		if(Index[i]>n) Index[i-1]++;			// if the current position's value overpasses n "send one" to the position before
		else break;								// if not then there is no need to keep ascending
	}
	if(Index[i]>n) Index[i] = 0;				// if Index[i]>n it means i==0 so no way to "send one" as there is no position before, so cycle back to 0
	if(Index[++i]>Index[i-1])
		for(i; i<Index.size(); i++)
			Index[i] = Index[i-1];
}

void GD::Classify(TNGraph &G, GD::GraphList *kSubgraphs, long unsigned int n) {
	kSubgraphs->cursor = kSubgraphs->ini->next;

    while(kSubgraphs->cursor) {
		graph g[n*n];
		graph canon[n*n];
		int lab[n],ptn[n],orbits[n];
		static DEFAULTOPTIONS_DIGRAPH(options);
		options.getcanon = TRUE;
		statsblk stats;
		set *gv;
		int m = SETWORDSNEEDED(n);
		nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
		EMPTYGRAPH(canon,m,n);
		EMPTYGRAPH(g,m,n);
		for(long unsigned int i=0; i<n; i++) {
			long unsigned int a = kSubgraphs->cursor->vertices.at(i);
			for(long unsigned int j=0; j<n; j++) {
				long unsigned int b = kSubgraphs->cursor->vertices.at(j);
			    if(G.GetNI(a).IsOutNId(b)) {
					ADDONEARC(g,i,j,m);
				}
			}
		}

		densenauty(g,lab,ptn,orbits,&options,&stats,m,n,canon);
		for(long unsigned int i=0; i<m*n; i++)
			kSubgraphs->cursor->label.push_back(canon[i]);
		#ifdef DEBUG
			if(DEBUG_LEVEL>=3) {
				for(long unsigned int i=0; i<kSubgraphs->cursor->label.size(); i++) {
					printf("%lu ", kSubgraphs->cursor->label.at(i));
				}
				printf("= ");
				for(long unsigned int i=0; i<n; i++) {
					printf("%lu ", kSubgraphs->cursor->vertices.at(i));
				}
				printf("\n");
			}
		#endif

	    kSubgraphs->cursor = kSubgraphs->cursor->next;
	}
}

TNGraph* GD::Randomize(TNGraph &G) {
	TNGraph *R = new TNGraph(G);
	TNGraph::TNodeI a, b, c, d;
	for(long unsigned int i=0; i<3; i++) {
		for (a = R->BegNI(); a < R->EndNI(); a++) {
		//for(long unsigned int j=0; j<R->GetNodes(); j++) {
			//a = R->GetNI(j);
			if(a.GetOutDeg()==0) continue;
			do {
				b = R->GetRndNI();
			} while(a.GetId()==b.GetId() || b.GetOutDeg()==0);
			std::default_random_engine generator;
			std::uniform_int_distribution<int> distributionA(0,a.GetOutDeg()-1);
			std::uniform_int_distribution<int> distributionB(0,b.GetOutDeg()-1);
			c = R->GetNI(a.GetOutNId(distributionA(generator)));
			d = R->GetNI(b.GetOutNId(distributionB(generator)));
			R->DelEdge(a.GetId(), c.GetId());
			R->DelEdge(b.GetId(), d.GetId());
			R->AddEdge(a.GetId(), d.GetId());
			R->AddEdge(b.GetId(), c.GetId());
		}
	}
	return R;
}

void GD::GetFrequencies(GD::GraphList* kSubgraphs, std::map<std::vector<long unsigned int>, long unsigned int> &Frequencies) {
	kSubgraphs->cursor = kSubgraphs->ini->next;
	while(kSubgraphs->cursor) {
		Frequencies[kSubgraphs->cursor->label]++;
		kSubgraphs->cursor = kSubgraphs->cursor->next;
	}
	std::map<std::vector<long unsigned int>, long unsigned int>::iterator it;
	#ifdef DEBUG
		if(DEBUG_LEVEL>=3)
			for(it = Frequencies.begin(); it != Frequencies.end(); it++) {
				for(long unsigned int i=0; i<it->first.size(); i++)
			    	printf("%lu ", (it->first).at(i));
			    printf(": %lu\n", it->second);
			}
	#endif
}

void GD::DiscoverMotifs(std::vector<std::map<std::vector<long unsigned int>, long unsigned int>> &FVector,
						std::vector<std::vector<long unsigned int>> &Motifs, std::vector<long unsigned int> &IDs,
						long unsigned int motif_size, std::string destination, TNGraph &G, GD::GraphList *kSubgraphs)
{
	std::map<std::vector<long unsigned int>, long unsigned int>::iterator i;
	std::stringstream ss;
	ss << destination << "statistic_measures.txt";

	std::string s = ss.str();
	const char* path = s.c_str();
	std::ofstream TXT;
	TXT.open(path);

	TXT << "\t\t\t\t\t\t\t[Original Network]\t[\t\t\t\t\t\tRandom Networks\t\t\t\t\t\t\t]\n"
		<< "ID\t\tAdjacency Matrix\tFrequency\t\tMean Frequency\t\tStandard Deviation\t\tZ-Score\t\t\tP-value\n";

	for(i = FVector[0].begin(); i != FVector[0].end(); i++) {
    	double MeanFrequency = 0.0;
    	double StandardDeviation = 0.0;
    	double Zscore = 0.0;
    	double Pvalue = 0.0;
    	for(int index=1; index<FVector.size(); index++) {
    		std::map<std::vector<long unsigned int>, long unsigned int>::iterator j = FVector[index].find(i->first);
    		if(j!=FVector[index].end()) {
	    		if(j->second > i->second) Pvalue++;
	    		MeanFrequency += j->second;
	    	}
    	}
    	MeanFrequency /= FVector.size()-1;
    	for(long unsigned int index=1; index<FVector.size(); index++) {
    		std::map<std::vector<long unsigned int>, long unsigned int>::iterator j = FVector[index].find(i->first);
    		if(j!=FVector[index].end()) StandardDeviation += (MeanFrequency-j->second)*(MeanFrequency-j->second);
    	}
    	StandardDeviation /= FVector.size()-1;
    	StandardDeviation = sqrt(StandardDeviation);
    	Zscore = (i->second - MeanFrequency)/StandardDeviation;
    	Pvalue /= FVector.size()-1;

    	std::vector<long unsigned int> *vertices;
    	kSubgraphs->cursor = kSubgraphs->ini->next;
		while(kSubgraphs->cursor) {
			if(kSubgraphs->cursor->label==i->first) {
				vertices = &kSubgraphs->cursor->vertices;
				break;
			}
			kSubgraphs->cursor = kSubgraphs->cursor->next;
		}

		long unsigned int ID = 0;
		long unsigned int maxpow = motif_size*motif_size - 1;
		bool adjMatrix[motif_size][motif_size];

		for(int j=0; j<motif_size; j++) {
    		for(int k=0; k<motif_size; k++) {
    			adjMatrix[j][k] = false;
	    		if(G.GetNI(vertices->at(motif_size-1-j)).IsOutNId(vertices->at(motif_size-1-k))) {
	    			adjMatrix[j][k] = true;
	    			ID += (long) (pow(2, (maxpow - (j*motif_size+k))));
	    		}
	    	}
    	}

    	if(i->second > 4 && Pvalue < 0.01 && Zscore > 1) {
    		Motifs.push_back(i->first);
    		IDs.push_back(ID);
    	}

    	/*for(int j=motif_size-1; j>=0; j--) {
    		if(G.GetNI(vertices->at(motif_size-1)).IsOutNId(vertices->at(j)))
    			TXT << "1 ";
    		else
    			TXT << "0 ";
    	}*/

    	TXT << ID << "\t";
    	if(ID<999)
    		TXT << "\t";

    	for(int j=0; j<motif_size; j++) {
    		if(adjMatrix[0][j])
    			TXT << "1 ";
    		else
    			TXT << "0 ";
    	}

    	TXT << std::fixed << std::setprecision(2) << "\t\t\t"
    		<< i->second << "\t\t\t\t" << MeanFrequency << "\t\t\t\t" << StandardDeviation << "\t\t\t\t\t"
    		<< Zscore << "\t\t\t\t" << Pvalue << "\n";

    	/*for(int j=motif_size-2; j>=0; j--) {
    		for(int k=motif_size-1; k>=0; k--) {
	    		if(G.GetNI(vertices->at(j)).IsOutNId(vertices->at(k)))
	    			TXT << "1 ";
	    		else
	    			TXT << "0 ";
	    	}
	    	TXT << "\n";
    	}*/

	    for(int j=1; j<motif_size; j++) {
	    	TXT << "\t\t";
	    	for(int k=0; k<motif_size; k++) {
	    		if(adjMatrix[j][k])
	    			TXT << "1 ";
	    		else
	    			TXT << "0 ";
	    	}
	    	TXT << "\n";
	    }

    	TXT << "\n";

    	#ifdef DEBUG
    		if(DEBUG_LEVEL>=0) {
    			for(long unsigned int j=0; j<i->first.size(); j++) {
    				printf("%lu ", i->first.at(j));
    			}
    			printf(": Original Frequency = %lu : Mean Frequency = %f : Standard Deviation = %f : Zscore = %f : Pvalue = %f\n", i->second, MeanFrequency, StandardDeviation, Zscore, Pvalue);
    		}
    	#endif
	}

    TXT.close();

	#ifdef DEBUG
        if(DEBUG_LEVEL>=0) {
            for(int i=0; i<Motifs.size(); i++) {
                printf("\nMotif %d label: ", i+1);
                for(int j=0; j<Motifs[i].size(); j++)
                	printf("%lu ", Motifs[i][j]);
                printf("\n");
            }
        }
    #endif
}

TNGraph* GD::ConcatMotifs(TNGraph &G, std::vector<std::vector<long unsigned int>> &Motifs, GD::GraphList *kSubgraphs)
{
	#ifdef DEBUG
		if(DEBUG_LEVEL>=1) {
			printf("Concatenation Started\nInserting motifs as nodes...\n");
		}
	#endif

	TNGraph *H = new TNGraph(G);
	for(int m=0; m<Motifs.size(); m++) {
		kSubgraphs->cursor = kSubgraphs->ini->next;
		while(kSubgraphs->cursor) {
			if(kSubgraphs->cursor->label==Motifs[m]) {
				long unsigned int mNodeId = H->AddNode(-1);
				#ifdef DEBUG
					if(DEBUG_LEVEL>=2) {
						printf("Node %lu = (", mNodeId);
						int j;
						for(j=0; j<Motifs[m].size()-1; j++)
							printf("%lu, ", kSubgraphs->cursor->vertices.at(j));
						printf("%lu) added\n", kSubgraphs->cursor->vertices.at(j));
					}
				#endif
				for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
					TNGraph::TNodeI u = H->GetNI(kSubgraphs->cursor->vertices.at(i));
					for(int j=0; j<u.GetOutDeg(); j++)
						if(mNodeId != u.GetOutNId(j))
							H->AddEdge(mNodeId, u.GetOutNId(j));
					for(int j=0; j<u.GetInDeg(); j++)
						if(mNodeId != u.GetInNId(j))
							H->AddEdge(u.GetInNId(j), mNodeId);
				}
				#ifdef DEBUG
					if(DEBUG_LEVEL>=2) {
						TNGraph::TNodeI mNodeNI = H->GetNI(mNodeId);
						printf("\tDestinations: ");
						for(int i=0; i<mNodeNI.GetOutDeg(); i++) {
							printf("%d ", mNodeNI.GetOutNId(i));
						}
						printf("\n\tIncomings: ");
						for(int i=0; i<mNodeNI.GetInDeg(); i++) {
							printf("%d ", mNodeNI.GetInNId(i));
						}
						printf("\n");
					}
				#endif
			}
			kSubgraphs->cursor = kSubgraphs->cursor->next;
		}
	}

	#ifdef DEBUG
		if(DEBUG_LEVEL>=1) {
			printf("Removing vertices from motifs...\n");
		}
	#endif

	for(int m=0; m<Motifs.size(); m++) {
		kSubgraphs->cursor = kSubgraphs->ini->next;
		while(kSubgraphs->cursor) {
			if(kSubgraphs->cursor->label==Motifs[m]) {
				for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
					if(H->IsNode(kSubgraphs->cursor->vertices.at(i))) {
						H->DelNode(kSubgraphs->cursor->vertices.at(i));
						#ifdef DEBUG
							if(DEBUG_LEVEL>=2) {
								printf("Node %lu removed\n", kSubgraphs->cursor->vertices.at(i));
							}
						#endif
					}

				}
			}
			kSubgraphs->cursor = kSubgraphs->cursor->next;
		}
	}

	H->Defrag();

	return H;
}

void GD::ExportGDF(TNGraph &G, std::vector<std::vector<long unsigned int>> *Motifs, std::vector<long unsigned int> *IDs,
					 GD::GraphList *kSubgraphs, std::string destination)
{
	int size;
	if(Motifs!=NULL)
		size = Motifs->size();
	else
		size = 1;
	for(int m=0; m<size; m++) {
		std::vector<bool> marker(G.GetNodes(),false);
		if(Motifs!=NULL) {
			kSubgraphs->cursor = kSubgraphs->ini->next;
			while(kSubgraphs->cursor) {
				if(kSubgraphs->cursor->label==Motifs->at(m)) {
					for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++)
						marker[kSubgraphs->cursor->vertices.at(i)] = true;
				}
				kSubgraphs->cursor = kSubgraphs->cursor->next;
			}
		}

		std::stringstream ss;
		if(Motifs!=NULL)
			ss << destination << "motif_ID" << IDs->at(m) << "_highlighted.gdf";
		else
			ss << destination << "concatenated_motifs.gdf";
		std::string s = ss.str();
		const char* path = s.c_str();
		std::ofstream GDF;
		GDF.open(path);
		GDF << "nodedef>name VARCHAR,color VARCHAR\n";
		/*for(long unsigned int i=0; i<marker.size(); i++) {
			GDF << i << ",";
			if(marker[i]) GDF << "'255,0,0'\n";
			else GDF << "'0,0,0'\n";
		}*/
		for (TNGraph::TNodeI NI = G.BegNI(); NI < G.EndNI(); NI++) {
			int id = NI.GetId();
		    GDF << id << ",";
			if(Motifs!=NULL && marker[id]) GDF << "'255,0,0'\n";
			else GDF << "'0,0,0'\n";
		}
		GDF << "edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,color VARCHAR\n";

		TNGraph T(G);
		if(Motifs!=NULL) {
			kSubgraphs->cursor = kSubgraphs->ini->next;
			while(kSubgraphs->cursor) {
				if(kSubgraphs->cursor->label==Motifs->at(m)) {
					for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
						for(int j=0; j<kSubgraphs->cursor->vertices.size(); j++) {
							if((T.GetNI(kSubgraphs->cursor->vertices.at(i))).IsOutNId(kSubgraphs->cursor->vertices.at(j))) {
								GDF << kSubgraphs->cursor->vertices.at(i) << "," << kSubgraphs->cursor->vertices.at(j) << ",true,'255,0,0'\n";
								T.DelEdge(kSubgraphs->cursor->vertices.at(i), kSubgraphs->cursor->vertices.at(j));
							}
						}
					}
				}
				kSubgraphs->cursor = kSubgraphs->cursor->next;
			}
		}

		for(TNGraph::TEdgeI EI = T.BegEI(); EI < T.EndEI(); EI++)
			GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'0,0,0'\n";

		GDF.close();
	}
}

void GD::PrintMotifs(TNGraph &G, std::vector<std::vector<long unsigned int>> &Motifs, std::vector<long unsigned int> &IDs,
					 GD::GraphList *kSubgraphs, std::string destination)
{
	for(int m=0; m<Motifs.size(); m++) {
		PNGraph motif = TNGraph::New();
		kSubgraphs->cursor = kSubgraphs->ini->next;
		unsigned int group = 0;
		while(kSubgraphs->cursor->label!=Motifs[m])
			kSubgraphs->cursor = kSubgraphs->cursor->next;

		for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++)
			motif->AddNode(kSubgraphs->cursor->vertices.at(i));

		for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++)
			for(int j=0; j<kSubgraphs->cursor->vertices.size(); j++)
				if((G.GetNI(kSubgraphs->cursor->vertices.at(i))).IsOutNId(kSubgraphs->cursor->vertices.at(j)))
					motif->AddEdge(kSubgraphs->cursor->vertices.at(i),kSubgraphs->cursor->vertices.at(j));

		std::stringstream ss;
		ss << destination << "motif_ID" << IDs[m] << ".png";
		std::string s = ss.str();
		const char* path = s.c_str();
		TSnap::DrawGViz<PNGraph>(motif, gvlDot, path, "", NULL);

	}
}

void GD::SaveResults(TNGraph &G, std::vector<std::vector<long unsigned int>> &Motifs, GD::GraphList *kSubgraphs) {
	/*std::stringstream ss;
	ss << "../results/kavosh_results.txt";
	std::string s = ss.str();
	const char* path = s.c_str();
	std::ofstream file;
	file.open (path);
	file << "Kavosh 0.9.2 by Rodrigo M. Oliveira\n"
			"-----------------------------------\n"
			"\n"
			"Results:\n"
			"\n";
	for(int m=0; m<Motifs.size(); m++) {
		std::vector<bool> marker(G.GetNodes(),false);
		kSubgraphs->cursor = kSubgraphs->ini->next;
		while(kSubgraphs->cursor) {
			if(kSubgraphs->cursor->label==Motifs[m]) {
				for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
					marker[kSubgraphs->cursor->vertices.at(i)] = true;
				}
			}
			kSubgraphs->cursor = kSubgraphs->cursor->next;
		}
		std::stringstream ss;
		ss << "../results/motif_" << m << "_highlighted.gdf";
		std::string s = ss.str();
		const char* path = s.c_str();
		std::ofstream file;
		file.open (path);
		file << "nodedef>name VARCHAR,color VARCHAR\n";
		for(int i=0; i<marker.size(); i++) {
			GDF << i << ",";
			if(marker[i]) GDF << "'255,0,0'\n";
			else GDF << "'0,0,0'\n";
		}
		GDF << "edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,color VARCHAR\n";
		for(TNGraph::TEdgeI EI = G.BegEI(); EI < G.EndEI(); EI++) {
			if(marker[EI.GetSrcNId()]>0 && marker[EI.GetDstNId()]) {
				kSubgraphs->cursor = kSubgraphs->ini;
				while(kSubgraphs->cursor->next || ((GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'0,0,0'\n") && printf("miracle\n") && 0==1) ) {
					kSubgraphs->cursor = kSubgraphs->cursor->next;
					int i;
					for(i=0; i<kSubgraphs->cursor->vertices.size(); i++) if(kSubgraphs->cursor->vertices.at(i)==EI.GetSrcNId()) break;
					if(i==kSubgraphs->cursor->vertices.size()) continue;
					for(i=0; i<kSubgraphs->cursor->vertices.size(); i++) if(kSubgraphs->cursor->vertices.at(i)==EI.GetDstNId()) break;
					if(i==kSubgraphs->cursor->vertices.size()) continue;
					GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'255,0,0'\n";
					break;
				}
				//if(!kSubgraphs->cursor->next) GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'0,0,0'\n";
			}
			else GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'0,0,0'\n";
		}
		GDF.close();
	}
	file.close();*/
}
