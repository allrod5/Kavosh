#include <GraphData.hpp>
#include <GetTimeMs64.hpp>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

struct GD::Cell {
    std::vector<int> vertices;
    std::multiset<unsigned long> label;
    Cell* next = NULL;
};

GD::GraphList::GraphList() {
    ini = new Cell;
    cursor = ini;
}

GD::GraphList::~GraphList() {
    cursor = ini;
    while (cursor) {
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

GD::MetaObject::MetaObject(PNGraph &G) {
    this->G = &*G;
    metaMap = new std::map<int, std::vector<int>*>;
    for (TNGraph::TNodeI NI = this->G->BegNI(); NI < this->G->EndNI(); NI++) {
        std::vector<int> *v = new std::vector<int>;
        v->push_back(NI.GetId());
        metaMap->emplace(NI.GetId(), v);
    }
}

GD::MetaObject::~MetaObject() {

}

GD::MetaObject::MetaObject(GD::MetaObject &metaObj) {
    this->G = new TNGraph(*metaObj.G);
    this->metaMap = metaObj.metaMap;
}

#ifdef DEBUG
    int LEVEL = 0;
    int COUNTER = 0;
#endif

void GD::Enumerate(TNGraph &G, unsigned long k, GD::GraphList *kSubgraphs) {
    uint64 t0 = GetTimeMs64();
    std::vector<std::vector<int>> Neighbors((size_t) G.GetMxNId());
    uint64 t1 = GetTimeMs64();
    GD::mapNeighbors(G, Neighbors);  // map neighbors for all vertices to avoid repeatly discovering the neighborhood of a vertex
    uint64 t2 = GetTimeMs64();
    #ifdef DEBUG
        if (DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1) {
            printf("GetNodes: %lu ms\nMap Neighbors: %lu ms\n", t1-t0, t2-t1);
        }
    #endif
    //for(unsigned long i=0; G.IsNode(i); i++) {  // iterate through all vertices of G
    for (TNGraph::TNodeI NI = G.BegNI(); NI < G.EndNI(); NI++) {
        int i = NI.GetId();
        #ifdef DEBUG
            if (DEBUG_LEVEL>=2) {
                printf("Root: %d\n", i);
            }
        #endif
        std::vector<bool> Visited((size_t) G.GetMxNId());  // list to keep track of visited vertices
        Visited[i] = true;  // mark the current root as visited
        std::vector<int> S(1);  // vertices selected in current level
        S[0] = i;  // only the root is selected in the first level
        std::vector<int> subgraph(S);  // vector to store the vertices of the subgraph
        GD::Explore(G, Neighbors, i, Visited, S, subgraph, kSubgraphs, k-1);  // Explore all possible combinations that generete different trees
        Visited[i] = false;  // unmark current root as visited at the end
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

void GD::Explore(
        TNGraph &G, std::vector<std::vector<int>> &Neighbors, int root, std::vector<bool> &Visited, std::vector<int> &S,
        std::vector<int> &subgraph, GD::GraphList *kSubgraphs, unsigned long Remainder)
{
    #ifdef DEBUG
        LEVEL++;
    #endif
    uint64 t0 = GetTimeMs64();
    if (Remainder==0) {  // if there are no more vertices to select the k-subgraph is formed (trivial case of the recursion)
        kSubgraphs->AddGraph();  // then add it to my
        kSubgraphs->cursor->vertices = subgraph;  // my data structure
        #ifdef DEBUG
            COUNTER++;  // count one more k-subgraph discovered
            if (DEBUG_LEVEL>=2) {
                for (int p=0; p<LEVEL; p++) {
                    printf("\t");
                }
                printf("Subgraph: ");
                for (int p=0; p<subgraph.size(); p++) {
                    printf("%lu ", subgraph[p]);
                }
                printf("\n");
            }
            LEVEL--;
        #endif
        return;
    }
    std::vector<int> ValList;  // create a vector to store valid children at this level of the tree
    GD::Validate(Neighbors, Visited, S, root, ValList);  // fill the vector with the actual children
    uint64 t1 = GetTimeMs64();
    unsigned long n = std::min(ValList.size(), Remainder);  // get the number of vertices to select at this level
    #ifdef DEBUG
        if (DEBUG_LEVEL>=4 || DEBUG_LEVEL==-1) {
            printf("\tExplore > Validade: %lu ms\n", t1-t0);
        }
        if (DEBUG_LEVEL>=2) {
            for (unsigned long p=0; p<LEVEL; p++) {
                printf("\t");
            }
            printf("Valid Children: ");
            for (unsigned long p=0; p<ValList.size(); p++) {
                printf("%lu ", ValList[p]);
            }
            printf("\n");
        }
    #endif
    for (unsigned long i=1; i<=n; i++) {  // iterate from selecting a single vertex to selecting n vertices
        std::vector<int> Index(ValList.size()-i, 0);  // create a vector to keep track of which vertices to delete from ValList in order to generate a combination
        std::vector<int> C = GD::genComb(i, ValList, Index, 0);  // generate the initial combination
        do {  // repeat
            #ifdef DEBUG
                if(DEBUG_LEVEL>=2) {
                    for (unsigned long p=0; p<LEVEL; p++) {
                        printf("\t");
                    }
                    printf("Selection: ");
                    for (unsigned long p=0; p<C.size(); p++) {
                        printf("%lu ", C[p]);
                    }
                    printf("\n");
                }
            #endif
            subgraph.insert(subgraph.end(), C.begin(), C.end());  // compose subgraph with the current combination
            GD::Explore(G, Neighbors, root, Visited, C, subgraph, kSubgraphs, Remainder-i);  // recursive call to advance to next level of the implicit tree
            uint64 t2 = GetTimeMs64();
            GD::updateIndex(Index, i);  // update Index to generate the next combination
            C = GD::genComb(i, ValList, Index, 0);  // generate next combination
            uint64 t3 = GetTimeMs64();
            #ifdef DEBUG
                if (DEBUG_LEVEL>=4 || DEBUG_LEVEL==-1) {
                    printf("Get next combination time: %lu ms\n", t3-t2);
                }
            #endif
            subgraph.erase(subgraph.end() - C.size(), subgraph.end());  // clear old combination from the subgraph "constructor"
        } while (!Index.empty() && Index.back()!=0);
    }
    for(size_t i=0; i<ValList.size(); i++) {
        Visited.at((size_t) ValList[i]) = false;  // unmark selected vertices at current level
    }
    #ifdef DEBUG
        LEVEL--;
    #endif
}

void GD::mapNeighbors(TNGraph &G, std::vector<std::vector<int>> &Neighbors) {
    #ifdef DEBUG
        if(DEBUG_LEVEL>=2) {
            printf("Mapping neighborhood...\n");
        }
    #endif
    for (TNGraph::TNodeI NI1 = G.BegNI(); NI1 < G.EndNI(); NI1++) {
        for (TNGraph::TNodeI NI2 = G.BegNI(); NI2 < G.EndNI(); NI2++) {
            if(NI1.IsNbrNId(NI2.GetId())) {
                Neighbors[NI1.GetId()].push_back(NI2.GetId());
            }
        }
    }
    /*for(unsigned long i=0; i<G.GetNodes(); i++)
        for(unsigned long j=0; j<G.GetNodes(); j++)
            if((G.GetNI(i)).IsNbrNId(j))
                Neighbors[i].push_back(j);*/
    #ifdef DEBUG
        if(DEBUG_LEVEL>=2) {
            printf("Neighborhood mapped.\n");
        }
    #endif
}

void GD::Validate(
        std::vector<std::vector<int>>& Neighbors, std::vector<bool> &Visited, std::vector<int> &S, int root,
        std::vector<int> &ValList)
{
    for (unsigned long i=0; i<S.size(); i++) {  // S is the selected children of the last level (or root)
        for (unsigned long j=0; j<Neighbors[S[i]].size(); j++) {  // j is the label of possible children, which can't be greater than root
            if (Neighbors[S[i]][j]>root && !Visited[Neighbors[S[i]][j]]) {  // if not visited and it's in the neighborhood of S[i]
                Visited[Neighbors[S[i]][j]] = true;  // mark as visited
                ValList.push_back(Neighbors[S[i]][j]);  // include as valid children
            }
        }
    }
}

std::vector<int> GD::genComb(unsigned long i, std::vector<int> ValList, std::vector<int> &Index, unsigned long level) {
    if (!Index.empty()) {
        while (ValList.size()!=i) {
            ValList.erase(ValList.begin() + Index[level++]);
        }
    }
    return ValList;
}

void GD::updateIndex(std::vector<int> &Index, unsigned long n) {
    if (Index.empty()) {
        return;  // nothing to update if Index is empty
    }
    unsigned long i = Index.size() - 1;  // go to last position of Index
    Index[i]++;  // increase it
    for (i; i>0; i--) {  // ascend Index
        if (Index[i]>n) {
            Index[i-1]++;  // if the current position's value overpasses n "send one" to the position before
        } else {
            break;  // if not then there is no need to keep ascending
        }
    }
    if (Index[i]>n) {
        Index[i] = 0;  // if Index[i]>n it means i==0 so no way to "send one" as there is no position before, so cycle back to 0
    }
    if (Index[++i]>Index[i-1]) {
        for (i; i<Index.size(); i++) {
            Index[i] = Index[i-1];
        }
    }
}

void GD::Classify(TNGraph &G, GD::GraphList *kSubgraphs, int n, optionblk &options, int m, set *dnwork) {
    kSubgraphs->cursor = kSubgraphs->ini->next;

    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(graph, canon, canon_sz);
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    statsblk stats;

    DYNALLOC2(graph, g, g_sz, m, n, "malloc");
    DYNALLOC2(graph, canon, canon_sz, m, n, "malloc");
    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n, "malloc");

    while (kSubgraphs->cursor) {
        EMPTYGRAPH(g,m,n);
        for (int i=0; i<n; i++) {
            int a = kSubgraphs->cursor->vertices.at((size_t) i);
            for (int j=0; j<n; j++) {
                int b = kSubgraphs->cursor->vertices.at((size_t) j);
                if (G.GetNI(a).IsOutNId(b)) {
                    ADDONEARC(g,i,j,m);
                }
            }
        }

        //densenauty(g,lab,ptn,orbits,&options,&stats,m,n,canon);
        nauty(g, lab, ptn, NULL, orbits, &options, &stats, dnwork, 2*60*m, m, n, canon);
        for (int i=0; i<m*n; i++) {
            kSubgraphs->cursor->label.insert(canon[i]);
        }
        #ifdef DEBUG
            if (DEBUG_LEVEL>=3) {
                for (auto&& it : kSubgraphs->cursor->label) {
                    printf("%lu ", it);
                }
                printf("= ");
                for (unsigned long i=0; i<n; i++) {
                    printf("%d ", kSubgraphs->cursor->vertices.at(i));
                }
                printf("\n");
            }
        #endif

        kSubgraphs->cursor = kSubgraphs->cursor->next;
    }

    DYNFREE(g,g_sz);
    DYNFREE(canon,canon_sz);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    //nauty_freedyn();
    //nautil_freedyn();
    //naugraph_freedyn();
}

TNGraph* GD::Randomize(TNGraph &G) {
    TNGraph *R = new TNGraph(G);
    TNGraph::TNodeI a, b, c, d;
    for (int i=0; i<3; i++) {
        for (a = R->BegNI(); a < R->EndNI(); a++) {
        //for(unsigned long j=0; j<R->GetNodes(); j++) {
            //a = R->GetNI(j);
            if (a.GetOutDeg()==0) {
                continue;
            }
            do {
                b = R->GetRndNI();
            } while (a.GetId()==b.GetId() || b.GetOutDeg()==0);
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

void GD::GetFrequencies(GD::GraphList* kSubgraphs, std::map<std::multiset<unsigned long>, int> &Frequencies) {
    kSubgraphs->cursor = kSubgraphs->ini->next;
    while (kSubgraphs->cursor) {
        Frequencies[kSubgraphs->cursor->label]++;
        kSubgraphs->cursor = kSubgraphs->cursor->next;
    }

    #ifdef DEBUG
        if (DEBUG_LEVEL>=3) {
            std::map<std::multiset<unsigned long>, int>::iterator it;
            for (it = Frequencies.begin(); it != Frequencies.end(); it++) {
                for (auto&& el : it->first) {
                    printf("%lu ", el);
                }
                printf(": %d\n", it->second);
            }
        }
    #endif
}

void GD::DiscoverMotifs(
        std::vector<std::map<std::multiset<unsigned long>, int>> &FVector,
        std::vector<std::multiset<unsigned long>> &Motifs, std::vector<unsigned long> &IDs, int motif_size,
        std::string destination, TNGraph &G, GD::GraphList *kSubgraphs, GD::MetaObject *metaObj, bool expand)
{
    std::map<std::multiset<unsigned long>, int>::iterator i;
    std::stringstream ss;
    ss << destination << "statistic_measures.txt";

    std::string s = ss.str();
    const char* path = s.c_str();
    std::ofstream TXT;
    TXT.open(path);

    std::ofstream ETXT;
    if (expand) {
        std::stringstream ess;
        ess << destination << "expanded_adjacency_matrixes.txt";
        std::string es = ess.str();
        const char* epath = es.c_str();
        ETXT.open(epath);
    }

    TXT << "\t\t\t\t\t\t\t[Original Network]\t[\t\t\t\t\t\tRandom Networks\t\t\t\t\t\t\t]\n"
        << "ID\t\tAdjacency Matrix\tFrequency\t\tMean Frequency\t\tStandard Deviation\t\tZ-Score\t\t\tP-value\n";
    
    for (i = FVector[0].begin(); i != FVector[0].end(); i++) {
        double MeanFrequency = 0.0;
        double StandardDeviation = 0.0;
        double Zscore;
        double Pvalue = 0.0;
        for (int index=1; index<FVector.size(); index++) {
            std::map<std::multiset<unsigned long>, int>::iterator j = FVector[index].find(i->first);
            if (j!=FVector[index].end()) {
                if (j->second > i->second) {
                    Pvalue++;
                }
                MeanFrequency += j->second;
            }
        }
        MeanFrequency /= FVector.size()-1;
        for (size_t index=1; index<FVector.size(); index++) {
            std::map<std::multiset<unsigned long>, int>::iterator j = FVector[index].find(i->first);
            if (j!=FVector[index].end()) {
                StandardDeviation += (MeanFrequency-j->second)*(MeanFrequency-j->second);
            }
        }
        StandardDeviation /= FVector.size()-1;
        StandardDeviation = sqrt(StandardDeviation);
        Zscore = (i->second - MeanFrequency)/StandardDeviation;
        Pvalue /= FVector.size()-1;

        std::vector<int> *vertices = nullptr;
        kSubgraphs->cursor = kSubgraphs->ini->next;
        while (kSubgraphs->cursor) {
            if (kSubgraphs->cursor->label==i->first) {
                vertices = &kSubgraphs->cursor->vertices;
                break;
            }
            kSubgraphs->cursor = kSubgraphs->cursor->next;
        }

        if (vertices == nullptr) {
            return;
        }

        unsigned long ID = 0;
        unsigned long maxpow = (unsigned long) motif_size * (unsigned long) motif_size - 1;
        bool adjMatrix[motif_size][motif_size];
        std::vector<bool> expandedAdjMatrix[motif_size];

        for (int j=0; j<motif_size; j++) {
            for (int k=0; k<motif_size; k++) {
                adjMatrix[j][k] = false;
                if (G.GetNI(vertices->at((size_t) motif_size-1-j)).IsOutNId(vertices->at((size_t) motif_size-1-k))) {
                    adjMatrix[j][k] = true;
                    ID += (long) (pow(2, (maxpow - (j*motif_size+k))));
                }
            }
            if (expand) {
                std::vector<int> aux(*(metaObj->metaMap->find(vertices->at((size_t) j))->second));
                for (std::vector<int>::iterator it0 = aux.begin(); it0 != aux.end(); ++it0) {
                    for (std::vector<int>::iterator it1 = aux.begin(); it1 != aux.end(); ++it1) {
                        if (*it0 != *it1 && metaObj->G->GetNI(*it0).IsOutNId(*it1)) {
                            expandedAdjMatrix[j].push_back(1);
                        } else {
                            expandedAdjMatrix[j].push_back(0);
                        }
                    }
                }
            }
        }

        if (i->second > 4 && Pvalue < 0.01 && Zscore > 1) {
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
        if (ID<999) {
            TXT << "\t";
        }
        
        for (int j=0; j<motif_size; j++) {
            if (adjMatrix[0][j]) {
                TXT << "1 ";
            } else {
                TXT << "0 ";
            }
        }

        if (expand) {
            ETXT << "ID " << ID << "\n";

            for (int h=0; h<motif_size; h++) {
                int expandedMatrixSize = (int) sqrt(expandedAdjMatrix[h].size());

                for (int j = 0; j < expandedMatrixSize; j++) {
                    ETXT << "\t\t";
                    for (int k = 0; k < expandedMatrixSize; k++) {
                        if (expandedAdjMatrix[h][j * expandedMatrixSize + k]) {
                            ETXT << "1 ";
                        } else {
                            ETXT << "0 ";
                        }
                    }
                    ETXT << "\n";
                }
                
                ETXT << "\n";
            }
            
            ETXT << "\n";
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

        for (int j=1; j<motif_size; j++) {
            TXT << "\t\t";
            for (int k=0; k<motif_size; k++) {
                if (adjMatrix[j][k]) {
                    TXT << "1 ";
                } else {
                    TXT << "0 ";
                }
            }
            TXT << "\n";
        }

        TXT << "\n";

        #ifdef DEBUG
            if (DEBUG_LEVEL>=0) {
                for (auto&& j : i->first) {
                    printf("%lu ", j);
                }
                printf(
                        ": Original Frequency = %d : Mean Frequency = %f : Standard Deviation = %f : Zscore = %f"
                                " : Pvalue = %f\n", i->second, MeanFrequency, StandardDeviation, Zscore, Pvalue);
            }
        #endif
    }

    TXT.close();

    if (expand) {
        ETXT.close();
    }
    
    #ifdef DEBUG
        if (DEBUG_LEVEL>=0) {
            for (int i=0; i<Motifs.size(); i++) {
                printf("\nMotif %d label: ", i+1);
                for (auto&& j : Motifs[i]) {
                    printf("%lu ", j);
                }
                printf("\n");
            }
        }
    #endif
}

TNGraph* GD::ConcatMotifs(
        TNGraph &G, std::vector<std::multiset<unsigned long>> &Motifs, GD::GraphList *kSubgraphs,
        GD::MetaObject *metaObj)
{
    #ifdef DEBUG
        if (DEBUG_LEVEL>=1) {
            printf("Concatenation Started\nInserting motifs as nodes...\n");
        }
    #endif

    std::map<int, std::vector<int>*> *metaMap = new std::map<int, std::vector<int>*>;

    TNGraph *H = new TNGraph(G);
    for (int m=0; m<Motifs.size(); m++) {
        kSubgraphs->cursor = kSubgraphs->ini->next;
        while (kSubgraphs->cursor) {
            if (kSubgraphs->cursor->label==Motifs[m]) {
                int mNodeId = H->AddNode(-1);
                std::vector<int> *aux = new std::vector<int>;
                for (int i=0; i<Motifs[m].size(); i++) {
                    std::vector<int> *tmp;
                    tmp = metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) i))->second;
                    aux->insert(aux->end(), tmp->begin(), tmp->end());
                }
                metaMap->emplace(mNodeId, aux);
                #ifdef DEBUG
                    if (DEBUG_LEVEL>=2) {
                        printf("Node %lu = (", mNodeId);
                        int j;
                        for (j=0; j<Motifs[m].size()-1; j++) {
                            printf("%lu, ", kSubgraphs->cursor->vertices.at(j));
                        }
                        printf("%lu) added\n", kSubgraphs->cursor->vertices.at(j));
                    }
                #endif
                for (int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
                    TNGraph::TNodeI u = H->GetNI(kSubgraphs->cursor->vertices.at((size_t) i));
                    for (int j=0; j<u.GetOutDeg(); j++) {
                        if (mNodeId != u.GetOutNId(j)) {
                            H->AddEdge(mNodeId, u.GetOutNId(j));
                        }
                    }
                    for (int j=0; j<u.GetInDeg(); j++) {
                        if (mNodeId != u.GetInNId(j)) {
                            H->AddEdge(u.GetInNId(j), mNodeId);
                        }
                    }
                }
                #ifdef DEBUG
                    if (DEBUG_LEVEL>=2) {
                        TNGraph::TNodeI mNodeNI = H->GetNI(mNodeId);
                        printf("\tDestinations: ");
                        for (int i=0; i<mNodeNI.GetOutDeg(); i++) {
                            printf("%d ", mNodeNI.GetOutNId(i));
                        }
                        printf("\n\tIncomings: ");
                        for (int i=0; i<mNodeNI.GetInDeg(); i++) {
                            printf("%d ", mNodeNI.GetInNId(i));
                        }
                        printf("\n");
                    }
                #endif
            }
            kSubgraphs->cursor = kSubgraphs->cursor->next;
        }
    }

    // FIXME Memory leak!?
    metaObj->metaMap = metaMap;

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1) {
            printf("Removing vertices from motifs...\n");
        }
    #endif

    for (int m=0; m<Motifs.size(); m++) {
        kSubgraphs->cursor = kSubgraphs->ini->next;
        while (kSubgraphs->cursor) {
            if (kSubgraphs->cursor->label==Motifs[m]) {
                for (int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
                    if (H->IsNode(kSubgraphs->cursor->vertices.at((size_t) i))) {
                        H->DelNode(kSubgraphs->cursor->vertices.at((size_t) i));
                        #ifdef DEBUG
                            if (DEBUG_LEVEL>=2) {
                                printf("Node %lu removed\n", kSubgraphs->cursor->vertices.at(i));
                            }
                        #endif
                    }

                }
            }
            kSubgraphs->cursor = kSubgraphs->cursor->next;
        }
    }

    #ifdef DEBUG
        if (DEBUG_LEVEL>=1) {
            printf("Defragmenting motifs graph...\n");
        }
    #endif
    H->Defrag();

    return H;
}

void GD::ExportGDF(
        TNGraph &G, std::vector<std::multiset<unsigned long>> *Motifs, std::vector<unsigned long> *IDs,
        GD::GraphList *kSubgraphs, std::string destination, GD::MetaObject *metaObj)
{
    int size;
    if(Motifs!=NULL) {
        size = (int) Motifs->size();
    } else {
        size = 1;
    }

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1) {
            printf("Exporting to GDF...\n");
        }
    #endif

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1) {
            printf("Marking vertices to be colored as motifs...\n");
        }
    #endif
    for (int m=0; m<size; m++) {
        std::vector<bool> marker((size_t) metaObj->G->GetNodes(), false);  // vector to mark which nodes will be colored as being part of a motif
        if (Motifs!=NULL) {
            kSubgraphs->cursor = kSubgraphs->ini->next;
            while (kSubgraphs->cursor) {
                if (kSubgraphs->cursor->label==Motifs->at((size_t) m)) {
                    for (int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
                        std::vector<int> *tmp = metaObj->metaMap->find(
                                kSubgraphs->cursor->vertices.at((size_t) i))->second;
                        for (std::vector<int>::iterator it = tmp->begin(); it != tmp->end(); ++it)
                            marker[*it] = true;
                    }
                }
                kSubgraphs->cursor = kSubgraphs->cursor->next;
            }
        }
        std::stringstream ss;
        if(Motifs!=NULL) {
            ss << destination << "motif_ID" << IDs->at((size_t) m) << "_highlighted.gdf";
        } else {
            ss << destination << "concatenated_motifs.gdf";
        }
        std::string s = ss.str();
        const char* path = s.c_str();
        std::ofstream GDF;
        GDF.open(path);
        GDF << "nodedef>name VARCHAR,color VARCHAR\n";

        /*for(unsigned long i=0; i<marker.size(); i++) {
            GDF << i << ",";
            if(marker[i]) GDF << "'255,0,0'\n";
            else GDF << "'0,0,0'\n";
        }*/
        for (TNGraph::TNodeI NI = metaObj->G->BegNI(); NI < metaObj->G->EndNI(); NI++) {
            int id = NI.GetId();
            GDF << id << ",";
            if (Motifs!=NULL && marker[id]) {
                GDF << "'255,0,0'\n";
            } else {
                GDF << "'0,0,0'\n";
            }
        }
        GDF << "edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,color VARCHAR\n";

        #ifdef DEBUG
            if (DEBUG_LEVEL>=1) {
                printf("Marking edges to be colored as motifs...\n");
            }
        #endif
        TNGraph T(*(metaObj->G));
        if (Motifs!=NULL) {
            kSubgraphs->cursor = kSubgraphs->ini->next;
            while (kSubgraphs->cursor) {
                //std::cerr << "\n";
                if (kSubgraphs->cursor->label==Motifs->at((size_t) m)) {
                    for (int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
                        for (int j=0; j<kSubgraphs->cursor->vertices.size(); j++) {
                            //std::cerr << "···";
                            if ((G.GetNI(kSubgraphs->cursor->vertices.at((size_t) i)))
                                    .IsOutNId(kSubgraphs->cursor->vertices.at((size_t) j)))
                            {
                                std::vector<int> tmp(*(metaObj->metaMap->find(
                                        kSubgraphs->cursor->vertices.at((size_t) i))->second));
                                std::vector<int> aux(*(metaObj->metaMap->find(
                                        kSubgraphs->cursor->vertices.at((size_t) j))->second));
                                aux.insert(aux.end(), tmp.begin(), tmp.end());
                                for (std::vector<int>::iterator it0 = aux.begin(); it0 != aux.end(); ++it0) {
                                    for (std::vector<int>::iterator it1 = aux.begin(); it1 != aux.end(); ++it1) {
                                        if (*it0!=*it1 && T.GetNI(*it0).IsOutNId(*it1)) {
                                            GDF << *it0 << "," << *it1 << ",true,'255,0,0'\n";
                                            T.DelEdge(*it0, *it1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                //std::cerr << "\n";
                kSubgraphs->cursor = kSubgraphs->cursor->next;
            }
            //printf("\n\nHMM\n\n");
        }

        #ifdef DEBUG
            if (DEBUG_LEVEL>=1) {
                printf("Inserting remaining edges...\n");
            }
        #endif
        /*for(TNGraph::TEdgeI EI = T.BegEI(); EI < T.EndEI(); EI++) {
            if(EI.GetSrcNId()==EI.GetSrcNId()) continue;
            std::cerr << EI.GetSrcNId() << " -> " << EI.GetSrcNId() << "\n";
            std::vector<unsigned long> aux(*(metaObj->metaMap->find(EI.GetSrcNId())->second));
            std::vector<unsigned long> tmp(*(metaObj->metaMap->find(EI.GetDstNId())->second));
            aux.insert(aux.end(), tmp.begin(), tmp.end());
            for(std::vector<unsigned long>::iterator it0 = aux.begin(); it0 != aux.end(); ++it0) {
                for(std::vector<unsigned long>::iterator it1 = aux.begin(); it1 != aux.end(); ++it1) {
                    if(*it0!=*it1 && T.GetNI(*it0).IsOutNId(*it1)) {
                        GDF << *it0 << "," << *it1 << ",true,'0,0,0'\n";
                    }
                }
            }
        }*/
        for (TNGraph::TEdgeI EI = T.BegEI(); EI < T.EndEI(); EI++)
            GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'0,0,0'\n";

        #ifdef DEBUG
            if (DEBUG_LEVEL>=1) {
                printf("Saving GDF...\n");
            }
        #endif
        GDF.close();
    }
    //printf("metaObj->G nodes: %d\n", metaObj->G->GetNodes());
}

void GD::PrintMotifs(
        TNGraph &G, std::vector<std::multiset<unsigned long>> &Motifs, std::vector<unsigned long> &IDs,
        GD::GraphList *kSubgraphs, std::string destination, GD::MetaObject *metaObj, bool expand)
{
    for (int m=0; m<Motifs.size(); m++) {

        PNGraph motif = TNGraph::New();
        
        PNGraph expanded;
        if (expand) {
            expanded = TNGraph::New();
        }
        
        kSubgraphs->cursor = kSubgraphs->ini->next;
        while (kSubgraphs->cursor->label!=Motifs[m]) {
            kSubgraphs->cursor = kSubgraphs->cursor->next;
        }

        std::vector<int> *v = &kSubgraphs->cursor->vertices;
        
        for (int i=0; i<v->size(); i++) {
            motif->AddNode(v->at((size_t) i));
            
            if (expand) {
                std::vector<int> *aux = metaObj->metaMap->find(v->at((size_t) i))->second;
                for (std::vector<int>::iterator it = aux->begin(); it != aux->end(); ++it) {
                    if (expanded->IsNode(*it)) {
                        continue;
                    }
                    expanded->AddNode(*it);
                }
            }
        }
        
        for (int i=0; i<v->size(); i++) {
            for (int j=0; j<v->size(); j++) {
                
                if (i!=j && G.GetNI(v->at((size_t) i)).IsOutNId(v->at((size_t) j))) {
                    motif->AddEdge(v->at((size_t) i), v->at((size_t) j));
                }
                
                if (expand) {
                    std::vector<int> aux(*(metaObj->metaMap->find(v->at((size_t) i))->second));
                    std::vector<int> tmp(*(metaObj->metaMap->find(v->at((size_t) j))->second));
                    aux.insert(aux.end(), tmp.begin(), tmp.end());
                    for (std::vector<int>::iterator it0 = aux.begin(); it0 != aux.end(); ++it0) {
                        for (std::vector<int>::iterator it1 = aux.begin(); it1 != aux.end(); ++it1) {
                            if (*it0 != *it1 && metaObj->G->GetNI(*it0).IsOutNId(*it1)) {
                                expanded->AddEdge(*it0, *it1);
                            }
                        }
                    }
                }
            }
        }
        
        std::stringstream ss;
        ss << destination << "motif_ID" << IDs[m] << ".png";
        std::string s = ss.str();
        const char* path = s.c_str();
        TSnap::DrawGViz<PNGraph>(motif, gvlDot, path, "", NULL);
        
        if (expand) {
            std::stringstream ssE;
            ssE << destination << "motif_ID" << IDs[m] << "-expanded.png";
            std::string sE = ssE.str();
            const char* pathE = sE.c_str();
            TSnap::DrawGViz<PNGraph>(expanded, gvlDot, pathE, "", NULL);
        }
    }
}
