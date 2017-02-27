#include <GraphData.hpp>
#include <GetTimeMs64.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

GD::MetaObject::MetaObject(PNGraph &G) {
    this->G = &*G;
    metaMap = new std::map<int, std::vector<int>*>;
    for (TNGraph::TNodeI NI = this->G->BegNI(); NI < this->G->EndNI(); NI++) {
        std::vector<int> *v = new std::vector<int>;
        v->push_back(NI.GetId());
        metaMap->emplace(NI.GetId(),v);
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

void GD::Enumerate(TNGraph &G, int k, std::shared_ptr<graphmap> &kSubgraphs, optionblk &options, int m, set *dnwork,
                   bool original, int pos)
{
    uint64 t0 = GetTimeMs64();
    std::vector<std::vector<int>> Neighbors((size_t) G.GetMxNId());
    uint64 t1 = GetTimeMs64();
    // map neighbors for all vertices to avoid repeatly discovering the neighborhood of a vertex
    GD::mapNeighbors(G, Neighbors);
    uint64 t2 = GetTimeMs64();
    #ifdef DEBUG
        if(DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1) printf("GetNodes: %lu ms\nMap Neighbors: %lu ms\n", t1-t0, t2-t1);
    #endif
    // iterate through all vertices of G
    //for(unsigned long i=0; G.IsNode(i); i++) {
    for (TNGraph::TNodeI NI = G.BegNI(); NI < G.EndNI(); NI++) {
        int i = NI.GetId();
        #ifdef DEBUG
            if(DEBUG_LEVEL>=2) printf("Root: %d\n", i);
        #endif
        // list to keep track of visited vertices
        std::vector<bool> Visited((size_t) G.GetMxNId());
        // mark the current root as visited
        Visited[i] = true;
        // vertices selected in current level
        std::vector<int> S(1);
        // only the root is selected in the first level
        S[0] = i;
        // vector to store the vertices of the subgraph
        std::vector<int> subgraph(S);
        // Explore all possible combinations that generete different trees
        GD::Explore(G, Neighbors, i, Visited, S, subgraph, kSubgraphs, k-1, k, options, m, dnwork, original, pos);
        // unmark current root as visited at the end
        Visited[i] = false;
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

void GD::Explore(TNGraph &G, std::vector<std::vector<int>> &Neighbors, int root, std::vector<bool> &Visited,
                 std::vector<int> &S, std::vector<int> &subgraph, std::shared_ptr<graphmap> &kSubgraphs, int Remainder,
                 int motif_size, optionblk &options, int m, set *dnwork, bool original, int pos)
{
    #ifdef DEBUG
        LEVEL++;
    #endif
    uint64 t0 = GetTimeMs64();
    // if there are no more vertices to select the k-subgraph is formed
    if(Remainder==0) {
        int *v = (*kSubgraphs)[Classify(G, subgraph, motif_size, options, m, dnwork)];

        if (v == nullptr) {
            if (!original) {
                return;
            }

            v = new int[pos]();
        }

        original ? v[0]++ : v[pos]++;

        #ifdef DEBUG
            COUNTER++;                                    // count one more k-subgraph discovered
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
    // create a vector to store valid children at this level of the tree
    std::vector<int> ValList;
    // fill the vector with the actual children
    GD::Validate(Neighbors, Visited, S, root, ValList);
    uint64 t1 = GetTimeMs64();
    // get the number of vertices to select at this level
    int n = (int) std::min(ValList.size(), (const unsigned long &) Remainder);
    #ifdef DEBUG
        if(DEBUG_LEVEL>=4 || DEBUG_LEVEL==-1) printf("\tExplore > Validade: %lu ms\n", t1-t0);
        if(DEBUG_LEVEL>=2) {
            for(unsigned long p=0; p<LEVEL; p++) printf("\t");
            printf("Valid Children: ");
            for(unsigned long p=0; p<ValList.size(); p++) printf("%lu ", ValList[p]);
            printf("\n");
        }
    #endif
    // iterate from selecting a single vertex to selecting n vertices
    for(int i=1; i<=n; i++) {
        // create a vector to keep track of which vertices to delete from ValList in order to generate a combination
        std::vector<int> Index(ValList.size()-i, 0);
        // generate the initial combination
        std::vector<int> C = GD::genComb(i, ValList, Index, 0);
        do {
            #ifdef DEBUG
                if(DEBUG_LEVEL>=2) {
                    for(unsigned long p=0; p<LEVEL; p++) printf("\t");
                    printf("Selection: ");
                    for(unsigned long p=0; p<C.size(); p++) printf("%lu ", C[p]);
                    printf("\n");
                }
            #endif
            // compose subgraph with the current combination
            subgraph.insert(subgraph.end(), C.begin(), C.end());
            // recursive call to advance to next level of the implicit tree
            GD::Explore(G, Neighbors, root, Visited, C, subgraph, kSubgraphs, Remainder-i, motif_size, options,
                        m, dnwork, original, pos);
            uint64 t2 = GetTimeMs64();
            // update Index to generate the next combination
            GD::updateIndex(Index, i);
            // generate next combination
            C = GD::genComb(i, ValList, Index, 0);
            uint64 t3 = GetTimeMs64();
            #ifdef DEBUG
                if(DEBUG_LEVEL>=4 || DEBUG_LEVEL==-1) printf("Get next combination time: %lu ms\n", t3-t2);
            #endif
            // clear old combination from the subgraph "constructor"
            subgraph.erase(subgraph.end() - C.size(), subgraph.end());
        } while(!Index.empty() && Index.back()!=0);
    }
    // unmark selected vertices at current level
    for(size_t i=0; i<ValList.size(); i++) Visited.at((size_t) ValList[i]) = false;
    #ifdef DEBUG
        LEVEL--;
    #endif
}


std::multiset<unsigned long>& GD::Classify(TNGraph &G, std::vector<int> &subgraph, int n, optionblk &options,
                                           int m, set *dnwork)
{
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(graph,canon,canon_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);

    statsblk stats;

    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,canon,canon_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

    EMPTYGRAPH(g,m,n);
    for(int i=0; i<n; i++) {
        int a = subgraph.at((size_t) i);
        for(int j=0; j<n; j++) {
            int b = subgraph.at((size_t) j);
            if(G.GetNI(a).IsOutNId(b)) {
                ADDONEARC(g,i,j,m);
            }
        }
    }

    //densenauty(g,lab,ptn,orbits,&options,&stats,m,n,canon);
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,dnwork,2*60*m,m,n,canon);

    std::multiset<unsigned long> *label = new std::multiset<unsigned long>;

    for(int i=0; i<m*n; i++)
        label->insert(canon[i]);
#ifdef DEBUG
        if(DEBUG_LEVEL>=3) {
                for(unsigned long i=0; i<kSubgraphs->cursor->label.size(); i++) {
                    printf("%lu ", kSubgraphs->cursor->label.at(i));
                }
                printf("= ");
                for(unsigned long i=0; i<n; i++) {
                    printf("%lu ", kSubgraphs->cursor->vertices.at(i));
                }
                printf("\n");
            }
#endif


    DYNFREE(g,g_sz);
    DYNFREE(canon,canon_sz);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    //nauty_freedyn();
    //nautil_freedyn();
    //naugraph_freedyn();

    return *label;
}

void GD::mapNeighbors(TNGraph &G, std::vector<std::vector<int>> &Neighbors) {
    #ifdef DEBUG
        if(DEBUG_LEVEL>=2) printf("Mapping neighborhood...\n");
    #endif
    for (TNGraph::TNodeI NI1 = G.BegNI(); NI1 < G.EndNI(); NI1++)
        for (TNGraph::TNodeI NI2 = G.BegNI(); NI2 < G.EndNI(); NI2++)
            if(NI1.IsNbrNId(NI2.GetId()))
                Neighbors[NI1.GetId()].push_back(NI2.GetId());
    /*for(unsigned long i=0; i<G.GetNodes(); i++)
        for(unsigned long j=0; j<G.GetNodes(); j++)
            if((G.GetNI(i)).IsNbrNId(j))
                Neighbors[i].push_back(j);*/
    #ifdef DEBUG
        if(DEBUG_LEVEL>=2) printf("Neighborhood mapped.\n");
    #endif
}

void GD::Validate(std::vector<std::vector<int>>& Neighbors, std::vector<bool> &Visited, std::vector<int> &S,
                  int root, std::vector<int> &ValList)
{
    // S is the selected children of the last level (or root)
    for(unsigned long i=0; i<S.size(); i++) {
        // j is the label of possible children, which can't be greater than root
        for(unsigned long j=0; j<Neighbors[S[i]].size(); j++) {
            // if not visited and it's in the neighborhood of S[i]
            if(Neighbors[S[i]][j]>root && !Visited[Neighbors[S[i]][j]]) {
                // mark as visited
                Visited[Neighbors[S[i]][j]] = true;
                // include as valid children
                ValList.push_back(Neighbors[S[i]][j]);
            }
        }
    }
}

std::vector<int> GD::genComb(int i, std::vector<int> ValList, std::vector<int> &Index, int level) {
    if(!Index.empty())
        while(ValList.size()!=i)
            ValList.erase(ValList.begin() + Index[level++]);
    return ValList;
}

void GD::updateIndex(std::vector<int> &Index, int n) {
    if(Index.empty()) return;

    int i = (int) (Index.size() - 1);
    Index[i]++;

    for(i; i>0; i--) {
        // if the current position's value overpasses n "send one" to the position before
        if(Index[i]>n) Index[i-1]++;
        // if not then there is no need to keep ascending
        else break;
    }

    // if Index[i]>n it means i==0 so no way to "send one" as there is no position before, so cycle back to 0
    if(Index[i]>n) Index[i] = 0;

    if(Index[++i]>Index[i-1])
        for(i; i<Index.size(); i++)
            Index[i] = Index[i-1];
}

TNGraph* GD::Randomize(TNGraph &G) {
    TNGraph *R = new TNGraph(G);
    TNGraph::TNodeI a, b, c, d;
    for(int i=0; i<3; i++) {
        for (a = R->BegNI(); a < R->EndNI(); a++) {
        //for(unsigned long j=0; j<R->GetNodes(); j++) {
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

void GD::DiscoverMotifs(std::vector<std::multiset<unsigned long>> &Motifs, std::vector<unsigned long> &IDs,
                        int motif_size, std::string destination, TNGraph &G,
                        std::shared_ptr<graphmap> &kSubgraphs, int r)
{
    //std::map<std::multiset<unsigned long>, int>::iterator i;
    std::stringstream ss;
    ss << destination << "statistic_measures.txt";

    std::string s = ss.str();
    const char* path = s.c_str();
    std::ofstream TXT;
    TXT.open(path);

    TXT << "\t\t\t\t\t\t\t[Original Network]\t[\t\t\t\t\t\tRandom Networks\t\t\t\t\t\t\t]\n"
        << "ID\t\tAdjacency Matrix\tFrequency\t\tMean Frequency\t\tStandard Deviation\t\tZ-Score\t\t\tP-value\n";

    double MeanFrequency, StandardDeviation, Zscore, Pvalue;

    for (auto&& it : *kSubgraphs) {
        MeanFrequency = 0.0;
        StandardDeviation = 0.0;
        Zscore = 0.0;
        Pvalue = 0.0;

        for (int i = 1; i <= r; i++) {
            if (it.second[0] > it.second[i]) {
                Pvalue++;
            }
            MeanFrequency += it.second[i];
        }

        Pvalue /= r;
        MeanFrequency /= r;

        for (int i = 1; i <= r; i++) {
            double diff = MeanFrequency - it.second[i];
            StandardDeviation += diff * diff;
        }

        StandardDeviation = sqrt(StandardDeviation / r);
        Zscore = (it.second[0] - MeanFrequency) / StandardDeviation;

        std::vector<int> *vertices = nullptr;
    }
        /*if (vertices == nullptr)
            return;

        unsigned long ID = 0;
        unsigned long maxpow = (unsigned long) motif_size * (unsigned long) motif_size - 1;
        bool adjMatrix[motif_size][motif_size];

        for(int j=0; j<motif_size; j++) {
            for(int k=0; k<motif_size; k++) {
                adjMatrix[j][k] = false;
                if(G.GetNI(vertices->at((size_t) motif_size-1-j)).IsOutNId(vertices->at((size_t) motif_size-1-k))) {
                    adjMatrix[j][k] = true;
                    ID += (long) (pow(2, (maxpow - (j*motif_size+k))));
                }
            }
        }

        if(i->second > 4 && Pvalue < 0.01 && Zscore > 1) {
            Motifs.push_back(i->first);
            IDs.push_back(ID);
        }

        *//*for(int j=motif_size-1; j>=0; j--) {
            if(G.GetNI(vertices->at(motif_size-1)).IsOutNId(vertices->at(j)))
                TXT << "1 ";
            else
                TXT << "0 ";
        }*//*

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

        *//*for(int j=motif_size-2; j>=0; j--) {
            for(int k=motif_size-1; k>=0; k--) {
                if(G.GetNI(vertices->at(j)).IsOutNId(vertices->at(k)))
                    TXT << "1 ";
                else
                    TXT << "0 ";
            }
            TXT << "\n";
        }*//*

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
                for(unsigned long j=0; j<i->first.size(); j++) {
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
    #endif*/
}

TNGraph* GD::ConcatMotifs(TNGraph &G, std::vector<std::multiset<unsigned long>> &Motifs,
                          std::unordered_map<std::multiset<unsigned long>, int> &kSubgraphs, GD::MetaObject *metaObj)
{
    /*#ifdef DEBUG
        if(DEBUG_LEVEL>=1) {
            printf("Concatenation Started\nInserting motifs as nodes...\n");
        }
    #endif

    std::map<int, std::vector<int>*> *metaMap = new std::map<int, std::vector<int>*>;

    TNGraph *H = new TNGraph(G);
    for(int m=0; m<Motifs.size(); m++) {
        kSubgraphs->cursor = kSubgraphs->ini->next;
        while(kSubgraphs->cursor) {
            if(kSubgraphs->cursor->label==Motifs[m]) {
                int mNodeId = H->AddNode(-1);
                std::vector<int> *aux = new std::vector<int>;
                for(int i=0; i<Motifs[m].size(); i++) {
                    std::vector<int> *tmp;
                    tmp = metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) i))->second;
                    aux->insert(aux->end(), tmp->begin(), tmp->end());
                }
                metaMap->emplace(mNodeId, aux);
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
                    TNGraph::TNodeI u = H->GetNI(kSubgraphs->cursor->vertices.at((size_t) i));
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

    // FIXME Memory leak!?
    metaObj->metaMap = metaMap;

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
                    if(H->IsNode(kSubgraphs->cursor->vertices.at((size_t) i))) {
                        H->DelNode(kSubgraphs->cursor->vertices.at((size_t) i));
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

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1) {
            printf("Defragmenting motifs graph...\n");
        }
    #endif
    H->Defrag();

    return H;*/
}

void GD::ExportGDF(TNGraph &G, std::vector<std::multiset<unsigned long>> *Motifs, std::vector<unsigned long> *IDs,
                   std::unordered_map<std::multiset<unsigned long>, int> *kSubgraphs, std::string destination,
                   GD::MetaObject *metaObj)
{
    /*int size;
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
    for(int m=0; m<size; m++) {
        std::vector<bool> marker((size_t) metaObj->G->GetNodes(), false);    // vector to mark which nodes will be colored as being part of a motif
        if(Motifs!=NULL) {
            kSubgraphs->cursor = kSubgraphs->ini->next;
            while(kSubgraphs->cursor) {
                if(kSubgraphs->cursor->label==Motifs->at((size_t) m)) {
                    for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
                        std::vector<int> *tmp = metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) i))->second;
                        for(std::vector<int>::iterator it = tmp->begin(); it != tmp->end(); ++it)
                            marker[*it] = true;
                    }
                }
                kSubgraphs->cursor = kSubgraphs->cursor->next;
            }
        }
        std::stringstream ss;
        if(Motifs!=NULL)
            ss << destination << "motif_ID" << IDs->at((size_t) m) << "_highlighted.gdf";
        else
            ss << destination << "concatenated_motifs.gdf";
        std::string s = ss.str();
        const char* path = s.c_str();
        std::ofstream GDF;
        GDF.open(path);
        GDF << "nodedef>name VARCHAR,color VARCHAR\n";

        *//*for(unsigned long i=0; i<marker.size(); i++) {
            GDF << i << ",";
            if(marker[i]) GDF << "'255,0,0'\n";
            else GDF << "'0,0,0'\n";
        }*//*
        for (TNGraph::TNodeI NI = metaObj->G->BegNI(); NI < metaObj->G->EndNI(); NI++) {
            int id = NI.GetId();
            GDF << id << ",";
            if(Motifs!=NULL && marker[id]) GDF << "'255,0,0'\n";
            else GDF << "'0,0,0'\n";
        }
        GDF << "edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN,color VARCHAR\n";

        #ifdef DEBUG
            if(DEBUG_LEVEL>=1) {
                printf("Marking edges to be colored as motifs...\n");
            }
        #endif
        TNGraph T(*(metaObj->G));
        if(Motifs!=NULL) {
            kSubgraphs->cursor = kSubgraphs->ini->next;
            while(kSubgraphs->cursor) {
                //std::cerr << "\n";
                if(kSubgraphs->cursor->label==Motifs->at((size_t) m)) {
                    for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
                        for(int j=0; j<kSubgraphs->cursor->vertices.size(); j++) {
                            //std::cerr << "···";
                            if((G.GetNI(kSubgraphs->cursor->vertices.at((size_t) i))).IsOutNId(kSubgraphs->cursor->vertices.at((size_t) j))) {
                                std::vector<int> tmp(*(metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) i))->second));
                                std::vector<int> aux(*(metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) j))->second));
                                aux.insert(aux.end(), tmp.begin(), tmp.end());
                                for(std::vector<int>::iterator it0 = aux.begin(); it0 != aux.end(); ++it0) {
                                    for(std::vector<int>::iterator it1 = aux.begin(); it1 != aux.end(); ++it1) {
                                        if(*it0!=*it1 && T.GetNI(*it0).IsOutNId(*it1)) {
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
            printf("\n\nHMM\n\n");
        }

        #ifdef DEBUG
            if(DEBUG_LEVEL>=1) {
                printf("Inserting remaining edges...\n");
            }
        #endif
        *//*for(TNGraph::TEdgeI EI = T.BegEI(); EI < T.EndEI(); EI++) {
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
        }*//*
        for(TNGraph::TEdgeI EI = T.BegEI(); EI < T.EndEI(); EI++)
            GDF << EI.GetSrcNId() << "," << EI.GetDstNId() << ",true,'0,0,0'\n";

        #ifdef DEBUG
            if(DEBUG_LEVEL>=1) {
                printf("Saving GDF...\n");
            }
        #endif
        GDF.close();
    }
    printf("metaObj->G nodes: %d\n", metaObj->G->GetNodes());*/
}

void GD::PrintMotifs(TNGraph &G, std::vector<std::multiset<unsigned long>> &Motifs, std::vector<unsigned long> &IDs,
                     std::unordered_map<std::multiset<unsigned long>, int> &kSubgraphs, std::string destination,
                     GD::MetaObject *metaObj)
{
    /*for(int m=0; m<Motifs.size(); m++) {
        PNGraph motif = TNGraph::New();
        kSubgraphs->cursor = kSubgraphs->ini->next;
        while(kSubgraphs->cursor->label!=Motifs[m])
            kSubgraphs->cursor = kSubgraphs->cursor->next;

        for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
            std::vector<int> *aux = metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) i))->second;
            for(std::vector<int>::iterator it = aux->begin(); it != aux->end(); ++it) {
                if(motif->IsNode(*it)) continue;
                motif->AddNode(*it);
            }
        }
        for(int i=0; i<kSubgraphs->cursor->vertices.size(); i++) {
            for(int j=0; j<kSubgraphs->cursor->vertices.size(); j++) {
                std::vector<int> aux(*(metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) i))->second));
                std::vector<int> tmp(*(metaObj->metaMap->find(kSubgraphs->cursor->vertices.at((size_t) j))->second));
                aux.insert(aux.end(), tmp.begin(), tmp.end());
                for(std::vector<int>::iterator it0 = aux.begin(); it0 != aux.end(); ++it0) {
                    for(std::vector<int>::iterator it1 = aux.begin(); it1 != aux.end(); ++it1) {
                        if(*it0!=*it1 && metaObj->G->GetNI(*it0).IsOutNId(*it1)) {
                            motif->AddEdge(*it0,*it1);
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
    }*/
}
