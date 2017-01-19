#ifndef GRAPHDATA_HPP_INCLUDED
#define GRAPHDATA_HPP_INCLUDED

#include <Snap.h>
//#define MAXN 100
#undef max
#undef min
#include <vector>
#include <map>
#include <set>
#include "gtools.h"
#include "nauty.h"
//#include "nausparse.h"
#include "nautinv.h"

#ifdef DEBUG
    extern int DEBUG_LEVEL;
#endif

namespace GD {
    struct Cell;
    struct Class;
    class GraphList {
    public:
        GraphList();
        ~GraphList();
        void AddGraph();
        Cell* ini;
        Cell* cursor;
    };

    class MetaObject {
    public:
        MetaObject(MetaObject&);
        MetaObject(PNGraph &);
        ~MetaObject();
        TNGraph *G;
        std::map<long unsigned int, std::vector<long unsigned int>*> *metaMap;
    };

    void PrintThat(GD::GraphList*);

    void Enumerate(TNGraph&, long unsigned int, GD::GraphList*);

    void Explore(TNGraph&, std::vector<std::vector<long unsigned int>>&, long unsigned int, std::vector<bool>&,
                    std::vector<long unsigned int>&, std::vector<long unsigned int>&, GD::GraphList*,
                 long unsigned int);

    void mapNeighbors(TNGraph&, std::vector<std::vector<long unsigned int>>&);

    void Validate(std::vector<std::vector<long unsigned int>>&, std::vector<bool>&, std::vector<long unsigned int>&,
                    long unsigned int, std::vector<long unsigned int>&);

    std::vector<long unsigned int> genComb(long unsigned int, std::vector<long unsigned int>, long unsigned int,
                    std::vector<long unsigned int>&, long unsigned int level = 0);

    void updateIndex(std::vector<long unsigned int>&, long unsigned int);

    void Classify(TNGraph&, GD::GraphList*, long unsigned int, optionblk&, int, set*);

    TNGraph* Randomize(TNGraph&);

    void GetFrequencies(GD::GraphList*, std::map<std::multiset<long unsigned int>, long unsigned int>&);

    void DiscoverMotifs(std::vector<std::map<std::multiset<long unsigned int>, long unsigned int>>&,
                        std::vector<std::multiset<long unsigned int>>&, std::vector<long unsigned int>&,
                        long unsigned int, std::string, TNGraph&, GD::GraphList*);

    TNGraph* ConcatMotifs(TNGraph&, std::vector<std::multiset<long unsigned int>>&, GD::GraphList*, GD::MetaObject*);

    void ExportGDF(TNGraph&, std::vector<std::multiset<long unsigned int>>*, std::vector<long unsigned int>*,
                    GD::GraphList*, std::string, GD::MetaObject*);

    void PrintMotifs(TNGraph&, std::vector<std::multiset<long unsigned int>>&, std::vector<long unsigned int>&,
                    GD::GraphList*, std::string, GD::MetaObject*);

    void SaveResults(TNGraph&, std::vector<std::vector<long unsigned int>>&, GD::GraphList*);

}

//#include "GraphData.tpp"

#endif // GRAPHDATA_HPP_INCLUDED
