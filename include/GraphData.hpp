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
        std::map<int, std::vector<int>*> *metaMap;
    };

    void Enumerate(TNGraph&, unsigned long, GD::GraphList*);

    void Explore(TNGraph&, std::vector<std::vector<int>>&, int, std::vector<bool>&,
                    std::vector<int>&, std::vector<int>&, GD::GraphList*,
                 unsigned long);

    void mapNeighbors(TNGraph&, std::vector<std::vector<int>>&);

    void Validate(std::vector<std::vector<int>>&, std::vector<bool>&, std::vector<int>&,
                    int, std::vector<int>&);

    std::vector<int> genComb(unsigned long, std::vector<int>, std::vector<int> &, unsigned long level);

    void updateIndex(std::vector<int>&, unsigned long);

    void Classify(TNGraph&, GD::GraphList*, int, optionblk&, int, set*);

    TNGraph* Randomize(TNGraph&);

    void GetFrequencies(GD::GraphList*, std::map<std::multiset<unsigned long>, int>&);

    void DiscoverMotifs(std::vector<std::map<std::multiset<unsigned long>, int>>&,
                        std::vector<std::multiset<unsigned long>>&, std::vector<unsigned long>&,
                        int, std::string, TNGraph&, GD::GraphList*);

    TNGraph* ConcatMotifs(TNGraph&, std::vector<std::multiset<unsigned long>>&, GD::GraphList*, GD::MetaObject*);

    void ExportGDF(TNGraph&, std::vector<std::multiset<unsigned long>>*, std::vector<unsigned long>*,
                    GD::GraphList*, std::string, GD::MetaObject*);

    void PrintMotifs(TNGraph&, std::vector<std::multiset<unsigned long>>&, std::vector<unsigned long>&,
                    GD::GraphList*, std::string, GD::MetaObject*);

}

//#include "GraphData.tpp"

#endif // GRAPHDATA_HPP_INCLUDED
