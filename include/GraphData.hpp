#ifndef GRAPHDATA_HPP_INCLUDED
#define GRAPHDATA_HPP_INCLUDED

#include <Snap.h>
//#define MAXN 100
#undef max
#undef min
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include "gtools.h"
#include "nauty.h"
//#include "nausparse.h"
#include "nautinv.h"

#ifdef DEBUG
    extern int DEBUG_LEVEL;
#endif

namespace GD {

    class MetaObject {
    public:
        MetaObject(MetaObject&);
        MetaObject(PNGraph &);
        ~MetaObject();
        TNGraph *G;
        std::map<int, std::vector<int>*> *metaMap;
    };

    void Enumerate(TNGraph&, int, std::unordered_map<std::multiset<unsigned long>, int>&, optionblk&, int, set*);

    void Explore(TNGraph&, std::vector<std::vector<int>>&, int, std::vector<bool>&, std::vector<int>&,
                 std::vector<int>&, std::unordered_map<std::multiset<unsigned long>, int>&, int, int, optionblk&, int, set*);

    std::multiset<unsigned long>& Classify(TNGraph&, std::vector<int>&, int, optionblk&, int, set*);

    void mapNeighbors(TNGraph&, std::vector<std::vector<int>>&);

    void Validate(std::vector<std::vector<int>>&, std::vector<bool>&, std::vector<int>&,
                    int, std::vector<int>&);

    std::vector<int> genComb(unsigned long, std::vector<int>, std::vector<int> &, unsigned long level);

    void updateIndex(std::vector<int>&, unsigned long);

    TNGraph* Randomize(TNGraph&);

    void DiscoverMotifs(std::vector<std::map<std::multiset<unsigned long>, int>>&,
                        std::vector<std::multiset<unsigned long>>&, std::vector<unsigned long>&,
                        int, std::string, TNGraph&, std::unordered_map<std::multiset<unsigned long>, int>&);

    TNGraph* ConcatMotifs(TNGraph&, std::vector<std::multiset<unsigned long>>&,
                          std::unordered_map<std::multiset<unsigned long>, int>&, GD::MetaObject*);

    void ExportGDF(TNGraph&, std::vector<std::multiset<unsigned long>>*, std::vector<unsigned long>*,
                   std::unordered_map<std::multiset<unsigned long>, int>*, std::string, GD::MetaObject*);

    void PrintMotifs(TNGraph&, std::vector<std::multiset<unsigned long>>&, std::vector<unsigned long>&,
                     std::unordered_map<std::multiset<unsigned long>, int>&, std::string, GD::MetaObject*);

}

//#include "GraphData.tpp"

#endif // GRAPHDATA_HPP_INCLUDED
