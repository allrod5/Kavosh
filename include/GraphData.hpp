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
#include <atomic>
#include <memory>
#include <functional>
#include "gtools.h"
#include "nauty.h"
//#include "nausparse.h"
#include "nautinv.h"

#ifdef DEBUG
    extern int DEBUG_LEVEL;
#endif

namespace std {
    template <>
    struct hash<std::multiset<unsigned long>>
    {
        std::size_t operator()(const std::multiset<unsigned long>& k) const
        {
            std::size_t r = 0;

            bool shift = false;
            for (auto&& it : k) {
                r = (r >> !shift) ^ (std::hash<unsigned long>()(it) << shift);
                shift = !shift;
            }

            return r;
        }
    };
}

namespace GD {

    typedef std::unordered_map<std::multiset<unsigned long>, int*> graphmap;

    class MetaObject {
    public:
        MetaObject(MetaObject&);
        MetaObject(PNGraph &);
        ~MetaObject();
        TNGraph *G;
        std::map<int, std::vector<int>*> *metaMap;
    };

    void Enumerate(TNGraph&, int, std::shared_ptr<graphmap>&, optionblk&, int, set*, bool, int);

    void Explore(TNGraph&, std::vector<std::vector<int>>&, int, std::vector<bool>&, std::vector<int>&,
                 std::vector<int>&, std::shared_ptr<graphmap>&, int, int, optionblk&, int, set*, bool, int);

    std::multiset<unsigned long>& Classify(TNGraph&, std::vector<int>&, int, optionblk&, int, set*);

    void mapNeighbors(TNGraph&, std::vector<std::vector<int>>&);

    void Validate(std::vector<std::vector<int>>&, std::vector<bool>&, std::vector<int>&,
                    int, std::vector<int>&);

    std::vector<int> genComb(int, std::vector<int>, std::vector<int> &, int level);

    void updateIndex(std::vector<int>&, int);

    TNGraph* Randomize(TNGraph&);

    void DiscoverMotifs(std::vector<std::multiset<unsigned long>>&, std::vector<unsigned long>&, int, std::string,
                        TNGraph&, std::shared_ptr<graphmap>&, int);

    TNGraph* ConcatMotifs(TNGraph&, std::vector<std::multiset<unsigned long>>&,
                          std::unordered_map<std::multiset<unsigned long>, int>&, GD::MetaObject*);

    void ExportGDF(TNGraph&, std::vector<std::multiset<unsigned long>>*, std::vector<unsigned long>*,
                   std::unordered_map<std::multiset<unsigned long>, int>*, std::string, GD::MetaObject*);

    void PrintMotifs(TNGraph&, std::vector<std::multiset<unsigned long>>&, std::vector<unsigned long>&,
                     std::unordered_map<std::multiset<unsigned long>, int>&, std::string, GD::MetaObject*);

}

//#include "GraphData.tpp"

#endif // GRAPHDATA_HPP_INCLUDED
