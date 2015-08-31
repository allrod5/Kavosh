/*
Compile with "g++ -Ofast -o main main.cpp GraphData.cpp GetTimeMs64.cpp ../Snap-2.4/snap-core/Snap.o -std=gnu++11 -I../include -I../Snap-2.4/snap-core -I../Snap-2.4/glib-core -I../nauty25r9 -L../nauty25r9/ -l nauty"
*/
#include <iostream>
#include <Snap.h>
#include <nauty.h>
#include <GraphData.hpp>
#include <GetTimeMs64.hpp>
#include <string>
#include <iostream>
#include <sstream>

#ifdef DEBUG
    int DEBUG_LEVEL;
#endif

int main(int argc, char* argv[]) {

    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << "-i INPUT_GRAPH -o OUTPUT_FOLDER -s MOTIF_SIZE -r NUMBER_OF_RANDOM_GRAPHS" << std::endl;
        return 1;
    }
    std::string source;
    std::string destination;
    int motif_size;
    int num_null_models;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-i") {
            if (i + 1 < argc) {                                         // Make sure we aren't at the end of argv!
                source = argv[++i];                                     // Increment 'i' so we don't get the argument as the next argv[i].
            } else {                                                    // Unsufficient arguments
                std::cerr << "Option -i requires one argument!" << std::endl;
                return 1;
            }
        } else if(std::string(argv[i]) == "-o") {
            if (i + 1 < argc) {                                         // Make sure we aren't at the end of argv!
                destination = argv[++i];
            } else {                                                    // Unsufficient arguments
                std::cerr << "Option -o requires one argument!" << std::endl;
                return 1;
            }
        } else if(std::string(argv[i]) == "-s") {
            if (i + 1 < argc) {                                         // Make sure we aren't at the end of argv!
                std::stringstream value;
                value << argv[++i];
                value >> motif_size;
            } else {                                                    // Unsufficient arguments
                std::cerr << "Option -s requires one argument!" << std::endl;
                return 1;
            }
        } else if(std::string(argv[i]) == "-r") {
            if (i + 1 < argc) {                                         // Make sure we aren't at the end of argv!
                std::stringstream value;
                value << argv[++i];
                value >> num_null_models;
            } else {                                                    // Unsufficient arguments
                std::cerr << "Option -r requires one argument!" << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Can't understand options!" << std::endl;
            return 1;
        }
    }

    #ifdef DEBUG
        printf("    !DEBUG MODE ON!\n  Choose debugging level:\n -1 - full benchmark\n  0 - minimalistic, only for progress observation\n  1 - superficial, for fast debug\n  2 - detailed, gives better overview\n  3 - in depth, as much details as possible\n>>> ");
        std::cin >> DEBUG_LEVEL;
        std::cout << "\n::::::::::::::::START::::::::::::::::\n";
    #endif

    uint64 start_time = GetTimeMs64();

    uint64 enum_time = 0;
    uint64 class_time = 0;
    uint64 freq_time = 0;

	char const* file_path = source.c_str();;

	GD::GraphList* kSubgraphs = new GD::GraphList;
    std::vector<long unsigned int> IDs;
    std::vector<std::vector<long unsigned int>> Motifs;
    std::vector<std::map<std::vector<long unsigned int>, long unsigned int>> FVector(num_null_models+1);
    //read input
    PNGraph G = TSnap::LoadEdgeList<PNGraph>(file_path, 0, 1);

    #ifdef DEBUG
        printf("\n::::::::::ENUMERATION AND CLASSIFICATION ON ORIGNIAL GRAPH:::::::::::\n");
    #endif

    uint64 t0 = GetTimeMs64();
    GD::Enumerate(*G, motif_size, kSubgraphs);
    uint64 t1 = GetTimeMs64();
    GD::Classify(*G, kSubgraphs, motif_size);
    uint64 t2 = GetTimeMs64();
    GD::GetFrequencies(kSubgraphs, FVector[0]);
    uint64 t3 = GetTimeMs64();

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1) printf("\nEnumeration: %lu ms\nClassification: %lu ms\nGetting Frequencies: %lu ms\n", t1-t0, t2-t1, t3-t2);
        printf("\n::::::::::ENUMERATION AND CLASSIFICATION ON RANDOM GRAPHS::::::::::::\n");
    #endif
    enum_time += t1-t0;
    class_time += t2-t1;
    freq_time += t3-t2;

    for(int i=0; i<num_null_models; i++) {
        GD::GraphList* kRSubgraphs = new GD::GraphList;
        uint64 t4 = GetTimeMs64();
        TNGraph *R = GD::Randomize(*G);
        uint64 t5 = GetTimeMs64();
        #ifdef DEBUG
            if(DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1) printf("\n------\nRandomization: %lu ms\n", t5-t4);
            if(DEBUG_LEVEL>=3) {
                printf("\nNew random graph edges: ");
                long unsigned int j=0;
                while(R->IsNode(j)) {
                    for(long unsigned int e=0; e<(R->GetNI(j)).GetOutDeg(); e++) {
                        printf("%lu-%d ", j, (R->GetNI(j)).GetOutNId(e));
                    }
                    j++;
                }
                printf("\n");
            }
        #endif
        uint64 t6 = GetTimeMs64();
        GD::Enumerate(*R, motif_size, kRSubgraphs);
        uint64 t7 = GetTimeMs64();
        GD::Classify(*R, kRSubgraphs, motif_size);
        uint64 t8 = GetTimeMs64();
        GD::GetFrequencies(kRSubgraphs, FVector[i+1]);
        uint64 t9 = GetTimeMs64();
        #ifdef DEBUG
            if(DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1)
                printf("\nEnumeration: %lu ms\nClassification: %lu ms\nGetting Frequencies: %lu ms\n", t7-t6, t8-t7, t9-t8);
            else if (DEBUG_LEVEL<=0) 
                printf("%d null models processed - Time elapsed: %lu ms\n", i+1, t9-t3);
        #endif
        enum_time += t7-t6;
        class_time += t8-t7;
        freq_time += t9-t8;
        delete R;
        delete kRSubgraphs;
    }
    uint64 t10 = GetTimeMs64();
    GD::DiscoverMotifs(FVector, Motifs, IDs, motif_size, destination, *G, kSubgraphs);
    uint64 t11 = GetTimeMs64();

    uint64 finish_time = GetTimeMs64();

    printf("\n%lu motifs found.\nTime elapsed: %lu ms\n Enumeration took %lu ms\n Classification took %lu ms\n Getting frequencies took %lu ms\n Discovering motifs took %lu ms\n", Motifs.size(), finish_time - start_time, enum_time, class_time, freq_time, t11-t10);

    GD::ExportGDF(*G, Motifs, IDs, kSubgraphs, destination);
    GD::PrintMotifs(*G, Motifs, IDs, kSubgraphs, destination);

    return 0;
}