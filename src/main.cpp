/*
Compile with "g++ -Ofast -DDEBUG -march=native -flto -fwhole-program -o main main.cpp GraphData.cpp GetTimeMs64.cpp ../Snap-2.4/snap-core/Snap.o -std=gnu++11 -I../include -I../Snap-2.4/snap-core -I../Snap-2.4/glib-core -I../nauty25r9 -L../nauty25r9/ -l nauty"
*/
#include <iostream>
#include <Snap.h>
#include <GraphData.hpp>
#include <GetTimeMs64.hpp>
#include <sstream>
#include <atomic>
#include <thread>

#ifdef DEBUG
    int DEBUG_LEVEL;
#endif

struct KavoshData;
void Kavosh(std::string, std::string, int, int, int, int t);
KavoshData Kavosh(PNGraph &, std::string, int, int, int, GD::MetaObject &, int t);
void proccessRandomNetwork(PNGraph G, int motif_size, std::shared_ptr<std::atomic_int> counter, int num_null_models,
                           std::shared_ptr<std::vector<std::map<std::multiset<
                                   long unsigned int>, long unsigned int>>> FVector, optionblk options, int m, set *dnwork);

int main(int argc, char* argv[]) {

#if !HAVE_TLS
    fprintf(stderr,">E This program needs to be linked with a version\n");
    fprintf(stderr,"  of nauty successfully configured with --enable-tls.\n");
    exit(1);
#endif

    if (argc == 2) {
        if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
            std::cerr   << "\nKavosh algorithm implemented by Rodrigo Martins de Oliveira" << std::endl
                        << "This software takes as input a directed network and identifies the motifs "
                        << "in it. The algorithm was extended to permit identifying metamotifs if "
                        << "desired." << std::endl
                        << "\n Usage: ./Kavosh -ARGUMENTS [-OPTIONS]" << std::endl
                        << "\nARGUMENTS:" << std::endl
                        << "-i\t Input network. Accepts relative path. Must be an edgelist." << std::endl
                        << "-o\t Output folder. Accepts relative path. Must to be created already." << std::endl
                        << "-s\t Motif size. The size of the motifs to look for. The higher, the slower." << std::endl
                        << "-r\t Number of random networks. The more random networks, the more "
                        << "precise is the motif identification." << std::endl
                        << "\nOPTIONS:" << std::endl
                        << "-t\t Number of threads. Use an appropriate number of threads for your machine." << std::endl
                        << "--metamotifs\t Perform metamotifs identification. For this "
                        << "option to work correctly a folder named 'metamotifs' must exist inside "
                        << "the output folder informed as argument.\n" << std::endl;
            return 1;
        }
    }
    if (argc < 8) {
        std::cerr << "\nUnsufficient arguments!" << std::endl;
        std::cerr << "Usage: " << argv[0] << " -i INPUT_GRAPH -o OUTPUT_FOLDER -s MOTIF_SIZE -r "
                  << "NUMBER_OF_RANDOM_GRAPHS [-t NUMBER_OF_THREADS] [--metamotifs NUMBER_OF_REPETITIONS]" << std::endl;
        std::cerr << "For more help run: " << argv[0] << " -h\n" << std::endl;
        return 1;
    }
    std::string source;
    std::string destination;
    int metamotifs = 0;
    int motif_size;
    int num_null_models;
    int t = 1;
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
        } else if(std::string(argv[i]) == "-t") {
            if (i + 1 < argc) {                                         // Make sure we aren't at the end of argv!
                std::stringstream value;
                value << argv[++i];
                value >> t;
            } else {                                                    // Unsufficient arguments
                std::cerr << "Option -t requires one argument!" << std::endl;
                return 1;
            }
        }  else if(std::string(argv[i]) == "--metamotifs") {
            if (i + 1 < argc) {                                         // Make sure we aren't at the end of argv!
                std::stringstream value;
                value << argv[++i];
                value >> metamotifs;
            } else {                                                    // Unsufficient arguments
                std::cerr << "Option --metamotifs requires one argument!" << std::endl;
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

    Kavosh(source, destination, metamotifs, motif_size, num_null_models, t);

    return 0;
}

struct KavoshData {
    PNGraph G;
    std::string destination;
    int metamotifs;
    int motif_size;
    int num_null_models;
    GD::MetaObject *metaObj;
    int t;

    KavoshData(PNGraph G, std::string destination, int metamotifs, int motif_size, int num_null_models,
               GD::MetaObject metaObj, int t) {
        this->G = G;
        this->destination = destination;
        this->metamotifs = metamotifs;
        this->motif_size = motif_size;
        this->num_null_models = num_null_models;
        this->metaObj = new GD::MetaObject(metaObj);
        this->t = t;
    }
};

void Kavosh(std::string source, std::string destination, int metamotifs, int motif_size, int num_null_models, int t) {
    char const* file_path = source.c_str();
    PNGraph G = TSnap::LoadEdgeList<PNGraph>(file_path, 0, 1);
    GD::MetaObject metaObj(G);
    KavoshData kd(G, destination, metamotifs, motif_size, num_null_models, metaObj, t);
    while (kd.metamotifs>=0) {
        kd = Kavosh(kd.G, kd.destination, kd.metamotifs, kd.motif_size, kd.num_null_models, *kd.metaObj, kd.t);
    }
}

KavoshData Kavosh(PNGraph &G, std::string destination, int metamotifs, int motif_size, int num_null_models,
            GD::MetaObject &metaObj, int t)
{

    // FIXME: Cross-platform
    // Make sure output directories exist
    /*std::string command = "mkdir " + destination;
    system(command.c_str());
    if (metamotifs > 0) {
        command =  "mkdir " + destination;
        for (int i=0; i<metamotifs; i++) {
            command += "/metamotifs";
            system(command.c_str());
        }
    }*/

    uint64 start_time = GetTimeMs64();

    uint64 enum_time = 0;
    uint64 class_time = 0;
    uint64 freq_time = 0;

    //char const* file_path = source.c_str();

    GD::GraphList* kSubgraphs = new GD::GraphList;
    std::vector<long unsigned int> IDs;
    std::vector<std::multiset<long unsigned int>> Motifs;
    //std::vector<std::map<std::vector<long unsigned int>, long unsigned int>> FVector(num_null_models+1);
    std::shared_ptr<std::vector<std::map<std::multiset<long unsigned int>, long unsigned int>>> FVectorPtr
            (new std::vector<std::map<std::multiset<long unsigned int>, long unsigned int>>(num_null_models+1));
    std::vector<std::thread> threads;
    std::shared_ptr<std::atomic_int> shared_counter (new std::atomic_int (0));
    int m = SETWORDSNEEDED(motif_size);
    nauty_check(WORDSIZE,m,motif_size,NAUTYVERSIONID);
    DYNALLSTAT(set,dnwork,dnwork_sz);
    DYNALLOC1(set,dnwork,dnwork_sz,2*60*m,"densenauty malloc");
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    //read input
    //PNGraph G = TSnap::LoadEdgeList<PNGraph>(file_path, 0, 1);

    #ifdef DEBUG
        printf("\n::::::::::ENUMERATION AND CLASSIFICATION ON ORIGNIAL GRAPH:::::::::::\n");
    #endif

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1)
            printf("Starting enumeration...\n");
    #endif
    uint64 t0 = GetTimeMs64();
    GD::Enumerate(*G, motif_size, kSubgraphs);
    #ifdef DEBUG
        if(DEBUG_LEVEL>=1)
            printf("Enumeration done!\nStarting classification...\n");
    #endif
    uint64 t1 = GetTimeMs64();
    GD::Classify(*G, kSubgraphs, motif_size, options, m, dnwork);
    #ifdef DEBUG
        if(DEBUG_LEVEL>=1)
            printf("Classification done!\nGetting frequencies...\n");
    #endif
    uint64 t2 = GetTimeMs64();
    GD::GetFrequencies(kSubgraphs, FVectorPtr->at(shared_counter->fetch_add(1)));
    uint64 t3 = GetTimeMs64();
    #ifdef DEBUG
        if(DEBUG_LEVEL>=1)
            printf("Frequencies obtained!\n");
    #endif

    #ifdef DEBUG
        if(DEBUG_LEVEL>=1) printf("\nEnumeration: %lu ms\nClassification: %lu ms\nGetting Frequencies: %lu ms\n", t1-t0, t2-t1, t3-t2);
        printf("\n::::::::::ENUMERATION AND CLASSIFICATION ON RANDOM GRAPHS::::::::::::\n");
    #endif
    enum_time += t1-t0;
    class_time += t2-t1;
    freq_time += t3-t2;

    for(int i=0; i<t; i++) threads.push_back(
                std::thread(proccessRandomNetwork, G, motif_size, shared_counter, num_null_models, FVectorPtr, options, m, dnwork));

    for(auto& th : threads) th.join();

    /*if(motif_size==3) {
        GD::PrintThat(kSubgraphs);
    }*/

    uint64 t10 = GetTimeMs64();
    GD::DiscoverMotifs(*FVectorPtr, Motifs, IDs, motif_size, destination, *G, kSubgraphs);
    uint64 t11 = GetTimeMs64();

    uint64 finish_time = GetTimeMs64();

    printf("\n%lu motifs found.\nTime elapsed: %u ms\n Enumeration took %u ms\n Classification took %u ms\n"
                   "Getting frequencies took %u ms\n Discovering motifs took %u ms\n", Motifs.size(),
           finish_time - start_time, enum_time, class_time, freq_time, t11-t10);

    GD::ExportGDF(*G, &Motifs, &IDs, kSubgraphs, destination, &metaObj);
    GD::PrintMotifs(*G, Motifs, IDs, kSubgraphs, destination, &metaObj);

    TNGraph *H = GD::ConcatMotifs(*G, Motifs, kSubgraphs, &metaObj);
    GD::ExportGDF(*G, NULL, NULL, NULL, destination, &metaObj);
    PNGraph H2 = H;
    TSnap::SaveEdgeList(H2, "concatenated_motifs.txt");
    //Kavosh(H2, destination + "metamotifs/", metamotifs - 1, motif_size, num_null_models, metaObj, 0);

    return KavoshData(H2, destination + "metamotifs/", metamotifs - 1, motif_size, num_null_models, metaObj, t);
}

void proccessRandomNetwork(PNGraph G, int motif_size, std::shared_ptr<std::atomic_int> counter, int num_null_models,
                           std::shared_ptr<std::vector<std::map<std::multiset
                                   <long unsigned int>, long unsigned int>>> FVector, optionblk options, int m, set *dnwork)
{
    uint64 t3 = GetTimeMs64();
    for(int pos = counter->fetch_add(1); pos <= num_null_models; pos = counter->fetch_add(1)) {
        if (pos%100==0)
            std::cerr << pos << std::endl;

        GD::GraphList *kRSubgraphs = new GD::GraphList;
        uint64 t4 = GetTimeMs64();
        #ifdef DEBUG
                if(DEBUG_LEVEL>=2)
                            printf("Starting randomization...\n");
        #endif
        TNGraph *R = GD::Randomize(*G);
        #ifdef DEBUG
                if(DEBUG_LEVEL>=2)
                            printf("Randomization done!\n");
        #endif
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
        GD::Classify(*R, kRSubgraphs, motif_size, options, m, dnwork);
        uint64 t8 = GetTimeMs64();
        GD::GetFrequencies(kRSubgraphs, FVector->at(pos));
        uint64 t9 = GetTimeMs64();
        #ifdef DEBUG
                if(DEBUG_LEVEL>=1 || DEBUG_LEVEL==-1)
                            printf("\nEnumeration: %lu ms\nClassification: %lu ms\nGetting Frequencies: %lu ms\n", t7-t6, t8-t7, t9-t8);
                        else if (DEBUG_LEVEL<=0)
                            printf("%d. null model processed - Time consumed: %lu ms\n", pos, t9-t3);
        #endif
        delete R;
        delete kRSubgraphs;
    }
}