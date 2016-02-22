#include <Snap.h>

int main() {
	// create a directed random graph on 100 nodes and 1k edges
	printf("Creating new random graph with 100 nodes and 1000 edges...\n");
	PNGraph Graph = TSnap::GenRndGnm<PNGraph>(100, 1000);

	printf("EndNI's ID: %d\n", Graph->EndNI().GetId());

	printf("Deleting node 0...\n");
	Graph->DelNode(0);
	Graph->DelNode(42);
	Graph->DelNode(43);
	Graph->DelNode(44);
	Graph->DelNode(45);
	Graph->DelNode(46);
	Graph->DelNode(47);
	Graph->DelNode(48);
	Graph->DelNode(52);
	Graph->DelNode(73);
	Graph->DelNode(74);
	Graph->DelNode(75);

	printf("Deleted! EndNI's ID: %d\n", Graph->EndNI().GetId());
	// traverse the nodes
	for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
	    printf("node id %d with out-degree %d and in-degree %d\n",
		NI.GetId(), NI.GetOutDeg(), NI.GetInDeg());
	}

	printf("\nFixing...\n");

	for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		printf("%d\n", NI.GetId());
	}

	int i=0;
	TNGraph::TNodeI aux;
	TNGraph::TNodeI NI = Graph->BegNI();
	do {
	//for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
	    if(NI.GetId()!=i) {
			int newNodeId = Graph->AddNode(i);
			printf("New node: %d, i: %d, NId to be deleted: %d, EndNI's ID: %d\n", newNodeId, i, NI.GetId(), Graph->EndNI().GetId());
			for(int j=0; j<NI.GetOutDeg(); j++)
				Graph->AddEdge(newNodeId, NI.GetOutNId(j));
			for(int j=0; j<NI.GetInDeg(); j++)
				Graph->AddEdge(NI.GetInNId(j), newNodeId);
			aux = NI;
			aux--;
			Graph->DelNode(NI.GetId());
			NI = aux++;
		} else {
			NI++;
		}
		i++;
	} while(NI < Graph->EndNI());

	// traverse the nodes
	for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
	    printf("node id %d with out-degree %d and in-degree %d\n",
		NI.GetId(), NI.GetOutDeg(), NI.GetInDeg());
	}

	printf("First node id is %d\n", Graph->BegNI().GetId());
	TNGraph::TNodeI zero = Graph->GetNI(0);
	printf("node id 0 with out-degree %d and in-degree %d\n", zero.GetOutDeg(), zero.GetInDeg());

}
