#include "random.h"
using namespace std;


int main(int argc, char **argv) {
	if (argc <= 3) {
		cout << "error in number of arguments" << endl;
	}
	string exec_type = argv[1];
	int left = stoi(argv[2]);
	int right = stoi(argv[3]);
	double edges = stod(argv[4]);
	
	cout << "start generation."  << endl;
	// total_random_graph(left,right,edges,exec_type);
	generate_Barabasi_Albert_random_graph(left,right,edges,exec_type);
	// generate_Erdos_Renyi_random_graph(left,right,edges,exec_type);
	return 0;
}