#include "random.h"
#include "static.h"
#include "dynamic.h"
#include "gephi.h"
#include "preprocessing.h"
#include "incremental_test.h"
#include <ctime>
#include "experiment.h"
#include "paper.h"
#include "mischievous.h"
#include "Multithread.h"
using namespace std;

vector<string> DATASET = { "wiki-en-cat_", "flickr-groupmemberships_",
"epinions-rating_", "edit-dewiki_", "reuters_", "gottron-trec_",
"delicious-ui_", "livejournal-groupmemberships_", "trackers_", "orkut-groupmemberships_"
};
vector<string> ALIAS = { "WC", "FG",
"EP", "DE", "RE", "TR",
"DUI", "LJ", "WT", "OG"
};



int main(int argc, char **argv) {
	if (argc == 1) {
		cout << "error in number of arguments" << endl;
	}
	string exec_type = argv[1];
	if (exec_type == "-BasicDecom") {
		cout << "start BasicDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		auto start = chrono::system_clock::now();
		coreIndexBasic(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "run time: " << time.count() << endl;
	}
	else if (exec_type == "-ComShrDecom") {
		cout << "start ComShrDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "run time: " << time.count() << endl;
	}
	else if (exec_type == "-ParallelDecom") {
		cout << "start ParallelDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		multithread_index_construction(g, stoi(argv[3]));
	}
	else if (exec_type == "-Query") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		vector<bool> left; vector<bool> right;
		// all the vertices in query result are set as true
		left.resize(g.num_v1, false); right.resize(g.num_v2, false);
		auto start = chrono::system_clock::now();
		retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left, right, stoi(argv[3]), stoi(argv[4]));
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "query time: " << time.count() << endl;
	}
	else if (exec_type == "-BiCore-Index-Ins") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), true);
		cout << "BiCore-Index-Ins running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Rem") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), false);
		cout << "BiCore-Index-Rem running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Ins*") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), true);
		cout << "BiCore-Index-Ins* running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Rem*") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), false);
		cout << "BiCore-Index-Rem* running time: " << time << endl;
	}
	else if (exec_type == "-ParallelIns") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), true, stoi(argv[5]));
		cout << "ParallelIns running time: " << time << endl;
	}
	else if (exec_type == "-ParallelRem") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), false, stoi(argv[5]));
		cout << "ParallelRem running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Batch") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		batch_set I, R;
		string RAND_MAX_V1_string = argv[5];
        string RAND_MAX_V2_string = argv[4];
        int RAND_MAX_V1 = stoi(RAND_MAX_V1_string);
        int RAND_MAX_V2 = stoi(RAND_MAX_V2_string);
	    int insert_count = 0;
	    vector<int> insert_first; vector<int> insert_second;
        string num_string = argv[3];
        int num = stoi(num_string);
        while (insert_count <= num){
            int t1 = rand()%RAND_MAX_V1; int t2 = rand()%RAND_MAX_V2;
            int i = 0;
            for (i; i < insert_first.size(); i++){
                if (insert_first[i] == t1 && insert_second[i] == t2) break;
            }
            if (i == insert_first.size()){
                insert_first.push_back(t1); insert_second.push_back(t2);
                insert_count += 1;
            }
        };
		bool flag = false;
        vector <vector<int>> insert_vector;
        insert_vector.push_back(insert_first);
		insert_vector.push_back(insert_second);
        for (int i=0; i<insert_count; i++){
            vector<vid_t>& tmp_neigh_ = g.neighbor_v1[insert_vector[0][i]];
            int ss = tmp_neigh_.size();
            for (int j = 0; j < ss; j++) {
                if (tmp_neigh_[j] == insert_vector[1][i]){
                    R.insert(make_pair(insert_vector[0][i], insert_vector[1][i]));
                    flag = true;
                    continue;
                }
            }
            if(flag){
                flag = false;
                continue;
            }
            I.insert(make_pair(insert_vector[0][i], insert_vector[1][i]));
        }
		auto time = Dyn_rebuild::update_bicore_index_batch(g, bicore_index_u, bicore_index_v, I, R);
		cout << "BiCore-Index-Batch running average time: " << time << endl;
	}
	else if (exec_type == "-Multi-BiCore-Index-Ins*") {
		vector<double> time_vector;
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        char* dir = argv[3];
        int u, v;
	    int r;
        string edgeFile = dir;
        bool flag = false;
        int i = 0;
        FILE * edgenum = fopen(edgeFile.c_str(), "r");
        while ((r = fscanf(edgenum, "%d %d", &u, &v)) != EOF)
        {
            i++;
            if (r != 2)
            {
                fprintf(stderr, "Bad edges format: u v incorrect\n");
                exit(1);
            }
            vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
            int ss = tmp_neigh_.size();
            for (int j = 0; j < ss; j++) {
                if (tmp_neigh_[j] == v){
                    flag = true;
                    continue;
                }
            }
            if(flag){
                flag = false;
                continue;
            }
			auto time = Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, true);
            time_vector.push_back(time);
        }
        fclose(edgenum);
        double sum = std::accumulate(std::begin(time_vector), std::end(time_vector), 0.0);
        double mean = sum/time_vector.size();
    
		cout << "Multi-BiCore-Index-Ins* running average time: " << mean << endl;
	}
	else if (exec_type == "-Multi-BiCore-Index-Rem*") {
		vector<double> time_vector;
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        char* dir = argv[3];
        int u, v;
	    int r;
        string edgeFile = dir;
        FILE * edgenum = fopen(edgeFile.c_str(), "r");
        while ((r = fscanf(edgenum, "%d %d", &u, &v)) != EOF)
        {
            if (r != 2)
            {
                fprintf(stderr, "Bad edges format: u v incorrect\n");
                exit(1);
            }
			auto time = Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, false);
            time_vector.push_back(time);
        }
        fclose(edgenum);

        double sum = std::accumulate(std::begin(time_vector), std::end(time_vector), 0.0);
        double mean = sum/time_vector.size();
		cout << "Multi-BiCore-Index-Rem* running average time: " << mean << endl;
	}
	else {
		cout << "illegal arguments" << endl;
	}
	return 0;
}