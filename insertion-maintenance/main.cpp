#include "random.h"
#include "gephi.h"
#include "preprocessing.h"
#include <ctime>
#include "paper.h"
#include <numeric>

using namespace std;



int main(int argc, char **argv) {
	if (argc == 1) {
		cout << "error in number of arguments" << endl;
	}
	string exec_type = argv[1];
	if (exec_type == "-ComShrDecom") {
		cout << "start ComShrDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "run time: " << time.count() << endl;
	}
    else if (exec_type == "-Build-BiCore"){
        char* a = argv[2];

        auto start = chrono::system_clock::now();
        BiGraph g(a);
        coreIndexKCore(g);
        vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        auto end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        
        cout << "BiCore-Core build: " << elapsed_seconds.count() << endl;
	}
    else if (exec_type == "-Build-BiCore-Number") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		vector<skyline_block*> skyline_index_u; vector<skyline_block*> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);

        //string txt_name = "path";
        //ofstream output_stream;  //(txt_name)
        //output_stream.open(txt_name,ofstream::app);

        //for (vid_t u = 0; u < g.num_v1; u++) {
        //    for(int i = 0; i < skyline_index_u[u]->vector_set.size()/2; i++){
        //        output_stream << "u ";
        //        output_stream << u;
         //       output_stream << " ";
         //       output_stream << skyline_index_u[u]->vector_set[i*2];
        //        output_stream << " ";
        //        output_stream << skyline_index_u[u]->vector_set[i*2+1];
        //        output_stream << endl;
         //   }
        //}
        //for (vid_t v = 0; v < g.num_v2; v++) {
        //    for(int i = 0; i < skyline_index_v[v]->vector_set.size()/2; i++){
        //        output_stream << "v ";
        //        output_stream << v;
        //        output_stream << " ";
        //        output_stream << skyline_index_v[v]->vector_set[i*2];
        //        output_stream << " ";
        //        output_stream << skyline_index_v[v]->vector_set[i*2+1];
        //        output_stream << endl;
        //    }
        //}
        //output_stream.close();
	}
	else if (exec_type == "-Edge-Insert") {
        BiGraph g(argv[2]);
        //BiGraph g_check = g;
        coreIndexKCore(g);
        vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block*> skyline_index_u; vector<skyline_block*> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
        auto time = Dyn_rebuild::update_skyline_index_swap_with_bfs(g, skyline_index_u, skyline_index_v, stoi(argv[3]), stoi(argv[4]));
		cout << "Skyline-Index-Ins running time: " << time << endl;

        //coreIndexKCore(g_check);
		//vector<vector<bicore_index_block*>> bicore_index_u_check; vector<vector<bicore_index_block*>> bicore_index_v_check;
		//build_bicore_index(g_check, bicore_index_u_check, bicore_index_v_check);
        //auto time_check = Dyn_rebuild::update_bicore_index_with_limit_swap(g_check, bicore_index_u_check, bicore_index_v_check, stoi(argv[3]), stoi(argv[4]), true);
		//cout << "BiCore-Index-Ins* running time: " << time_check << endl;
        //for (int i = 0; i < g_check.left_index.size(); i++) {
		//	for (int j = 1; j < g_check.left_index[i].size(); j++) {
		//		if(g_check.left_index[i][j] != g.left_index[i][j]) cout << "u " << i << " " << j << " " << g.left_index[i][j] << " " << g_check.left_index[i][j] << endl;
		//	}
		//}

		//for (int i = 0; i < g_check.right_index.size(); i++) {
		//	for (int j = 1; j < g_check.right_index[i].size(); j++) {
		//		if(g_check.right_index[i][j] != g.right_index[i][j]) cout << "v " << i << " " << j << " " << g.right_index[i][j] << " " << g_check.right_index[i][j] << endl;
		//	}
		//}
	}
    else if (exec_type == "-Batch-Edge-Insert") {
		vector<double> time_vector;
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block*> skyline_index_u; vector<skyline_block*> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
        char* dir = argv[3];
        int u, v;
	    int r;
        string edgeFile = dir;
        bool flag = false;
        FILE * edgenum = fopen(edgeFile.c_str(), "r");
        while ((r = fscanf(edgenum, "%d %d", &u, &v)) != EOF)
        {
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

            auto time = Dyn_rebuild::update_skyline_index_swap_with_bfs(g, skyline_index_u, skyline_index_v, u, v);
            time_vector.push_back(time);
        }
        fclose(edgenum);
        double sum = std::accumulate(std::begin(time_vector), std::end(time_vector), 0.0);
        double mean = sum/time_vector.size();
		cout << "Skyline-Index-Ins running average time: " << mean << endl;
	}
	else {
		cout << "illegal arguments" << endl;
	}
	return 0;
}
