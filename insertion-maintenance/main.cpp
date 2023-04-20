#include "random.h"
#include "gephi.h"
#include "preprocessing.h"
#include <ctime>
#include "paper.h"
#include <cmath>
#include <numeric>

using namespace std;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "error in number of arguments" << endl;
    }
    string exec_type = argv[1];
    if (exec_type == "-Query-Index-BiCore")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<bool> left;
        vector<bool> right;
        // all the vertices in query result are set as true
        left.resize(g.num_v1, false);
        right.resize(g.num_v2, false);
        auto start = chrono::system_clock::now();
        retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left, right, stoi(argv[3]), stoi(argv[4]));
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "bicore query time: " << time.count() << endl;

        // int i_right = 0; int i_left = 0;
        // int j_right = 0; int j_left = 0;
        // for (auto node: right) {
        //     if (node) i_right++;
        //     else j_right++;
        // }
        // for (auto node: left) {
        //     if (node) i_left++;
        //     else j_left++;
        // }
        // cout << i_left << " " << i_right << " " << j_left << " " << j_right << endl;
    }
    else if (exec_type == "-Query-BiCore")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block *> skyline_index_u;
        vector<skyline_block *> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
        vector<bool> left;
        vector<bool> right;
        // all the vertices in query result are set as true
        left.resize(g.num_v1, false);
        right.resize(g.num_v2, false);
        auto start = chrono::system_clock::now();
        retrieve_bicore_via_bicore_number(g, skyline_index_u, skyline_index_v, left, right, stoi(argv[3]), stoi(argv[4]));
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "bicore query time: " << time.count() << endl;

        // int i_right = 0; int i_left = 0;
        // int j_right = 0; int j_left = 0;
        // for (auto node: right) {
        //     if (node) i_right++;
        //     else j_right++;
        // }
        // for (auto node: left) {
        //     if (node) i_left++;
        //     else j_left++;
        // }
        // cout << i_left << " " << i_right << " " << j_left << " " << j_right << endl;
    }
    else if (exec_type == "-Query-Index-node")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        batch_set pair;
        auto start = chrono::system_clock::now();
        retrieve_node_via_bicore_index(g, bicore_index_u, bicore_index_v, pair, stoi(argv[3]), stoi(argv[4]));
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "Offset query time: " << time.count() << endl;
        // cout << "pair size: " << pair.size() << endl;
    }
    else if (exec_type == "-Query-node")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block *> skyline_index_u;
        vector<skyline_block *> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
        batch_set pair;
        auto start = chrono::system_clock::now();
        retrieve_node_via_index_skyline(g, skyline_index_u, skyline_index_v, pair, stoi(argv[3]), stoi(argv[4]));
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "Offset query time: " << time.count() << endl;
        // cout << "pair size: " << pair.size() << endl;
    }
    else if (exec_type == "-Query-Index-community")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        // node
        int node = -1;
        vector<bool> left_node;
        vector<bool> right_node;
        // all the vertices in query result are set as true
        left_node.resize(g.num_v1, false);
        right_node.resize(g.num_v2, false);
        retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left_node, right_node, stoi(argv[3]), stoi(argv[4]));
        for (int i = 0; i < left_node.size(); i++)
        {
            if (left_node[i])
            {
                node = i;
                break;
            }
        }
        if (node == -1)
        {
            cout << "No node" << endl;
        }
        else
        {
            vector<bool> left;
            vector<bool> right;
            // all the vertices in query result are set as true
            auto start = chrono::system_clock::now();
            set<int> left_set, right_set;
            vector<int> left_visit;
            vector<int> right_visit;
            left.resize(g.num_v1, false);
            right.resize(g.num_v2, false);
            left_visit.resize(g.num_v1);
            right_visit.resize(g.num_v2);
            fill_n(left_visit.begin(), left_visit.size(), 0);
            fill_n(right_visit.begin(), right_visit.size(), 0);
            retrieve_hierarchy_via_bicore_index(g, bicore_index_u, bicore_index_v, left, right, left_set, right_set, left_visit, right_visit, stoi(argv[5]), stoi(argv[6]), stoi(argv[3]), stoi(argv[4]));
            auto end = chrono::system_clock::now();
            chrono::duration<double> time = end - start;
            cout << "Community query time: " << time.count() << endl;
            // cout << left_set.size() << " " << right_set.size() << endl;
        }
    }
    else if (exec_type == "-Query-community")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block *> skyline_index_u;
        vector<skyline_block *> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);

        int node = -1;
        vector<bool> left_node;
        vector<bool> right_node;
        // all the vertices in query result are set as true
        left_node.resize(g.num_v1, false);
        right_node.resize(g.num_v2, false);
        retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left_node, right_node, stoi(argv[3]), stoi(argv[4]));
        for (int i = 0; i < left_node.size(); i++)
        {
            if (left_node[i])
            {
                node = i;
                break;
            }
        }
        if (node == -1)
        {
            cout << "No node" << endl;
        }
        else
        {
            vector<int> left;
            vector<int> right;
            // all the vertices in query result are set as true
            left.resize(g.num_v1, 0);
            right.resize(g.num_v2, 0);
            auto start = chrono::system_clock::now();
            queue<int> node_queue;
            queue<bool> is_left_queue;
            // node_queue.push(node); is_left_queue.push(true);
            node_queue.push(stoi(argv[5]));
            is_left_queue.push(stoi(argv[6]));
            retrieve_hierarchy_via_index_skyline(g, skyline_index_u, skyline_index_v, left, right, node_queue, is_left_queue, stoi(argv[3]), stoi(argv[4]));
            set<int> left_set, right_set;
            for (int i = 0; i < right.size(); i++)
            {
                if (right[i] == 1)
                    right_set.insert(i);
            }
            for (int i = 0; i < left.size(); i++)
            {
                if (left[i] == 1)
                    left_set.insert(i);
            }

            auto end = chrono::system_clock::now();
            chrono::duration<double> time = end - start;
            cout << "Community query time: " << time.count() << endl;
            // cout << left_set.size() << " " << right_set.size() << endl;
        }
    }
    else if (exec_type == "-ComShrDecom")
    {
        cout << "start ComShrDecom for " << argv[2] << endl;
        BiGraph g(argv[2]);
        auto start = chrono::system_clock::now();
        coreIndexKCore(g);
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "run time: " << time.count() << endl;
    }
    else if (exec_type == "-Build-BiCore")
    {
        char *a = argv[2];

        BiGraph g(a);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);

        cout << "BiCore-index build finished!" << endl;
    }
    else if (exec_type == "-Build-BiCore-Number")
    {
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block *> skyline_index_u;
        vector<skyline_block *> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);

        cout << "BiCore Number build finished!" << endl;

        // string txt_name = "path";
        // ofstream output_stream;  //(txt_name)
        // output_stream.open(txt_name,ofstream::app);

        // for (vid_t u = 0; u < g.num_v1; u++) {
        //     for(int i = 0; i < skyline_index_u[u]->vector_set.size()/2; i++){
        //         output_stream << "u ";
        //         output_stream << u;
        //        output_stream << " ";
        //        output_stream << skyline_index_u[u]->vector_set[i*2];
        //        output_stream << " ";
        //        output_stream << skyline_index_u[u]->vector_set[i*2+1];
        //        output_stream << endl;
        //   }
        //}
        // for (vid_t v = 0; v < g.num_v2; v++) {
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
        // output_stream.close();
    }
    else if (exec_type == "-Edge-Insert")
    {
        BiGraph g(argv[2]);
        // BiGraph g_check = g;
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block *> skyline_index_u;
        vector<skyline_block *> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
        auto time = Dyn_rebuild::update_skyline_index_swap_with_bfs(g, skyline_index_u, skyline_index_v, stoi(argv[3]), stoi(argv[4]));
        cout << "Edge-Insert running time: " << time << endl;

        // coreIndexKCore(g_check);
        // vector<vector<bicore_index_block*>> bicore_index_u_check; vector<vector<bicore_index_block*>> bicore_index_v_check;
        // build_bicore_index(g_check, bicore_index_u_check, bicore_index_v_check);
        // auto time_check = Dyn_rebuild::update_bicore_index_with_limit_swap(g_check, bicore_index_u_check, bicore_index_v_check, stoi(argv[3]), stoi(argv[4]), true);
        // cout << "BiCore-Index-Ins* running time: " << time_check << endl;
        // for (int i = 0; i < g_check.left_index.size(); i++) {
        //	for (int j = 1; j < g_check.left_index[i].size(); j++) {
        //		if(g_check.left_index[i][j] != g.left_index[i][j]) cout << "u " << i << " " << j << " " << g.left_index[i][j] << " " << g_check.left_index[i][j] << endl;
        //	}
        // }

        // for (int i = 0; i < g_check.right_index.size(); i++) {
        //	for (int j = 1; j < g_check.right_index[i].size(); j++) {
        //		if(g_check.right_index[i][j] != g.right_index[i][j]) cout << "v " << i << " " << j << " " << g.right_index[i][j] << " " << g_check.right_index[i][j] << endl;
        //	}
        // }
    }
    else if (exec_type == "-Multi-Edges-Insert")
    {
        vector<double> time_vector;
        BiGraph g(argv[2]);
        coreIndexKCore(g);
        vector<vector<bicore_index_block *>> bicore_index_u;
        vector<vector<bicore_index_block *>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        vector<skyline_block *> skyline_index_u;
        vector<skyline_block *> skyline_index_v;
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
        char *dir = argv[3];
        int u, v;
        int r;
        string edgeFile = dir;
        bool flag = false;
        int i = 0;
        FILE *edgenum = fopen(edgeFile.c_str(), "r");
        while ((r = fscanf(edgenum, "%d %d", &u, &v)) != EOF)
        {
            i++;
            if (r != 2)
            {
                fprintf(stderr, "Bad edges format: u v incorrect\n");
                exit(1);
            }
            vector<vid_t> &tmp_neigh_ = g.neighbor_v1[u];
            int ss = tmp_neigh_.size();
            for (int j = 0; j < ss; j++)
            {
                if (tmp_neigh_[j] == v)
                {
                    flag = true;
                    continue;
                }
            }
            if (flag)
            {
                flag = false;
                continue;
            }
            auto time = Dyn_rebuild::update_skyline_index_swap_with_bfs(g, skyline_index_u, skyline_index_v, u, v);
            time_vector.push_back(time);
        }
        fclose(edgenum);
        double sum = std::accumulate(std::begin(time_vector), std::end(time_vector), 0.0);
        double mean = sum / time_vector.size();

        cout << "Multi-Edges-Insert running average time: " << mean << endl;
    }
    else if (exec_type == "-Recompute")
    {
        BiGraph g(argv[2]);
        string u_ = argv[3];
        int u = stoi(u_);
        string v_ = argv[4];
        int v = stoi(v_);
        g.addEdge(u, v);

        cout << "start ComShrDecom for " << argv[2] << endl;
        auto start = chrono::system_clock::now();
        coreIndexKCore(g);
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "run time: " << time.count() << endl;
    }
    else if (exec_type == "-Multi-Recompute")
    {
        int count = 0;
        BiGraph g(argv[2]);
        char *dir = argv[3];
        int u, v;
        int r;
        string edgeFile = dir;
        FILE *edgenum = fopen(edgeFile.c_str(), "r");
        while ((r = fscanf(edgenum, "%d %d", &u, &v)) != EOF)
        {
            if (r != 2)
            {
                fprintf(stderr, "Bad edges format: u v incorrect\n");
                exit(1);
            }
            g.addEdge(u, v);
            count++;
        }
        fclose(edgenum);

        cout << "start ComShrDecom for " << argv[2] << endl;
        auto start = chrono::system_clock::now();
        coreIndexKCore(g);
        auto end = chrono::system_clock::now();
        chrono::duration<double> time = end - start;
        cout << "run time: " << time.count() / count << endl;
    }
    else
    {
        cout << "illegal arguments" << endl;
    }
    return 0;
}