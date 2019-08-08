//  Heaader files from Lemon
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/network_simplex.h>
#include <lemon/capacity_scaling.h>
#include <lemon/cost_scaling.h>
#include <lemon/cycle_canceling.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/heap.h>
#include <lemon/concept_check.h>
// other Headers
#include <vector>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace lemon;

typedef ListDigraph::Node Node;
typedef ListDigraph::Arc Arc;
typedef CapacityScaling<ListDigraph, int, double> CSAlg;


////////////////////////    Utility functions
    //  nchoosk
double nchoosek(int x, int y, int ncount = 0, int result = 1)
{
    if (y > x)        return -1;
    if (x == y)       return 1;
    if (y == 0)       return 1;
    if (y == 1)       return x;
    return x * nchoosek(x-1, y-1) / y;
}
    //  generating arrival times according to Poisson process and deadlines according to uniform dist between min and max deadline
void gen_arrival_deadline(const int K, const double lambda,
                          int MinDeadline, int MaxDeadline,
                          vector<int >& arrival, vector<int>& deadline){
    
    default_random_engine generator_exp, generator_uni;
    exponential_distribution<double> distribution_exp(lambda);
    uniform_int_distribution<int> distribution_uni(MinDeadline, MaxDeadline);
    
    int Ti = 0;
    for (int ni = 0; ni < K;  ni++) {
        double dT = distribution_exp(generator_exp);
        int dD = distribution_uni(generator_uni);
        arrival.push_back(Ti);
        deadline.push_back(Ti + dD);
        Ti = floor(Ti + dT);
    }
    
}
    //  converting a vector of int to string to generate a key
string vec_int_2_str(vector<int>& vec_int, string prefix){
    string to_str = prefix;
    for (int ni = 0; ni < vec_int.size(); ni++) {
        to_str += to_string(vec_int[ni]);
        if(ni < vec_int.size()-1) to_str += ",";
    }
    return to_str;
}
//  between all combination of length t subset of {0,1,...,K-1} generates next subset of cur_set;
void nextCombin(vector<int>& cur_set, int K, int t){
    if(t == 0) return;
    int count = t - 1;
    int Val = K-1;
    while ((count >= 0) && (cur_set[count] == Val)){
        count--;
        Val--;
    }
    
    if (!(count < 0)){
        cur_set[count]++;
        for (int ni = count + 1; ni < cur_set.size(); ni++)
            cur_set[ni] = cur_set[count] + ni - count;
    }
}

template <typename T>
T min(T x, T y) {return (x < y) ? x : y;}
//  Ascending sorting with keep tracking of indeces
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
    
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    
    return idx;
}
/////////////////////////////////////////////
//      Classes
////////////////////////////////////////////
//  A Class to construct graph and add node or arc
class DGraph{
    friend class MCF_Solver;
    friend class ProblemInstance;
private:
    ListDigraph G;
    ListDigraph::NodeMap<int> supply;   // suply of each node
    ListDigraph::NodeMap<string> node_label;   // label of nodes
    ListDigraph::ArcMap<int> up_bound;  //  upper bound of each arc
    ListDigraph::ArcMap<double> costs;  //  cost of each arc
    ListDigraph::ArcMap<string> arc_label;  //  label of arcs
    ListDigraph::ArcMap<double> gamma;      //  dual variable
    ListDigraph::ArcMap<double> flow;       //  primal solution x^(i)_u(l)
    
public:
    DGraph() : supply(G), up_bound(G), costs(G), node_label(G), arc_label(G), gamma(G), flow(G) {}
    ~DGraph() {}
    Node AddNode(int sup) {
        Node v = G.addNode();
        supply[v] = sup;
        return v;
    }
    Arc AddArc(Node v_1, Node v_2, int ub_bnd, double cst){
        Arc e = G.addArc(v_1, v_2);
        up_bound[e] = ub_bnd;
        costs[e] = cst;
        return e;
    }
    void update_cost(Arc e, double new_cst){
        costs[e] = new_cst;
    }
};
////////////////////////////////////////////
// A class for solving MCF
class MCF_Solver{
private:
    CSAlg CapScaling;
    CSAlg::ProblemType IsFeasible;
    
public:
    MCF_Solver(DGraph& dg) : CapScaling(dg.G) {
        CapScaling.upperMap(dg.up_bound);
        CapScaling.costMap(dg.costs);
        CapScaling.supplyMap(dg.supply);
    }
    void run_cost(ListDigraph::ArcMap<double>& costMap){
        CapScaling.costMap(costMap);
    }
    bool run(){
        IsFeasible = CapScaling.run();
        if (IsFeasible) return true;
        return false;
    }
    double get_flow(Arc e){
        return CapScaling.flow(e);
    }
    
    
};

////////////////////////////////////////////
/*      Creating Problem Instance and returning a graph */
class ProblemInstance{
private:
    
    vector<int> arrival_times;
    vector<int> deadlines;
    vector<vector<bool> > usr_time;
    vector<int> len_time;
    unordered_map<string, vector<int> > subfile_id;
    vector<double> zeta;
    unordered_set<string> all_usr_grp;
    
    

    void gen_usr_vs_time(vector<int>& arrival, vector<int>& deadline, vector<vector<bool> >& usr_time, vector<int>& len_time);
    vector<vector<int> > get_possible_ug(int l /* l: time interval*/);  // finding all possible user groups at time interval l
    void gen_subfile_id();          //  generates subfile id where file id of user k is
    void gen_source_sink_nodes();   //  generates source and sink nodes
    void gen_subfile_nodes();       //  generates nodes corresponding to missing subfiles and edges between source node and missing subfiles
    void gen_time_interval_nodes(); //  generates nodes corresponding to time intervals
    void gen_usr_grp_nodes();       //  generates nodes corresponding to user groups and edges between user group to time intervals, subfile to user groups
    vector<string> subfile_carried_by_usr_grp(vector<int>& usr_grp, int nk);    //  missing subfiles of user "nk" that can be carried by user group "usr_grp"
//    void update_cost_gamma();     //  update cost and gamma of edges
//    void update_zeta();     //  update zeta value, dual cofficient
    
    
    
public:
    DGraph G;
    const int K;
    const int t;
    const int r;
    int beta;
    double primal_solution;
    double dual_solution;
    vector<Node> source_nodes;
    vector<Node> sink_nodes;
    vector<vector<Node> > subfile_nodes;
    vector<vector<Node> > usr_grp_nodes;
    vector<vector<Node> > time_interval_nodes;
    vector<vector<Arc> > usr_time_arcs;                     //  which edges are connected to the same user groups at a given time interval
    unordered_map<string, vector<Node> > usr_group_nodes;//  key: user group 
    vector<int> time_interval_usr_grp;                      //  time interval of the user groups in "usr_time_arcs"
    vector<pair<int, int> > begining_end_time_interval;     //  begining and end time interval of each user
//    vector<double> primal_values;                           //  primal variables, i.e. x_u(l)'s
    
    ProblemInstance(vector<int>& arrival, vector<int>& deadline, const int num_users, const int cach_size, const int time_per_subfile) : arrival_times(arrival), deadlines(deadline), K(num_users), t(cach_size), r(time_per_subfile) {    }
    void ConstructProblem();
    void PrintGraph();           //  prints all nodes and edges of the constructed graph with their labels
    void SolveProblem(int Tmaxitr, double alpha, int offItrs);   //  solve the problem with Tmaxitr iteration and step size n^{-alpha}
    void project_gamma(vector<double>& tmp_gamma);
    void PrintNodes();
    void PrintArcs();
    void PrintOutput();
    
};


void ProblemInstance::gen_usr_vs_time(vector<int>& arrival, vector<int>& deadline,
                     vector<vector<bool> >& usr_time, vector<int>& len_time){
    vector<int> all_times = arrival;
    all_times.insert(all_times.end(), deadline.begin(), deadline.end());
    sort(all_times.begin(), all_times.end());
    vector<int>::iterator it = unique(all_times.begin(), all_times.end());
    all_times.resize(distance(all_times.begin(), it));
    for (int ni =0 ; ni < all_times.size()-1; ni++) {
        len_time.push_back(all_times[ni+1] - all_times[ni]);
//        cout << len_time.back() << endl;
    }
    for(int nl = 0; nl < all_times.size()-1; nl++){
        vector<bool> tmp;
        usr_time.push_back(tmp);
        int l_time = all_times[nl];
        int r_time = all_times[nl+1];
        for (int nk = 0; nk < arrival.size(); nk++) {
            begining_end_time_interval.push_back({0,0});
            if(arrival[nk] <= l_time && deadline[nk] >= r_time) usr_time[nl].push_back(true);
            else usr_time[nl].push_back(false);
            if(arrival[nk] == l_time) begining_end_time_interval[nk].first = nl;
            if(deadline[nk] == r_time) begining_end_time_interval[nk].second = nl;
            //            cout << usr_time[nl][nk] ;
        }
        //        cout << "" << endl;
    }
    
}

//      generating possible user groups at time interval "l"
vector<vector<int> > ProblemInstance::get_possible_ug(int l){
    vector<vector<int> > possible_ug;
    vector<int> active_usr;
    for (int nk = 0; nk < K; nk++) {
        if(usr_time[l][nk]) {active_usr.push_back(nk);}
    }
    if(active_usr.empty()) return possible_ug;
    // if(active_usr.size() == t+1) {possible_ug.push_back(active_usr); return possible_ug;}
    int t_max = (t+1 < active_usr.size()) ? t+1 : active_usr.size();
    for (int t0 = 1; t0 <= t_max; t0++) {
        vector<int> idx;
        for(int ni = 0; ni < t0; ni++) idx.push_back(ni);
        int num_ug_size_t0 = nchoosek(active_usr.size(), t0);
        for(int ni_ug = 0; ni_ug < num_ug_size_t0; ni_ug++){
            vector<int> usr_grp;
            for(int nj = 0; nj < t0; nj++) usr_grp.push_back(active_usr[idx[nj]]);
            nextCombin(idx, active_usr.size(), t0);
            possible_ug.push_back(usr_grp);
        }
    }
    
    return possible_ug;
}
//  generating an id for each missing subfile so that subfile with id "id" in user "k" is subfile_id[id][nk]-th missing subfile of user "k"
void ProblemInstance::gen_subfile_id(){
    int num_sub_files = nchoosek(K, t);
    vector<int> idx;
    vector<int> count(K, 0);
    for(int ni = 0; ni < t; ni++) idx.push_back(ni);
    for(int nf = 0; nf < num_sub_files; nf++){
        vector<int> count_nf(K);
        string label = vec_int_2_str(idx, "");
        unordered_set<int> usr_has_subfile(idx.begin(), idx.end());
        for (int nk = 0; nk < K; nk++) {
            if(usr_has_subfile.count(nk) != 0) count_nf[nk] = -1;
            else count_nf[nk] = count[nk]++;
        }
        subfile_id[label] = count_nf;
        nextCombin(idx, K, t);
    }
}

void ProblemInstance::gen_source_sink_nodes(){
    int num_mis_subfile = nchoosek(K-1, t);
    for (int nk = 0; nk < K; nk++) {
        Node s = G.G.addNode();
        Node t = G.G.addNode();
        G.supply[s] = r * num_mis_subfile;
        G.supply[t] = -r * num_mis_subfile;
        source_nodes.push_back(s);
        sink_nodes.push_back(t);
        G.node_label[s] = "source " + to_string(nk);
        G.node_label[t] = "sink " + to_string(nk);
    }
}
//  generates nodes corresponding to missing subfiles and edges between source and the missing subiles
void ProblemInstance::gen_subfile_nodes(){
    for(int nk = 0; nk < K; nk++) {vector<Node> tmp; subfile_nodes.push_back(tmp);}
    int num_subfiles = nchoosek(K,t);
    for(auto it = subfile_id.begin(); it != subfile_id.end(); ++it){
        vector<int> count_nf(it->second.begin(), it->second.end());
        string key = it->first;
        for(int nk = 0; nk < K; nk++){
            if(count_nf[nk] >= 0){
                Node v = G.G.addNode(); // node corresponding to "count_nf[nk]" missing subfile
                G.node_label[v] = "subfile {" + key + "} of user " + to_string(nk);
                G.supply[v] = 0;
                subfile_nodes[nk].push_back(v);
                Arc e = G.G.addArc(source_nodes[nk], v);    // arc between source node and node "v"
                G.arc_label[e] = "user " + to_string(nk) + " source to subfile {" + key + "}";
                G.costs[e] = 0;
                G.up_bound[e] = r;
                G.gamma[e] = 0;
                G.flow[e] = 0;
            }
        }
    }
}
// generates nodes corresponding to time intervals and edges from time intervals to sink
void ProblemInstance::gen_time_interval_nodes(){
    for(int nk = 0; nk < K; nk++){
        vector<Node> tmp_nk;
        for(int nl = begining_end_time_interval[nk].first; nl <= begining_end_time_interval[nk].second; nl++){
            Node u = G.G.addNode();
            G.node_label[u] = "time interval " + to_string(nl) + " user " + to_string(nk);
            G.supply[u] = 0;
            tmp_nk.push_back(u);
            Arc e = G.G.addArc(u, sink_nodes[nk]);
            G.arc_label[e] = "user " + to_string(nk) + " time " + to_string(nl) + " to sink";
            G.up_bound[e] = len_time[nl];
            G.costs[e] = 0;
        }
        time_interval_nodes.push_back(tmp_nk);
    }
}
vector<string> ProblemInstance::subfile_carried_by_usr_grp(vector<int>& usr_grp, int nk){
    vector<int> fix_idx;
    for(int ni = 0; ni < usr_grp.size(); ni++) if(usr_grp[ni] != nk) fix_idx.push_back(usr_grp[ni]);
    vector<int> var_idx;
    unordered_set<int> fix_set(usr_grp.begin(), usr_grp.end());
    for (int ni = 0; ni < K; ni++)
        if(fix_set.count(ni) == 0) var_idx.push_back(ni);
    vector<string> miss_subfiles;
    int num_fix = fix_idx.size();
    int num_var = K - 1 - num_fix;
    int num_subfiles = nchoosek(num_var, t-num_fix);
    vector<int> idx;
    for(int ni = 0; ni < t-num_fix; ni++) idx.push_back(ni);
    for(int nf = 0; nf < num_subfiles; nf++){
        vector<int> cur_subfile(fix_idx.begin(), fix_idx.end());
        for(int ni = 0; ni < t-num_fix; ni++) cur_subfile.push_back(var_idx[idx[ni]]);
        sort(cur_subfile.begin(), cur_subfile.end());
        miss_subfiles.push_back(vec_int_2_str(cur_subfile, ""));
        if(!idx.empty()) nextCombin(idx, num_var, t-num_fix);
    }
    return miss_subfiles;
}

// generating nodes corresponding to user groups and edges coming and leaving these nodes
void ProblemInstance::gen_usr_grp_nodes(){
    for (int nl = 0; nl < beta; nl++) {
        zeta.push_back(0);
        vector<vector<int> > usr_grp_cur_time = get_possible_ug(nl);
        int num_cur_ug = usr_grp_cur_time.size();
        // cout << "+++++++++++++++++++++++" << num_cur_ug << "++++++++++++++++++" << endl;
        for(int nu = 0; nu < num_cur_ug; nu++){
            string key_ug = vec_int_2_str(usr_grp_cur_time[nu], "usr_grp {") + "}";
            // cout << key_ug << endl;
            vector<Arc> arcs_cur_usr_grp;   // edges going from current user group to time intervals
            time_interval_usr_grp.push_back(nl);
            if(usr_group_nodes.count(key_ug) > 0){  // the user group already exsited
                for(int nk = 0; nk < usr_grp_cur_time[nu].size(); nk++){
                    int cur_usr = usr_grp_cur_time[nu][nk];
                    Node u = usr_group_nodes[key_ug][nk];
                    Arc e = G.G.addArc(u, time_interval_nodes[cur_usr][nl-begining_end_time_interval[cur_usr].first]);   // edge between user group and time interval
                    G.up_bound[e] = len_time[nl];
                    G.costs[e] = 1 / (double) usr_grp_cur_time[nu].size();
                    G.arc_label[e] = "user " + to_string(cur_usr) + " user group " + key_ug + " to time " + to_string(nl);
                    G.gamma[e] = G.costs[e];
                    G.flow[e] = 0;
                    arcs_cur_usr_grp.push_back(e);    
                }
            }else{
                vector<Node> tmp_node_usr_grp;
                for(int nk = 0; nk < usr_grp_cur_time[nu].size(); nk++){
                    int cur_usr = usr_grp_cur_time[nu][nk];
                    Node u = G.G.addNode();     //  node corresponding to this user group and user "cur_usr"
                    tmp_node_usr_grp.push_back(u);
                    G.node_label[u] =  key_ug + " user " + to_string(cur_usr);
                    G.supply[u] = 0;
                    Arc e = G.G.addArc(u, time_interval_nodes[cur_usr][nl-begining_end_time_interval[cur_usr].first]);   // edge between user group and time interval
                    G.up_bound[e] = len_time[nl];
                    G.costs[e] = 1 / (double) usr_grp_cur_time[nu].size();
                    G.arc_label[e] = "user " + to_string(cur_usr) + " user group " + key_ug + " to time " + to_string(nl);
                    G.gamma[e] = G.costs[e];
                    G.flow[e] = 0;
                    arcs_cur_usr_grp.push_back(e);
                    vector<string> miss_subfile = subfile_carried_by_usr_grp(usr_grp_cur_time[nu], cur_usr);
                    for(int nf = 0; nf < miss_subfile.size(); nf++){
                        int nf_id = subfile_id[miss_subfile[nf]][cur_usr];
                        Arc e_f = G.G.addArc(subfile_nodes[cur_usr][nf_id], u); // edge from missing subfile to current user group
                        G.arc_label[e_f] = "user " + to_string(cur_usr) + " subfile {" + miss_subfile[nf] + "} " +  " user group " + key_ug;
                        G.costs[e_f] = 0;
                        G.up_bound[e_f] = r;
                        G.gamma[e_f] = 0;
                        G.flow[e_f] = 0;
                }
                usr_group_nodes[key_ug] = tmp_node_usr_grp;
            }

            }
            usr_time_arcs.push_back(arcs_cur_usr_grp);
        }
    }
}

void ProblemInstance::ConstructProblem(){
    
    //  generate users vs interval time matrix
    gen_usr_vs_time(arrival_times, deadlines, usr_time, len_time);
    beta = len_time.size(); //  number of time intervals
    //  generate subfile id
    gen_subfile_id();
    //  generate source and sinck nodes
    gen_source_sink_nodes();
    //  generate subfile nodes
    gen_subfile_nodes();
    //  generate time interval nodes
    gen_time_interval_nodes();
    //  generate user group nodes and all edges connected to them
    gen_usr_grp_nodes();

}

//  print information about nodes and edges of the constructed graph, shall be run after constructing the graph
void ProblemInstance::PrintNodes(){
    cout << "===========\t Printing Nodes \t===============" << endl;
    cout << "Node id\t" << "label\t\t" << "supply\t" << endl;
    for(ListDigraph::NodeIt v(G.G); v != INVALID; ++v){
        cout << G.G.id(v)  << " \t" <<  G.node_label[v] << "\t" << G.supply[v] << endl;
    }
}
void ProblemInstance::PrintArcs(){
    cout << "===========\t Printing Edges \t===============" << endl;
    cout << "Edge id\t" << "label\t\t\t" << "upper bound\t" << "cost\t" << "gamma\t" << "flow\t" << endl;
    for(ListDigraph::ArcIt e(G.G); e != INVALID; ++e){
        cout << G.G.id(e) << "\t" << G.arc_label[e] << "\t" << G.up_bound[e] << "\t" << G.costs[e] << "\t" << G.gamma[e] << "\t" << G.flow[e] << endl;
    }
    
}
void ProblemInstance::PrintOutput(){
    for(int nu = 0; nu < usr_time_arcs.size(); nu++){
        for(int nk = 0; nk < usr_time_arcs[nu].size(); nk++){
            Arc e = usr_time_arcs[nu][nk];
            if(G.flow[e] > 0.05) cout << "flow of edge " << G.arc_label[e] << " is " << G.flow[e] << endl;
        }
    }
}
void ProblemInstance::PrintGraph(){
    cout << "===========\t Printing Nodes \t===============" << endl;
    cout << "Node id\t" << "label\t\t" << "supply\t" << endl;
    for(ListDigraph::NodeIt v(G.G); v != INVALID; ++v){
        cout << G.G.id(v)  << " \t" <<  G.node_label[v] << "\t" << G.supply[v] << endl;
    }
    cout << "===========\t Printing Edges \t===============" << endl;
    cout << "Edge id\t" << "label\t\t\t" << "upper bound\t" << "cost\t" << "gamma\t" << "flow\t" << endl;
    for(ListDigraph::ArcIt e(G.G); e != INVALID; ++e){
        cout << G.G.id(e) << "\t" << G.arc_label[e] << "\t" << G.up_bound[e] << "\t" << G.costs[e] << "\t" << G.gamma[e] << "\t" << G.flow[e] << endl;
    }
    cout << "===========\t Printing User Groups \t===============" << endl;
    for(int nu = 0; nu < usr_time_arcs.size(); nu++){
        for(int ni = 0; ni < usr_time_arcs[nu].size(); ni++){
            cout << G.arc_label[usr_time_arcs[nu][ni]] << endl;
        }
        cout << "*********************************************" << endl;
    }
}

//  Project gamma
void ProblemInstance::project_gamma(vector<double>& tmp_gamma){
    vector<size_t> sorted_idx = sort_indexes(tmp_gamma);
    int k_hat = 1;
    double a_u = tmp_gamma[sorted_idx[0]];
    while (k_hat < tmp_gamma.size() && (1 - a_u > -k_hat * tmp_gamma[sorted_idx[k_hat]])) {
//        cout << a_u << "\t" << k_hat << endl;
        a_u += tmp_gamma[sorted_idx[k_hat++]];
    }
//    cout << "k_hat:\t" << k_hat << "\t a_u:\t" << a_u << endl;
    for(int nk = 0; nk < tmp_gamma.size(); nk++)
    {
//        cout << sorted_idx[nk] << endl;
        if(nk < k_hat) tmp_gamma[sorted_idx[nk]] += (1 - a_u) / k_hat;
        else tmp_gamma[sorted_idx[nk]] = 0;
    }
    
}

//  Reading the arrivale and deadlines from file
void readFromFile(vector<int>& arriv, vector<int>& dead){
    ifstream file("data.txt");
    string str;
    int count = 0;
    while (getline(file, str)) {
        int tmp = stoi(str, 0, 10);
        if(count % 2 == 0) arriv.push_back(tmp);
        else dead.push_back(tmp);
        count++;
    }
    file.close();
}
//  Solving the constructed problem instance
void ProblemInstance::SolveProblem(int Tmaxitr, double alpha, int offItrs){
    ofstream out_file;
    out_file.open("out_K20_t2.txt");
    out_file << "Itr\tRateP\t\tTime\t\tRc\n";
    int Rc = 1969; //nchoosek(K, t+1);
    int exe_time = 0;
    MCF_Solver mcfSolver(G);
    double step_size;
    for(int na = 0; na < Tmaxitr; na++){
        auto start_mcf = high_resolution_clock::now();
        //  solve MCF
        bool isFeasible = mcfSolver.run();
        if(! isFeasible) cout << "Problem is not feasible" << endl;
        //  update costs and gamma
        primal_solution = 0;
        dual_solution = 0;
        double primal_non_xonvex_sol = 0;
        step_size = pow(na+1, -alpha);
        vector<double> gradient_zeta(beta, 0);
        for(int nu = 0; nu < usr_time_arcs.size(); nu++){
            double tmp_primal_value = 0;
            double tmp_primal_value_non_convex = 0;
            int cur_time_interval = time_interval_usr_grp[nu];
            vector<double> unprojected_gamma;
            for(int nk = 0; nk < usr_time_arcs[nu].size(); nk++){
                Arc e = usr_time_arcs[nu][nk];
                double tmp_cost = (1+zeta[cur_time_interval]) * G.gamma[e];
                double tmp_flow = mcfSolver.get_flow(e);
                G.costs[e] = tmp_cost;
                dual_solution += tmp_cost * tmp_flow;
                G.gamma[e] = G.gamma[e] + step_size * (1+zeta[cur_time_interval]) * tmp_flow;
                unprojected_gamma.push_back(G.gamma[e]);
                if(na >= offItrs) G.flow[e] = ((na-offItrs) * G.flow[e] + tmp_flow) / ((na-offItrs)+ 1) ;
                //if(tmp_flow > 0) cout << "flow of edge " << G.arc_label[e] << " is " << tmp_flow << " and primal " << G.flow[e] << endl;
                tmp_primal_value = (tmp_primal_value < G.flow[e]) ? G.flow[e] : tmp_primal_value;
//                tmp_primal_value_non_convex = (tmp_primal_value_non_convex < tmp_flow) ? tmp_flow : tmp_primal_value_non_convex;
            }
            //  project gamma
            project_gamma(unprojected_gamma);
            for(int nk = 0; nk < usr_time_arcs[nu].size(); nk++){
                Arc e = usr_time_arcs[nu][nk];
                G.gamma[e] = unprojected_gamma[nk];
                gradient_zeta[cur_time_interval] += G.gamma[e] * mcfSolver.get_flow(e);
            }
//            primal_non_xonvex_sol += tmp_primal_value_non_convex;
            primal_solution += tmp_primal_value;
        }
        //  update zeta
        for(int nl = 0; nl < beta; nl++){
            gradient_zeta[nl] -= len_time[nl];
            zeta[nl] = zeta[nl] + step_size * gradient_zeta[nl];
            zeta[nl] = (zeta[nl] > 0) ? zeta[nl] : 0 ;
            dual_solution -= zeta[nl] * len_time[nl];
        }
        mcfSolver.run_cost(G.costs);
        auto stop_mcf = high_resolution_clock::now();
        auto duration_mcf = duration_cast<microseconds>(stop_mcf - start_mcf);
        exe_time += duration_mcf.count();
        if(na % 10 == 0){
            out_file << na-offItrs << "\t" << primal_solution << "\t" << exe_time  << "\t" << Rc << "\n";
            cout << "***********************************" << endl;
            cout << "iteration #: " << na-offItrs << "\tPrimal Solution:\t " << primal_solution << "\tDual Solution:\t " << dual_solution << endl;
        }
        
        // PrintGraph();
        
    }
    
}
