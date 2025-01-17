#include <random>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <list>
#include <cmath>
#include <bits/stdc++.h>
#include <fstream>
#include <typeinfo>
#include <string>



using namespace std;

// Global constants

// ---- PARAMETERS -------------------------------------------------------------
const int V = 100;              // Total number of nodes
const int v = V-2;              // Nodes that can fail (source and sink don't fail)
const int m1 = 28;              // Number of nodes in the first group/class of nodes
const int m2 = 17;
const int m3 = 15;
const int m4 = 18;              // Number of nodes in the second group/class of nodes
const int m5 = v-(m1+m2+m3+m4);      // Number of nodes in the third group/class of nodes
const int M = 100;             // Number of MC replications

// -----------------------------------------------------------------------------

// -------- Data Structures ----------------------------------------------------
double survSig[m1+1][m2+1][m3+1][m4+1][m5+1];
bool survSig_filing[m1+1][m2+1][m3+1][m4+1][m5+1];

list<array<int, 8>> ndpoints[V];
int Nindex[V] = { 0 };
int d[V][5] = { 0 };
bool InH[V] = { false };

#define N 8

int c[V][5];
int cij[V][5];
int cji[V][5];


// ---- CLASSES AND FUNCTIONS -----------------------------------------------

class Graph
{
public:
    int V;

    list<int> *adj;
    list<int> *pred;

    Graph(int V);

    void addEdge(int v, int w);

    bool BFS(int s, vector<int> inactive, int v);

    list<int> succ_list (int x);

    list<int> pred_list (int x);

    void deleteLists();

    ~Graph();
};

class MaxHeap{
private:
    int _size{};
    vector<array<int, N>> vect = {{-1,0,0,0,0,0,0,0}};;

    int p(int i) {return i >> 1;};
    int l(int i) {return i << 1;};
    int r(int i) {return (i << 1) + 1;};

public:

    bool isEmpty() const {return _size == 0;};

    array<int, N> find_max() const {return vect.at(1);};

    void shiftUp (int i);

    void shiftDown (int i);

    void insert (array<int, N> lab);

    void delete_max ();

    void increase_key (array<int, N> lab );

    void print();

    void print_lable(array<int, N> x);

    ~MaxHeap();

};


void node_labeling(vector<int>& group1, vector<int>& group2, vector<int>& group3, vector<int>& group4, vector<int>& group5, Graph g);

bool NewCandLabel (Graph g, int i, array<int, N> l_star, int r, array<int, 8> &l_new, MaxHeap* Heap);

void RelaxProcess (Graph g, int i, array<int, N> l_star, MaxHeap* Heap);

vector<array<int, 7>> lable;

void print_signature(Graph g, int m1, int m2, int m3, int m4, int m5);

vector<int> sample_WO (vector<int> m, int sizeComb, int j);

void BO_maxCapPath (vector<int> perm1, vector<int> perm2, vector<int> perm3, vector<int> perm4, vector<int> perm5, Graph g, MaxHeap* Heap);

void signature_estimation();


void graphGen_RGG(Graph & g, int V, int m1, int m2, int m3, int m4, int m5);

void power2000(Graph & g, int V, int m1, int m2, int m3, int m4, int m5);


// ***********************************************************************************


clock_t s1, s2, s3, s4, s5;

// ------------------ MAIN ----------------------------------------------------------
int main()
{
    s1 = clock();

    lable.reserve(v);
    MaxHeap* Heap = new MaxHeap();

    string net = "RGG 10";
    Graph g(V);


// Either construct a graph by adding edges (in case of a small network) or
// run a function to add edges.

// function graphGen_RGG generates a RGG with m1 nodes in the first class, m2 nodes in the second class and m3 nodes in the third class;

    //graphGen_RGG(g, V, m1, m2, m3);


// function power2000 generates the Power 2000 bus electric grid network. Make sure V=4000 since in this function
// all edges are hard-coded.

    power2000(g, V, m1, m2, m3, m4, m5);


// Either assign the nodes to each class "manually" or following the piece of code below

    std::vector<int> class_1 = {1, 2, 3, 4, 5, 6, 7, 11, 16, 17, 36, 37, 38, 39, 43, 44, 45, 53, 54, 55, 56, 57, 58, 59, 69, 73, 74, 75};
    std::vector<int> class_2 = {8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 23, 24, 25, 26, 27, 76, 83};
    std::vector<int> class_3 = {90, 93, 94, 95, 96, 97, 98, 63, 64, 65, 66, 67, 68, 84, 85};
    std::vector<int> class_4 = {21, 41, 51, 61, 71, 91, 46, 47, 48, 49, 86, 87, 88, 89, 18, 19, 28, 35};
    std::vector<int> class_5 = {22, 32, 42, 52, 62, 72, 82, 92, 29, 31, 33, 34, 77, 78, 79, 81, 12, 13, 14, 15}; 





// Either assign the nodes to each class "manually" or following the piece of code below

    vector<int> group_1;
    vector<int> group_2;
    vector<int> group_3;
    vector<int> group_4;
    vector<int> group_5;


// Assigns the first m1 even nodes to the first class and the remaining m2 nodes to the second class


    int ct1 = 0, ct2 = 0, ct3 = 0, ct4 = 0, ct5 = 0;

    for (int i = 0; i < class_1.size() && ct1 < m1; ++i) {
        group_1.push_back(class_1[i]);
        ++ct1;
    }

    for (int i = 0; i < class_2.size() && ct2 < m2; ++i) {
        group_2.push_back(class_2[i]);
        ++ct2;
    }

    for (int i = 0; i < class_3.size() && ct3 < m3; ++i) {
        group_3.push_back(class_3[i]);
        ++ct3;
    }
    for (int i = 0; i < class_4.size() && ct4 < m4; ++i) {
        group_4.push_back(class_4[i]);
        ++ct4;
    }
    for (int i = 0; i < class_5.size() && ct5 < m5; ++i) {
        group_5.push_back(class_5[i]);
        ++ct5;
    }

    vector<int> perm1;
    vector<int> perm2;
    vector<int> perm3;
    vector<int> perm4;
    vector<int> perm5;

    s2 = clock();


    for (int j=1; j <= M; j++)
    {
       

        perm1 = sample_WO(group_1, group_1.size(), j);
		//perm1.clear();
		//perm1.push_back(7);
		//perm1.push_back(4);
		//perm1.push_back(1);
		//perm1.push_back(10);

        perm2 = sample_WO(group_2, group_2.size(), j);
		//perm2.clear();
		//perm2.push_back(2);
		//perm2.push_back(5);
		//perm2.push_back(8);
		//perm2.push_back(11);
		
        perm3 = sample_WO(group_3, group_3.size(), j);
        //perm3.clear();
		//perm3.push_back(12);
		//perm3.push_back(6);
		//perm3.push_back(9);
		//perm3.push_back(3);
        
        perm4 = sample_WO(group_4, group_4.size(), j);

        perm5 = sample_WO(group_5, group_5.size(), j);
        
		

        

        // solves the bi-objective max capacity path problem
        BO_maxCapPath (perm1, perm2, perm3, perm4, perm5, g, Heap);




        // estimates thesi survival signature
        signature_estimation();

    }

    s3 = clock();

    print_signature(g, m1, m2, m3, m4, m5);

    s4 = clock();

    // Reliability estimation present numerical issues for large networks due to the very large and very small
    // binomial coefficients, result in "NaN" values.

    // reliability_function (m1, m2, m3);

    s5 = clock();

    cout << endl;
    cout << endl <<"MO Algorithm - Run-time Performance:" <<endl;
    cout << "Network system: " << net << ", m1 = " << m1 << ", m2 = " << m2 << ", m3 = " << m3 << ", m4 = " << m4 << ", M = " << M <<endl<<endl;

    cout << "Initialization: " << double(s2 - s1) / double(CLOCKS_PER_SEC) << setprecision(6) << endl << endl;
    cout << "Survival Signature estimation: " << double(s3 - s2) / double(CLOCKS_PER_SEC) << setprecision(6) << endl <<endl;
    cout << "Survival Signature rate: " << (double(s3 - s2) / double(CLOCKS_PER_SEC))/(double(M)) << setprecision(6) << endl <<endl;
    cout << "Printing Signature: " << (double(s4 - s3) / double(CLOCKS_PER_SEC)) << setprecision(6) << endl <<endl;
//    cout << "Reliability estimation: " << double(s4 - s3) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;
//    cout << "Total Time: " << double(s4 - s1) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;


    g.deleteLists();
    delete Heap;
    return 0;
}
// *************************************************************************





/*  ~~~~~~~ Class Graph  ~~~~~~~~~~~ */


Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
    pred = new list<int>[V];
}

void Graph::addEdge(int u, int w)
{
    adj[u].push_back(w);
    pred[w].push_back(u);
}


bool Graph::BFS(int s, vector<int> inactive, int v)
{

    int nodes = 0;

    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
        visited[i] = false;

    list<int> queue;

    visited[s] = true;
    queue.push_back(s);

    list<int>::iterator i;

    while(!queue.empty())
    {

        s = queue.front();

        if(find(inactive.begin(), inactive.end(), s) != inactive.end()){

                queue.pop_front();
                continue;
        }

//      cout << s << " ";
        if (s == V-1)
        {
            nodes = true;
        }
        queue.pop_front();

        for (i = adj[s].begin(); i != adj[s].end(); ++i)
        {
            if (!visited[*i])
            {
                visited[*i] = true;
                queue.push_back(*i);
            }
        }
    }

    delete[] visited;
    return nodes;
}


list<int> Graph::succ_list (int x)
{
    return adj[x];
}

list<int> Graph::pred_list (int x)
{
    return pred[x];
}

void Graph:: deleteLists ()
{
    delete [] adj;
    delete [] pred;
}

Graph::~Graph()
{

}



void MaxHeap::shiftUp(int i){
    if (i > _size) return;
    if (i == 1) return;
    if (vect.at(i)[1] > vect.at(p(i))[1] || (vect.at(i)[1] == vect.at(p(i))[1] && vect.at(i)[2] > vect.at(p(i))[2]) || (vect.at(i)[1] == vect.at(p(i))[1] && vect.at(i)[2] == vect.at(p(i))[2] && vect.at(i)[3] > vect.at(p(i))[3]) || (vect.at(i)[1] == vect.at(p(i))[1] && vect.at(i)[2] == vect.at(p(i))[2] && vect.at(i)[3] == vect.at(p(i))[3] && vect.at(i)[4] > vect.at(p(i))[4]) || (vect.at(i)[1] == vect.at(p(i))[1] && vect.at(i)[2] == vect.at(p(i))[2] && vect.at(i)[3] == vect.at(p(i))[3] && vect.at(i)[4] == vect.at(p(i))[4] && vect.at(i)[5] > vect.at(p(i))[5])){
        std::swap(vect.at(p(i)), vect.at(i));
    }
    shiftUp(p(i));
}

void MaxHeap::insert(array<int, N> lab){

    if (_size + 1 >= vect.size()){
        vect.push_back({0, 0, 0, 0, 0, 0, 0, 0});
    }
    vect[++_size] = lab;
    shiftUp(_size);
    return;

}

void MaxHeap::shiftDown(int i) {
    
    if (i > _size) return;

    int swapId = i;

    if (l(i) <= _size && (vect.at(i)[1] < vect.at(l(i))[1] || (vect.at(i)[1] == vect.at(l(i))[1] && vect.at(i)[2] < vect.at(l(i))[2]) || (vect.at(i)[1] == vect.at(l(i))[1] && vect.at(i)[2] == vect.at(l(i))[2] && vect.at(i)[3] < vect.at(l(i))[3]) || (vect.at(i)[1] == vect.at(l(i))[1] && vect.at(i)[2] == vect.at(l(i))[2] && vect.at(i)[3] == vect.at(l(i))[3] && vect.at(i)[4] < vect.at(l(i))[4]) || (vect.at(i)[1] == vect.at(l(i))[1] && vect.at(i)[2] == vect.at(l(i))[2] && vect.at(i)[3] == vect.at(l(i))[3] && vect.at(i)[4] == vect.at(l(i))[4] && vect.at(i)[5] < vect.at(l(i))[5])) ) {
        swapId = l(i);
    }

    if (r(i) <= _size && (vect.at(swapId)[1] < vect.at(r(i))[1] || (vect.at(swapId)[1] == vect.at(r(i))[1] && vect.at(swapId)[2] < vect.at(r(i))[2]) || (vect.at(swapId)[1] == vect.at(r(i))[1] && vect.at(swapId)[2] == vect.at(r(i))[2] && vect.at(swapId)[3] < vect.at(r(i))[3]) || (vect.at(swapId)[1] == vect.at(r(i))[1] && vect.at(swapId)[2] == vect.at(r(i))[2] && vect.at(swapId)[3] == vect.at(r(i))[3] && vect.at(swapId)[4] < vect.at(r(i))[4]) || (vect.at(swapId)[1] == vect.at(r(i))[1] && vect.at(swapId)[2] == vect.at(r(i))[2] && vect.at(swapId)[3] == vect.at(r(i))[3] && vect.at(swapId)[4] == vect.at(r(i))[4] && vect.at(swapId)[5] < vect.at(r(i))[5])) ) {
        swapId = r(i);
    }

    if (swapId != i) {
        std::swap(vect.at(i), vect.at(swapId));
        shiftDown(swapId);
    }

    return;

}

void MaxHeap::delete_max(){
    array<int, N> maxNum = vect.at(1);
    std::swap(vect.at(1), vect.at(_size--));
    shiftDown(1);
    vect.pop_back();
    return;
}

void MaxHeap::increase_key(array<int, N> lab){
    int s;
    for (int i=0; i < vect.size(); ++i)
    {
        if (vect.at(i)[0] == lab[0])
        {
            vect.at(i)[1] = lab[1];
            vect.at(i)[2] = lab[2];
            vect.at(i)[3] = lab[3];
            vect.at(i)[4] = lab[4];
            vect.at(i)[5] = lab[5];
            vect.at(i)[6] = lab[6];
            vect.at(i)[7] = lab[7];
            s = i;
        }
    }

    shiftUp(s);
}




MaxHeap::~MaxHeap(){

}


bool NewCandLabel (Graph g, int i, array<int, N> l_star, int r, array<int, N> & l_new, MaxHeap* Heap)
{
    bool key = false;

    array<int, 5> f;
    array<int, N> l;
    array<int, N> cl;

    list<int> pred;

    array<int, 5> d;

    d[0] = 0;
    d[1] = 0;
    d[2] = 0;
    d[3] = 0;
    d[4] = 0;

    pred = g.pred_list(i);

    for (list<int>::iterator j = pred.begin(); j != pred.end(); j++)
    {

        for (list<array<int, N>>::iterator lj=ndpoints[*j].begin(); lj != ndpoints[*j].end(); lj++)
        {

                f[0] = min(c[*j][0], c[i][0]);
                f[1] = min(c[*j][1], c[i][1]);
                f[2] = min(c[*j][2], c[i][2]);
                f[3] = min(c[*j][3], c[i][3]);
                f[4] = min(c[*j][4], c[i][4]);

                l=*lj;

                f[0] = min(f[0], l[1]);
                f[1] = min(f[1], l[2]);
                f[2] = min(f[2], l[3]);
                f[3] = min(f[3], l[4]);
                f[4] = min(f[4], l[5]);

                
                if (f[0] > d[0] || (f[0] == d[0] && f[1] > d[1]) || (f[0] == d[0] && f[1] == d[1] && f[2] > d[2]) || (f[0] == d[0] && f[1] == d[1] && f[2] == d[2] && f[3] > d[3])  || (f[0] == d[0] && f[1] == d[1] && f[2] == d[2] && f[3] == d[3] && f[4] > d[4])){
                   
                    bool dominated = false;
   
                    for (list<array<int, N>>::iterator iter=ndpoints[i].begin(); iter != ndpoints[i].end(); iter++){
                        cl = *iter;
                        int first_count = 0;
                        int second_count = 0;
                        int third_count = 0;


                        for (int k = 0; k < 5; k++) {

                            if (f[k] > cl[k+1]){
                                first_count += 1;
                            }
                            else if (f[k] < cl[k+1]){
                                second_count += 1;
                            }
                            else if (f[k] == cl[k+1]){
                                third_count += 1;
                            }
                        }
                        if (first_count == 0){
                            dominated = true;
                            break;
                        }
                    }
                    if (dominated == false){
                        d[0] = f[0];
                        d[1] = f[1];
                        d[2] = f[2];
                        d[3] = f[3];
                        d[4] = f[4];
                        l_new[0] = i;
                        l_new[1] = d[0];
                        l_new[2] = d[1];
                        l_new[3] = d[2];
                        l_new[4] = d[3];
                        l_new[5] = d[4];
                        l_new[6] = *j;
                        l_new[7] = Nindex[*j];
                        key = true;
                    }
                }
        }
    }

    return key;
}

void node_labeling (vector<int>& perm1, vector<int>& perm2, vector<int>& perm3, vector<int>& perm4, vector<int>& perm5, Graph g)
{

    vector<int>::iterator it;

    for(int a = 0; a < V; ++a)
    {
        if (a == 0)
        {
            c[0][0] = V;
            c[0][1] = V;
            c[0][2] = V;
            c[0][3] = V;
            c[0][4] = V;
        }

        else if(find(perm1.begin(), perm1.end(), a) != perm1.end())
        {
            it = find(perm1.begin(), perm1.end(), a);
            c[a][0] = (it - perm1.begin()) + 1;
            c[a][1] = V;
            c[a][2] = V;
            c[a][3] = V;
            c[a][4] = V;
        }

        else if (find(perm2.begin(), perm2.end(), a) != perm2.end())
        {
            it = find(perm2.begin(), perm2.end(), a);
            c[a][0] = V;
            c[a][1] = (it - perm2.begin()) + 1;
            c[a][2] = V;
            c[a][3] = V;
            c[a][4] = V;
        }

        else if (find(perm3.begin(), perm3.end(), a) != perm3.end())
        {
            it = find(perm3.begin(), perm3.end(), a);
            c[a][0] = V;
            c[a][1] = V;
            c[a][2] = (it - perm3.begin()) + 1;
            c[a][3] = V;
            c[a][4] = V;
        }

        else if (find(perm4.begin(), perm4.end(), a) != perm4.end())
        {
            it = find(perm4.begin(), perm4.end(), a);
            c[a][0] = V;
            c[a][1] = V;
            c[a][2] = V;
            c[a][3] = (it - perm4.begin()) + 1;
            c[a][4] = V;
        }


        else if (find(perm4.begin(), perm4.end(), a) != perm4.end())
        {
            it = find(perm4.begin(), perm4.end(), a);
            c[a][0] = V;
            c[a][1] = V;
            c[a][2] = V;
            c[a][3] = V;
            c[a][4] = (it - perm5.begin()) + 1;
        }

        else
        {
            c[a][0] = V;
            c[a][1] = V;
            c[a][2] = V;
            c[a][3] = V;
            c[a][4] = V;
        }
    }

}


void RelaxProcess (Graph g, int i, array<int, N> l_star, MaxHeap* Heap)
{
    array<int, 5> f;
    //array<int, N> l;
    array<int, N> l_new;
    array<int, N> cl;

    list<int> succ;

    succ = g.succ_list(i);

    for(list<int>::iterator j = succ.begin(); j != succ.end(); j++)
    {

        f[0] = min(c[*j][0], c[i][0]);
        f[1] = min(c[*j][1], c[i][1]);
        f[2] = min(c[*j][2], c[i][2]);
        f[3] = min(c[*j][3], c[i][3]);
        f[4] = min(c[*j][4], c[i][4]);

        f[0] = min(f[0], l_star[1]);
        f[1] = min(f[1], l_star[2]);
        f[2] = min(f[2], l_star[3]);
        f[3] = min(f[3], l_star[4]);
        f[4] = min(f[4], l_star[5]);


        if (f[0] > d[*j][0] || (f[0] == d[*j][0] && f[1] > d[*j][1]) || (f[0] == d[*j][0] && f[1] == d[*j][1] && f[2] > d[*j][2]) || (f[0] == d[*j][0] && f[1] == d[*j][1] && f[2] == d[*j][2] && f[3] > d[*j][3]) || (f[0] == d[*j][0] && f[1] == d[*j][1] && f[2] == d[*j][2] && f[3] == d[*j][3] && f[4] > d[*j][4]))
        {
            //l = ndpoints[*j].back();

            bool dominated = false;

            for (list<array<int, N>>::iterator iter=ndpoints[*j].begin(); iter != ndpoints[*j].end(); iter++){

                cl = *iter;
                int rp_counter_1 = 0;
                int rp_counter_2 = 0;
                int rp_counter_3 = 0;

                for  (int k = 0; k < 5; k++){
                    if (f[k] > cl[k+1]){
                        rp_counter_1 += 1;
                    }
                    if (f[k] < cl[k+1]){
                        rp_counter_2 += 1;
                    }
                    if (f[k] == cl[k+1]){
                        rp_counter_3 += 1;
                    }
                }

                if (rp_counter_1 == 0 ){
                    dominated = true;
                    break;
                }
            }

            if (Nindex[*j] == 0 || dominated == false)
            {
                d[*j][0] = f[0];
                d[*j][1] = f[1];
                d[*j][2] = f[2];
                d[*j][3] = f[3];
                d[*j][4] = f[4];

                l_new[0] = *j;
                l_new[1] = d[*j][0];
                l_new[2] = d[*j][1];
                l_new[3] = d[*j][2];
                l_new[4] = d[*j][3];
                l_new[5] = d[*j][4];
                l_new[6] = i;
                l_new[7] = Nindex[i];

                if (InH[*j] == false)
                {
                    Heap->insert(l_new);
                    InH[*j] = true;
                }
                else
                {
                    Heap->increase_key(l_new);
                }
            }
        }
    }
}



vector<int> sample_WO (vector<int> m, int sizeComb, int j)
{
    random_device myRandomDevice;
    unsigned seed = j; //myRandomDevice(); // 2*j; (when generating the same random numbers)
	// unsigned seed = chrono::system_clock::now().time_since_epoch().count(); \\(Windows)
    default_random_engine myRandomEngine(j);

    shuffle(m.begin(), m.end(), myRandomEngine);
    vector<int> result(m.begin(), m.begin() + sizeComb);

    return result;

}


void BO_maxCapPath(vector<int> perm1, vector<int> perm2, vector<int> perm3, vector<int> perm4, vector<int> perm5, Graph g, MaxHeap* Heap)
{
    array<int, N> l_star;
    array<int, 8> l_new ;

    for(int i=0; i<V; i++)
        ndpoints[i].clear();

    int i;

    fill(&d[0][0], &d[0][0] + V*5, 0);
    fill(&Nindex[0], &Nindex[0] + V, 0);
    fill(&InH[0], &InH[0] + V, false);

    node_labeling (perm1, perm2, perm3, perm4, perm5, g);

    d[0][0] = V;
    d[0][1] = V;
    d[0][2] = V;
    d[0][3] = V;
    d[0][4] = V;

    array<int, N> l = {{0, d[0][0], d[0][1], d[0][2], d[0][3], d[0][4], 0, 0}};

    Heap->insert(l);

    InH[0] = true;

    while (!Heap->isEmpty())
    {

        l_star = Heap->find_max();
        Heap->delete_max();

        i = l_star[0];

        d[i][0] = 0;
        d[i][1] = 0;
        d[i][2] = 0;
        d[i][3] = 0;
        d[i][4] = 0;

        Nindex[i] = Nindex[i] + 1;

        ndpoints[i].push_back(l_star);

        InH[i] = false;

        if(NewCandLabel(g, i, l_star, Nindex[i], l_new, Heap))
        {

            Heap->insert(l_new);
            InH[i] = true;
            d[i][0] = l_new[1];
            d[i][1] = l_new[2];
            d[i][2] = l_new[3];
            d[i][3] = l_new[4];
            d[i][4] = l_new[5];
        }
        

        RelaxProcess (g, i, l_star, Heap);
    }
}

void signature_estimation()
{
    array<int, N> l;
    int l1, l2, l3, l4, l5;


    fill(&survSig_filing[0][0][0][0][0], &survSig_filing[0][0][0][0][0] + (m1+1)*(m2+1)*(m3+1)*(m4+1)*(m5+1), false);

   for (list<array<int, N>>::iterator lj=ndpoints[V-1].begin(); lj != ndpoints[V-1].end(); lj++)
   {
        l = *lj;
        if (l[1]<=m1)
           l1 = m1 - l[1];
        else
           l1 = -1;

        if (l[2]<= m2)
           l2 = m2 - l[2];
        else
           l2 = -1;

        if (l[3]<= m3)
           l3 = m3 - l[3];
        else
           l3 = -1;

        if (l[4]<= m4)
           l4 = m4 - l[4];
        else
           l4 = -1;

        if (l[5]<= m5)
           l5 = m5 - l[5];
        else
           l5 = -1;

        for (int i = 0; i <= m1; i++)
        {
          for (int j = 0; j <= m2; j++)
          {
             for (int k = 0; k <= m3; k++)
             {
                for (int y = 0; y <= m4; y++)
                {
                    for (int x = 0; x <= m5; x++)
                    {
                        if (i > l1 && j > l2 && k > l3 && y > l4)
                            survSig_filing[i][j][k][y][x] = true;
                    }
                }
             }
          }
        }
    }


  for(int i=0; i<=m1; i++)
  {
    for(int j=0; j<=m2; j++)
    {
        for(int k=0; k<=m3; k++)
        {
            for(int y=0; y<=m4; y++)
            {
                for (int x = 0; x <= m5; x++)
                {
                    if (survSig_filing[i][j][k][y][x] == true)
                        survSig[i][j][k][y][x] = survSig[i][j][k][y][x] + 1;
                }
            }
        }
    }
   }
}


void print_signature(Graph g, int m1, int m2, int m3, int m4, int m5)
{
    ofstream myfile;
	myfile.open("signature.txt");

    cout <<endl<<endl;
    for(int l1=0; l1<=m1; l1++)
    {
        for(int l2=0; l2<=m2; l2++)
        {
            for(int l3=0; l3<=m3; l3++)
            {
                for(int l4=0; l4<=m4; l4++)
                {
                    for(int l5=0; l5<=m5; l5++)
                    {
                        survSig[l1][l2][l3][l4][l5] = survSig[l1][l2][l3][l4][l5]/(double) M;
                        cout << setprecision(6) << survSig[l1][l2][l3][l4][l5] << " ";
                        myfile << setprecision(6) << survSig[l1][l2][l3][l4][l5] << " ";
                    }
                }
                    cout << endl;
                    myfile << endl;
            }
            cout << "\n";  
        }
        cout << "\n";
    }
    cout <<endl<<endl;
    myfile.close();
}





void power2000(Graph & g, int V, int m1, int m2, int m3, int m4, int m5)
{   
    // Adding edges
g.addEdge(0, 1);
g.addEdge(0, 10);
g.addEdge(1, 2);
g.addEdge(1, 11);
g.addEdge(2, 3);
g.addEdge(2, 12);
g.addEdge(3, 4);
g.addEdge(3, 13);
g.addEdge(4, 5);
g.addEdge(4, 14);
g.addEdge(5, 6);
g.addEdge(5, 15);
g.addEdge(6, 7);
g.addEdge(6, 16);
g.addEdge(7, 8);
g.addEdge(7, 17);
g.addEdge(8, 9);
g.addEdge(8, 18);
g.addEdge(9, 19);
g.addEdge(10, 11);
g.addEdge(10, 20);
g.addEdge(11, 12);
g.addEdge(11, 21);
g.addEdge(12, 13);
g.addEdge(12, 22);
g.addEdge(13, 14);
g.addEdge(13, 23);
g.addEdge(14, 15);
g.addEdge(14, 24);
g.addEdge(15, 16);
g.addEdge(15, 25);
g.addEdge(16, 17);
g.addEdge(16, 26);
g.addEdge(17, 18);
g.addEdge(17, 27);
g.addEdge(18, 19);
g.addEdge(18, 28);
g.addEdge(19, 29);
g.addEdge(20, 21);
g.addEdge(20, 30);
g.addEdge(21, 22);
g.addEdge(21, 31);
g.addEdge(22, 23);
g.addEdge(22, 32);
g.addEdge(23, 24);
g.addEdge(23, 33);
g.addEdge(24, 25);
g.addEdge(24, 34);
g.addEdge(25, 26);
g.addEdge(25, 35);
g.addEdge(26, 27);
g.addEdge(26, 36);
g.addEdge(27, 28);
g.addEdge(27, 37);
g.addEdge(28, 29);
g.addEdge(28, 38);
g.addEdge(29, 39);
g.addEdge(30, 31);
g.addEdge(30, 40);
g.addEdge(31, 32);
g.addEdge(31, 41);
g.addEdge(32, 33);
g.addEdge(32, 42);
g.addEdge(33, 34);
g.addEdge(33, 43);
g.addEdge(34, 35);
g.addEdge(34, 44);
g.addEdge(35, 36);
g.addEdge(35, 45);
g.addEdge(36, 37);
g.addEdge(36, 46);
g.addEdge(37, 38);
g.addEdge(37, 47);
g.addEdge(38, 39);
g.addEdge(38, 48);
g.addEdge(39, 49);
g.addEdge(40, 41);
g.addEdge(40, 50);
g.addEdge(41, 42);
g.addEdge(41, 51);
g.addEdge(42, 43);
g.addEdge(42, 52);
g.addEdge(43, 44);
g.addEdge(43, 53);
g.addEdge(44, 45);
g.addEdge(44, 54);
g.addEdge(45, 46);
g.addEdge(45, 55);
g.addEdge(46, 47);
g.addEdge(46, 56);
g.addEdge(47, 48);
g.addEdge(47, 57);
g.addEdge(48, 49);
g.addEdge(48, 58);
g.addEdge(49, 59);
g.addEdge(50, 51);
g.addEdge(50, 60);
g.addEdge(51, 52);
g.addEdge(51, 61);
g.addEdge(52, 53);
g.addEdge(52, 62);
g.addEdge(53, 54);
g.addEdge(53, 63);
g.addEdge(54, 55);
g.addEdge(54, 64);
g.addEdge(55, 56);
g.addEdge(55, 65);
g.addEdge(56, 57);
g.addEdge(56, 66);
g.addEdge(57, 58);
g.addEdge(57, 67);
g.addEdge(58, 59);
g.addEdge(58, 68);
g.addEdge(59, 69);
g.addEdge(60, 61);
g.addEdge(60, 70);
g.addEdge(61, 62);
g.addEdge(61, 71);
g.addEdge(62, 63);
g.addEdge(62, 72);
g.addEdge(63, 64);
g.addEdge(63, 73);
g.addEdge(64, 65);
g.addEdge(64, 74);
g.addEdge(65, 66);
g.addEdge(65, 75);
g.addEdge(66, 67);
g.addEdge(66, 76);
g.addEdge(67, 68);
g.addEdge(67, 77);
g.addEdge(68, 69);
g.addEdge(68, 78);
g.addEdge(69, 79);
g.addEdge(70, 71);
g.addEdge(70, 80);
g.addEdge(71, 72);
g.addEdge(71, 81);
g.addEdge(72, 73);
g.addEdge(72, 82);
g.addEdge(73, 74);
g.addEdge(73, 83);
g.addEdge(74, 75);
g.addEdge(74, 84);
g.addEdge(75, 76);
g.addEdge(75, 85);
g.addEdge(76, 77);
g.addEdge(76, 86);
g.addEdge(77, 78);
g.addEdge(77, 87);
g.addEdge(78, 79);
g.addEdge(78, 88);
g.addEdge(79, 89);
g.addEdge(80, 81);
g.addEdge(80, 90);
g.addEdge(81, 82);
g.addEdge(81, 91);
g.addEdge(82, 83);
g.addEdge(82, 92);
g.addEdge(83, 84);
g.addEdge(83, 93);
g.addEdge(84, 85);
g.addEdge(84, 94);
g.addEdge(85, 86);
g.addEdge(85, 95);
g.addEdge(86, 87);
g.addEdge(86, 96);
g.addEdge(87, 88);
g.addEdge(87, 97);
g.addEdge(88, 89);
g.addEdge(88, 98);
g.addEdge(89, 99);
g.addEdge(90, 91);
g.addEdge(91, 92);
g.addEdge(92, 93);
g.addEdge(93, 94);
g.addEdge(94, 95);
g.addEdge(95, 96);
g.addEdge(96, 97);
g.addEdge(97, 98);
g.addEdge(98, 99);


}
/*
void power2000(Graph & g, int V, int m1, int m2, int m3) {
    // Assuming V >= 14 for this example

    for (int i = 0; i < V; ++i) {
        for (int j = i + 1; j < V; ++j) {
            if (i != j) {
                // Add edges to connect nodes
                g.addEdge(i, j);
            }
        }
    }
}
*/