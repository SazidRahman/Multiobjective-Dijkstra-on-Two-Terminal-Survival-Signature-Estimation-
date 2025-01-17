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
const int V = 12;              // Total number of nodes
const int v = V-2;              // Nodes that can fail (source and sink don't fail)
const int m1 = 3;              // Number of nodes in the first group/class of nodes
const int m2 = 4;            // Number of nodes in the second group/class of nodes
const int m3 = v-(m1+m2);      // Number of nodes in the third group/class of nodes
const int M = 1000;              // Number of MC replications

// -----------------------------------------------------------------------------

// -------- Data Structures ----------------------------------------------------
double survSig[m1+1][m2+1][m3+1];
bool survSig_filing[m1+1][m2+1][m3+1];

list<array<int, 5>> ndpoints[V];
int Nindex[V] = { 0 };
int d[V][2] = { 0 };
bool InH[V] = { false };
bool Active[V] = { true };

#define N 6
#define N_BO 5

int c[V][2];
int cij[V][2];
int cji[V][2];


// ---- CLASSES AND FUNCTIONS -----------------------------------------------

class Graph
{
public:
    int V;

    list<int> *adj;
    list<int> *pred;
    bool *InGraph;
    bool *Active;

    Graph(int V);

    void addEdge(int v, int w);

    bool BFS(int s, vector<int> inactive, int v);

    list<int> succ_list (int x);

    list<int> pred_list (int x);

    void deleteLists();

    void printGraph();


    ~Graph();
};

class MaxHeap{
private:
    int _size{};
    vector<array<int, N_BO>> vect = {{-1,0,0,0,0}};;

    int p(int i) {return i >> 1;};
    int l(int i) {return i << 1;};
    int r(int i) {return (i << 1) + 1;};

public:

    bool isEmpty() const {return _size == 0;};

    array<int, N_BO> find_max() const {return vect.at(1);};

    void shiftUp (int i);

    void shiftDown (int i);

    void insert (array<int, N_BO> lab);

    void delete_max ();

    void increase_key (array<int, N_BO> lab );

    void print();

    void print_lable(array<int, N_BO> x);

    ~MaxHeap();

};




void node_labeling(vector<int>& group2, vector<int>& group3, Graph g);

bool NewCandLabel (Graph g, int i, array<int, N> l_star, int r, array<int, 5> &l_new, MaxHeap* Heap);

void RelaxProcess (Graph g, int i, array<int, N> l_star, MaxHeap* Heap);

vector<array<int, 5>> lable;

void print_signature(Graph g, int m1, int m2, int m3);

vector<int> sample_WO (vector<int> m, int sizeComb, int j);

void BO_maxCapPath (vector<int> perm2, vector<int> perm3, Graph g, MaxHeap* Heap);

void signature_estimation(int l1);


void graphGen_RGG(Graph & g, int V, int m1, int m2, int m3);

void power2000(Graph & g, int V, int m1, int m2, int m3);


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

    power2000(g, V, m1, m2, m3);


// Either assign the nodes to each class "manually" or following the piece of code below

    vector<int> group_1;
    vector<int> group_2;
    vector<int> group_3;


// Assigns the first m1 even nodes to the first class and the remaining m2 nodes to the second class

    int ct1 = 0, ct2 = 0, ct3 = 0;

    for (int i = 1; i <= m1+m2+m3; ++i)
    {
        if((i == 1 || i == 3 || i == 9) && (ct1 < m1)){
            group_1.push_back(i);
            ++ct1;}
        if(i  > 3 && i < 8 && ct2 < m2){
            group_2.push_back(i);
            ++ct2;}
        if((i == 2 || i == 8 || i == 10) && ct3 < m3){
            group_3.push_back(i);
            ++ct3;}
    }




    vector<int> perm1;
    vector<int> perm2;
    vector<int> perm3;




    s2 = clock();


            for (int i = 0; i < V; i++) {
            Active[i] = true;
        }      





    for (int j=1; j <= M; j++)
    {
        
        perm1.clear();
        perm2.clear();
        perm3.clear();


        vector<int> surviving_nodes;
        surviving_nodes.clear();



        perm1 = sample_WO(group_1, group_1.size(), j);

		
		//perm1.push_back(7);
		//perm1.push_back(4);
		//perm1.push_back(1);
		//perm1.push_back(10);


        perm2 = sample_WO(group_2, group_2.size(), j);

		//perm2.push_back(2);
		//perm2.push_back(5);
		//perm2.push_back(8);
		//perm2.push_back(11);
		

        perm3 = sample_WO(group_3, group_3.size(), j);

		//perm3.push_back(12);
		//perm3.push_back(6);
		//perm3.push_back(9);
		//perm3.push_back(3);

        surviving_nodes = perm1;






        
        
        for(int k=m1; k>=0; k--)
        {
            
            
            // Run BO_maxCapPath to find non-dominated pouints. Arguments are perm2 and perm 3.
            BO_maxCapPath (perm2, perm3, g, Heap);

            signature_estimation(k);

            

            int node = surviving_nodes.front(); // Get the first element in perm1
            Active[node] = false; // Set Active status to false for the node
            if (!surviving_nodes.empty()){
                surviving_nodes.erase(surviving_nodes.begin()); // Remove the first element from perm1
            }

        }



        for (int node : perm1)
        {
            Active[node] = true;
        }

    }
    


    


    s3 = clock();

    print_signature(g, m1, m2, m3);

    s4 = clock();

    // Reliability estimation present numerical issues for large networks due to the very large and very small
    // binomial coefficients, result in "NaN" values.

    // reliability_function (m1, m2, m3);

    s5 = clock();

    cout << endl;
    cout << endl <<"MO Algorithm - Run-time Performance:" <<endl;
    cout << "Network system: " << net << ", m1 = " << m1 << ", m2 = " << m2 << ", m3 = " << m3 << ", M = " << M <<endl<<endl;

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
    InGraph = new bool[V];
    Active = new bool[V];
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
    if (vect.at(i)[1] > vect.at(p(i))[1] || (vect.at(i)[1] == vect.at(p(i))[1] && vect.at(i)[2] > vect.at(p(i))[2])){
        std::swap(vect.at(p(i)), vect.at(i));
    }
    shiftUp(p(i));
}

void MaxHeap::insert(array<int, N_BO> lab){

    if (_size + 1 >= vect.size()){
        vect.push_back({0, 0, 0, 0, 0});
    }
    vect[++_size] = lab;
    shiftUp(_size);
    return;

}

void MaxHeap::shiftDown(int i) {
    
    if (i > _size) return;

    int swapId = i;

    if (l(i) <= _size && (vect.at(i)[1] < vect.at(l(i))[1] || (vect.at(i)[1] == vect.at(l(i))[1] && vect.at(i)[2] < vect.at(l(i))[2])) ) {
        swapId = l(i);
    }

    if (r(i) <= _size && (vect.at(swapId)[1] < vect.at(r(i))[1] || (vect.at(swapId)[1] == vect.at(r(i))[1] && vect.at(swapId)[2] < vect.at(r(i))[2])) ) {
        swapId = r(i);
    }

    if (swapId != i) {
        std::swap(vect.at(i), vect.at(swapId));
        shiftDown(swapId);
    }

    return;

}

void MaxHeap::delete_max(){
    array<int, N_BO> maxNum = vect.at(1);
    std::swap(vect.at(1), vect.at(_size--));
    shiftDown(1);
    vect.pop_back();
    return;
}

void MaxHeap::increase_key(array<int, N_BO> lab){
    int s;
    for (int i=0; i < vect.size(); ++i)
    {
        if (vect.at(i)[0] == lab[0])
        {
            vect.at(i)[1] = lab[1];
            vect.at(i)[2] = lab[2];
            vect.at(i)[3] = lab[3];
            vect.at(i)[4] = lab[4];
            s = i;
        }
    }

    shiftUp(s);
}




MaxHeap::~MaxHeap(){

}


bool NewCandLabel (Graph g, int i, array<int, N_BO> l_star, int r, array<int, N_BO> & l_new, MaxHeap* Heap)
{
    bool key = false;

    int f1, f2;
    array<int, N_BO> l;

    list<int> pred;

    int d1 = 0;
    int d2 = 0;

    pred = g.pred_list(i);

    for (list<int>::iterator j = pred.begin(); j != pred.end(); j++)
    {

        for (list<array<int, N_BO>>::iterator lj=ndpoints[*j].begin(); lj != ndpoints[*j].end(); lj++)
        {
                f1 = min(c[*j][0], c[i][0]);
                f2 = min(c[*j][1], c[i][1]);

                l=*lj;

                f1 = min(f1, l[1]);
                f2 = min(f2, l[2]);

                if(f1 > d1 || (f1==d1 && f2 > d2))
                {

                    if(f1 < l_star[1] && f2 > l_star[2])
                    {
                        d1 = f1;
                        d2 = f2;
                        l_new[0] = i;
                        l_new[1] = d1;
                        l_new[2] = d2;
                        l_new[3] = *j;
                        l_new[4] = Nindex[*j];
                        key = true;
                    }
                }
        }
    }

    return key;
}

void node_labeling (vector<int>& perm2, vector<int>& perm3, Graph g)
{

    vector<int>::iterator it;

    for(int a = 0; a < V; ++a)
    {
        if (a == 0)
        {
            c[0][0] = V;
            c[0][1] = V;
        }

        else if(find(perm2.begin(), perm2.end(), a) != perm2.end())
        {
            it = find(perm2.begin(), perm2.end(), a);
            c[a][0] = (it - perm2.begin()) + 1;
            c[a][1] = V;
        }

        else if (find(perm3.begin(), perm3.end(), a) != perm3.end())
        {
            it = find(perm3.begin(), perm3.end(), a);
            c[a][0] = V;
            c[a][1] = (it - perm3.begin()) + 1;
        }


        else
        {
            c[a][0] = V;
            c[a][1] = V;
        }


    }
}


void RelaxProcess (Graph g, int i, array<int, N_BO> l_star, MaxHeap* Heap)
{
    array<int, 2> f;
    //array<int, N> l;
    array<int, N_BO> l_new;
    array<int, N_BO> cl;

    list<int> succ;

    succ = g.succ_list(i);



    for(list<int>::iterator j = succ.begin(); j != succ.end(); j++)
    {
        if (Active[*j] == true){

            f[0] = min(c[*j][0], c[i][0]);
            f[1] = min(c[*j][1], c[i][1]);


            f[0] = min(f[0], l_star[1]);
            f[1] = min(f[1], l_star[2]);





            if (f[0] > d[*j][0] || (f[0] == d[*j][0] && f[1] > d[*j][1]))
            {
                //l = ndpoints[*j].back();

                bool dominated = false;

                for (list<array<int, N_BO>>::iterator iter=ndpoints[*j].begin(); iter != ndpoints[*j].end(); iter++){

                    cl = *iter;
                    int rp_counter_1 = 0;
                    int rp_counter_2 = 0;
                    int rp_counter_3 = 0;

                    for  (int k = 0; k < 3; k++){
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


                    l_new[0] = *j;
                    l_new[1] = d[*j][0];
                    l_new[2] = d[*j][1];
                    l_new[3] = i;
                    l_new[4] = Nindex[i];

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


void BO_maxCapPath(vector<int> perm2, vector<int> perm3, Graph g, MaxHeap* Heap)
{
    array<int, N_BO> l_star;
    array<int, 5> l_new ;

    for(int i=0; i<V; i++)
        ndpoints[i].clear();

    int i;



    fill(&d[0][0], &d[0][0] + V*2, 0);
    fill(&Nindex[0], &Nindex[0] + V, 0);
    fill(&InH[0], &InH[0] + V, false);

    node_labeling (perm2, perm3, g);

    d[0][0] = V;
    d[0][1] = V;

    array<int, N_BO> l = {{0, d[0][0], d[0][1], 0, 0}};

    Heap->insert(l);

    InH[0] = true;

    while (!Heap->isEmpty())
    {

        l_star = Heap->find_max();
        Heap->delete_max();

        i = l_star[0];

        d[i][0] = 0;
        d[i][1] = 0;


        Nindex[i] = Nindex[i] + 1;

        ndpoints[i].push_back(l_star);

        InH[i] = false;

        if(NewCandLabel(g, i, l_star, Nindex[i], l_new, Heap))
        {  

            Heap->insert(l_new);
            InH[i] = true;
            d[i][0] = l_new[1];
            d[i][1] = l_new[2];
        }
        
        RelaxProcess (g, i, l_star, Heap);
    }
}

void signature_estimation(int l1)
{
    array<int, N_BO> l;
    int l2, l3;


    


    fill(&survSig_filing[0][0][0], &survSig_filing[0][0][0] + (m1+1)*(m2+1)*(m3+1), false);

   for (list<array<int, N_BO>>::iterator lj=ndpoints[V-1].begin(); lj != ndpoints[V-1].end(); lj++)
   {
        l = *lj;
        if (l[1]<=m2)
           l2 = m2 - l[1];
        else
           l2 = -1;
        if (l[2]<= m3)
           l3 = m3 - l[2];
        else
           l3 = -1;
        for (int i = 0; i <= m2; i++)
        {
          for (int j = 0; j <= m3; j++)
          {
            if (i > l2 && j > l3)
                survSig_filing[l1][i][j] = true;
            
          }
        }
    }


  for(int i=0; i<=m2; i++)
  {
    for(int j=0; j<=m3; j++)
    {
        if (survSig_filing[l1][i][j] == true)
            survSig[l1][i][j] = survSig[l1][i][j] + 1;
    }
   }
}



void print_signature(Graph g, int m1, int m2, int m3)
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
                survSig[l1][l2][l3] = survSig[l1][l2][l3]/(double) M;
                cout << setprecision(6) << survSig[l1][l2][l3] << " ";
                myfile << setprecision(6) << survSig[l1][l2][l3] << " ";
            }
                cout << endl;
                myfile << endl;
        }
        cout << "\n";
    }
    cout <<endl<<endl;
    myfile.close();
}





void power2000(Graph & g, int V, int m1, int m2, int m3)
{   
    // Adding edges
g.addEdge(0, 1);
g.addEdge(0, 4);
g.addEdge(1, 2);
g.addEdge(1, 5);
g.addEdge(2, 3);
g.addEdge(2, 6);
g.addEdge(3, 7);
g.addEdge(4, 5);
g.addEdge(5, 6);
g.addEdge(6, 7);
g.addEdge(4, 8);
g.addEdge(5, 9);
g.addEdge(6, 10);
g.addEdge(7, 11);
g.addEdge(8, 9);
g.addEdge(9, 11); 
g.addEdge(9, 3);
g.addEdge(9, 10);
g.addEdge(10, 11);


}
