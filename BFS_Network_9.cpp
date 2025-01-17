#include <random>
#include <random>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <list>
#include <cmath>
#include <bits/stdc++.h>
#include <fstream>
#include <string>


using namespace std;

// ---- PARAMETERS -------------------------------------------------------------
const int V = 81;              // Total number of nodes
const int v = V-2;              // Nodes that can fail (source and sink don't fail)
const int m1 = 30;              // Number of nodes in the first group/class of nodes
const int m2 = 30;            // Number of nodes in the second group/class of nodes
const int m3 = v-(m1+m2);      // Number of nodes in the third group/class of nodes
const int M = 10;             // Number of MC replications

// -----------------------------------------------------------------------------


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


clock_t s1, s2, s3, s4, s5;

double survSig[m1+1][m2+1][m3+1]  = {};


// --------  CLASSES AND FUNCTIONS --------------------------------------------------------

class Graph
{
public:
	int V;

	list<int> *adj;
	list<int> *pred;

	Graph(int V);

	void addEdge(int v, int w);

	bool BFS(int s, vector<int> inactive);

	list<int> succ_list (int x);

    list<int> pred_list (int x);

	//void genGrid(vector<int> grid[], vector <int> &group_1, vector <int> & group_2, vector <int> & group_3, int m1, int m2, int m3, int v);

	void deleteAdj();

    ~Graph();
};

vector<int> sample_WO (vector<int> m, int sizeComb, int j);



void signature_estimation(vector<int> group_1, vector<int> group_2, vector<int> group_3, Graph g, int m1, int m2, int m3, int v, int V);

void print_signature(Graph g, int m1, int m2, int m3);



//void graphGen_RGG(Graph & g, int V, int m1, int m2, int m3);

void power2000(Graph & g, int V, int m1, int m2, int m3);


// ***********************************************************************************


// ------------------ MAIN ----------------------------------------------------------
int main()
{

    s1 = clock();

    cout << "Time at the start : " << double(s1) / double(CLOCKS_PER_SEC) << setprecision(5) << endl;

    string net = "RGG 350";
    Graph g(V);

// Either construct a graph by adding edges (in case of a small network) or
// run a function to add edges.

// function graphGen_RGG generates a RGG with m1 nodes in the first class and m2 nodes in the second class;

    //graphGen_RGG(g, V, m1, m2, m3);


// function power2000 generates the Power 2000 bus electric grid network. Make sure V=4000 since in this function
// all edges are hard-coded.

    power2000(g, V, m1, m2, m3);



    std::vector<int> class_1;
    std::vector<int> class_2;
    std::vector<int> class_3;

    // Fill class_1 with numbers from 0 to 29
    for (int i = 0; i <= 29; ++i) {
        class_1.push_back(i);
    }

    // Fill class_2 with numbers from 30 to 59
    for (int i = 30; i <= 59; ++i) {
        class_2.push_back(i);
    }

    // Fill class_3 with numbers from 60 to 80
    for (int i = 60; i <= 80; ++i) {
        class_3.push_back(i);
    }






// Either assign the nodes to each class "manually" or following the piece of code below

    vector<int> group_1;
    vector<int> group_2;
    vector<int> group_3;


// Assigns the first m1 even nodes to the first class and the remaining m2 nodes to the second class


    int ct1 = 0, ct2 = 0, ct3 = 0;

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

    cout << "Group 1: ";
    for (int node : group_1) {
        cout << node << " ";
    }
    cout << endl;

    cout << "Group 2: ";
    for (int node : group_2) {
        cout << node << " ";
    }
    cout << endl;

    cout << "Group 3: ";
    for (int node : group_3) {
        cout << node << " ";
    }
    cout << endl;

    s2 = clock();

    cout << "Time at the initialization : " << double(s2) / double(CLOCKS_PER_SEC) << setprecision(5) << endl;

    signature_estimation(group_1, group_2, group_3, g, m1, m2, m3, v, V);

    s3 = clock();

    print_signature(g, m1, m2, m3);

    s4 = clock();

    // Reliability estimation present numerical issues for large networks due to the very large and very small
    // binomial coefficients, result in "NaN" values.

    //  reliability_function (m1, m2);

    s5 = clock();

    cout << endl;
    cout << endl <<"Naive Algorithm - Run-time Performance:" <<endl;
    cout << "Network system: " << net << ", m1 = " << m1 << ", m2 = " << m2 <<", m3 = " << m3 <<", M = " << M <<endl<<endl;

    cout << "Initialization: " << double(s2 - s1) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;
    cout << "Survival Signature estimation: " << double(s3 - s2) / double(CLOCKS_PER_SEC) << setprecision(5) << endl <<endl;
    cout << "Survival Signature rate: " << (double(s3 - s2) / double(CLOCKS_PER_SEC))/(double(M)) << setprecision(5) << endl <<endl;
    cout << "Printing Signature: " << (double(s4 - s3) / double(CLOCKS_PER_SEC)) << setprecision(5) << endl <<endl;
//    cout << "Reliability estimation: " << double(s4 - s3) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;
//    cout << "Total Time: " << double(s4 - s1) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;


    g.deleteAdj();
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


bool Graph::BFS(int s, vector<int> inactive)
{

	int nodes = 0;
    bool key = false;

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

        if (s == V-1)
		{
            key = true;

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
	return key;
}

list<int> Graph::succ_list (int x)
{
    return adj[x];
}

list<int> Graph::pred_list (int x)
{
    return pred[x];
}


void Graph:: deleteAdj ()
{
    delete [] adj;
    delete [] pred;
}


Graph::~Graph()
{

}



// ~~~~~~~~~~~~~~~ Sampling Without Replacement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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



// ~~~~~~~~~~~~~~ survival signature ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void signature_estimation(vector<int> group_1, vector<int> group_2, vector<int> group_3, Graph g, int m1, int m2, int m3, int v, int V)
{
    vector<int> perm1;
    vector<int> perm2;
    vector<int> perm3;

    vector<int> failed_nodes;
    vector<int> failed_nodes_2;
    vector<int> sampling1;
    vector<int> sampling2;
    vector<int> sampling3;


    for (int j=1; j <= M; j++)
    {
        perm1.clear();
        perm2.clear();
        perm3.clear();

        perm1 = sample_WO(group_1, group_1.size(), j);
        perm2 = sample_WO(group_2, group_2.size(), j);
        perm3 = sample_WO(group_3, group_3.size(), j);

        failed_nodes=perm1;
        failed_nodes_2=perm2;
        failed_nodes.insert(failed_nodes.end(), perm2.begin(), perm2.end());
        failed_nodes.insert(failed_nodes.end(), perm3.begin(), perm3.end());

        for(int l1=0; l1<=m1; l1++)
        {
            for(int l2=0; l2<=m2; l2++)
            {
                for(int l3=0; l3<=m3; l3++)
                {
                    if(g.BFS(0, failed_nodes))
                    {
                        survSig[l1][l2][l3]=survSig[l1][l2][l3] + 1;
                    }

                    if(!failed_nodes.empty())
                        failed_nodes.pop_back();
                }

                
                if(!failed_nodes_2.empty())
                    failed_nodes_2.pop_back();

                failed_nodes.clear();
                failed_nodes = perm1;
                failed_nodes.insert(failed_nodes.end(), failed_nodes_2.begin(), failed_nodes_2.end());
                failed_nodes.insert(failed_nodes.end(), perm3.begin(), perm3.end());
                

            }

            if(!perm1.empty())
                perm1.pop_back();

            failed_nodes.clear();
            failed_nodes = perm1;
            failed_nodes.insert(failed_nodes.end(), perm2.begin(), perm2.end());
            failed_nodes.insert(failed_nodes.end(), perm3.begin(), perm3.end());

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
    // Create the network with exactly 4 neighbors for each node
g.addEdge(24, 29);
g.addEdge(16, 24);
g.addEdge(16, 39);
g.addEdge(36, 66);
g.addEdge(36, 74);
g.addEdge(34, 45);
g.addEdge(46, 19);
g.addEdge(46, 68);
g.addEdge(46, 6);
g.addEdge(10, 19);
g.addEdge(10, 55);
g.addEdge(31, 19);
g.addEdge(68, 8);
g.addEdge(68, 11);
g.addEdge(68, 37);
g.addEdge(71, 76);
g.addEdge(72, 71);
g.addEdge(76, 8);
g.addEdge(76, 25);
g.addEdge(76, 33);
g.addEdge(67, 71);
g.addEdge(67, 30);
g.addEdge(67, 26);
g.addEdge(48, 63);
g.addEdge(48, 5);
g.addEdge(58, 70);
g.addEdge(3, 80);
g.addEdge(63, 17);
g.addEdge(74, 63);
g.addEdge(41, 51);
g.addEdge(50, 72);
g.addEdge(30, 64);
g.addEdge(73, 30);
g.addEdge(37, 64);
g.addEdge(12, 64);
g.addEdge(20, 29);
g.addEdge(20, 34);
g.addEdge(20, 14);
g.addEdge(20, 54);
g.addEdge(14, 1);
g.addEdge(14, 66);
g.addEdge(54, 18);
g.addEdge(15, 29);
g.addEdge(15, 32);
g.addEdge(44, 38);
g.addEdge(44, 69);
g.addEdge(44, 34);
g.addEdge(38, 60);
g.addEdge(40, 38);
g.addEdge(40, 36);
g.addEdge(40, 28);
g.addEdge(79, 38);
g.addEdge(60, 78);
g.addEdge(60, 33);
g.addEdge(0, 44);
g.addEdge(0, 16);
g.addEdge(0, 59);
g.addEdge(0, 70);
g.addEdge(62, 73);
g.addEdge(62, 15);
g.addEdge(62, 27);
g.addEdge(25, 73);
g.addEdge(25, 77);
g.addEdge(25, 31);
g.addEdge(2, 73);
g.addEdge(2, 43);
g.addEdge(56, 47);
g.addEdge(47, 7);
g.addEdge(77, 56);
g.addEdge(33, 9);
g.addEdge(22, 43);
g.addEdge(22, 5);
g.addEdge(43, 16);
g.addEdge(5, 26);
g.addEdge(49, 22);
g.addEdge(49, 23);
g.addEdge(75, 65);
g.addEdge(65, 50);
g.addEdge(39, 32);
g.addEdge(39, 13);
g.addEdge(32, 66);
g.addEdge(32, 42);
g.addEdge(42, 10);
g.addEdge(35, 13);
g.addEdge(35, 37);
g.addEdge(35, 10);
g.addEdge(35, 5);
g.addEdge(51, 79);
g.addEdge(51, 69);
g.addEdge(27, 15);
g.addEdge(27, 21);
g.addEdge(26, 47);
g.addEdge(80, 45);
g.addEdge(80, 74);
g.addEdge(59, 54);
g.addEdge(70, 3);
g.addEdge(70, 9);
g.addEdge(21, 40);
g.addEdge(8, 1);
g.addEdge(8, 42);
g.addEdge(1, 57);
g.addEdge(1, 6);
g.addEdge(78, 58);
g.addEdge(78, 31);
g.addEdge(4, 78);
g.addEdge(4, 58);
g.addEdge(4, 28);
g.addEdge(45, 59);
g.addEdge(18, 55);
g.addEdge(18, 14);
g.addEdge(28, 79);
g.addEdge(55, 12);
g.addEdge(16, 30);
g.addEdge(30, 55);
g.addEdge(55, 12);
g.addEdge(12, 28);
g.addEdge(28, 80);

}