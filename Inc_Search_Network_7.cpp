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

// ---- PARAMETERS -------------------------------------------------------------

const int V = 25;                      // Total number of nodes
const int v = V-2;                     // Nodes that can fail (source and sink don't fail)
const int m1 = 11;
const int m2 = 9;                      // Number of nodes in the first group/class of nodes
const int m3 = v-(m1+m2);              // Number of nodes in the second group/class of nodes
const int M = 100;                      // Number of MC replications

// -----------------------------------------------------------------------------


// ~~~~~~~~~~  class Graph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




void power2000(Graph & g, int V, int m1, int m2, int m3);

double survSig[m1+1][m2+1][m3+1]  = {};


vector<int> sample_WO (vector<int> m, int sizeComb, int j);



void print_signature(Graph g, int m1, int m2, int m3);



bool initial_BFS (int s, vector<bool> & active_nodes, vector<bool> &marked, Graph g, list<int> &search_list);

bool update_BFS(int k, vector<bool> & active_nodes, vector<bool>& marked, Graph g, list<int> &search_list);

// ***********************************************************************************

clock_t s1, s2, s3, s4, s5;

// ------------------ MAIN ----------------------------------------------------------
int main()
{

    s1 = clock();

    string net = "RGG 350";
    Graph g(V);

// Either construct a graph by adding edges (in case of a small network) or
// run a function to add edges.

// function graphGen_RGG generates a RGG with m1 nodes in the first class and m2 nodes in the second class;



// function power2000 generates the Power 2000 bus electric grid network. Make sure V=4000 since in this function
// all edges are hard-coded.

    power2000(g, V, m1, m2, m3);


// Either assign the nodes to each class "manually" or following the piece of code below

    std::vector<int> class_1 = {1, 2, 3, 4, 10, 11, 12, 21, 19, 23, 17};
    std::vector<int> class_2 = {5, 6, 7, 13, 20, 22, 18, 15, 16};
    std::vector<int> class_3 = {8, 9, 14};





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


    vector<int> perm1;
    vector<int> perm2;
	vector<int> perm3;

    vector<bool> active_nodes(V,false);
    vector<bool> marked(V,false);
    list<int> search_list;


    int pos1;
    int pos2;
	int pos3;

    int node;

    s2 = clock();

    for (int j=1; j <= M; j++)
    {
        perm1.clear();
        perm2.clear();
		perm3.clear();

        perm1 = sample_WO(group_1, group_1.size(), j);
        perm2 = sample_WO(group_2, group_2.size(), j);
		perm3 = sample_WO(group_3, group_3.size(), j);

        search_list.clear();

        for(int i=0; i < active_nodes.size(); i++)
        {
            if (i == 0 || i == V-1)
                active_nodes[i] =  true;
            else
                active_nodes[i] =  false;
        }

        pos1 = perm1.size()-1;
        for(int l1=0; l1<=m1; l1++)
        {
            pos2 = perm2.size()-1;

            for(int l2=0; l2<=m2; l2++)
            {
				
				pos3 = perm3.size()-1;

				for(int l3=0; l3<=m3; l3++)
				{
                	if(l3 == 0)
                	{
                    	if(initial_BFS(0,active_nodes, marked, g, search_list))
                        	survSig[l1][l2][l3]=survSig[l1][l2][l3] + 1;
                	}
                	else
                	{
                        
                    	if(update_BFS(node, active_nodes, marked, g, search_list))
                        {   	
                            survSig[l1][l2][l3]=survSig[l1][l2][l3] + 1;
                            
                        }
                	}


					if(pos3 >= 0)
                	{
                    	node = perm3[pos3];
                    	active_nodes[node] = true;
                    	pos3 = pos3-1;
                	}
                	else
                	{
                    	for (int k : perm3)
                    	{
                        	active_nodes[k] = false;
                    	}
                	}

				

				}

            	if(pos2 >= 0)
				{
                	node = perm2[pos2];
                	active_nodes[node] = true;
                	pos2 = pos2-1;
            	}
                else
                {
                    for (int z: perm2)
                    {
                        active_nodes[z] = false;
                    }
                }

            }

            if(pos1 >= 0){
                node = perm1[pos1];
                active_nodes[node] = true;
                pos1 = pos1-1;
            }
        }
    }


    s3 = clock();

    print_signature(g, m1, m2, m3);

    s4 = clock();

    // Reliability estimation present numerical issues for large networks due to the very large and very small
    // binomial coefficients, result in "NaN" values.

    // reliability_function (m1, m2);

    s5 = clock();

    cout << endl;
    cout << endl <<"IncSearch Algorithm - Run-time Performance:" <<endl;
    cout << "Network system: " << net << ", m1 = " << m1 << ", m2 = " << m2 <<", m3 = " << m3 <<", M = " << M <<endl<<endl;

    cout << "Initialization: " << double(s2 - s1) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;
    cout << "Survival Signature estimation: " << double(s3 - s2) / double(CLOCKS_PER_SEC) << setprecision(5) << endl <<endl;
    cout << "Survival Signature rate: " << (double(s3 - s2) / double(CLOCKS_PER_SEC))/(double(M)) << setprecision(5) << endl <<endl;
    cout << "Printing Signature: " << (double(s4 - s3) / double(CLOCKS_PER_SEC)) << setprecision(5) << endl <<endl;
//    cout << "Reliability estimation: " << double(s4 - s3) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;
//    cout << "Total Time: " << double(s4 - s1) / double(CLOCKS_PER_SEC) << setprecision(5) << endl << endl;

g.deleteLists();
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

void Graph:: deleteLists ()
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


// ~~~~~~~~~~~~~~~~~~~~~~~ BFS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool initial_BFS(int s, vector<bool> & active_nodes, vector<bool> & marked, Graph g, list<int> & search_list)
{

    int j;

   	for(int i=0; i<marked.size();i++)
    {
        marked[i] = false;
    }

   	marked[s] = true;
   	search_list.push_back(s);

	while(!search_list.empty())
	{
        j = search_list.front();
        search_list.pop_front();


		for (list<int>::iterator i = g.adj[j].begin(); i != g.adj[j].end(); i++)
		{
			if (!marked[*i] && (active_nodes[*i]))
			{
                marked[*i] = true;
                
                search_list.push_back(*i);
			}
		}
	}


    
	return marked[V-1];

}


bool update_BFS(int k, vector<bool> & active_nodes, vector<bool> &marked, Graph g, list<int> & search_list)
{

    int j;

    for(list<int>::iterator i = g.pred[k].begin();i != g.pred[k].end(); ++i){
        if(marked[*i] == true)
        {
            marked[k] = true;
            search_list.push_back(k);
        }
    }

	while(!search_list.empty())
	{

        j = search_list.front();
        search_list.pop_front();


		for (list<int>::iterator i = g.adj[j].begin(); i != g.adj[j].end(); ++i)
		{
			if (active_nodes[*i] && !marked[*i])
			{
                marked[*i] = true;
				search_list.push_back(*i);
			}
		}
	}

	return marked[V-1];

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
g.addEdge(6, 15);
g.addEdge(6, 2);
g.addEdge(15, 24);
g.addEdge(15, 14);
g.addEdge(15, 13);
g.addEdge(5, 6);
g.addEdge(16, 6);
g.addEdge(16, 4);
g.addEdge(2, 12);
g.addEdge(2, 5);
g.addEdge(14, 1);
g.addEdge(14, 21);
g.addEdge(13, 23);
g.addEdge(8, 24);
g.addEdge(8, 11);
g.addEdge(8, 5);
g.addEdge(7, 24);
g.addEdge(7, 2);
g.addEdge(18, 24);
g.addEdge(18, 17);
g.addEdge(18, 16);
g.addEdge(3, 7);
g.addEdge(3, 22);
g.addEdge(3, 1);
g.addEdge(3, 20);
g.addEdge(22, 7);
g.addEdge(22, 5);
g.addEdge(1, 23);
g.addEdge(1, 13);
g.addEdge(20, 12);
g.addEdge(20, 17);
g.addEdge(12, 4);
g.addEdge(12, 9);
g.addEdge(11, 4);
g.addEdge(17, 4);
g.addEdge(17, 19);
g.addEdge(9, 11);
g.addEdge(9, 21);
g.addEdge(10, 22);
g.addEdge(10, 11);
g.addEdge(10, 23);
g.addEdge(0, 9);
g.addEdge(0, 20);
g.addEdge(0, 16);
g.addEdge(0, 13);
g.addEdge(21, 19);
g.addEdge(21, 10);
g.addEdge(19, 18);
g.addEdge(19, 8);
g.addEdge(23, 14);


}