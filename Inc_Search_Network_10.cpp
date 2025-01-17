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

const int V = 100;              // Total number of nodes
const int v = V-2;              // Nodes that can fail (source and sink don't fail)
const int m1 = 20;              // Number of nodes in the first group/class of nodes
const int m2 = 20;
const int m3 = 20;
const int m4 = 20;              // Number of nodes in the second group/class of nodes
const int m5 = v-(m1+m2+m3+m4);      // Number of nodes in the third group/class of nodes
const int M = 10;             // Number of MC replications

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




void power2000(Graph & g, int V, int m1, int m2, int m3, int m4, int m5);

double survSig[m1+1][m2+1][m3+1][m4+1][m5+1]  = {};


vector<int> sample_WO (vector<int> m, int sizeComb, int j);



void print_signature(Graph g, int m1, int m2, int m3, int m4, int m5);



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

    power2000(g, V, m1, m2, m3, m4, m5);



    std::vector<int> class_1 = {3, 18, 11, 15, 12, 9, 13, 20, 2, 16, 6, 7, 10, 19, 4, 17, 14, 1, 8, 5};
    std::vector<int> class_2 = {21, 31, 30, 39, 37, 34, 40, 24, 35, 27, 22, 32, 23, 26, 28, 25, 36, 38, 29, 33};
    std::vector<int> class_3 = {57, 53, 56, 49, 60, 50, 52, 43, 44, 41, 45, 58, 51, 42, 46, 59, 55, 47, 54, 48};
    std::vector<int> class_4 = {69, 70, 63, 78, 80, 79, 65, 77, 67, 76, 71, 61, 66, 68, 62, 73, 64, 74, 72, 75};
    std::vector<int> class_5 = {88, 82, 87, 81, 85, 90, 86, 95, 83, 98, 96, 92, 91, 94, 84, 93, 89, 97}; 




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

    vector<bool> active_nodes(V,false);
    vector<bool> marked(V,false);
    list<int> search_list;


    int pos1;
    int pos2;
	int pos3;
    int pos4;
    int pos5;

    int node;

    s2 = clock();

    for (int j=1; j <= M; j++)
    {
        perm1.clear();
        perm2.clear();
		perm3.clear();
        perm4.clear();
        perm5.clear();

        perm1 = sample_WO(group_1, group_1.size(), j);
        perm2 = sample_WO(group_2, group_2.size(), j);
		perm3 = sample_WO(group_3, group_3.size(), j);
        perm4 = sample_WO(group_4, group_4.size(), j);
        perm5 = sample_WO(group_5, group_5.size(), j);

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
                    pos4 = perm4.size()-1;

                    for(int l4=0; l4<=m4; l4++)
                    {

                        pos5 = perm5.size()-1;

                        for(int l5=0; l5<=m5; l5++)
                        {
                            if(l5 == 0)
                            {
                                if(initial_BFS(0,active_nodes, marked, g, search_list))
                                    survSig[l1][l2][l3][l4][l5]=survSig[l1][l2][l3][l4][l5] + 1;
                            }
                            else
                            {
                                if(update_BFS(node, active_nodes, marked, g, search_list))
                                {   	
                                    survSig[l1][l2][l3][l4][l5]=survSig[l1][l2][l3][l4][l5] + 1;
                                }
                            }
                            if(pos5 >= 0)
                            {
                                node = perm5[pos5];
                                active_nodes[node] = true;
                                pos5 = pos5-1;
                            }
                            else
                            {
                                for (int k : perm3)
                                {
                                    active_nodes[k] = false;
                                }
                            }



                        }


                        if(pos4 >= 0)
                        {
                            node = perm4[pos4];
                            active_nodes[node] = true;
                            pos4 = pos4-1;
                        }
                        else
                        {
                            for (int k : perm3)
                            {
                                active_nodes[k] = false;
                            }
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

    print_signature(g, m1, m2, m3, m4, m5);

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

            }
                cout << endl;
                myfile << endl;
        }
        cout << "\n";
    }
    cout <<endl<<endl;
    myfile.close();
}





void power2000(Graph & g, int V, int m1, int m2, int m3, int m4, int m5)
{
g.addEdge(73, 72);
g.addEdge(72, 83);
g.addEdge(72, 28);
g.addEdge(12, 72);
g.addEdge(12, 60);
g.addEdge(28, 44);
g.addEdge(28, 45);
g.addEdge(81, 73);
g.addEdge(81, 33);
g.addEdge(81, 66);
g.addEdge(90, 73);
g.addEdge(90, 87);
g.addEdge(35, 73);
g.addEdge(39, 81);
g.addEdge(39, 92);
g.addEdge(39, 54);
g.addEdge(33, 22);
g.addEdge(66, 64);
g.addEdge(24, 42);
g.addEdge(24, 93);
g.addEdge(24, 38);
g.addEdge(42, 9);
g.addEdge(63, 24);
g.addEdge(63, 47);
g.addEdge(63, 86);
g.addEdge(63, 78);
g.addEdge(93, 70);
g.addEdge(93, 52);
g.addEdge(38, 5);
g.addEdge(38, 17);
g.addEdge(38, 41);
g.addEdge(82, 42);
g.addEdge(82, 53);
g.addEdge(82, 94);
g.addEdge(82, 2);
g.addEdge(9, 70);
g.addEdge(9, 99);
g.addEdge(9, 10);
g.addEdge(15, 42);
g.addEdge(15, 50);
g.addEdge(79, 97);
g.addEdge(79, 76);
g.addEdge(79, 87);
g.addEdge(99, 65);
g.addEdge(99, 77);
g.addEdge(76, 48);
g.addEdge(76, 94);
g.addEdge(76, 74);
g.addEdge(87, 89);
g.addEdge(87, 18);
g.addEdge(7, 79);
g.addEdge(7, 54);
g.addEdge(69, 97);
g.addEdge(69, 68);
g.addEdge(65, 89);
g.addEdge(65, 41);
g.addEdge(65, 33);
g.addEdge(77, 86);
g.addEdge(77, 48);
g.addEdge(26, 7);
g.addEdge(26, 66);
g.addEdge(26, 61);
g.addEdge(26, 25);
g.addEdge(41, 7);
g.addEdge(61, 70);
g.addEdge(61, 22);
g.addEdge(25, 85);
g.addEdge(25, 15);
g.addEdge(19, 55);
g.addEdge(19, 50);
g.addEdge(19, 71);
g.addEdge(55, 85);
g.addEdge(55, 27);
g.addEdge(50, 28);
g.addEdge(71, 83);
g.addEdge(89, 19);
g.addEdge(47, 55);
g.addEdge(47, 99);
g.addEdge(47, 27);
g.addEdge(27, 85);
g.addEdge(62, 90);
g.addEdge(62, 93);
g.addEdge(18, 62);
g.addEdge(18, 45);
g.addEdge(18, 59);
g.addEdge(6, 62);
g.addEdge(6, 52);
g.addEdge(6, 83);
g.addEdge(31, 90);
g.addEdge(31, 15);
g.addEdge(31, 50);
g.addEdge(70, 68);
g.addEdge(22, 96);
g.addEdge(22, 49);
g.addEdge(68, 61);
g.addEdge(68, 35);
g.addEdge(0, 51);
g.addEdge(0, 37);
g.addEdge(0, 57);
g.addEdge(0, 95);
g.addEdge(51, 11);
g.addEdge(51, 64);
g.addEdge(51, 49);
g.addEdge(37, 17);
g.addEdge(37, 44);
g.addEdge(37, 84);
g.addEdge(57, 4);
g.addEdge(57, 43);
g.addEdge(95, 34);
g.addEdge(11, 64);
g.addEdge(11, 45);
g.addEdge(11, 49);
g.addEdge(49, 13);
g.addEdge(91, 89);
g.addEdge(91, 21);
g.addEdge(91, 75);
g.addEdge(29, 20);
g.addEdge(29, 75);
g.addEdge(20, 52);
g.addEdge(20, 36);
g.addEdge(58, 20);
g.addEdge(58, 95);
g.addEdge(58, 10);
g.addEdge(58, 6);
g.addEdge(52, 98);
g.addEdge(36, 3);
g.addEdge(84, 29);
g.addEdge(84, 60);
g.addEdge(40, 29);
g.addEdge(40, 30);
g.addEdge(40, 99);
g.addEdge(40, 71);
g.addEdge(75, 53);
g.addEdge(75, 98);
g.addEdge(34, 92);
g.addEdge(34, 23);
g.addEdge(34, 31);
g.addEdge(92, 27);
g.addEdge(92, 83);
g.addEdge(23, 10);
g.addEdge(10, 2);
g.addEdge(53, 23);
g.addEdge(53, 74);
g.addEdge(8, 23);
g.addEdge(8, 77);
g.addEdge(30, 4);
g.addEdge(30, 8);
g.addEdge(30, 71);
g.addEdge(4, 74);
g.addEdge(13, 4);
g.addEdge(13, 35);
g.addEdge(74, 44);
g.addEdge(46, 21);
g.addEdge(46, 2);
g.addEdge(43, 21);
g.addEdge(43, 91);
g.addEdge(78, 21);
g.addEdge(78, 67);
g.addEdge(67, 46);
g.addEdge(67, 14);
g.addEdge(80, 46);
g.addEdge(80, 57);
g.addEdge(80, 59);
g.addEdge(80, 48);
g.addEdge(48, 96);
g.addEdge(94, 25);
g.addEdge(94, 2);
g.addEdge(16, 95);
g.addEdge(16, 5);
g.addEdge(16, 66);
g.addEdge(16, 12);
g.addEdge(32, 85);
g.addEdge(32, 8);
g.addEdge(32, 59);
g.addEdge(32, 88);
g.addEdge(59, 13);
g.addEdge(88, 54);
g.addEdge(88, 3);
g.addEdge(88, 60);
g.addEdge(14, 99);
g.addEdge(14, 69);
g.addEdge(14, 78);
g.addEdge(17, 35);
g.addEdge(17, 69);
g.addEdge(44, 39);
g.addEdge(5, 67);
g.addEdge(5, 56);
g.addEdge(56, 54);
g.addEdge(56, 1);
g.addEdge(56, 36);
g.addEdge(3, 41);
g.addEdge(60, 33);
g.addEdge(1, 64);
g.addEdge(1, 12);
g.addEdge(1, 36);
g.addEdge(98, 3);
g.addEdge(98, 84);
g.addEdge(96, 45);
g.addEdge(96, 86);
g.addEdge(86, 43);
g.addEdge(51, 86);
g.addEdge(57, 43);
g.addEdge(43, 65);
g.addEdge(65, 86);
g.addEdge(86, 99);


}