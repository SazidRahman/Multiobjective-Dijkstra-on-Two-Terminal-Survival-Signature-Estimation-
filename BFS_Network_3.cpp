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
const int V = 54;              // Total number of nodes
const int v = V-2;              // Nodes that can fail (source and sink don't fail)
const int m1 = 24;              // Number of nodes in the first group/class of nodes
const int m2 = 16;            // Number of nodes in the second group/class of nodes
const int m3 = v-(m1+m2);      // Number of nodes in the third group/class of nodes
const int M = 1000;             // Number of MC replications

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

double binomial(int n, int k);

void signature_estimation(vector<int> group_1, vector<int> group_2, vector<int> group_3, Graph g, int m1, int m2, int m3, int v, int V);

void print_signature(Graph g, int m1, int m2, int m3);

double Weibull(double t, double k, double lambda);

vector<double> reliability_function (int m1, int m2, int m3);

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



    std::vector<int> class_1 = {1, 2, 3, 4, 10, 11, 12, 13, 19, 20, 21, 22, 28, 29, 30, 31, 37, 38, 39, 40, 46, 47, 48, 49};
    std::vector<int> class_2 = {7, 8, 9, 16, 17, 18, 25, 26, 27, 34, 35, 36, 43, 44, 45, 52};
    std::vector<int> class_3 = {5, 6, 14, 15, 23, 24, 32, 33, 41, 42, 50, 51};





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


//~~~~~~~~~~~~~ Binomial Coefficient ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double binomial(int n, int k)
{
	if (n < k)
		{
			cout << "Error: n < k !";
			return 0;
		}

	double res = 1;

	for (int i=1; i<=k; ++i)
	{
		res = res * ((n + 1 - i)/ (double) i);
	}

	return res;
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


// -------------- Mission Times Function ---------------
double Weibull (double t, double k, double lambda)
{
	if (t >= 0)
		return 1 - exp(-pow(t/lambda, k));
	else
		return 0;
}


// ---------------- survival function --------------------

vector<double> reliability_function (int m1, int m2, int m3)
{
	vector<double> t = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200};
	vector<double> R;

    ofstream myfile;
	myfile.open("Reliability.txt");

	myfile << "Reliability estimation" << endl << endl;
	myfile << "t <-c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200)" << endl << endl;
	myfile << "R <-c(";

	double F1;
	double F1_c;
	double F2;
	double F2_c;
    double F3;
	double F3_c;
	double prod;
	double T;

	for (int v = 0; v < t.size(); ++v)
	{
		F1 = Weibull(t.at(v), 2.5, 100);
		F1_c = 1-F1;
		F2 = Weibull(t.at(v), 3.0, 100);
		F2_c = 1-F2;
        F3 = Weibull(t.at(v), 2.0, 100);
		F3_c = 1-F3;

		prod = 0.0;
		T = 0.0;

		for(int l1=0; l1<=m1; l1++)
		{
			for(int l2=0; l2<=m2; l2++)
			{
                for(int l3=0; l3<=m3; l3++)
                {
                    if (survSig[l1][l2][l3] != 0)
                    {
                        prod = survSig[l1][l2][l3] * binomial(m1, l1) * pow(F1, (m1-l1)) * pow(F1_c, l1) * binomial(m2, l2) * pow(F2, (m2-l2)) * pow(F2_c, l2)*binomial(m3, l3) * pow(F3, (m3-l3)) * pow(F3_c, l3);
                        T = T + prod;
                    }
                }
                
			}
            
		}
		cout <<endl << "P(T > " << t.at(v) <<") = " << T << endl;
		myfile << T <<", ";
		R.push_back(T);

	}

    myfile << ")" <<endl;
	myfile.close();
    return R;


}


void power2000(Graph & g, int V, int m1, int m2, int m3)
{
    // Create the network with exactly 4 neighbors for each node
g.addEdge(0, 1);
g.addEdge(0, 9);
g.addEdge(1, 2);
g.addEdge(1, 10);
g.addEdge(2, 3);
g.addEdge(2, 11);
g.addEdge(3, 4);
g.addEdge(3, 12);
g.addEdge(4, 5);
g.addEdge(4, 13);
g.addEdge(5, 6);
g.addEdge(5, 14);
g.addEdge(6, 7);
g.addEdge(6, 15);
g.addEdge(7, 8);
g.addEdge(7, 16);
g.addEdge(9, 10);
g.addEdge(9, 18);
g.addEdge(10, 11);
g.addEdge(10, 19);
g.addEdge(11, 12);
g.addEdge(11, 20);
g.addEdge(12, 13);
g.addEdge(12, 21);
g.addEdge(13, 14);
g.addEdge(13, 22);
g.addEdge(14, 15);
g.addEdge(14, 23);
g.addEdge(15, 16);
g.addEdge(15, 24);
g.addEdge(16, 25);
g.addEdge(18, 19);
g.addEdge(18, 27);
g.addEdge(19, 20);
g.addEdge(19, 28);
g.addEdge(20, 21);
g.addEdge(20, 29);
g.addEdge(21, 22);
g.addEdge(21, 30);
g.addEdge(22, 23);
g.addEdge(22, 31);
g.addEdge(23, 24);
g.addEdge(23, 32);
g.addEdge(24, 25);
g.addEdge(24, 33);
g.addEdge(25, 34);
g.addEdge(27, 28);
g.addEdge(27, 36);
g.addEdge(28, 29);
g.addEdge(28, 37);
g.addEdge(29, 30);
g.addEdge(29, 38);
g.addEdge(30, 31);
g.addEdge(30, 39);
g.addEdge(31, 32);
g.addEdge(31, 40);
g.addEdge(32, 33);
g.addEdge(32, 41);
g.addEdge(33, 34);
g.addEdge(33, 42);
g.addEdge(34, 43);
g.addEdge(36, 37);
g.addEdge(36, 45);
g.addEdge(37, 38);
g.addEdge(37, 46);
g.addEdge(38, 39);
g.addEdge(38, 47);
g.addEdge(39, 40);
g.addEdge(39, 48);
g.addEdge(40, 41);
g.addEdge(40, 49);
g.addEdge(41, 42);
g.addEdge(41, 50);
g.addEdge(42, 43);
g.addEdge(42, 51);
g.addEdge(43, 52);
g.addEdge(45, 46);
g.addEdge(46, 47);
g.addEdge(47, 48);
g.addEdge(48, 49);
g.addEdge(49, 50);
g.addEdge(50, 51);
g.addEdge(51, 52);
g.addEdge(33, 44);
g.addEdge(8, 17);
g.addEdge(17, 26);
g.addEdge(42, 35);
g.addEdge(44, 53);
g.addEdge(52, 53);
g.addEdge(27, 34);
g.addEdge(42, 53);
g.addEdge(26, 35);
g.addEdge(43, 44);
g.addEdge(35, 34);
g.addEdge(26, 25);


}