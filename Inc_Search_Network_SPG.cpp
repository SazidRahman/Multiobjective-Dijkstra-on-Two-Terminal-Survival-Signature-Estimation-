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

const int V = 238;              // Total number of nodes
const int v = V-2;              // Nodes that can fail (source and sink don't fail)
const int m1 = 128;              // Number of nodes in the first group/class of nodes
const int m2 = 75;            // Number of nodes in the second group/class of nodes
const int m3 = v-(m1+m2);      // Number of nodes in the third group/class of nodes
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

    std::vector<int> class_1;
    std::vector<int> class_2;
    std::vector<int> class_3;

    // Fill class_1 with numbers from 0 to 29
    for (int i = 1; i <= 128; ++i) {
        class_1.push_back(i);
    }

    // Fill class_2 with numbers from 30 to 59
    for (int i = 129; i <= 203; ++i) {
        class_2.push_back(i);
    }

    // Fill class_3 with numbers from 60 to 80
    for (int i = 204; i <= 236; ++i) {
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
g.addEdge(0, 1);
g.addEdge(1, 2);
g.addEdge(2, 3);
g.addEdge(3, 4);
g.addEdge(4, 5);
g.addEdge(5, 6);
g.addEdge(6, 7);
g.addEdge(7, 8);
g.addEdge(8, 9);
g.addEdge(9, 10);
g.addEdge(10, 11);
g.addEdge(11, 12);
g.addEdge(12, 13);
g.addEdge(13, 14);
g.addEdge(14, 15);
g.addEdge(15, 16);
g.addEdge(16, 17);
g.addEdge(17, 18);
g.addEdge(18, 19);
g.addEdge(19, 20);
g.addEdge(20, 21);
g.addEdge(21, 22);
g.addEdge(22, 23);
g.addEdge(23, 24);
g.addEdge(24, 25);
g.addEdge(25, 26);
g.addEdge(26, 27);
g.addEdge(27, 28);
g.addEdge(28, 29);
g.addEdge(29, 30);
g.addEdge(30, 31);
g.addEdge(31, 32);
g.addEdge(32, 33);
g.addEdge(33, 34);
g.addEdge(34, 35);
g.addEdge(35, 36);
g.addEdge(36, 37);
g.addEdge(37, 38);
g.addEdge(38, 39);
g.addEdge(39, 40);
g.addEdge(40, 41);
g.addEdge(41, 42);
g.addEdge(42, 43);
g.addEdge(43, 44);
g.addEdge(44, 45);
g.addEdge(45, 46);
g.addEdge(46, 47);
g.addEdge(47, 48);
g.addEdge(48, 49);
g.addEdge(49, 50);
g.addEdge(50, 51);
g.addEdge(51, 52);
g.addEdge(52, 53);
g.addEdge(53, 54);
g.addEdge(54, 55);
g.addEdge(55, 56);
g.addEdge(56, 57);
g.addEdge(57, 58);
g.addEdge(58, 59);
g.addEdge(59, 60);
g.addEdge(60, 61);
g.addEdge(61, 62);
g.addEdge(62, 63);
g.addEdge(63, 64);
g.addEdge(64, 65);
g.addEdge(65, 66);
g.addEdge(66, 67);
g.addEdge(67, 68);
g.addEdge(68, 69);
g.addEdge(69, 70);
g.addEdge(70, 71);
g.addEdge(71, 72);
g.addEdge(72, 73);
g.addEdge(73, 74);
g.addEdge(74, 75);
g.addEdge(75, 76);
g.addEdge(76, 77);
g.addEdge(77, 78);
g.addEdge(78, 79);
g.addEdge(79, 80);
g.addEdge(80, 81);
g.addEdge(81, 82);
g.addEdge(82, 83);
g.addEdge(83, 84);
g.addEdge(84, 85);
g.addEdge(85, 86);
g.addEdge(86, 87);
g.addEdge(87, 88);
g.addEdge(88, 89);
g.addEdge(89, 90);
g.addEdge(90, 91);
g.addEdge(91, 92);
g.addEdge(92, 93);
g.addEdge(93, 94);
g.addEdge(94, 95);
g.addEdge(95, 96);
g.addEdge(96, 97);
g.addEdge(97, 98);
g.addEdge(98, 99);
g.addEdge(99, 100);
g.addEdge(100, 101);
g.addEdge(101, 102);
g.addEdge(102, 103);
g.addEdge(103, 104);
g.addEdge(104, 105);
g.addEdge(105, 106);
g.addEdge(106, 107);
g.addEdge(107, 108);
g.addEdge(108, 109);
g.addEdge(109, 110);
g.addEdge(110, 111);
g.addEdge(111, 112);
g.addEdge(112, 113);
g.addEdge(113, 114);
g.addEdge(114, 115);
g.addEdge(115, 116);
g.addEdge(116, 117);
g.addEdge(117, 118);
g.addEdge(118, 119);
g.addEdge(119, 120);
g.addEdge(120, 121);
g.addEdge(121, 122);
g.addEdge(122, 123);
g.addEdge(123, 124);
g.addEdge(124, 125);
g.addEdge(125, 126);
g.addEdge(126, 127);
g.addEdge(127, 128);
g.addEdge(128, 129);
g.addEdge(129, 130);
g.addEdge(130, 131);
g.addEdge(131, 132);
g.addEdge(132, 133);
g.addEdge(133, 134);
g.addEdge(134, 135);
g.addEdge(135, 136);
g.addEdge(136, 137);
g.addEdge(137, 138);
g.addEdge(138, 139);
g.addEdge(139, 140);
g.addEdge(140, 141);
g.addEdge(141, 142);
g.addEdge(142, 143);
g.addEdge(143, 144);
g.addEdge(144, 145);
g.addEdge(145, 146);
g.addEdge(146, 147);
g.addEdge(147, 148);
g.addEdge(148, 149);
g.addEdge(149, 150);
g.addEdge(150, 151);
g.addEdge(151, 152);
g.addEdge(152, 153);
g.addEdge(153, 154);
g.addEdge(154, 155);
g.addEdge(155, 156);
g.addEdge(156, 157);
g.addEdge(157, 158);
g.addEdge(158, 159);
g.addEdge(159, 160);
g.addEdge(160, 161);
g.addEdge(161, 162);
g.addEdge(162, 163);
g.addEdge(163, 164);
g.addEdge(164, 165);
g.addEdge(165, 166);
g.addEdge(166, 167);
g.addEdge(167, 168);
g.addEdge(168, 169);
g.addEdge(169, 170);
g.addEdge(170, 171);
g.addEdge(171, 172);
g.addEdge(172, 173);
g.addEdge(173, 174);
g.addEdge(174, 175);
g.addEdge(175, 176);
g.addEdge(176, 177);
g.addEdge(177, 178);
g.addEdge(178, 179);
g.addEdge(179, 180);
g.addEdge(180, 181);
g.addEdge(181, 182);
g.addEdge(182, 183);
g.addEdge(183, 184);
g.addEdge(184, 185);
g.addEdge(185, 186);
g.addEdge(186, 187);
g.addEdge(187, 188);
g.addEdge(188, 189);
g.addEdge(189, 190);
g.addEdge(190, 191);
g.addEdge(191, 192);
g.addEdge(192, 193);
g.addEdge(193, 194);
g.addEdge(194, 195);
g.addEdge(195, 196);
g.addEdge(196, 197);
g.addEdge(197, 198);
g.addEdge(198, 199);
g.addEdge(199, 200);
g.addEdge(200, 201);
g.addEdge(201, 202);
g.addEdge(202, 203);
g.addEdge(203, 204);
g.addEdge(204, 205);
g.addEdge(205, 206);
g.addEdge(206, 207);
g.addEdge(207, 208);
g.addEdge(208, 209);
g.addEdge(209, 210);
g.addEdge(210, 211);
g.addEdge(211, 212);
g.addEdge(212, 213);
g.addEdge(213, 214);
g.addEdge(214, 215);
g.addEdge(215, 216);
g.addEdge(216, 217);
g.addEdge(217, 218);
g.addEdge(218, 219);
g.addEdge(219, 220);
g.addEdge(220, 221);
g.addEdge(221, 222);
g.addEdge(222, 223);
g.addEdge(223, 224);
g.addEdge(224, 225);
g.addEdge(225, 226);
g.addEdge(226, 227);
g.addEdge(227, 228);
g.addEdge(228, 229);
g.addEdge(229, 230);
g.addEdge(230, 231);
g.addEdge(231, 232);
g.addEdge(232, 233);
g.addEdge(233, 234);
g.addEdge(234, 235);
g.addEdge(235, 236);
g.addEdge(236, 237);
g.addEdge(237, 18);
g.addEdge(1, 3);
g.addEdge(4, 7);
g.addEdge(8, 13);
g.addEdge(14, 18);
g.addEdge(22, 25);
g.addEdge(27, 30);
g.addEdge(32, 36);
g.addEdge(40, 45);
g.addEdge(48, 52);
g.addEdge(55, 60);
g.addEdge(62, 67);
g.addEdge(70, 74);
g.addEdge(77, 81);
g.addEdge(84, 89);
g.addEdge(91, 95);
g.addEdge(98, 102);
g.addEdge(105, 109);
g.addEdge(112, 116);
g.addEdge(119, 123);
g.addEdge(126, 130);
g.addEdge(133, 137);
g.addEdge(140, 144);
g.addEdge(147, 151);
g.addEdge(154, 158);
g.addEdge(161, 165);
g.addEdge(168, 172);
g.addEdge(175, 179);
g.addEdge(182, 186);
g.addEdge(189, 193);
g.addEdge(196, 200);
g.addEdge(203, 207);
g.addEdge(210, 214);
g.addEdge(217, 221);
g.addEdge(224, 228);
g.addEdge(231, 235);
g.addEdge(237, 1);
g.addEdge(5, 12);
g.addEdge(15, 21);
g.addEdge(23, 29);
g.addEdge(31, 37);
g.addEdge(39, 43);
g.addEdge(46, 50);
g.addEdge(53, 58);
g.addEdge(60, 65);
g.addEdge(68, 73);
g.addEdge(76, 80);
g.addEdge(83, 87);
g.addEdge(90, 94);
g.addEdge(97, 101);
g.addEdge(104, 108);
g.addEdge(111, 115);
g.addEdge(118, 122);
g.addEdge(125, 129);
g.addEdge(132, 136);
g.addEdge(139, 143);
g.addEdge(146, 150);
g.addEdge(153, 157);
g.addEdge(160, 164);
g.addEdge(167, 171);
g.addEdge(174, 178);
g.addEdge(181, 185);
g.addEdge(188, 192);
g.addEdge(195, 199);
g.addEdge(202, 206);
g.addEdge(209, 213);
g.addEdge(216, 220);
g.addEdge(223, 227);
g.addEdge(230, 234);
g.addEdge(237, 2);



}