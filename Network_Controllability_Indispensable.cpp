#include <iostream>
#include<bits/stdc++.h> 
#include <unordered_map>


#pragma GCC diagnostic ignored "-Wc++11-extensions"

using namespace std;


#define NIL 0 
#define INF INT_MAX


enum regulationType {
	Node, Promoter, Enhancer, Both
};
static std::map<std::string, regulationType> regType_map;

void regTypeInitialize(){
	regType_map["Node"] = Node; //Not implemented here, different code exists for this
	regType_map["Promoter"] = Promoter;
	regType_map["Enhancer"] = Enhancer;
	regType_map["Both"] = Both;
}

class Edge{
  public:
    int index1, index2, type;
    //type: 2 = Promoter, 3 = Enhancer, 1 = Both 
    vector<int> neighbors(int i);

  Edge(int index1, int index2, int type) : index1(index1), index2(index2) , type(type) {}
};


class biPartiteGraph{
  int leftSize, rightSize;
  list<int> *adj; 
  int *pairU, *pairV, *dist; 

public:
  biPartiteGraph(int l, int r);
  ~biPartiteGraph();
  void addEdge(int u, int v);
  bool bfs(int u_r, int v_r);
  bool dfs(int u, int u_r, int v_r);
  int hopcroftKarp(unordered_map<string, int> nodes, unordered_map<string, int> V, int n, int u, int v); 

};

//Inspired from https://www.geeksforgeeks.org/hopcroft-karp-algorithm-for-maximum-matching-set-2-implementation/
int biPartiteGraph::hopcroftKarp(unordered_map<string, int> nodes, unordered_map<string, int> V, int n_r, int u_r, int v_r) { 
	pairU = new int[leftSize+1]; //stores pairs of matching u on left side
	pairV = new int[rightSize+1]; //stores pairs of matching v 
	dist = new int[leftSize+1]; 

	for (int u=0; u<=leftSize; u++) {
		pairU[u] = NIL; 
	}
		
	for (int v=0; v<=rightSize; v++) {
		pairV[v] = NIL; 
	}

	int result = 0; 

	while (bfs(u_r, v_r)) { 
		//find free vertex
		for (int u=1; u<=leftSize; u++) {
			// if free and there is augmenting path
			if (u_r != u && pairU[u]==NIL && dfs(u, u_r, v_r)) {
				result++; 
			}
		}
	} 
	int dNodes = 0;

	for(int v=1; v<=rightSize; v++){
		if( v_r != v && pairV[v]==NIL){
			dNodes++;
		}
	}

	for(auto n : nodes){
		if(n.second != n_r){
			string v1 = "V_" + n.first;
			if(V[v1] == 0 && V[v1] != v_r){
				dNodes++;
			}
		}
	}

	//cout << "Matching Size: " << result << endl;
	//cout << "Driver Nodes Size: " << dNodes << endl;
	return dNodes; 
} 

//returns if augmenting path exists
bool biPartiteGraph::bfs(int u_r, int v_r) { 
	queue<int> Q; 
	for (int u=1; u<=leftSize; u++) {
		if(u != u_r){
			if (pairU[u]==NIL) { 
				dist[u] = 0; 
				Q.push(u); 
			} 
			else dist[u] = INF; 
		} 		
	} 

	dist[NIL] = INF; 
	while (!Q.empty()) { 
		int u = Q.front(); 
		Q.pop(); 
		if (dist[u] < dist[NIL]) { 
			list<int>::iterator i; 
			for (i=adj[u].begin(); i!=adj[u].end(); ++i) { 
				int v = *i; 
				if (v != v_r && dist[pairV[v]] == INF) { 
					dist[pairV[v]] = dist[u] + 1; 
					Q.push(pairV[v]); 
				} 
			} 
		} 
	} 

	return (dist[NIL] != INF); 
} 

//returns if augmenting path existing from free vertex u
bool biPartiteGraph::dfs(int u, int u_r, int v_r) { 
	if (u != u_r && u != NIL) { 
		list<int>::iterator i; 
		for (i=adj[u].begin(); i!=adj[u].end(); ++i) { 
			int v = *i; 

			if(v != v_r){
				if (dist[pairV[v]] == dist[u]+1) { 

					if (dfs(pairV[v], u_r, v_r) == true) { 
						pairV[v] = u; 
						pairU[u] = v; 
						return true; 
					} 
				} 
			}			
		} 
		dist[u] = INF; 
		return false; 
	} 
	return true; 
} 

biPartiteGraph::~biPartiteGraph() { 

   delete[] adj; 
   adj = NULL;
   delete[] pairU;
   delete[] pairV;
   delete[] dist;
 
} 

biPartiteGraph::biPartiteGraph(int l, int r) { 
    this->leftSize = l; 
    this->rightSize = r; 
    adj = new list<int>[l+1]; 
} 
  
void biPartiteGraph::addEdge(int u, int v) { 
    adj[u].push_back(v); 
} 


void readEdges(unordered_map<string, int>& nodes, vector<Edge>& edges, string filename, unordered_map<string, int>& U, unordered_map<string, int>& V, vector<Edge>& bedges){
  ifstream file(filename);
  string line;
  int row, col, u_i, v_i, type;
  int node_count = 0, u_count=1, v_count=1;

  while(getline(file, line)){
    istringstream iss(line);

    int j = 0;
    
    
    do{
      string s, u1, v1;
      iss >> s;
      if (s != ""){
        if(j%3 == 0){
          u1 = "U_" + s;

          unordered_map<string, int>::const_iterator iter = nodes.find(s);
          
          if(iter != nodes.end()){
            row = nodes[s];

            unordered_map<string, int>::const_iterator it = U.find(u1);
            if(it == U.end()){
            	std::pair<std::string,int> new_uNode (u1,u_count);
	            u_count++;
	            U.insert(new_uNode);
            }

          }else{
            std::pair<std::string,int> newNode (s,node_count);
            row = node_count;
            node_count++;
            nodes.insert(newNode);

            std::pair<std::string,int> new_uNode (u1,u_count);
            u_count++;
            U.insert(new_uNode);
          }
          j++;
          u_i = U[u1];
        }else if(j%3 == 1){
          v1 = "V_" + s;
          unordered_map<string, int>::const_iterator iter = nodes.find(s);
          if(iter != nodes.end())
          {
            col = nodes[s];
            unordered_map<string, int>::const_iterator it = V.find(v1);
            if(it == V.end()){
            	std::pair<std::string,int> new_vNode (v1,v_count);
	            v_count++;
	            V.insert(new_vNode);
            }
          }else{
            std::pair<std::string,int> newNode (s,node_count);
            col = node_count;
            node_count++;
            nodes.insert(newNode);

			std::pair<std::string,int> new_vNode (v1,v_count);
            v_count++;
            V.insert(new_vNode);
            
          }
          j++;
          v_i = V[v1];
        }else{
          	switch (regType_map[s]){          	
          		case Promoter:
          			type = 1;
          			break;
          		case Enhancer:
          			type = 2;
          			break;
          		case Both:
          			type = 3;
          			break;
          		default:
          			type = 1;
          			break;
          		}

          }
      }
    }while(iss);
    //cout << row <<  "  " << col << " " << type << endl;
    edges.push_back(Edge(row, col, type));    
    bedges.push_back(Edge(u_i, v_i, type));
  }
}

string findIndex(unordered_map<string, int> nodes,int i){
	for (auto it = nodes.begin(); it != nodes.end(); ++it) {
  		if (it->second == i) return it->first;
	} 
	return 0;
}


int main(int argc, char** argv){

	std::string edgeListFile(argv[1]);
	std::string regType(argv[2]);

	regTypeInitialize();

	//string edgeListFile = "../py/Toy/Human/human_toy_unique.txt";
	//string edgeListFile = "../py/Data/tmp/01_neurons_fetal_brain.ncol";
	//string edgeListFile = "../py/LineGraph/HighLevel/LG_01_neurons_fetal_brain.ncol";
	
	unordered_map<string, int> nodes;
	unordered_map<string, int> U, V;

	vector<Edge> edges, bip_edges; 
  	//cout << "Reading Edges\n";
  
  	
	readEdges(nodes, edges, edgeListFile, U, V, bip_edges);

	int nodeSize = nodes.size();
	biPartiteGraph bg(U.size(), V.size()); 

	int count = 0;

	vector<Edge>::iterator it; 
	for(it = bip_edges.begin(); it != bip_edges.end(); it++) {
	    bg.addEdge(it->index1, it->index2);
        count++;
	}
	cout << count << " Bipartite Edges" << endl;

	//cout << "Calculating Matching\n";
	int originalDriverSize = bg.hopcroftKarp(nodes,V, -1, -1, -1);

	cout << "OrigialDriverSize: " << originalDriverSize << endl;
	cout << "Running " << edgeListFile << " Regulatory Type removal" << regType << endl;
    
    
    //return 0;

	for(auto n : nodes){	

		int n_r, u_r = -1, v_r = -1;
		int b = 0;

		string v = "V_" + n.first;
		string u = "U_" + n.first;
		n_r = n.second;

		if(U.count(u) == 1){
			u_r = U.at(u);
		}

		if(V.count(v) == 1){
			v_r = V.at(v);
		}

		biPartiteGraph new_bg(U.size(), V.size()); 

		vector<Edge>::iterator it; 
		for(it = bip_edges.begin(); it != bip_edges.end(); it++) {
			if(!(it->index1 == u_r && it->type == regType_map[regType])){
				new_bg.addEdge(it->index1, it->index2);
				b++;
			}    
		}

		
		int n_ds = new_bg.hopcroftKarp(nodes,V, -1, -1, -1);

		if(originalDriverSize < n_ds){
			cout << n.first << " indispensable " << originalDriverSize << "  " << n_ds <<  endl; 
		}else if(originalDriverSize == n_ds){
			cout << n.first << " neutral " << originalDriverSize << "  " << n_ds << endl; 
		}else{
			cout << n.first << " dispensable " << originalDriverSize << "  " << n_ds << endl; 
		}
		
		
	}
		
    return 0; 
  
}
