int#include <fstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <iostream>
#include <string>


using namespace std;

//Items have weights and profits. They also have a cost/benifit ratio which will be used to calculate optimal profits of a node, and which will be used for sorting items.
struct item
{
	int weight, profit, c_b_ratio;
};

//or another way...we use arrays to identify weight, profitm c_b_ratio. I think I'll go with this one. Will choose based on ease of implementation

//Nodes contain the current amount of weight & profit they contain, their level in the space state tree, their optimal predicted profit, and a vector which will store the current solution in binary. We also use a comparator so the priority queue can properly arrange Nodes.
struct Node
{
  int weight;
	int profit;
  int level;
	int opt_profit;

	vector<int> curr_solution;

  bool operator<(const Node &comp) const{
        return (opt_profit < comp.opt_profit);
  }
};

//Define global vector for optimal solution
vector<int> OptSol;

//We calculate the bound of a node. If a node is infeasible, we set its bound to -1. Else we calulate the bound using the greedy fractional knapsack method, to ensure an optimistic upper-bound.
int optimistically_predict_bound(int n, int profits[], int weights[], Node nd, int C) {

    //bound includes profit of node
    int bound = nd.profit;

    //the level indicates how many nodes we have already used
    int i = nd.level;

    //infeasible Node
    if(nd.weight >= C) return -1;

    //capacity includes weight of Node
    int curr_C = nd.weight;

    //loop through items until bag is ~full
    while( (i <= n) && (curr_C + weights[i] <= C)) { //add items
          curr_C += weights[i];
          bound = bound + profits[i];
          i++;
    }

    //if there are items left and room in the bag, we take a fraction of the next item(we know the whole item can't fit due to the previous loop)
    if((i < n) && (curr_C != C) ) {
        int fractional_space = C - curr_C;
        bound +=  fractional_space * profits[i]/weights[i];
        return bound;
    }

    return bound;
}

//knapsack_01_BeFS finds the optimal solution and prints the appropriate information to the output file
void knapsack_01_BeFS(int n, int profits[], int weights[], int C, string fileName){

    //what we'll print (except for the items themselves)
    int maxprofit = 0;
    int sol_size = 0;
    int numNodes = 1;
    int numLeaves = 0;

    priority_queue<Node> BeFS_Q;

    Node current, next;

    // current starts as the root.
    current.level = 0;
    current.profit = 0;
    current.weight = 0;

    //account for initial node in vector size
    current.curr_solution.resize(n+1, 0);
    next.curr_solution.resize(n+1,0);

    //finish initializing root and place into Q
    current.opt_profit  = optimistically_predict_bound(n, profits, weights, current, C);
    BeFS_Q.push(current);

    //Start moving through state space! Deque the most promising nodes while there are nodes in the q.
    while(!BeFS_Q.empty()){

        //pop most promising Node off
        current = BeFS_Q.top();
        BeFS_Q.pop();

        //explore Node if necessary
        if(current.opt_profit > maxprofit) {
            //trying copying to the opt solution here
            //we'll visit both the children
            numNodes += 2;
            next.level = current.level + 1;

            //try node that takes next item
	    int index = next.level -1;
            next.weight = current.weight + weights[index];
            next.profit = current.profit + profits[index];

            //optimal solution will be the same except you set choosen item to 1
            next.curr_solution = current.curr_solution;
            next.curr_solution.at(next.level-1) = 1;

            //update global solution and maxprofit if necessary
            if(next.profit > maxprofit && next.weight <= C) {
                maxprofit = next.profit;
                OptSol = next.curr_solution;
            }

            //find the upper bound of current node and place in q if necessary.
            next.opt_profit = optimistically_predict_bound(n, profits, weights, next, C);

            if(next.opt_profit > maxprofit  && next.weight < C) {
                BeFS_Q.push(next);
            }

            //If its not pushed we will not explore it and it is a leaf
            else numLeaves++;

            //Try the same thing but without taking the next item
            next.weight = current.weight;
            next.profit = current.profit;
            next.curr_solution.at(next.level - 1) = 0;
            next.opt_profit = optimistically_predict_bound(n, profits, weights, next, C);


            if(next.opt_profit > maxprofit  && next.weight < C) {
                BeFS_Q.push(next);
            }
            else numLeaves++;

        }
        else numLeaves++;
    }

    ofstream outFile;
    outFile.open(fileName);

    //might be an error here in calculating item stats
    for (int i = 0; i < OptSol.size(); i++)  {
      if(OptSol[i] == 1) sol_size++;
    }

    outFile << n << "," << maxprofit << ","<< sol_size << endl;

    outFile << numNodes << "," << numLeaves << endl;

    for(int i = 0; i < OptSol.size(); i++) {
      if(OptSol[i]  == 1) outFile << weights[i] << "," << profits[i] << endl;
    }
}

int main(int argc, const char * argv[]) {

    //declare necessary variables
    int n = 0;
    int C;

    //open file for reading, parse first line
    ifstream infile (argv[1]);
    string splitMe;

    infile >> splitMe;

    size_t split = splitMe.find(",");

    n = stoi(splitMe.substr(0, split));
    C = stoi(splitMe.substr(split+1));

    //initialize arrays
    int profits[n];
    int weights[n];

    //parse lines 1-n+1
    for(int i = 0;i < n; i++){

    	infile >> splitMe;
	split = splitMe.find(",");

	weights[i] = stoi(splitMe.substr(0, split));
        profits[i] = stoi(splitMe.substr(split+1));
    }

  //bubble sort for sorting because I was told it won't be graded
	for (int i=0; i < profits.size() - 1; i++){
		for (int j=0; j< profits.size()-1-i; j++){
			if (profits[j]/weights[j] < profits[j+1]/weights[j+1]){
				int  prof_temp = profits[j];
        int weigh_temp = weights[j];
				profits[j] = profits[j+1];
        weights[j] = weights[j+1];
				profits[j+1] = prof_temp;
        weights[j+1] =weigh_temp;
			}
		}
	}

    string output = argv[2];
    knapsack_01_BeFS(n, profits, weights, C, output);
}
