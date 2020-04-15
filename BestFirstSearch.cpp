#include <fstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <iostream>
#include <string>


using namespace std;

//Items have weights and profits. They also have a cost/benifit ratio which will be used to calculate optimal profits of a node, and which will be used for sorting items.
struct item
{
	double wieght, profit, c_b_ratio;
};

//or another way...we use arrays to identify weight, profitm c_b_ratio. I think I'll go with this one. Will choose based on ease of implementation

//Nodes contain the current amount of weight & profit they contain, their level in the space state tree, their optimal predicted profit, and a vector which will store the current solution in binary. We also use a comparator so the priority queue can properly arrange Nodes.
struct Node
{
  double weight;
	double profit;
  double level;
	double opt_profit;

	vector<double> curr_solution;

  bool operator<(const Node &comp) const{
        return (opt_profit < comp.opt_profit);
  }
};

//Define global vector for optimal solution
vector<double> OptSol;

//We calculate the bound of a node. If a node is infeasible, we set its bound to -1. Else we calulate the bound using the greedy fractional knapsack method, to ensure an optimistic upper-bound.
double optimistically_predict_bound(double n, double profits[], double weights[], Node nd, double C) {

    //bound includes profit of node
    double bound = nd.profit;

    //the level indicates how many nodes we have already used
    double i = nd.level;

    //infeasible Node
    if(nd.weight >= C) return -1;

    //capacity includes wieght of Node
    double curr_C = nd.weight;

    //loop through items until bag is ~full
    while( (i <= n) && (curr_C + weights[i] <= C)) { //add items
          curr_C += weights[j];
          bound = bound + p[j];
          i++;
    }

    //if there are items left and room in the bag, we take a fraction of the next item(we know the whole item can't fit due to the previous loop)
    if((i < n) && (curr_C != C) ) {
        double fractional_space = C - curr_C;
        bound +=  fractional_space * profits[i]/wieghts[i];
        return bound;
    }

    return bound;
}

//knapsack_01_BeFS finds the optimal solution and prints the appropriate information to the output file
void knapsack_01_BeFS(double n, double profits[], double wieghts[], double maxprofit, double C, string fileName){

    //what we'll print (except for the items themselves)
    int max_profit = 0;
    int sol_size = 0;
    int numNodes = 1;
    int numLeaves = 0;

    priority_queue<Node> BeFS_Q;

    Nodes current, next;

    // current starts as the root.
    current.level = 0;
    current.profit = 0;
    current.weight = 0;

    //account for initial node in vector size
    current.curr_solution.resize(n+1, 0);
    next.curr_solution.resize(n+1,0);

    //finish initializing root and place into Q
    current.bound = optimistically_predict_bound(n, profits, wieghts, current, C);
    BeFS_Q.push(v);

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
            next.weight = current.weight + wieghts[next.level-1];
            next.profit = current.profit + profits[neXt.level-1];

            //optimal solution will be the same except you set choosen item to 1
            next.curr_solution = current.curr_solution;
            next.curr_solution.at(next.level-1) = 1;

            //update global solution and maxprofit if necessary
            if(next.profit > maxprofit && next.weight <= C) {
                maxprofit = next.profit;
                OptSol = next.curr_solution;
            }

            //find the upper bound of current node and place in q if necessary.
            next.opt_profit = optimistically_predict_bound(n, profit, weights, C, next);

            if(next.opt_profit > maxprofit  && next.weight < C) {
                BeFS_Q.push(next);
            }

            //If its not pushed we will not explore it and it is a leaf
            else numLeaves++;

            //Try the same thing but without taking the next item
            next.weight = current.weight;
            next.profit = current.profit;
            next.curr_solution.at(u.level - 1) = 0;
            next.opt_profit = optimistically_predict_bound(n, profits, weights, C, next);


            if(next.opt_profit > maxprofit  && next.weight < C) {
                BeFS_Q.push(u);
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

    outFile << n << "," << max_profit << ","<< sol_size << endl;

    outFile << numNodes << "," << numLeaves << endl;

    for(int i = 0; i < OptSol.size(); i++) {
      if(OptSol[i]  == 1) outFile << weights[i] << "," << profits[i] << endl;
    }
}

int main(int argc, const char * argv[]) {

    //declare necessary variables
    int n = 0;
    int C;

    //open file for reading
    ifstream infile (argv[1]);
    string line;
    string buf;


    //parse first line
    if (getline(infile, line)) {
        stringstream s (line);
        if (getline(s, buf, ',')) {
            stringstream ss (buf);
            ss >> n;
        }
        if (getline(s, buf, ',')) {
            stringstream ss (buf);
            ss >> C;
        }
    }

    //initialize arrays
    int profits[n];
    int wieghts[n];

    int count = 0;
    while (getline (infile, line)) {
        stringstream s (line);
        if (getline(s,buf, ',')) {
            stringstream ss (buf);
            ss >> in;
            wieghts[count] = in;
        }
        if (getline(s, buf, ',')) {
            int in;
            stringstream ss (buf);
            ss >> in;
            profits[count] = in;
        }
        count ++;
    }

  //bubble sort for sorting because I was told it won't be graded
	for (int i=0; i < profits.size() - 1; i++){
		for (int j=0; j< profits.size()-1-i; j++){
			if (profits[j]/weight[j] < profits[j+1]/weights[j+1]){
				double  prof_temp = profits[j];
        double weigh_temp = weights[j];
				profits[j] = profits[j+1];
        wieghts[j] = wieghts[j+1]
				profits[j+1] = prof_temp;
        weights[j+1] =wiegh_temp;
			}
		}
	}

    output = argv[2];
    knapsack_01_BeFS(n, pSorted, wSorted, C, output);
}
