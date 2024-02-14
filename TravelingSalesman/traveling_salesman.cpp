#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <ctime>   // For time()
#include <climits>
#include <thread> // Trying to multithread
#include <fstream> // For reading in the file!
#include <sstream> // For parsing input

const static int CEIL = INT_MAX;

// Print 2D Vector
template<typename T>
void printVector(const std::vector<std::vector<T>>& vec) {
    for (const auto& row : vec) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

class Map{
    public:
        /**
         * @brief Construct a new Map object, for now the map size will be 4x5 (20 total cities)
         * 
         * @param maxSize Maximum "weight" of each city
         *  
         */
        Map(int maxSize,int rows,int cols){
            std::vector<std::vector<int>> tempMap;
            for (int i=0;i<rows;i++){
                std::vector<int> curr;
                for (int j=0;j<=i;j++){
                    if (j==i){
                        curr.push_back(0);
                    } else {
                        // Generate Random Number
                        int randomNum = (std::rand()%maxSize) + 1; // S hould generate a random number?
                        curr.push_back(randomNum);
                    }
                }
                tempMap.push_back(curr); // Append this level of the map to the map at large
            }
            map = tempMap; // Time to find out if I know how memory works!
        }

        Map(std::string file){
            std::ifstream graph; // Declare input file stream
            graph.open(file); // Open the stream
            std::string line;
            while (getline(graph,line)){
                std::istringstream iss(line);
                std::vector<int> nums;
                // Now we need to add each number to nums
                int num;
                // So this moves each number from the stream into the int var num
                while (iss>>num){
                    nums.push_back(num);
                }
                map.push_back(nums);
            }
            graph.close();
        }
        
        
        std::vector<std::vector<int>> getMap(){
            return map;
        }

        friend std::ostream& operator<<(std::ostream& os, const Map& m);
    private:
        std::vector<std::vector<int>> map;
};

// This may be useless with my new printVector helper
std::ostream& operator<<(std::ostream& os, const Map& m){
    std::string output = "";
    for (int i=0;i<m.map.size();i++){
        for (int j=0;j<m.map[i].size();j++){
            os << m.map[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}

/**
 * @brief Calculate the TSP problem using the Nearest Neighbor logic
 * 
 * @param graph 
 * @param startCity 
 * @param currLevel 
 * @param visited 
 * @param cost 
 * @return int 
 */
int tspNearestNeighbor(std::vector<std::vector<int>> graph,int startCity, int currLevel, std::vector<int> &visited,int cost){
    int mini = CEIL;
    // Mark the the current city as visited
    visited[currLevel] = 1;
    int neighborWeight = CEIL;
    int neighborCity = -1;
    for (int i=0;i<graph.size();i++){
        // Look for the smallest city
        if (i!=currLevel && visited[i] != 1){
            if (i<currLevel){
                if(graph[currLevel][i] < neighborWeight){
                    neighborWeight = graph[currLevel][i];
                    neighborCity = i;
                }
            } else{
                if(graph[i][currLevel] < neighborWeight){
                    neighborWeight = graph[i][currLevel];
                    neighborCity = i;
                }
            }
        }   
    }
    // Now determine if we can actually go to the city
    if (neighborCity == -1){
        // Add the cost to return to start and then return cost
        int toAdd = 0;
        if (startCity < currLevel){
            toAdd = graph[currLevel][startCity];
        } else {
            toAdd = graph[startCity][currLevel];
        }
        return (cost+toAdd);
    } else {
        return tspNearestNeighbor(graph,startCity,neighborCity,visited,(cost+neighborWeight));
    }
}

/**
 * @brief Brute forces the graph to attempt to solve the traveling salesman problem
 * 
 */
int tspBruteForce(std::vector<std::vector<int>> graph,int startCity, int currLevel,std::vector<int> &visited,int cost){
    // std::cout << "We're in city " << currLevel+1 << " the cost so far is " << cost << std::endl;
    // Create a "visited" graph?
    int mini = 9999999;
    bool flag = false;
    // std::cout << "Current cost is: " << cost << std::endl;
    for (int i=0;i<graph.size();i++){
        if (i != currLevel && visited[i] != 1){
            visited[i] = 1;
            flag = true;
            int toAdd = 0;
            if (i < currLevel){
                toAdd = graph[currLevel][i];
            }else{
                toAdd = graph[i][currLevel];
            }
            mini = std::min(mini,tspBruteForce(graph,startCity,i,visited,cost+toAdd));
            visited[i] = 0;
        }

    }
    if (!flag){
        // We need to add how much it costs to return to start
        cost = cost + (graph[currLevel][startCity]);
        // std::cout << "The cost when ending at " << currLevel+1 << " is: " << cost << std::endl;
        return cost;
    } else{
        // std::cout << "Returning from city " << currLevel+1 << " with a min of: " << mini << std::endl;
        return mini;
    }
}

// Ant Colony Optimization NOW!
// This function will have a single ant going through the graph one node at a time
// constrained to move in a cycle
void traverseGraph(std::vector<std::vector<int>> graph, int sourceNode, std::vector<std::vector<int>> pheromones){
    float ALPHA = 0.9;
    float BETA = 1.5;
    std::vector<int> visited(graph.size(),0);
    visited[sourceNode] = 1;

    std::vector<int> cycle = {sourceNode};
    int steps = 0;
    int current = sourceNode;
    int total_length = 0;
    while (steps < graph.size()){
        std::vector<int> jumpNeighbors;
        std::vector<int> jumpValues;

        for (int i=0;i<graph.size();i++){
            if (visited[i] == 0){
                float pheromoneLevel = std::max(1)
            }
        }
    }

}
int main(int argc, char** argv){
    // Seed the random number generator!
    // This allows for the humbers to be "true" random (not really)
    std::srand(std::time(nullptr));
    // Generate boundaries for random graph generation
    std::string filePath = argv[1];
    // Create the map object
    Map m = Map(filePath);
    // std::cout << m;
    // Blank!

    // Now lets try to bruteforce TSP!
    // A lot of websites are saying to use "next_permutation" but I want to do it on my own
    std::vector<std::vector<int>> static_cities = { { 0 },
                                                    { 10, 0, },
                                                    { 15, 35, 0, },
                                                    { 20, 25, 30, 0 } };
    // We can use the above deep copy to attempt to run a brute force on this
    std::cout << "Brute forcing the static city layout!" << std::endl;
    std::vector<int> visited(static_cities.size(), 0);
    visited[0] = 1; // Starting at city 0
    std::cout << "Starting at city: " << 0 << std::endl;
    int final1 = tspBruteForce(static_cities,0,0,visited,0);
    std::cout << "The best path is: " << final1 << std::endl;

    // Nearest Neighbor Time
    // Copy map
    std::vector<std::vector<int>> mapCopy = m.getMap();
    std::cout << "\nNow we're trying the nearest neighbor method on our input generated map!" << std::endl;
    // Declare mini
    int mini = CEIL;
    for(int i=0;i<mapCopy.size();i++){
        std::vector<int> neighbors(mapCopy.size(),0);
        // visited[i] = 1; Not needed here
        std::cout << "Stariting at city: " << i << std::endl;
        int final = tspNearestNeighbor(mapCopy,i,i,neighbors,0);
        std::cout << "The best route at this level is: " << final << std::endl;
        mini = std::min(mini,final);
    }
    std::cout << "The best route found is: " << mini << std::endl;

}