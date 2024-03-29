#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <ctime>   // For time()
#include <climits> // Ceiling
#include <cmath> // for pow()
#include <numeric> // The weird probability function
#include <fstream> // For reading in the file!
#include <sstream> // For parsing input
#include <algorithm> // For sort

const static int CEIL = 2147483646;

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

template<typename T>
void printVector1D(const std::vector<T>& vec) {
    std::cout << "Vector: { ";
    for (int i=0;i<vec.size();i++){
        std::cout << vec[i] << " ";   
    }
    std::cout << "}" << std::endl;

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

..
            }
            map = tempMap; // Time to find out if I know how memory works
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

int makeRandomPick(const std::vector<double>& weights) {
    // Calculate cumulative probabilities
    std::vector<double> cumulativeProbabilities(weights.size());
    std::partial_sum(weights.begin(), weights.end(), cumulativeProbabilities.begin());

    // Generate a random number between 0 and the sum of weights
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, cumulativeProbabilities.back());
    double randomValue = dis(gen);

    // Find the index corresponding to the first cumulative probability that is greater than the random number
    int selectedIndex = 0;
    while (selectedIndex < cumulativeProbabilities.size() && cumulativeProbabilities[selectedIndex] < randomValue) {
        selectedIndex++;
    }
    
    return selectedIndex;
}

// Ant Colony Optimization NOW!
// This function will have a single ant going through the graph one node at a time
// constrained to move in a cycle
std::pair<std::vector<int>,int> traverseGraph(std::vector<std::vector<int>> graph, int sourceNode, std::vector<std::vector<double>> &pheromones, int graphSize,std::vector<int> visited){
    double ALPHA = 0.9;
    double BETA = 1.5;

    
    visited[sourceNode] = 1; // Mark curr City to visited

    // Initialize The Path Taken
    std::vector<int> cycle;
    cycle.push_back(sourceNode);

    int cycleLength = 0;

    int steps = 0;
    int currentCity = sourceNode;
    int total_length = 0;
    std::vector<int> jumpNeighbors;
    std::vector<double> jumpValues;
    while (steps < graphSize-1){
        std::cout << "The path length so far is : " << total_length << std::endl;
        double cumProbabilities = 0.0; // Cumulative probablities, used for generating weighted jump value
        for (int i=0;i<graph.size();i++){
            if (visited[i] == 0){
                // We need to make sure we are reading the graph right due to the triangular format
                if (currentCity < i){
                    // If the current city is less than i, then we can just use graph[i][currentCity]
                    double pheromoneLevel = std::max(pheromones[i][currentCity],(0.00001));
                    double v = (pow(pheromoneLevel,ALPHA)) / (pow(double(graph[i][currentCity]),BETA));
                    jumpNeighbors.push_back(i);
                    jumpValues.push_back(v);

                } else {
                    // If the current city is GREATER than i, then we can just use graph[currentCity][i]
                    double pheromoneLevel = std::max(pheromones[currentCity][i],(0.00001));
                    double v = (pow(pheromoneLevel,ALPHA)) / (pow(double(graph[currentCity][i]),BETA));
                    jumpNeighbors.push_back(i);
                    jumpValues.push_back(v);

                }
            }
        }
        // Now we need to choose the next node
        int index = makeRandomPick(jumpValues);
        // std::cout << "Selected 'RANDOM INDEX' of : " << index << std::endl;
        int jump = jumpNeighbors[index];
        visited[jump] = 1;
        total_length += (currentCity < jump) ? graph[currentCity][jump] : graph[jump][currentCity];
        currentCity = jump;
        cycle.push_back(currentCity);
        steps++;
        jumpNeighbors.clear();
        jumpValues.clear();
    }
    // Adjust path to home
    if (sourceNode < currentCity){
        total_length += graph[currentCity][sourceNode];
    } else {
        total_length += graph[sourceNode][currentCity];
    }
    visited.clear();
    // std::cout << " Path len " << cycle.size() << std::endl;
    // std::cout << "total length:" << total_length << std::endl;
    return std::make_pair(cycle,total_length);
}

void degradePheromones(std::vector<std::vector<double>> &pheromones,double degradation){
    for(int i=0;i<pheromones.size();i++){
        for(int j=0;j<pheromones.size();j++){
            pheromones[i][j] *= degradation;
        }
    }
}

void deleteCycles(std::vector<int> &cycle){

}

int antColonyOptimization(std::vector<std::vector<int>> graph, std::vector<std::vector<double>> &pheromones, int iterations, int ants, double q=10.0, double degradation=0.95){
    std::vector<int> bestCycle;
    int bestLength = CEIL; // Max integer for the 'best cycle'
    std::vector<std::pair<std::vector<int>,int>> cycles;
    std::vector<int> visited(graph.size(),0);
    for (int i=0; i < iterations; i++){
        // This is ugly, probably will need helper function to sort!
        std::cout << "Iteration: " << i << std::endl;
        // printVector(pheromones);
        for (int j=0; j<ants; j++){
            int startingCity = (std::rand()%(graph.size()-1)+1);
            std::pair<std::vector<int>,int> currentPair = traverseGraph(graph,startingCity,pheromones,graph.size(),visited);
            if (currentPair.second > 0){
                cycles.push_back(currentPair);
            }
        }

        // So once we have generated our 50 random paths, we need to sort for the best fit ones

        std::sort(cycles.begin(),cycles.end(), [](const auto& a, const auto& b){
            return a.second > b.second;
        });
        int length = cycles.size();
        for (int k=0;k<length;k++){ // Update pheromones
            std::vector<int> cycle = cycles[k].first;
            int currLength = cycles[k].second;
            std::cout << "Curr length: " << currLength << std::endl;
            int dontClear = 0;
            if (currLength<bestLength){
                bestLength = currLength;
                bestCycle = cycle;
                dontClear = 1;
            }

            double delta = q / currLength;
            // std::cout << "Delta: " << delta << std::endl;
            int j = 0;
            while (j < cycle.size()-1 && j+1 < cycle.size()){
                int currCity = cycle[j];
                int nextCity = cycle[j+1];
                if (currCity < nextCity){
                    pheromones[nextCity][currCity] += delta;
                } else{
                    pheromones[currCity][nextCity] += delta;
                }
                j += 1;
            }
            // Now adjust pheromones for the route back to home
            if (cycle[0] < cycle[j]){
                pheromones[cycle[j]][cycle[0]] += delta;
            } else {
                pheromones[cycle[0]][cycle[j]] += delta;
            }
            degradePheromones(pheromones,degradation);
            cycle.clear();
        }
        cycles.clear();
    }
    std::cout << "The best path is: " << bestLength << std::endl;
    pheromones.clear();
    return bestLength;
}


// Generate base pheromone vector
std::vector<std::vector<double>> generatePheromones(std::vector<std::vector<int>> graph, double defaultValue){
    std::vector<std::vector<double>> pheromones;
    for (int i=0;i<graph.size();i++){
        std::vector<double> currentLevel;
        for (int j=0;j<graph[i].size();j++){
            if (i != j){
                currentLevel.push_back(defaultValue);
            } else{
                currentLevel.push_back(0);
            }
        }
        pheromones.push_back(currentLevel);
    }

    return pheromones;
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
    std::vector<std::vector<int>> g = m.getMap();
    // Blank!

    // Now lets try to bruteforce TSP!
    // A lot of websites are saying to use "next_permutation" but I want to do it on my own
    // std::vector<std::vector<int>> static_cities = { { 0 },
    //                                                 { 10, 0, },
    //                                                 { 15, 35, 0, },
    //                                                 { 20, 25, 30, 0 } };
    // std::cout << "The size of the static cities is " << static_cities.size() << std::endl;
    // // We can use the above deep copy to attempt to run a brute force on this
    // std::cout << "Brute forcing the static city layout!" << std::endl;
    // std::vector<int> visited(static_cities.size(), 0);
    // visited[0] = 1; // Starting at city 0
    // std::cout << "Starting at city: " << 0 << std::endl;
    // int final1 = tspBruteForce(static_cities,0,0,visited,0);
    // std::cout << "The best path is: " << final1 << std::endl;

    // Nearest Neighbor Time
    // // Copy map
    // std::vector<std::vector<int>> mapCopy = m.getMap();
    // std::cout << "\nNow we're trying the nearest neighbor method on our input generated map!" << std::endl;
    // // Declare mini
    // int mini = CEIL;
    // for(int i=0;i<mapCopy.size();i++){
    //     std::vector<int> neighbors(mapCopy.size(),0);
    //     // visited[i] = 1; Not needed here
    //     std::cout << "Stariting at city: " << i << std::endl;
    //     int final = tspNearestNeighbor(mapCopy,i,i,neighbors,0);
    //     std::cout << "The best route at this level is: " << final << std::endl;
    //     mini = std::min(mini,final);
    // }
    // std::cout << "The best route found is: " << mini << std::endl;
    

   // Generate pheromone vector
    std::vector<std::vector<double>> pheromones = generatePheromones(g,10);
    printVector(pheromones);
    //std::pair<std::vector<int>,int> returns = traverseGraph(static_cities,start_city,pheromones,sz);
    int length = antColonyOptimization(g,pheromones,50,50);

    std::cout << "Total length: " << length << std::endl;
    

    return 0;
}