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
#include <unordered_set>

/* JT Adams - 12050538
* My attempt at a genetic algorithm to solve the TSP
*/
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
class TravelingSalesmanGA{
public:
    TravelingSalesmanGA(std::string file){
        generateGraph(file);
        // printVector(graph);
    }

    // Main function
    void runGeneticAlgorithm(int numGenerations, int populationSize, double mutationRate) {
        // Generate initial population
        population = generateInitialPopulation(populationSize, graph.size());
        // Main loop for iterating over generations
        for (int generation = 1; generation <= numGenerations; ++generation) {
            // Evaluate fitness of each tour in the population
            std::vector<int> fitnessValues;
            std::vector<std::vector<int>> fitnessPaths;
            std::cout << "New Gen!" << std::endl;
            for (int f = 0;f < population.size();f++) {
                // std::cout << "Finding the length of this path: " << std::endl;
                // printVector1D(population[f]);
                int pathLen = fitness(population[f]);
                fitnessPaths.push_back(population[f]);
                fitnessValues.push_back(pathLen);
                //std::cout << "The length of this path is " << pathLen << std::endl;
            }
            // printVector1D(fitnessValues);
            // std::cout << "Selecting" << std::endl;
            // Perform selection to choose parents for mating

            std::vector<std::vector<int>> selectedParents = tournamentSelection(population, fitnessValues, 10, populationSize);
            // std::cout << "The length of the selected Parents is: " << selectedParents.size() << std::endl;
            // Perform crossover to generate offspring
            // std::cout << "Crossover TIME" << std::endl;
            // std::sort(selectedParents.begin(),selectedParents.end(),[this](const std::vector<int> &a, const std::vector<int> &b){
            //     return fitness(a) < fitness(b);
            // });
            std::vector<std::vector<int>> offspring;
            for (int i = 1; i < populationSize / 2; i++) {
                auto parents = std::make_pair(selectedParents[i], selectedParents[i-1]);
                auto children = orderCrossover(parents.first, parents.second);
                offspring.push_back(children.first);
                offspring.push_back(children.second);

            }
            // std::cout << "Crossover Complete " << std::endl;
            // Perform mutation on the offspring
            for (auto& tour : offspring) {
                if (rand() / (RAND_MAX + 1.0) < mutationRate) {
                    swapMutation(tour);
                }
            }

            // Replace current generation with the offspring
            generationalReplacement(population, offspring);

            // Output information about the current generation, e.g., best fitness
            int min_index = 0;
            int max_path = CEIL;
            for(int q=0;q<fitnessValues.size();q++){
                if (fitnessValues[q] < max_path){
                    max_path = fitnessValues[q];
                    min_index = q;
                }
            }
            //int bestFitness = *std::min_element(fitnessValues.begin(), fitnessValues.end());
            std::cout << "Generation " << generation << ", Best Fitness: " << max_path << std::endl;
            printVector1D(fitnessPaths[min_index]);


            // Try to take care of memory
            fitnessValues.clear();
        }
    }

private:
    std::vector<std::vector<int>> graph;
    std::vector<int> bestPath;
    int lowestCost;

    // TODO: Add more stuff
    std::vector<std::vector<int>> population;
    std::vector<std::vector<int>> nextpop;

    // Swap was not working, make it more simple?
    std::pair<std::vector<int>,std::vector<int>> simpleSwap(const std::vector<int>& parent1,
                                                             const std::vector<int>& parent2){
                                                                
    }
    
    // Builds the adj matrix from input file
    void generateGraph(std::string file){
            std::ifstream g; // Declare input file stream
            g.open(file); // Open the stream
            std::string line;
            while (getline(g,line)){
                std::istringstream iss(line);
                std::vector<int> nums;
                // Now we need to add each number to nums
                int num;
                // So this moves each number from the stream into the int var num
                while (iss>>num){
                    nums.push_back(num);
                }
                graph.push_back(nums);
            }
            g.close();
    }

    // This generates a small amount of nearest neighbors to create initial population
    std::vector<int> generateNN(int start,std::vector<int>&visited,int currLevel,std::vector<int> &path){
        path.push_back(currLevel);
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
        if (neighborCity == -1){
            // Add the cost to return to start and then return cost
            // path.push_back(start);
            return path;
        } else {
            return generateNN(start,visited,neighborCity,path);
        }
    }

    // Generate a "random" initial population of tours
    std::vector<std::vector<int>> generateInitialPopulation(int populationSize, int numCities) {

        // To generate a "semi" efficient population, I will run nearest neighbor on the first 100 nodes
        std::vector<std::vector<int>> population;
        
        for(int i=0;i<100;i++){
            std::vector<int> visited(numCities,0);
            std::vector<int> path;
            path = generateNN(i,visited,i,path);
            for(int j=0;j<(populationSize/100);j++){
                population.push_back(path);
            }
            population.push_back(path);
        }

        std::random_device rd;
        auto rng = std::default_random_engine { rd() };
        std::shuffle(std::begin(population),std::end(population),rng);
        printVector(population);
        return population;
    } 

        // Calculate path cost
    int fitness(const std::vector<int>& path){
        int length = 0;
        int cities = path.size();

        // Subtract 1 from cities
        for(int i=1;i< cities; i++){
            int currCity = path[i];
            int prevCity = path[i-1];

            if (prevCity < currCity){
                prevCity = path[i];
                currCity = path[i-1];
            }
            //std::cout << "We are going from city " << prevCity << " to " << currCity << std::endl;
            //std::cout << "Adding this distance to length: " << graph[prevCity][currCity] << std::endl;
            length += graph[prevCity][currCity];
        }
        // Add from the last city back to the front
        if (path[0] < path[path.size()-1]){
            length += graph[path[path.size()-1]][path[0]];
        } else {
            length += graph[path[0]][path[path.size()-1]];
        }
        return length;
    }

    // Use Roullete Wheel selection
    // After testing this kinda sucks and does not do what I want at all
    std::vector<std::vector<int>> selection(std::vector<std::vector<int>> population,std::vector<int> fitnessValues,int numParents){

        
        std::vector<std::vector<int>> selectedParents;
        long totalFitness = 0;
        for (auto &n : fitnessValues){
            totalFitness += n;
        }
        // std::cout << "In the selection, with a total fitness of: " << totalFitness << std::endl;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<long> dis(0, totalFitness - 1);

        for (int i = 0; i < numParents; ++i) {
            long randNum = dis(gen);
            // std::cout << "For loopin" << std::endl;
            //std::cout << "Trying to generate long number" << std::endl;
            //std::cout << "Success!" << std::endl;
            // std::cout << "Gen Num" << std::endl;
            long cumulativeFitness = 0;
            int j = 0;
            while (cumulativeFitness <= randNum) {
                cumulativeFitness += fitnessValues[j];
                j++;
            }
            // std::cout << "This is the _ route pushed back: " << i << std::endl;
            selectedParents.push_back(population[j - 1]);
        }
        return selectedParents;
    }

    // Roulette wheel selection sucks, so lets try tournament selection
    std::vector<std::vector<int>> tournamentSelection(const std::vector<std::vector<int>>& population, const std::vector<int>& fitnessValues, int tournamentSize, int numParents) {
        std::vector<std::vector<int>> selectedParents;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, population.size() - 1);

        while (selectedParents.size() < numParents) {
            // Randomly select individuals for the tournament
            std::vector<int> tournamentParticipants;
            for (int i = 0; i < tournamentSize; ++i) {
                tournamentParticipants.push_back(dis(gen));
            }

            // Find the winner of the tournament (individual with lowest fitness)
            int winner = tournamentParticipants[0];
            for (int i = 1; i < tournamentSize; ++i) {
                if (fitnessValues[tournamentParticipants[i]] < fitnessValues[winner]) {
                    winner = tournamentParticipants[i];
                }
            }

            // Add the winner to the selected parents
            selectedParents.push_back(population[winner]);
        }

        return selectedParents;
    }

    // Function to perform swap mutation
    void swapMutation(std::vector<int>& tour) {
        int numCities = tour.size();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, numCities - 2); // Exclude the first and last cities

        // Choose two random cities to swap
        int city1 = dis(gen);
        int city2 = dis(gen);
        while (city1 == city2) {
            city2 = dis(gen); // Ensure city1 and city2 are distinct
        }

        // Perform the swap
        std::swap(tour[city1], tour[city2]);
    }

    // This is a main, attempts to implmenet order crossover because it is one of the best for GA's
    // However, I am not sure exactly how to do it so this is incredbily inefficient, but works!
    std::pair<std::vector<int>, std::vector<int>> orderCrossover(const std::vector<int>& parent1,
                                                             const std::vector<int>& parent2) {

        //std::cout << "Attempting crossover " << std::endl;
        // printVector1D(parent1);
        // printVector1D(parent2);
        int numCities = parent1.size();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, numCities - 2); // Exclude the first and last cities

        // Choose two random crossover points
        int point1 = dis(gen);
        int point2 = dis(gen);
        while (point1 == point2) {
            // std::cout << "Stuck" << std::endl;
            point2 = dis(gen); // Ensure point1 and point2 are distinct
        }
        if (point1 > point2) {
            std::swap(point1, point2); // Swap if point1 > point2
        }
        

        // Create empty offspring tours
        std::vector<int> offspring1(numCities, -1);
        std::vector<int> offspring2(numCities, -1);
        // Create city sets and populate them
        std::unordered_set<int> addedCities1;
        std::unordered_set<int> addedCities2;
        for(int i=0;i<numCities;i++){
            addedCities1.insert(i);
            addedCities2.insert(i);
        }

        for (int i=point1; i<=point2;i++){
            offspring1[i] = parent1[i];
            addedCities1.erase(parent1[i]);
        }

        // Now we need to copy the remaining NON DUPLICANT cities from parent2 TO offspring1
        // Duplicates will be dealt with later
        for(int i=0;i<offspring1.size();i++){
            if (offspring1[i] == -1){
                if (addedCities1.find(parent2[i])!=addedCities1.end()){
                    // This means the number has not been removed from the set and can be added!
                    offspring1[i] = parent2[i];
                    addedCities1.erase(parent2[i]);
                }
            }
        }
        
        // Now same for offspring 2
        for (int i=point1; i<=point2;i++){
            offspring2[i] = parent2[i];
            addedCities2.erase(parent2[i]);
        }

        // Now we need to copy the remaining NON-ZERO, and NON DUPLICANT cities from parent1 TO offspring2
        // Duplicates will be dealt with later
        for(int i=0;i<offspring2.size();i++){
            if (offspring2[i] == -1){
                if (addedCities2.find(parent1[i])!=addedCities2.end()){
                    // This means this number is not in the set and can be added
                    offspring2[i] = parent1[i];
                    addedCities2.erase(parent1[i]);
                }
            }
        }
        // std::cout << "After adding the window and adding non-duplicates from the parents, we have:" << std::endl;
        // printVector1D(offspring1);
        // printVector1D(offspring2);
        // Now lets fix everything up
        std::vector<int> remainingCities1(addedCities1.begin(),addedCities1.end());
        std::vector<int> remainingCities2(addedCities2.begin(),addedCities2.end());
        auto rng = std::default_random_engine { rd() };
        std::shuffle(std::begin(remainingCities1),std::end(remainingCities1),rng);
        std::shuffle(std::begin(remainingCities2),std::end(remainingCities2),rng);
        // std::cout << "Remaining, randomly shifted vals: " << std::endl;
        // printVector1D(remainingCities1);
        // printVector1D(remainingCities2);
        int c1 = 0;
        int c2 = 0;
        for(int w=0;w<numCities;w++){
            if (offspring1[w] == -1){
                offspring1[w] = remainingCities1[c1];
                c1++;
            }
        }

        for(int w=0;w<numCities;w++){
            if (offspring2[w] == -1){
                offspring2[w] = remainingCities2[c2];
                c2++;
            }
        }
        // std::cout << "Crossover Complete " << std::endl;
        // printVector1D(offspring1);
        // printVector1D(offspring2);
        return {offspring1, offspring2};
    }

    // Function to perform generational replacement
    void generationalReplacement(std::vector<std::vector<int>>& population,
                                std::vector<std::vector<int>>& offspring) {
        // Replace the current generation population with the offspring
        population = offspring;
        offspring.clear();
}
};

// Lets gooooooooo
int main(int argc, char** argv){
    std::srand(std::time(0));
    std::string filePath = argv[1];
    // Create the map object
    TravelingSalesmanGA t = TravelingSalesmanGA(filePath);
    t.runGeneticAlgorithm(100,1000,0.05);
}

