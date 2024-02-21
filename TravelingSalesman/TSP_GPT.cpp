#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

// Define constants
const int NUM_CITIES = 5;
const int POPULATION_SIZE = 10;
const int NUM_GENERATIONS = 1000;
const double MUTATION_RATE = 0.1;

// Define type for representing tours
using Tour = vector<int>;

// Function to generate a random lower triangular adjacency matrix
vector<vector<int>> generateRandomAdjacencyMatrix() {
    vector<vector<int>> adjacencyMatrix(NUM_CITIES, vector<int>(NUM_CITIES, 0));
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(1, 100); // Adjust range as needed

    for (int i = 0; i < NUM_CITIES; ++i) {
        for (int j = i + 1; j < NUM_CITIES; ++j) {
            adjacencyMatrix[i][j] = dist(gen);
        }
    }
    return adjacencyMatrix;
}

// Function to calculate tour length
int calculateTourLength(const Tour& tour, const vector<vector<int>>& adjacencyMatrix) {
    int length = 0;
    for (int i = 0; i < NUM_CITIES - 1; ++i) {
        length += adjacencyMatrix[tour[i]][tour[i + 1]];
    }
    length += adjacencyMatrix[tour[NUM_CITIES - 1]][tour[0]]; // Return to the starting city
    return length;
}

// Function to initialize a random tour
Tour generateRandomTour() {
    Tour tour(NUM_CITIES);
    for (int i = 0; i < NUM_CITIES; ++i) {
        tour[i] = i;
    }
    random_shuffle(tour.begin() + 1, tour.end());
    return tour;
}

// Genetic algorithm function
Tour solveTSPGenetic(const vector<vector<int>>& adjacencyMatrix) {
    // Initialize population
    vector<Tour> population;
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population.push_back(generateRandomTour());
    }

    // Main loop
    for (int gen = 0; gen < NUM_GENERATIONS; ++gen) {
        // Evaluate fitness
        vector<pair<int, Tour>> fitness;
        for (const auto& tour : population) {
            fitness.push_back({calculateTourLength(tour, adjacencyMatrix), tour});
        }
        sort(fitness.begin(), fitness.end());

        // Select parents (roulette wheel selection)
        vector<Tour> parents;
        for (int i = 0; i < POPULATION_SIZE / 2; ++i) {
            int idx1 = rand() % POPULATION_SIZE;
            int idx2 = rand() % POPULATION_SIZE;
            parents.push_back(fitness[idx1].second);
            parents.push_back(fitness[idx2].second);
        }

        // Crossover (partially mapped crossover)
        vector<Tour> offspring;
        for (int i = 0; i < POPULATION_SIZE / 2; ++i) {
            int crossoverPoint = rand() % (NUM_CITIES - 1) + 1;
            Tour parent1 = parents[i];
            Tour parent2 = parents[i + 1];
            Tour child1 = parent1;
            Tour child2 = parent2;
            for (int j = 0; j < crossoverPoint; ++j) {
                auto it1 = find(child1.begin(), child1.end(), parent2[j]);
                auto it2 = find(child2.begin(), child2.end(), parent1[j]);
                iter_swap(it1, it2);
            }
            offspring.push_back(child1);
            offspring.push_back(child2);
        }

        // Mutation
        for (auto& tour : offspring) {
            if ((double)rand() / RAND_MAX < MUTATION_RATE) {
                int idx1 = rand() % NUM_CITIES;
                int idx2 = rand() % NUM_CITIES;
                swap(tour[idx1], tour[idx2]);
            }
        }

        // Replace population with offspring
        population = offspring;
    }

    // Return the best tour
    return fitness[0].second;
}

int main() {
    // Generate a random lower triangular adjacency matrix
    vector<vector<int>> adjacencyMatrix = generateRandomAdjacencyMatrix();

    // Solve TSP using genetic algorithm
    Tour bestTour = solveTSPGenetic(adjacencyMatrix);

    // Output the best tour
    cout << "Best tour: ";
    for (int city : bestTour) {
        cout << city << " ";
    }
    cout << endl;

    return 0;
}