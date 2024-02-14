#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include <chrono>

using namespace std;

// Function to calculate the total tour distance
double tourDistance(const vector<int>& tour, const vector<vector<double>>& graph) {
    double totalDistance = 0.0;
    int numCities = tour.size();
    for (int i = 0; i < numCities; ++i) {
        int j = (i + 1) % numCities;
        totalDistance += graph[tour[i]][tour[j]];
    }
    return totalDistance;
}

// Function to initialize pheromone matrix
void initializePheromones(vector<vector<double>>& pheromones, double initialValue) {
    int numCities = pheromones.size();
    for (int i = 0; i < numCities; ++i) {
        for (int j = 0; j < numCities; ++j) {
            pheromones[i][j] = initialValue;
        }
    }
}

// Function to update pheromones based on ant's tour
void updatePheromones(vector<vector<double>>& pheromones, const vector<int>& tour, double evaporationRate, double Q) {
    double tourLength = tourDistance(tour, pheromones);
    int numCities = tour.size();

    for (int i = 0; i < numCities; ++i) {
        int j = (i + 1) % numCities;
        pheromones[tour[i]][tour[j]] += Q / tourLength;
        pheromones[tour[j]][tour[i]] += Q / tourLength;
    }

    for (int i = 0; i < numCities; ++i) {
        for (int j = 0; j < numCities; ++j) {
            pheromones[i][j] *= (1.0 - evaporationRate);
        }
    }
}

// Function to select the next city based on pheromone trails and heuristic information
int selectNextCity(int currentCity, const vector<vector<double>>& pheromones, const vector<vector<double>>& graph, const vector<bool>& visited, double alpha, double beta) {
    int numCities = graph.size();
    double totalProbability = 0.0;
    vector<double> probabilities(numCities, 0.0);

    for (int i = 0; i < numCities; ++i) {
        if (!visited[i]) {
            probabilities[i] = pow(pheromones[currentCity][i], alpha) * pow(1.0 / graph[currentCity][i], beta);
            totalProbability += probabilities[i];
        }
    }
    std::default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, totalProbability);
    double randomValue = distribution(generator);

    double cumulativeProbability = 0.0;
    for (int i = 0; i < numCities; ++i) {
        if (!visited[i]) {
            cumulativeProbability += probabilities[i];
            if (cumulativeProbability >= randomValue) {
                return i;
            }
        }
    }

    return -1; // Error: should not reach here
}

// Function to construct a solution (tour) using ant colony optimization
vector<int> constructSolution(const vector<vector<double>>& pheromones, const vector<vector<double>>& graph, double alpha, double beta) {
    int numCities = graph.size();
    vector<bool> visited(numCities, false);
    vector<int> tour(numCities, -1);
    std::default_random_engine generator;
    // Start from a random city
    uniform_int_distribution<int> distribution(0, numCities - 1);
    int initialCity = distribution(generator);
    tour[0] = initialCity;
    visited[initialCity] = true;

    // Construct the rest of the tour
    for (int i = 1; i < numCities; ++i) {
        int currentCity = tour[i - 1];
        int nextCity = selectNextCity(currentCity, pheromones, graph, visited, alpha, beta);
        tour[i] = nextCity;
        visited[nextCity] = true;
    }

    return tour;
}

// Function to perform ant colony optimization for the TSP
vector<int> solveACO(const vector<vector<double>>& graph, int numAnts, int numIterations, double evaporationRate, double alpha, double beta, double Q) {
    int numCities = graph.size();
    vector<vector<double>> pheromones(numCities, vector<double>(numCities, 1.0));
    initializePheromones(pheromones, 1.0);

    vector<int> bestTour;
    double bestTourLength = numeric_limits<double>::infinity();

    for (int iter = 0; iter < numIterations; ++iter) {
        for (int ant = 0; ant < numAnts; ++ant) {
            vector<int> tour = constructSolution(pheromones, graph, alpha, beta);
            double tourLength = tourDistance(tour, graph);
            if (tourLength < bestTourLength) {
                bestTourLength = tourLength;
                bestTour = tour;
            }
            updatePheromones(pheromones, tour, evaporationRate, Q);
        }
    }

    return bestTour;
}

int main() {
    // Example usage
    vector<vector<double>> graph = {{0, 10, 15, 20},
                                     {10, 0, 35, 25},
                                     {15, 35, 0, 30},
                                     {20, 25, 30, 0}};
    int numAnts = 10;
    int numIterations = 100;
    double evaporationRate = 0.1;
    double alpha = 1.0;
    double beta = 2.0;
    double Q = 1.0;

    auto start = chrono::steady_clock::now();
    vector<int> bestTour = solveACO(graph, numAnts, numIterations, evaporationRate, alpha, beta, Q);
    auto end = chrono::steady_clock::now();

    cout << "Best tour found: ";
    for (int city : bestTour) {
        cout << city << " ";
    }
    cout << endl;
    cout << "Length of the tour: " << tourDistance(bestTour, graph) << endl;

    cout << "Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " milliseconds" << endl;

    return 0;
}