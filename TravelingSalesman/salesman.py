import random


class Graph():
    def __init__(self, cities, l, pheromones):
        self.cities = cities
        self.adj_list = l
        self.pheromones = []
        # Build the pheromones
        for i in range(len(self.adj_list)):
            pheromone_level = []
            for j in range(len(self.adj_list[i])):
                pheromone_level.append(pheromones)
            self.pheromones.append(pheromone_level)
        
        # Now our graph is complete!


def ACO(graph, its, ants, q, degradation_factor = .9):
    lowest_journey = None
    cost = None
    print(graph.pheromones)
    for it in range(its):
        paths = [traverse_graph(graph,random.randint(0,graph.cities-1)) for _ in range(ants)]
        paths.sort(key= lambda x:x[1])
        paths = paths[:ants//2]

        if lowest_journey:
            paths.append([lowest_journey,cost])
        
        for p,c in paths:
            if not cost:
                lowest_journey = p
                cost = c
            else:
                if c < cost:
                    cost = c
                    lowest_journey = p 
            delta = q/c
            i = 0
            while i < len(p)-1 and i+1 < len(p):
                next_city = p[i+1]
                curr_city = p[i]
                if curr_city < next_city:
                    graph.pheromones[next_city][curr_city] += delta
                else:
                    graph.pheromones[curr_city][next_city] += delta
                i += 1
            # Now pheromone the end
            if p[i] < p[0]:
                graph.pheromones[p[0]][p[i]] += delta
            else:
                graph.pheromones[p[i]][p[0]] += delta
            
            for row in range(len(graph.pheromones)):
                for col in range(len(graph.pheromones[row])):
                    graph.pheromones[row][col] *= degradation_factor


    return [cost,lowest_journey]

# Send the ant through the map, start at a random point
def traverse_graph(g,s):
    ALPHA = 0.5
    BETA = 1.5
    visited = [0] * len(g.adj_list) # Empty visited set
    visited[s] = 1 # Set the first city as visited

    path = [s] # Begin path tracker
    steps = 0 # Count how far we've come
    journey = 0 # How long the entire journey takes
    city = s # Tracks where we're at
    while steps < len(g.adj_list)-1:
        prob_neighbors = []
        probibilities = []

        for i in range(len(g.adj_list)):
            if visited[i] == 0:
                if i < s:
                    # use g.adj_list[s][i]
                    this_pheromone = max(g.pheromones[s][i],1e-5)
                    v = (this_pheromone**ALPHA) / (g.adj_list[s][i]**BETA) # This does the P**alpha/D**beta
                    prob_neighbors.append(i) # Store the city we may go to 
                    probibilities.append(v) # Store the weight of the pheromone 
                else:
                    # use g.adj_list[i][s]
                    this_pheromone = max(g.pheromones[i][s],1e-5)
                    v = (this_pheromone**ALPHA) / (g.adj_list[i][s]**BETA) # This does the P**alpha/D**beta
                    prob_neighbors.append(i) # Store the city we may go to 
                    probibilities.append(v) # Store the weight of the pheromone 
        next_city = random.choices(prob_neighbors,weights=probibilities)[0]

        # Add the distance to journey length
        if next_city < city:
            # adj_list[city][next_city]
            journey += g.adj_list[city][next_city]
        else:
            journey += g.adj_list[next_city][city]
        visited[next_city] = 1
        city = next_city
        path.append(city)
        steps += 1
    # Return home
    if city < s:
        journey += g.adj_list[s][city]
    else:
        journey += g.adj_list[city][s]
    #print(f"The length of this path is {journey}")
    #print(f"{path}")
    return path,journey

# Get the input file
file = input()
graph = []
with open(file, 'r') as f:
    for line in f:
        l = [int(num) for num in line.split(' ')]
        graph.append(l)

inputGraph = Graph(len(graph),graph,50)

r = ACO(inputGraph,100,50,10)
print(f"The best path has length {r[0]} and is the path\n{r[1]}")
