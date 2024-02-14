# Learning the Traveling Salesman problem!
# Our nice graph of cities
# Of which, we need to find the shortest path to all
graph = [
    [5,2,1],
    [2,3,5],
    [9,3,6]
]

# Essentially backtracking no?
# Permutations
def find_shortest_path(start,cost,visited,route):
    print(f"Current cost: {cost}")
    shortest = float('inf')
    path = []
    visited[start[0]][start[1]] = 1
    all = 0
    for i in range(len(graph)):
        for j in range(len(graph[0])):
            if visited[i][j] == 0:
                all = 1
                print(f"Visiting {i},{j}")
                route.append((i,j))
                res = find_shortest_path([i,j],cost+(graph[start[0]][start[1]]+graph[i][j]),visited,route)
                if res[0] < shortest:
                    shortest = res[0]
                    path = res[1]
                visited[i][j] = 0
                route.pop()
    if all == 1:
        return [shortest,path]
    else:
        return [cost,route]

    
path = []
shortest_length = float('inf') # Store a large number as the shortest length to start!
for i in range(len(graph)):
    for j in range(len(graph[0])):
        # For each "city" in the graph, we're going to try and find the shortest possible route
        # Horrible time complexity in this form!

        # We do need to keep track of what cities we have visited
        visited = [[0] * len(graph[0])] * len(graph)
        print(f"Starting at {i},{j}")
        res = find_shortest_path([i,j],0,visited,[(i,j)])
        if res[0] < shortest_length:
            shortest_length = res[0]
            path = res[1]

print(shortest_length)
print(path)
        
