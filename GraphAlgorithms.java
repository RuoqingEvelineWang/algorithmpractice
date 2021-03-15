// Adjacency List Graph with ArrayList (Java)
// 1. iterative BFS
// 2. iterative DFS
// 3. recursive DFS
// 4. is cyclic BFS (directed)
// 5. is cyclic DFS (directed)
// 6. is cyclic (undirected)
// 7. is bipartite (undirected, BFS)
// 8. BFS find shortest path
// 9. bidirectional search find shortest path
// 10. Kruskal's minimum spanning tree (negatively weighted, undirected)
// 11. Prim's minimum spanning tree (negatively weighted, undirected)
// 12. Dijkstra's shortest path algorithm (weighted)
// 13. Dijkstra's shortest path algorithm with priority queue (weighted)
// 14. Bellman-Ford algorithm (negatively weighted)
// 15. Hamiltonian cycle
// 16. Travelling Salesman Problem


// adjacency list representation of graph with an ArrayList of ArrayLists, can be directed, undirected or weighted
public class Graph {
    public static int nV;
    public static ArrayList<ArrayList<Integer>> edges;
    public static ArrayList<ArrayList<Integer>> weights;
    
    Graph(int nV, boolean weighted) {
        this.nV = nV;
        this.edges = new ArrayList<ArrayList<Integer>>(nV);
        for (int i = 0; i < nV; i++)
            edges.add(new ArrayList<Integer>());
        if (weighted) {
            this.weights = new ArrayList<ArrayList<Integer>>(nV);
            for (int i = 0; i < nV; i++)
                weights.add(new ArrayList<Integer>());
        }
    }
    
    public void addEdge(int v1, int v2, boolean directed) {
        edges.get(v1).add(v2);
        if (!directed)
            edges.get(v2).add(v1);
    }
    
    public void addEdgeWeighted(int v1, int v2, int weight) {
        edges.get(v1).add(v2);
        edges.get(v2).add(v1);
        weights.get(v1).add(weight);
        weights.get(v2).add(weight);
    }
    
    public void print() {
        for (int i = 0; i < nV; i++) {
            System.out.print("vertex " + i + ":");
            for (int a : edges.get(i))
                System.out.print(" " + a);
            System.out.println();
        }
    }
    
    // 1. Breadth First Search Algorithm
    // add the source vertex in queue
    // dequeue one vertex at a time and enqueue all adjacent unvisited vertices until queue is empty
    public static void BFS_iterative(int v) {
        boolean[] visited = new boolean[nV];
        Queue<Integer> queue = new LinkedList<>();
        
        visited[v] = true;
        queue.add(v);
        
        while (queue.size() != 0) {
            int cur = queue.remove();
            System.out.println(cur);
            for (int i : edges.get(cur)) {
                if (!visited[i]) {
                    visited[i] = true;
                    queue.add(i);
                }
            }
        }
    }
    
    // 2. Depth First Search Algorithm
    // add the source vertex in stack
    // pop one vertex at a time and push all adjacent unvisited vertices until stack is empty
    public static void DFS_iterative(int v) {
        boolean[] visited = new boolean[nV];
        Stack<Integer> stack = new Stack<>();
        
        visited[v] = true;
        stack.push(v);
        
        while (!stack.isEmpty()) {
            int cur = stack.pop();
            System.out.println(cur);
            for (int i : edges.get(cur)) {
                if (!visited[i]) {
                    visited[i] = true;
                    stack.push(i);
                }
            }
        }
    }
    
    // 3. Depth First Search Recursive
    // in every recursive step, set current vertex to visited and recurse on each adjacent unvisited vertex
    // the recursive structure naturally forms a stack
    private static void DFS_recurse(int v, boolean[] visited) {
        visited[v] = true;
        System.out.println(v);
        for (int i : edges.get(v)) {
            if (!visited[i])
                DFS_recurse(i, visited);
        }
    }
    
    public static void DFS_recursive(int v) {
        boolean[] visited = new boolean[nV];
        DFS_recurse(v, visited);
    }
    
    // 4. BFS Cycle Detection
    // count the number of incoming edges on each vertex, 0 means the vertex is a root
    // add roots in queue, remove one at a time and detach that root from the graph by removing 1 incoming edge from all children, then add all new roots in queue
    // when queue is empty, if the number of roots detached is less than nV, it means some vertices will never become roots because there is a cycle
    public static boolean BFS_isCyclic() {
        int[] incoming_edges = new int[nV];
        Queue<Integer> queue = new LinkedList<>();
        int count = 0;
        for (ArrayList<Integer> i : edges) {
            for (int j : i)
                incoming_edges[j]++;
        }
        for (int i = 0; i < nV; i++) {
            if (incoming_edges[i] == 0)
                queue.add(i);
        }
        while (!queue.isEmpty()) {
            int cur = queue.remove();
            for (int i : edges.get(cur)) {
                incoming_edges[i]--;
                if (incoming_edges[i] == 0)
                    queue.add(i);
            }
            count++;
        }
        return count != nV;
    }
    
    // 5. DFS Cycle Detection (recursive)
    // for each unvisited vertex, set it as possibly part of a cycle, then recurse
    // if it comes out of the recursion, it is not part of a cycle
    // if it meets itself again in the recursion, it is part of a cycle
    public static boolean DFS_isCyclicHelper(boolean[] visited, boolean[] possible, int v) {
        if (possible[v])
            return true;
        if (visited[v])
            return false;
        visited[v] = true;
        possible[v] = true;
        for (int i : edges.get(v)) {
            if (DFS_isCyclicHelper(visited, possible, i))
                return true;
        }
        possible[v] = false;
        return false;
    }
    
    public static boolean DFS_isCyclic() {
        boolean[] visited = new boolean[nV];
        boolean[] possible = new boolean[nV];
        for (int i = 0; i < nV; i++) {
            if (DFS_isCyclicHelper(visited, possible, i))
                return true;
        }
        return false;
    }
    
    // 6. Cycle Detection in Undirected Graph (recursive)
    // for each unvisited vertex, run recursive step on all adjacent unvisited vertices, and check if there exists an adjacent visited vertex that is not the parent
    // if such vertex exist, it can be visited by two different ways and thus is part of a cycle
    private static boolean isCyclicUndirectedHelper(boolean[] visited, int cur, int parent) {
        visited[cur] = true;
        for (int i : edges.get(cur)) {
            if (!visited[i]) {
                if (isCyclicUndirectedHelper(visited, i, cur))
                    return true;
            }
            else if (parent != i)
                return true;
        }
        return false;
    }
    
    public static boolean isCyclicUndirected() {
        boolean[] visited = new boolean[nV];
        for (int i = 0; i < nV; i++) {
            if (!visited[i]) {
                if (isCyclicUndirectedHelper(visited, i, -1))
                    return true;
            }
        }
        return false;
    }
    
    // 7. Check if graph is bipartite with BFS
    // add the first vertex in queue
    // remove one vertex from queue at a time until empty, check if adjacent vertices are of opposite colour
    // same colour means not bipartite
    // if uncoloured, colour it opposite colour and add in queue
    public static boolean isBipartite() {
        Queue<Integer> queue = new LinkedList<>();
        int[] colour = new int[nV];
        colour[0] = 1;
        queue.add(0);
        while (!queue.isEmpty()) {
            int cur = queue.remove();
            for (int i : edges.get(cur)) {
                if (cur == i)
                    return false;
                if (colour[cur] == colour[i])
                    return false;
                if (colour[i] == 0) {
                    colour[i] = 3 - colour[cur];
                    queue.add(i);
                }
            }
        }
        return true;
    }
    
    // 8. BFS find shortest path
    // same as BFS except stops when reaching dest, and records parent to recover shortest path
    private static boolean BFS_shortest_path_helper(int[] prev, int v, int w) {
        boolean[] visited = new boolean[nV];
        Queue<Integer> queue = new LinkedList<>();
        for (int i = 0; i < nV; i++)
            prev[i] = -1;
        visited[v] = true;
        queue.add(v);
        
        while (!queue.isEmpty()) {
            int cur = queue.remove();
            if (cur == w)
                return true;

            for (int i : edges.get(cur)) {
                if (!visited[i]) {
                    visited[i] = true;
                    queue.add(i);
                    prev[i] = cur;
                }
            }
        }
        
        return false;
    }
    
    public static void BFS_shortest_path(int v, int w) {
        int[] prev = new int[nV];
        
        if (!BFS_shortest_path_helper(prev, v, w)) {
            System.out.println("No path");
            return;
        }
        
        int cur = w;
        LinkedList<Integer> path = new LinkedList<>();
        path.add(cur);
        
        while (cur != v) {
            cur = prev[cur];
            path.add(cur);
        }
        System.out.println("Shortest path length is " + (path.size() - 1));
        for (int i = path.size() - 1; i >= 0; i--)
            System.out.print(path.get(i) + " ");
    }
    
    // 9. Bidirectional Search (BFS starting from source and dest meeting in the middle)
    // maintain 2 visited arrays and queues
    
    private static void processOneFromQueue(int[] prev, Queue<Integer> queue, boolean[] visited) {
        int cur = queue.remove();
        for (int i : edges.get(cur)) {
            if (!visited[i]) {
                visited[i] = true;
                queue.add(i);
                prev[i] = cur;
            }
        }
    }
    
    public static void bidirectionalSearch(int v, int w) {
        boolean[] v_visited = new boolean[nV];
        boolean[] w_visited = new boolean[nV];
        Queue<Integer> v_queue = new LinkedList<>();
        Queue<Integer> w_queue = new LinkedList<>();
        int[] v_prev = new int[nV];
        int[] w_prev = new int[nV];
        v_visited[v] = true;
        v_queue.add(v);
        w_visited[w] = true;
        w_queue.add(w);
        while (!v_queue.isEmpty() && !w_queue.isEmpty()) {
            processOneFromQueue(v_prev, v_queue, v_visited);
            processOneFromQueue(w_prev, w_queue, w_visited);
            for (int i = 0; i < nV; i++) {
                if (v_visited[i] && w_visited[i]) {
                    LinkedList<Integer> path = new LinkedList<>();
                    int cur = i;
                    path.add(cur);
                    while (cur != v) {
                        cur = v_prev[cur];
                        path.add(cur);
                    }
                    Collections.reverse(path);
                    cur = i;
                    while (cur != w) {
                        cur = w_prev[cur];
                        path.add(cur);
                    }
                    System.out.println("Shortest path length is " + (path.size() - 1));
                    for (int j = 0; j < path.size(); j++)
                        System.out.print(path.get(j) + " ");
                    return;
                }
            }
        }
        System.out.println("No path");
    }
    
    // 10. Kruskal's Algorithm (find minimum spanning tree when the graph is sparse)
    // in order to sort the edges, store edges in array of edges representation instead
    static class Edge implements Comparable<Edge> {
        int v;
        int w;
        int weight;
        public int compareTo(Edge e) {
            return this.weight - e.weight;
        }
    }
    
    static class Subset {
        int parent;
        int rank;
    }
    
    private static int find_parent(Subset[] subsets, int v) {
        if (subsets[v].parent == v)
            return v;
        return find_parent(subsets, subsets[v].parent);
    }
    
    private static void unite(Subset[] subsets, int v, int w) {
        int vp = find_parent(subsets, v);
        int wp = find_parent(subsets, w);
        if (subsets[vp].rank > subsets[wp].rank)
            subsets[wp].parent = vp;
        else if (subsets[vp].rank < subsets[wp].rank)
            subsets[vp].parent = wp;
        else {
            subsets[vp].parent = wp;
            subsets[wp].rank++;
        }
    }
    
    public static int Kruskal_MST() {
        int n_edges = 0;
        for (int i = 0; i < edges.size(); i++)
            n_edges += edges.get(i).size();
        n_edges /= 2;
        Edge[] edge_array = new Edge[n_edges];
        int k = 0;
        for (int i = 0; i < edges.size(); i++) {
            for (int j = 0; j < edges.get(i).size(); j++) {
                if (i <= edges.get(i).get(j)) {
                    edge_array[k] = new Edge();
                    edge_array[k].v = i;
                    edge_array[k].w = edges.get(i).get(j);
                    edge_array[k].weight = weights.get(i).get(j);
                    k++;
                }
            }
        }
        Arrays.sort(edge_array);
        
        Subset[] subsets = new Subset[nV];
        for (int i = 0; i < nV; i++) {
            subsets[i] = new Subset();
            subsets[i].parent = i;
            subsets[i].rank = 0;
        }
        
        int chosen_edges = 0;
        int min_sum = 0;
        int cur = 0;
        while (chosen_edges < nV - 1) {
            int v = edge_array[cur].v;
            int w = edge_array[cur].w;
            int vp = find_parent(subsets, v);
            int wp = find_parent(subsets, w);
            if (vp != wp) {
                chosen_edges++;
                min_sum += edge_array[cur].weight;
                unite(subsets, vp, wp);
            }
            cur++;
        }
        return min_sum;
    }
    
    // 11. Prim's Algorithm (find minimum spanning tree when the graph is dense)
    // 
    public static int Prim_MST() {
        boolean[] in_tree = new boolean[nV];
        int[] length = new int[nV];
        for (int i = 0; i < nV; i++)
            length[i] = Integer.MAX_VALUE;
        length[0] = 0;
        
        int min_sum = 0;
        for (int i = 0; i < nV; i++) {
            int min = Integer.MAX_VALUE, min_idx = -1;
            for (int j = 0; j < nV; j++) {
                if (!in_tree[j] && length[j] < min) {
                    min = length[j];
                    min_idx = j;
                }
            }
            min_sum += min;
            in_tree[min_idx] = true;
            for (int k = 0; k < edges.get(min_idx).size(); k++) {
                int c = edges.get(min_idx).get(k);
                if (!in_tree[c])
                    length[c] = Math.min(length[c], weights.get(min_idx).get(k));
            }
        }
        return min_sum;
    }
    
    // 12. Dijkstra's Algorithm (complexity O(E + V^2), faster for dense graph)
    // BFS with edge weights (doesn't work on negative edge weights)
    // if a vertex is in tree, the dist array stores its shortest distance from the source
    // each time add the vertex with shortest distance from source that is not already in tree into the tree
    // then update all adjacent vertices of current vertex that are not in tree to the shortest distance, either original value or current value plus edge weight 
    public static int Dijkstra(int v, int w) {
        int[] dist = new int[nV];
        boolean[] in_tree = new boolean[nV];
        for (int i = 0; i < nV; i++)
            dist[i] = Integer.MAX_VALUE;
        dist[v] = 0;
        
        for (int i = 0; i < nV; i++) {
            int min = Integer.MAX_VALUE, min_idx = -1;
            for (int j = 0; j < nV; j++) {
                if (!in_tree[j] && dist[j] < min) {
                    min = dist[j];
                    min_idx = j;
                }
            }
            //System.out.println(min_idx);
            if (min_idx == -1)
                break;
            if (min_idx == w)
                return dist[w];
            in_tree[min_idx] = true;
            for (int k = 0; k < edges.get(min_idx).size(); k++) {
                int c = edges.get(min_idx).get(k);
                if (!in_tree[c])
                    dist[c] = Math.min(dist[c], min + weights.get(min_idx).get(k));
            }
        }
        return Integer.MAX_VALUE;
    }
    
    // 13. Dijkstra's Algorithm (complexity O(ElogV), faster for sparse graph)
    // change representation to an ArrayList of ArrayLists of VNodes: ArrayList i stores list of nodes, each node j containing an adjacent vertex to i and weight of edge (i, j) etc
    // initialise queue and distance from source
    // remove vertex from queue (with highest priority - shortest distance)
    // if distance stored in this vertex node is larger, discard it: this is faster than maintaining a set for settled nodes (add nodes in queue even if it's theoretically in set, and discard later here)
    // for each of that vertex's adjacent vertex where it's closer to come to from that vertex, reset distance and add in queue
    // loop V times, each time dequeue once and enqueue V times
    // thus complexity = O(V * (logV + VlogV)) = O(V^2 logV) = O(ElogV)
    static class VNode {
        int vertex;
        int distance;
        public VNode(int v, int d) {
            vertex = v;
            distance = d;
        };
    }
    
    static class compareVNode implements Comparator<VNode> {
        public int compare(VNode n1, VNode n2) {
            return n1.distance - n2.distance;
        }
    }
    
    public static int DijkstraPriorityQueue(int v, int w) {
        PriorityQueue<VNode> pq = new PriorityQueue<>(nV, new compareVNode());
        int[] dist = new int[nV];
        for (int i = 0; i < nV; i++)
            dist[i] = Integer.MAX_VALUE;
        dist[v] = 0;
        pq.add(new VNode(v, 0));
        
        ArrayList<ArrayList<VNode>> list_edges_nodes = new ArrayList<>();
        for (int i = 0; i < edges.size(); i++) {
            ArrayList<VNode> edge_nodes = new ArrayList<>();
            for (int j = 0; j < edges.get(i).size(); j++) {
                VNode node = new VNode(edges.get(i).get(j), weights.get(i).get(j));
                edge_nodes.add(node);
            }
            list_edges_nodes.add(edge_nodes);
        }
        
        while (!pq.isEmpty()) {
            VNode cur = pq.remove();
            int min_idx = cur.vertex, min_dist = cur.distance;
            if (min_dist > dist[min_idx])
                continue;
            if (min_idx == w)
                return min_dist;
            for (VNode i : list_edges_nodes.get(min_idx)) {
                int next = i.vertex, next_d = i.distance;
                if (dist[min_idx] + next_d < dist[next]) {
                    dist[next] = dist[min_idx] + next_d;
                    pq.add(new VNode(next, dist[next]));
                }
            }
        }
        return dist[w];
    }
    
    // 14. Bellman-Ford Algorithm
    // detects negative cycle and returns -1, otherwise return shortest distance
    // go through each edge updating the second vertex's distance, repeat nV - 1 times
    // if there is still a shorter path, then there is a negative weight cycle
    public static int BellmanFord(int v, int w) {
        int[] dist = new int[nV];
        for (int i = 0; i < nV; i++)
            dist[i] = Integer.MAX_VALUE;
        dist[v] = 0;
        
        for (int j = 1; j < nV; j++) {
            for (int v1 = 0; v1 < edges.size(); v1++) {
                for (int i = 0; i < edges.get(v1).size(); i++) {
                    int v2 = edges.get(v1).get(i);
                    int wei = weights.get(v1).get(i);
                    if (dist[v1] != Integer.MAX_VALUE && dist[v1] + wei < dist[v2])
                        dist[v2] = dist[v1] + wei;
                }
            }
        }
        
        for (int v1 = 0; v1 < edges.size(); v1++) {
            for (int i = 0; i < edges.get(v1).size(); i++) {
                int v2 = edges.get(v1).get(i);
                int wei = weights.get(v1).get(i);
                if (dist[v1] != Integer.MAX_VALUE && dist[v1] + wei < dist[v2])
                    return -1;
            }
        }
        
        return dist[w];
    }
    
    // 15. Hamiltonian cycle
    static int[] path;
    //check if v can be added to the constructed cycle of length n
    public static boolean canAdd(int v, int n) {
        if (!edges.get(path[n - 1]).contains(v))
            return false;
        for (int i = 0; i < n; i++) {
            if (path[i] == v)
                return false;
        }
        return true;
    }
    
    //n is the number of vertices in the cycle
    //recursively add one vertex at a time
    public static boolean isHamiltonianCycle(int n) {
        if (n == nV)
            return edges.get(path[n - 1]).contains(path[0]);
        for (int i = 1; i < nV; i++) {
            if (canAdd(i, n)) {
                path[n] = i;
                if (isHamiltonianCycle(n + 1))
                    return true;
                path[n] = -1;
            }
        }
        return false;
    }
    
    public static boolean hamiltonianCycle() {
        path = new int[nV];
        for (int i = 0; i < nV; i++)
            path[i] = -1;
        path[0] = 0;
        return isHamiltonianCycle(1);
    }
    
    // 16. Travelling Salesman Problem
    // there is an edge (weighted, undirected) between every pair of vertices
    // find the cost of the hamilton cycle with the lowest cost
    // at the start of each backtrack step, cur is the newly added vertex, count is the number of vertices in path, cost is the newly added weight (by adding cur), res is the minimum discovered up to now (before adding cur)
    // at the end of each backtrack step, res is updated to the minimum discovered after adding cur
    public static int TSP_backtrack(boolean[] visited, int cur, int count, int cost, int res) {
        //all vertices in path
        if (count == nV && edges.get(cur).contains(0)) {
            res = Math.min(res, cost + weights.get(cur).get(0));
            return res;
        }
        for (int j = 0; j < edges.get(cur).size(); j++) {
            int i = edges.get(cur).get(j);
            if (!visited[i]) {
                visited[i] = true;
                res = TSP_backtrack(visited, i, count + 1, cost + weights.get(cur).get(j), res);
                visited[i] = false;
            }
        }
        return res;
    }
    
    public static int TSP() {
        boolean[] visited = new boolean[nV];
        visited[0] = true;
        int res = Integer.MAX_VALUE;
        return TSP_backtrack(visited, 0, 1, 0, res);
    }
    
    public static void main(String[] args) {
        //example of an unweighted, undirected graph
        /*Graph g = new Graph(5, false);
        g.addEdge(0,1,false);
        g.addEdge(1,2,false);
        g.addEdge(0,3,false);
        g.addEdge(2,4,false);
        g.addEdge(1,3,false);
        g.addEdge(1,4,false);*/
        
        //example of an unweighted, directed graph
        /*Graph g = new Graph(4, false);
        g.addEdge(0,0,true);
        g.addEdge(1,2,true);
        g.addEdge(1,3,true);
        g.addEdge(2,3,true);
        g.addEdge(3,0,true);*/
        
        //example of a weighted, undirected graph
        Graph g = new Graph(4, true);
        g.addEdgeWeighted(0,0,0);
        g.addEdgeWeighted(0,1,10);
        g.addEdgeWeighted(0,2,15);
        g.addEdgeWeighted(0,3,20);
        g.addEdgeWeighted(1,0,10);
        g.addEdgeWeighted(1,1,0);
        g.addEdgeWeighted(1,2,35);
        g.addEdgeWeighted(1,3,25);
        g.addEdgeWeighted(2,0,15);
        g.addEdgeWeighted(2,1,35);
        g.addEdgeWeighted(2,2,0);
        g.addEdgeWeighted(2,3,30);
        g.addEdgeWeighted(3,0,20);
        g.addEdgeWeighted(3,1,25);
        g.addEdgeWeighted(3,2,30);
        g.addEdgeWeighted(3,3,0);
        
        System.out.println(TSP());
    }
}
