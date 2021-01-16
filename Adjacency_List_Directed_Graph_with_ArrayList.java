// 
import java.util.ArrayList;

public class Graph {
    public static int nV;
    public static ArrayList<ArrayList<Integer>> edges;
    
    Graph(int nV) {
        this.nV = nV;
        this.edges = new ArrayList<ArrayList<Integer>>(nV);
        for (int i = 0; i < nV; i++)
            edges.add(new ArrayList<Integer>());
    }
    
    public void addEdge(int v1, int v2) {
        if (!edges.get(v1).contains(v2))
            edges.get(v1).add(v2);
    }
    
    public void print() {
        for (int i = 0; i < nV; i++) {
            System.out.print("vertex " + i + ":");
            for (int a : edges.get(i))
                System.out.print(" " + a);
            System.out.println();
        }
    }
    
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
    
    public static void main(String[] args) {
        Graph g = new Graph(7);
        g.addEdge(0,1);
        g.addEdge(0,2);
        g.addEdge(1,3);
        g.addEdge(2,4);
        g.addEdge(3,4);
        g.addEdge(3,6);
        //g.print();
        DFS_iterative(0);
    }
}
