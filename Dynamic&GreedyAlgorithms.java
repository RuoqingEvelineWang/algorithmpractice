//dynamic programming and greedy problems
//1. egg dropping problem
//2. interval scheduling problem
//3. knapsack problem

public class Main {
    //1. egg dropping problem
    //n is the number of eggs and k is the number of floors, return the minimum number of drops to determine lowest floor from which egg will break
    public static int eggDropping(int n, int k) {
        int[][] save = new int[n + 1][k + 1];
        for (int i = 1; i <= n; i++) {
            save[i][0] = 0;
            save[i][1] = 1;
        }
        for (int i = 1; i <= k; i++)
            save[1][i] = i;
        
        for (int i = 2; i <= n; i++) {
            for (int j = 2; j <= k; j++) {
                save[i][j] = Integer.MAX_VALUE;
                for (int x = 1; x <= j; x++) {
                    int res = 1 + Math.max(save[i - 1][x - 1], save[i][j - x]);
                    if (res < save[i][j])
                        save[i][j] = res;
                }
            }
        }
        return save[n][k];
    }
    
    //2. interval scheduling problem (assume finish time is sorted ascendingly)
    //returns the maximum number of intervals not overlapping with each other
    public static int intervalScheduling(int[] start, int[] finish) {
        if (start.length <= 0)
            return 0;
        int count = 1, i = 0;
        for (int j = 1; j < start.length; j++) {
            if (start[j] >= finish[i]) {
                count++;
                i = j;
            }
        }
        return count;
    }
    
    //3. knapsack problem
    //n is number of items, c is weight capacity of knapsack
    public static int knapsack(int[] values, int[] weights, int c) {
        int n = values.length;
        int save[][] = new int[n + 1][c + 1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= c; j++) {
                if (i == 0 || j == 0)
                    save[i][j] = 0;
                else if (weights[i - 1] <= j)
                    save[i][j] = Math.max(values[i - 1] + save[i - 1][j - weights[i - 1]], save[i - 1][j]);
                else
                    save[i][j] = save[i - 1][j];
            }
        }
        return save[n][c];
    }
    
    public static void main(String[] args) {
        //System.out.println(eggDropping(2, 36));
        /*int[] s = {1,3,0,5,8,5};
        int[] e = {2,4,6,7,9,9};
        System.out.println(intervalScheduling(s, e));*/
        int[] v = {60, 100, 120};
        int[] weight = {10, 20, 30};
        System.out.println(knapsack(v, weight, 50));
    }
}
