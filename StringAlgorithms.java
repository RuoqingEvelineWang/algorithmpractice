//string algorithms
//1. distinct permutations of a string
//2. longest common subsequence of two strings
//3. edit distance
//4. KMP algorithm

public class Main {
    //1. distinct permutations of a string
    public static void permutation(String str, String res) {
        if (str.length() == 0)
            System.out.print(res + "\n");
        boolean[] used = new boolean[26];
        for (int i = 0; i < str.length(); i++) {
            char c = str.charAt(i);
            if (!used[c - 'a']) {
                permutation(str.substring(0, i) + str.substring(i + 1), res + c);
                used[c - 'a'] = true;
            }
        }
    }
    
    //2. find all longest common subsequences of two strings
    static int[][] results;
    static int length;
    
    //can also be done iteratively
    public static int recursive_lcs(String s1, String s2, int l1, int l2, int i, int j) {
        if (i == l1 || j == l2)
            return 0;
        if (results[i][j] != -1)
            return results[i][j];
        int r = 0;
        if (s1.charAt(i) == s2.charAt(j))
            r = recursive_lcs(s1, s2, l1, l2, i + 1, j + 1) + 1;
        else
            r = Math.max(recursive_lcs(s1, s2, l1, l2, i, j + 1), recursive_lcs(s1, s2, l1, l2, i + 1, j));
        results[i][j] = r;
        return r;
    }
    
    public static void print_lcs(String s1, String s2, int l1, int l2, int i, int j, int l, char[] lcs) {
        if (l == length) {
            lcs[length] = '\0';
            System.out.println(new String(lcs));
            return;
        }
        
        if (i == l1 || j == l2)
            return;
        
        for (char c = 'a'; c <= 'z'; c++) {
            boolean found = false;
            for (int k1 = i; k1 < l1; k1++) {
                if (s1.charAt(k1) == c) {
                    for (int k2 = j; k2 < l2; k2++) {
                        if (s2.charAt(k2) == c && results[k1][k2] == length - l) {
                            lcs[l] = c;
                            print_lcs(s1, s2, l1, l2, k1 + 1, k2 + 1, l + 1, lcs);
                            found = true;
                            break;
                        }
                    }
                }
                if (found)
                    break;
            }
        }
    }
    
    public static void lcs(String s1, String s2) {
        int l1 = s1.length(), l2 = s2.length();
        results = new int[l1][l2];
        for (int i = 0; i < l1; i++) {
            for (int j = 0; j < l2; j++)
                results[i][j] = -1;
        }
        length = recursive_lcs(s1, s2, l1, l2, 0, 0);
        print_lcs(s1, s2, l1, l2, 0, 0, 0, new char[length + 1]);
    }
    
    //3. find the edit distance (insert, remove, replace) between two strings
    public static int edit_distance(String s1, String s2, int l1, int l2) {
        if (l1 == 0)
            return l2;
        if (l2 == 0)
            return l1;
        if (s1.charAt(l1 - 1) == s2.charAt(l2 - 1))
            return edit_distance(s1, s2, l1 - 1, l2 - 1);
        return 1 + Math.min(edit_distance(s1, s2, l1 - 1, l2 - 1), Math.min(edit_distance(s1, s2, l1, l2 - 1), edit_distance(s1, s2, l1 - 1, l2)));
    }
    
    //4. KMP algorithm for pattern searching
    
    
    public static void main(String[] args) {
        //permutation("aab", "");
        //lcs("abcabcaa", "acbacba");
        //System.out.println(edit_distance("sunday", "saturday", 6, 8));
    }
}
