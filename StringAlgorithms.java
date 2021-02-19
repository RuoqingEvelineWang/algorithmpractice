//string algorithms
//1. distinct permutations of a string
//2. longest common subsequence of two strings
//3. edit distance
//4. KMP algorithm
//5. Rabin-Karp algorithm
//6. Finite Automata algorithm
//7. Boyer Moore algorithm
//8. if a string is a rotation of another string
//9. shortest substring containing all characters of another string
//10. longest palindromic substring

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
    //find the longest proper prefix that is also a suffix of the pattern
    public static void find_lps(String pattern, int m, int[] lps) {
        int k = 0;
        //j: upto which number is the loop currently calculating lps for
        //k: how many characters in sub pattern is matched
        for (int j = 1; j < m; j++) {
            while (k > 0 && pattern.charAt(k) != pattern.charAt(j))
                k = lps[k - 1];
            if (pattern.charAt(k) == pattern.charAt(j))
                k++;
            lps[j] = k;
        }
    }
    
    public static void KMP(String pattern, String text) {
        int m = pattern.length(), n = text.length();
        int[] lps = new int[m];
        find_lps(pattern, m, lps);
        int j = 0;
        //i: upto which character in text to compare
        //j: how many characters in pattern is matched
        for (int i = 0; i < n; i++) {
            while (j > 0 && pattern.charAt(j) != text.charAt(i))
                j = lps[j - 1];
            if (pattern.charAt(j) == text.charAt(i))
                j++;
            if (j == m) {
                System.out.println(i - j + 1);
                j = lps[j - 1];
            }
        }
    }
    
    //5. Rabin-Karp algorithm: calculate hash value and compare
    public static int d = 26; //number of characters in the alphabet of string
    public static int q = 101; //a prime number
    
    public static void RabinKarp(String pattern, String text) {
        int m = pattern.length(), n = text.length();
        int p = 0, t = 0, h = 1;
        
        for (int i = 0; i < m - 1; i++)
            h = (h * d) % q;
        
        for (int i = 0; i < m; i++) {
            p = (d * p + pattern.charAt(i)) % q;
            t = (d * t + text.charAt(i)) % q;
        }
        
        for (int i = 0; i <= n - m; i++) {
            if (p == t) {
                int j;
                for (j = 0; j < m; j++) {
                    if (pattern.charAt(j) != text.charAt(i + j))
                        break;
                }
                if (j == m)
                    System.out.println(i);
            }
            if (i < n - m) {
                t = (d * (t - text.charAt(i) * h) + (text.charAt(i + m))) % q;
                if (t < 0)
                    t += q;
            }
        }
    }
    
    //6. Finite Automata algorithm
    //meaning of the tf 2d array:
    //if (state) number of characters in the pattern are matched right now
    //how many characters will be matched after adding (x) at the end
    public static int getState(String pattern, int m, int state, int x) {
        if (state < m && x == pattern.charAt(state))
            return state + 1;
        for (int j = state; j > 0; j--) {
            if (pattern.charAt(j - 1) == x) {
                int i;
                for (i = 0; i < j - 1; i++) {
                    if (pattern.charAt(i) != pattern.charAt(state - j + 1 + i))
                        break;
                }
                if (i == j - 1)
                    return j;
            }
        }
        return 0;
    }
    
    public static void FiniteAutomata(String pattern, String text) {
        int m = pattern.length(), n = text.length();
        int[][] tf = new int[m + 1][256];
        
        for (int state = 0; state <= m; state++) {
            for (int x = 0; x < 256; x++) {
                tf[state][x] = getState(pattern, m, state, x);
            }
        }
        //state is how many characters match the pattern
        int state = 0;
        for (int i = 0; i < n; i++) {
            state = tf[state][text.charAt(i)];
            if (state == m)
                System.out.println(i - m + 1);
        }
    }
    
    //7. Boyer Moore algorithm
    public static void BoyerMoore(String pattern, String text) {
        int m = pattern.length(), n = text.length();
        //stores position of last occurrence of each character in the pattern string
        int last[] = new int[256];
        for (int i = 0; i < 256; i++)
            last[i] = -1;
        for (int i = 0; i < m; i++)
            last[pattern.charAt(i)] = i;
        //s is which character is being looked at (as start of pattern) in the text
        //once s is set, start comparing from the last character of the pattern
        int s = 0;
        while (s <= n - m) {
            int j = m - 1;
            while (j >= 0 && pattern.charAt(j) == text.charAt(s + j))
                j--;
            if (j < 0) {
                System.out.println(s);
                s += (s < n - m) ? m - last[text.charAt(s + m)] : 1;
            }
            else
                s += Math.max(1, j - last[text.charAt(s + j)]);
        }
    }
    
    //8. find if a string is a rotation of another string using lps from KMP
    public static boolean isRotation(String a, String b) {
        int m = a.length();
        if (m != b.length())
            return false;
        int[] lps = new int[m];
        int k = 0;
        //k is how many characters are matched for prefix/suffix
        //j is up to which character is the loop calculating lps for
        for (int j = 1; j < m; j++) {
            while (k > 0 && a.charAt(j) != b.charAt(k))
                k = lps[k - 1];
            if (a.charAt(j) == b.charAt(k))
                k++;
            lps[j] = k;
        }
        for (int i = lps[m - 1], j = 0; i < m; i++, j++) {
            if (a.charAt(j) != b.charAt(i))
                return false;
        }
        return true;
    }
    
    //9. find the shortest substring of a string that contains all characters of another string
    public static String shortestSubstring(String text, String pattern) {
        int n = text.length(), m = pattern.length();
        if (n < m)
            return "";
        int[] pattern_hash = new int[256];
        int[] text_hash = new int[256];
        for (int i = 0; i < m; i++)
            pattern_hash[pattern.charAt(i)]++;
        int count = 0, start = 0, min_idx = -1, min_len = Integer.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            int index = text.charAt(i);
            text_hash[index]++;
            //this is when the ith character of text is not extra, it is useful as part of pattern
            if (text_hash[index] <= pattern_hash[index])
                count++;
            if (count == m) {
                while (text_hash[text.charAt(start)] > pattern_hash[text.charAt(start)] || pattern_hash[text.charAt(start)] == 0) {
                    if (text_hash[text.charAt(start)] > pattern_hash[text.charAt(start)])
                        text_hash[text.charAt(start)]--;
                    start++;
                }
                if (i - start + 1 < min_len) {
                    min_len = i - start + 1;
                    min_idx = start;
                }
            }
        }
        if (min_idx == -1)
            return "";
        return text.substring(min_idx, min_idx + min_len);
    }
    
    //10. find the longest substring of a string that is a palindrome
    public static int longestPalindromeFromCentre(String s, int l, int c1, int c2) {
        int i = c1, j = c2;
        while (i >= 0 && j < l && s.charAt(i) == s.charAt(j)) {
            i--;
            j++;
        }
        return j - i - 1;
    }
    
    public static String longestPalindromicSubstring(String s) {
        int l = s.length(), start = 0, end = 0;
        for (int i = 0; i < l; i++) {
            int m1 = longestPalindromeFromCentre(s, l, i, i);
            int m2 = longestPalindromeFromCentre(s, l, i, i + 1);
            int m = Math.max(m1, m2);
            if (end - start + 1 < m) {
                start = i - (m - 1) / 2;
                end = i + m / 2;
            }
        }
        return s.substring(start, end + 1);
    }
    
    public static void main(String[] args) {
        //permutation("aab", "");
        //lcs("abcabcaa", "acbacba");
        //System.out.println(edit_distance("sunday", "saturday", 6, 8));
        //KMP("ABABCABAB", "ABABDABACDABABCABABABABCABABCABAB");
        //RabinKarp("ABABCABAB", "ABABDABACDABABCABABABABCABABCABAB");
        //System.out.println(isRotation("abcdefghijklmnopqz", "klmnopqyabcdefghij"));
        //FiniteAutomata("acacaga", "acacagacacaga");
        //BoyerMoore("acacaga", "acacagacacaga");
        //System.out.println(shortestSubstring("akjfkjeajkfhiewfehifejfahijfehfcccaddfdvdvvdccdb", "abccddd"));
        //System.out.println(longestPalindromicSubstring("aaabbbccccdccccbbbvvv"));
    }
}
