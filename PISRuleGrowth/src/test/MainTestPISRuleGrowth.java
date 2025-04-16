package test;

import ca.pfv.spmf.algorithms.sequential_rules.pisrulegrowth.AlgoPISRuleGrowth;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;

/**
 * Example of how to use the PISRuleGrowth Algorithm in source code.
 *
 * @author Nevena Milenkovic (Copyright 2025)
 */
public class MainTestPISRuleGrowth {

    public static void main(String[] arg) throws IOException {

        String input = fileToPath("contextPrefixSpan.txt");
        
        String output = ".//output.txt";  // the path for saving the frequent periodic sequential rules

        //  Applying SeqRuleGrowth algorithm with parameters
        double minconf = 0.2;
        double maxstd = 5;
        double minra = 0.01;

        double minsup_absolute = 0.05;  // it means 75 %
        AlgoPISRuleGrowth algo = new AlgoPISRuleGrowth();
        algo.runAlgorithm(input, output, minsup_absolute, minconf, maxstd, minra);

        // print statistics
        algo.printStats();
    }

    public static String fileToPath(String filename) throws UnsupportedEncodingException {
        URL url = MainTestPISRuleGrowth.class.getResource(filename);
        return java.net.URLDecoder.decode(url.getPath(), "UTF-8");
    }
}
