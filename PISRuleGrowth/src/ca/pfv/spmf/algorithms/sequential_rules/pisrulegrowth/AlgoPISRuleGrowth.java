package ca.pfv.spmf.algorithms.sequential_rules.pisrulegrowth;
/**
 *This code is © N. Milenkovic, 2025, and it is made available under the GPL 
 * license enclosed with the software. 
 */

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import ca.pfv.spmf.input.sequence_database_list_integers.Sequence;
import ca.pfv.spmf.input.sequence_database_list_integers.SequenceDatabase;
import ca.pfv.spmf.tools.MemoryLogger;
import java.util.Collections;

/**

/**
 * This is the original implementation of the PISRULEGROWTH algorithm for mining
 * periodic sequential rules common to several sequences where
 * antecedent and consequent are unordered itemsets. The PISRuleGrowth algorithm
 * is described in this paper:
 * <br/><br/>
 * Milenkovic,N. (2025). PISRuleGrowth: mining periodic intra-sequential rules 
 * common to multiple sequences. 
 * Knowledge and Information System, Springer
 * <br/><br/>
 * The main method of this algorithm is "runAlgorithm". It output the result to
 * a file.
 *
 * @see PISOccurence
 * @author Nevena Milenkovic
 */
public class AlgoPISRuleGrowth {
    //*** for statistics ***/

    /**
     * start time of latest execution
     */
    long timeStart = 0;

    /**
     * end time of latest execution
     */
    long timeEnd = 0;

    /**
     * number of rules generated
     */
    int ruleCount;

    //*** parameters ***/
    /**
     * minimum confidence
     */
    double minConf;

    /**
     * minimum support
     */
    double minSupp;

    /**
     * maximum standard deviation
     */
    double maxStd;

    /**
     * the minRa threshold
     */
    double minRa;

    /**
     * this is the sequence database
     */
    SequenceDatabase database;

    /**
     * * internal variables // This map contains for each item (key) a map of
     * occurences (value). // The map of occurences associates to sequence ID
     * (key), an occurence of the item (value).
     */
    Map<Integer, Map<Integer, PISOccurence>> mapItemsOccur;  // item, <tid, occurence>

    /**
     * object to write the output file
     */
    BufferedWriter writer = null;
    
    // number of sequences in the database
    int databaseSize;
    
    // lenght of each sequence in the database
    int[] seqLenghts;

    /**
     * Default constructor
     */
    public AlgoPISRuleGrowth() {
    }

    /**
     * The main method to run the algorithm
     *
     * @param input : an input file path of a sequence database
     * @param output : a file path for writing the output file containing the
     * seq. rules.
     * @param minSupport : the minimum support (percentage as a double value)
     * @param minConfidence : the minimum confidence threshold
     * @param maxStd : the maximum standard deviation of periods threshold
     * @param minRa : the minimum Ra threshold
     * @exception IOException if error reading/writing files
     */
    public void runAlgorithm(String input, String output, double minsupp, double minconf, double maxstd, double minra) throws IOException {

        try {
            // read the input database
            database = new SequenceDatabase();
            database.loadFile(input);
        } catch (Exception e) {
            e.printStackTrace();
        }

        // save the minimum support parameter
        this.minSupp = minsupp;
        // save the minimum confidence parameter
        this.minConf = minconf;
        // save the maximum standard deviation parameter
        this.maxStd = maxstd;
        // save the minimum Ra parameter
        this.minRa = minra;
        // reinitialize the number of rules found
        ruleCount = 0;

        // reset the stats for memory usage
        MemoryLogger.getInstance().reset();

        // prepare the object for writing the output file
        writer = new BufferedWriter(new FileWriter(output));

        // save the start time
        timeStart = System.currentTimeMillis(); // for stats

        // Remove infrequent items and count the support of each item 
        // and their occurrences from the database in one database scan.
        getOccurMapAndRemoveInfrequent(database);

        // save the number of sequences in the database
        databaseSize= database.size();
        //System.out.println("DATABASE SIZE "+databaseSize);
        
        // save the lenght of each sequence in the database
        seqLenghts = new int[database.size()];
        int suma=0;
        for(int i=0;i<database.size();i++){
        seqLenghts[i]=database.getSequences().get(i).size();
        suma+=database.getSequences().get(i).size();
        }
        double s=suma/database.size();
        //System.out.println("AVERAGE SEQUENCE LENGTH "+s);
        
        // after previous method, we don't need a reference to the database anymore.
        database = null;
         
        // We will now try to generate rules with one item in the
        // antecedent and one item in the consequent using
        // the frequent items.

        // For each pair of frequent items i  and j such that i != j
        for (int i : mapItemsOccur.keySet()) {
            // get the item I and its map of occurences
            Map<Integer, PISOccurence> occurencesI = mapItemsOccur.get(i);
            
            for (int j : mapItemsOccur.keySet()) {
                if((j==i) || (i>j)){continue;}
                // get the item j and its map of occurences
                Map<Integer, PISOccurence> occurencesJ = mapItemsOccur.get(j);
                //System.out.println(" run i "+i);
                //System.out.println(" run j "+j);
                // (1) We will now calculate the sidsets
                // of I -->J  and the rule J-->I.
                
                // create map of parameter values of itemset by sequence : 
                // key= sids of itemset, value = parameter values 
                // (support, confidence and standard deviation)
                // used when saving the rule
                Map<Integer, Double[]> parValIJ = new HashMap<Integer, Double[]>();
                Map<Integer, Double[]> parValJI = new HashMap<Integer, Double[]>();
                
                // Set of sequences of the rule that meet the support threshold,
                // but may not meet thresholds of confidence and standard deviation
                // so they are saved for expansion only
                Set<Integer> seqForExpIJ= new HashSet<Integer>();
                Set<Integer> seqForExpJI= new HashSet<Integer>();
                
                // for each sequence of an itemset
                for (Entry<Integer, PISOccurence> entryOccI : occurencesI.entrySet()) {
                    int lenghtOfSeq = seqLenghts[entryOccI.getKey()];
                    // get the occurence of J in the same sequence
                    PISOccurence occJ = occurencesJ.get(entryOccI.getKey());
                    // if J appears in that sequence
                    if (occJ != null) {
                        // if J appeared before I in that sequence,
                        // then we put this tid in the sidset of  J-->I
                        
                        // get occurence of single item rule (JI) in sequence
                        PISOccurence[] occurenceJI = getSingleItemRuleOccur(occJ.firstItemset, entryOccI.getValue().firstItemset, lenghtOfSeq);

                        if (occurenceJI != null) {
                            // calculate parameter values of rule JI in that sequence
                            Double[] parVal = isPisRule(occurenceJI[0], occurenceJI[1], occJ.firstItemset.size(), lenghtOfSeq);
                            if (parVal[0] != null) {
                                // add sequence of the rule to list for expansion
                                seqForExpJI.add(entryOccI.getKey());
                            }
                            if((parVal[0]!=null) && (parVal[1]!=null) && (parVal[2]!=null)){
                            // add values of frequent and periodic sequential rule
                                parValJI.put(entryOccI.getKey(), parVal);
                            }
                        }
                        // if I appeared before J in that sequence,
                        // then we put this sid in the tidset of  I-->J
                        
                        // get occurence of single item rule (IJ) in sequence
                        PISOccurence[] occurenceIJ = getSingleItemRuleOccur(entryOccI.getValue().firstItemset, occJ.firstItemset, lenghtOfSeq);
                        
                        if (occurenceIJ != null) {
                            // calculate parameter values of rule IJ in that sequence
                            Double[] parVal = isPisRule(occurenceIJ[0], occurenceIJ[1], entryOccI.getValue().firstItemset.size(), lenghtOfSeq);
                            if (parVal[0] != null) {
                                // add sequence of the rule to list for expansion
                                seqForExpIJ.add(entryOccI.getKey());
                            }
                            if((parVal[0]!=null) && (parVal[1]!=null) && (parVal[2]!=null)){
                            // add values of frequent and periodic sequential rule
                                parValIJ.put(entryOccI.getKey(), parVal);
                            }
                        }
                    }
                }

                // (2) check if the two itemsets have enough common sids
                // if not, we don't need to generate a rule for them.
                
                // create rule IJ
                
                // chech if rule IJ meets minRa threshold
                double raIJ = (double) parValIJ.keySet().size() / databaseSize;
                if ((!seqForExpIJ.isEmpty()) && (raIJ > minRa)) {
                    // create itemsets of the rule I ==> J
                    List<Integer> itemsetI = new ArrayList<Integer>();
                    itemsetI.add(i);
                    List<Integer> itemsetJ = new ArrayList<Integer>();
                    itemsetJ.add(j);
                    
                    // if rule meets the parameter thresholds
                     if(!parValIJ.isEmpty()){
                    // save the rule in a file
                    saveRule(raIJ, itemsetI, itemsetJ, parValIJ);}
                    
                    // recursive call to try to expand the rule on the left and
                    // right sides
                    
                    expandLeft(itemsetI, itemsetJ, seqForExpIJ);
                    
                    // create map of support of itemsetI in all sequences for expandRight
                    Map<Integer, Integer> occiCount = new HashMap<Integer, Integer>();
                    for (int c : occurencesI.keySet()) {
                        occiCount.put(c, occurencesI.get(c).firstItemset.size());
                    }
                    
                    expandRight(itemsetI, itemsetJ, seqForExpIJ, occiCount);
                  
                }

                // create rule IJ
                
                // chech if rule JI meets minRa threshold
                double raJI = (double) parValJI.keySet().size() / databaseSize;
                if ((!seqForExpJI.isEmpty()) && (raJI > minRa)) {
                    // create itemset of the rule J ==> I
                    List<Integer> itemsetI = new ArrayList<Integer>();
                    itemsetI.add(i);
                    List<Integer> itemsetJ = new ArrayList<Integer>();
                    itemsetJ.add(j);
                    
                    // if rule meets the parameter thresholds
                     if(!parValJI.isEmpty()){
                    // save the rule in a file, if all parameters are met
                    saveRule(raJI, itemsetJ, itemsetI, parValJI);}
                    
                    // recursive call to try to expand the rule on the left and
                    // right sides
                    
                    // create map of support of itemsetI in all sequences for expandRight
                    Map<Integer, Integer> occjCount = new HashMap<Integer, Integer>();
                    for (int c : occurencesJ.keySet()) {
                        occjCount.put(c, occurencesJ.get(c).firstItemset.size());
                    }

                    expandRight(itemsetJ, itemsetI, seqForExpJI, occjCount);

                    expandLeft(itemsetJ, itemsetI, seqForExpJI); 
                }
            }
        }

        // save end time
        timeEnd = System.currentTimeMillis();

        // close the file
        writer.close();

        
    }

    /**
     * Save a rule I ==> J to the output file
     *
     * @param ra         value of ra threshold of the rule IJ
     * @param itemsetI   items contained in the left part of the rule
     * @param itemsetJ   items contained in the right part of the rule
     * @param parSeqVal  a map where key = sid, value = array of parameter 
     *                   values of the rule IJ
     * @throws IOException exception if error writing the file
     */
    private void saveRule(double ra, List<Integer> itemsetI, List<Integer> itemsetJ, Map<Integer, Double[]> parSeqVal) throws IOException {
        // increase the number of rule found
        ruleCount++;
		
        //System.out.println("RC "+ruleCount);
        // create a string buffer
        StringBuilder buffer = new StringBuilder();

        // write itemset 1 (antecedent)
        for (int i = 0; i < itemsetI.size(); i++) {
            buffer.append(itemsetI.get(i));

            if (i != itemsetI.size() - 1) {
                buffer.append(",");
            }
        }
        
        // write separator
        buffer.append(" ==> ");

        // write itemset 2  (consequent)
        for (int i = 0; i < itemsetJ.size(); i++) {
            buffer.append(itemsetJ.get(i));

            if (i != itemsetJ.size() - 1) {
                buffer.append(",");
            }
        }
        // write ra value of the rule IJ
        buffer.append(" #RA: ");
        buffer.append(ra);

        // write support, confidence and standard deviation of the rule IJ
        // in each sequence
        for (Entry<Integer, Double[]> values : parSeqVal.entrySet()) {
            buffer.append("  #Seq:");
            buffer.append(values.getKey());
            buffer.append("  #SUP:");
            buffer.append(Math.floor(values.getValue()[0] * 100) / 100);
            buffer.append("  #CONF:");
            buffer.append(Math.floor(values.getValue()[1] * 100) / 100);
            buffer.append("  #STD:");
            buffer.append(Math.floor(values.getValue()[2] * 100) / 100);
        }
        writer.write(buffer.toString());
        writer.newLine();
    }

    /**
     * This method searches for any items c that can expand left side of a rule I --> J
     * This results in rules of the form IU{c} --> J. 
     * The method saves the rule and continues expansion of the rule only if 
     * it is periodic and frequent in sufficient number of sequences in the 
     * database (determined by the value of minRa parameter).
     *
     * @param itemsetI   items contained in the left part of the rule
     * @param itemsetJ   items contained in the right part of the rule
     * @param seqI       sids of sequences already containing previous rule IJ
     * @throws IOException
     */
    private void expandLeft(List<Integer> itemsetI, List<Integer> itemsetJ, Set<Integer> seqI) throws IOException {
        // Create a map that contains all possible items c that can extend rule IJ 
        // on the left side (IU{C}-->J). It contains only items c that have 
        // enough common sids with itemsetI, and are not already contained in 
        // itemsets I or J and are larger than any item in itemsetI.
        Map<Integer, List<Integer>> seqICAll = new HashMap<Integer, List<Integer>>();
        
        itemLoop: for (Entry<Integer, Map<Integer, PISOccurence>> entry : mapItemsOccur.entrySet()) {
            Integer itemC = entry.getKey();
            if ((itemsetI.contains(itemC)) || (itemsetJ.contains(itemC))) {
                continue;
            }
            for (Integer i : itemsetI) {

                if (i > itemC) {
                    continue itemLoop;
                }
            }
            List<Integer> seqIC = new ArrayList<Integer>();
            for (Integer i : seqI) {
                if (entry.getValue().keySet().contains(i)) {
                    seqIC.add(i);
                }
            }

            if (!seqIC.isEmpty() && (((double)seqIC.size() / databaseSize) > minRa)) {
                seqICAll.put(itemC, seqIC);
            }
        }
        
        // for every item in seqiCAll, try to expand the rule with it
        for (Map.Entry<Integer, List<Integer>> seqList : seqICAll.entrySet()) {
            // create map of occurences of itemsets IC and J contained in each 
            // rule IU{C}-->J in each sequence.
            Map<Integer, PISOccurence> occicAll = new HashMap<Integer, PISOccurence>();
            Map<Integer, PISOccurence> icoccjAll = new HashMap<Integer, PISOccurence>();
            // get value of itemC
            Integer itemC = seqList.getKey();
            // create itemset IC
            List<Integer> itemsetIC = new ArrayList<Integer>();
            itemsetIC.addAll(itemsetI);
            itemsetIC.add(itemC);
            
            //System.out.println("EXPAND LEFT ");
            //System.out.println("I je " +itemsetIC.toString());
            //System.out.println("J je " +itemsetJ.toString());
            // calculates number of occurences of IC in all common sequences
            // of items in itemsetI and item c.
            Map<Integer, Integer> icNum = countIC(itemsetIC, seqList.getValue());
            // check if number of occurences of an itemset IC meets the treshold of minRa
            if (icNum.isEmpty() || (((double)icNum.size() / databaseSize) <= minRa)) {
                continue;
            }

            // calculates occurence of ICJ in all sequence common to itemsetI and item C
            Map<Integer, PISOccurence[]> occicjAll = getRuleOccur(itemsetIC, itemsetJ, seqList.getValue());
            // check if rule ICJ meets the minRa threshold
            if ((occicjAll.isEmpty()) || (((double) occicjAll.keySet().size() / databaseSize)) <= minRa) {
                continue;
            }
            
            // for each sequence
            for (Integer seq : seqList.getValue()) {
                // check if rule ICJ exists in that sequence
                if (occicjAll.get(seq) == null) {
                    continue;
                }
                
                PISOccurence[] occicjSeq = occicjAll.get(seq);
                // check if itemsetI and itemsetJ of rule ICJ 
                // appear less than 2 times in that sequence
                if (occicjSeq[0] == null || occicjSeq[1] == null || occicjSeq[1].firstItemset.size() < 2) {
                    continue;
                }
                // add occurences of that sequence to the map of occurences of ICJ
                PISOccurence occic = occicjSeq[0];
                PISOccurence occj = occicjSeq[1];
                occicAll.put(seq, occic);
                icoccjAll.put(seq, occj);
            }
           
            // create map of parameter values of the rule ICJ by sequence : 
            // key= sids of the rule ICJ, value = parameter values 
            // (support, confidence and standard deviation)
            Map<Integer, Double[]> parValICJ = new HashMap<Integer, Double[]>();
            
            // Set of sequences of the rule that meet the support threshold,
            // but may not meet thresholds of confidence and standard deviation
            // so they are saved for expansion only
            Set<Integer> seqForExpICJ= new HashSet<Integer>();

           // for each sequence
            for (Integer s : occicAll.keySet()) {
                int lenghtOfSeq = seqLenghts[s];
                // calculate parameter values of rule ICJ in that sequence
                Double[] parVal = isPisRule(occicAll.get(s), icoccjAll.get(s), icNum.get(s), lenghtOfSeq);
                 if (parVal[0] != null) {
                    // add sequence of the rule to list for expansion
                    seqForExpICJ.add(s);
                      }
                    if((parVal[0]!=null) && (parVal[1]!=null) && (parVal[2]!=null)){
                        // add values of frequent and periodic sequential rule
                        parValICJ.put(s, parVal);
                            }
            }
            
            // (2) check if the rule ICJ has enough sids
            // if not, we don't need to generate a rule.
                
            // create rule ICJ
                
            // chech if rule ICJ meets minRa threshold
            double raICJ = (double)parValICJ.keySet().size() / databaseSize;
            if ((!seqForExpICJ.isEmpty()) && (raICJ > minRa)) {
                
                // if rule meets the parameter thresholds
                if(!parValICJ.isEmpty()){
                // save the rule in a file
                saveRule(raICJ, itemsetIC, itemsetJ, parValICJ);}
                
                // recursive call to try to expand the rule on the left side
                expandLeft(itemsetIC, itemsetJ, seqForExpICJ);
                
            }

        }

        // check the memory usage
        MemoryLogger.getInstance().checkMemory();
    }

    /**
     * This method searches for any items c that can expand right side of a rule I --> J. 
     * This results in rules of the form I --> J U�{c}. 
     * The method saves the rule and continues expansion of the rule only if 
     * it is periodic and frequent in sufficient number of sequences in the 
     * database (determined by the value of minRa parameter).
     *
     * @param itemsetI   items contained in the left part of the rule
     * @param itemsetJ   items contained in the right part of the rule
     * @param seqI       sids of sequences already containing previous rule IJ
     * @param occiCount  a map of occurences of itemsetI, where 
     *                   key= sids of sequences already containing itemset I, 
     *                   and value = number of occurences of the itemset I 
     *                   in that sequence
     * @throws IOException
     */
    private void expandRight(List<Integer> itemsetI, List<Integer> itemsetJ, Set<Integer> seqJ, Map<Integer, Integer> occiCount) throws IOException {
         // Create a map that contains all possible items c that can extend rule IJ 
        // on the right side (I-->JU{C}). It contains only items c that have 
        // enough common sids with itemsetJ, and are not already contained in 
        // itemsets I or J and are larger than any item in itemsetJ.
        Map<Integer, List<Integer>> seqJCAll = new HashMap<Integer, List<Integer>>();
        
 itemLoop: for (Entry<Integer, Map<Integer, PISOccurence>> entry : mapItemsOccur.entrySet()) {
            Integer itemC = entry.getKey();
            if ((itemsetI.contains(itemC)) || (itemsetJ.contains(itemC))) {
                continue;
            }
            for (Integer i : itemsetJ) {
                if (i > itemC) {
                    continue itemLoop;
                }
            }
            List<Integer> seqJC = new ArrayList<Integer>();
            for (Integer i : seqJ) {
                if (entry.getValue().keySet().contains(i)) {
                    seqJC.add(i);
                }
            }

            if (!seqJC.isEmpty() && (((double)seqJC.size() / databaseSize) > minRa)) {
                seqJCAll.put(itemC, seqJC);
            }
        }
 
        // for every item in seqiCAll, try to expand the rule with it
        for (Map.Entry<Integer, List<Integer>> seqList : seqJCAll.entrySet()) {
            // create map of occurences of itemsets I and JC contained in each 
            // rule I-->JU{C} in each sequence.
            Map<Integer, PISOccurence> jcocciAll = new HashMap<Integer, PISOccurence>();
            Map<Integer, PISOccurence> occjcAll = new HashMap<Integer, PISOccurence>();
            
            // get value of itemC
            Integer itemC = seqList.getKey();
            
            // create itemset JC
            List<Integer> itemsetJC = new ArrayList<Integer>();
            itemsetJC.addAll(itemsetJ);
            itemsetJC.add(itemC);
            
            //System.out.println("EXPAND RIGHT ");
            //System.out.println("I je " +itemsetI.toString());
            //System.out.println("J je " +itemsetJC.toString());
            
            
            // calculates occurence of IJC in all sequence common to itemsetI and item C
            Map<Integer, PISOccurence[]> occijcAll = getRuleOccur(itemsetI, itemsetJC, seqList.getValue());
            // check if the rule IJC meets the minRa threshold
            if ((occijcAll.isEmpty())|| (((double) occijcAll.keySet().size() / databaseSize ) <= minRa)) {
                continue;
            }
            
            // for each sequence
            for (Integer seq : seqList.getValue()) {
 
                // check if rule ICJ exists in that sequence
                if (occijcAll.get(seq) == null) {
                    continue;
                }
                PISOccurence[] occicjSeq = occijcAll.get(seq);
                // check if itemsetI and itemsetJ of rule IJC 
                // appear less than 2 times in that sequence
                if (occicjSeq[0] == null || occicjSeq[1] == null || occicjSeq[1].firstItemset.size() < 2) {
                    continue;
                }
                // add occurences of that sequence to the map of occurences of IJC
                PISOccurence occi = occicjSeq[0];
                PISOccurence occjc = occicjSeq[1];
                jcocciAll.put(seq, occi);
                occjcAll.put(seq, occjc);

            }

            // create map of parameter values of the rule IJC by sequence : 
            // key= sids of the rule IJC, value = parameter values 
            // (support, confidence and standard deviation)
            Map<Integer, Double[]> parValIJC = new HashMap<Integer, Double[]>();
            
            // Set of sequences of the rule that meet the support threshold,
            // but may not meet thresholds of confidence and standard deviation
            // so they are saved for expansion only
            Set<Integer> seqForExpIJC= new HashSet<Integer>();
            
            // for each sequence
            for (Integer s : occjcAll.keySet()) {
                int lenghtOfSeq = seqLenghts[s];
                // calculate parameter values of rule IJC in that sequence
                Double[] parVal = isPisRule(jcocciAll.get(s), occjcAll.get(s), occiCount.get(s), lenghtOfSeq);
                if (parVal[0] != null) {
                    // add sequence of the rule to list for expansion
                    seqForExpIJC.add(s);
                      }
                    if((parVal[0]!=null) && (parVal[1]!=null) && (parVal[2]!=null)){
                        // add values of frequent and periodic sequential rule
                        parValIJC.put(s, parVal);
                            }
            }
            
            // (2) check if the rule IJC has enough sids
            // if not, we don't need to generate a rule.
                
            // create rule IJC
                
            // chech if rule IJC meets minRa threshold
            double raIJC = (double)parValIJC.keySet().size() / databaseSize;
            if ((!seqForExpIJC.isEmpty()) && (raIJC > minRa)) {
                
                // if rule meets the parameter thresholds
                 if(!parValIJC.isEmpty()){
                // save the rule in a file
                saveRule(raIJC, itemsetI, itemsetJC, parValIJC);}
                
                // recursive call to try to expand the rule on the left and
                // right sides
                
                expandLeft(itemsetI, itemsetJC, seqForExpIJC);
                
                expandRight(itemsetI, itemsetJC, seqForExpIJC, occiCount);
      
            }

        }

        // check the memory usage
        MemoryLogger.getInstance().checkMemory();
    }

    /**
     * This method calculate the frequency of each item and removes items that
     * are not frequent in one database pass.
     *
     * @param database : a sequence database
     * @return A map such that key = item, value = a map, where a key = sid and
       a value = PISOccurence This map allows knowing the frequency of each item
       and all their occurences in each sequence.
     */
    private Map<Integer, Map<Integer, PISOccurence>> getOccurMapAndRemoveInfrequent(SequenceDatabase database) {
        // map of occurences of an item in database <item, Map<sid, occurence>>
        mapItemsOccur = new HashMap<Integer, Map<Integer, PISOccurence>>();

        // for each sequence in the database
        for (int k = 0; k < database.size(); k++) {
            // create set of non frequent items in sequence k
            HashSet<Integer> itemSeqRemove = new HashSet<Integer>();
            Sequence sequence = database.getSequences().get(k);
            // for each itemset in that sequence
            for (int j = 0; j < sequence.getItemsets().size(); j++) {
                List<Integer> itemset = sequence.get(j);
                // for each item in that itemset
                for (int i = 0; i < itemset.size(); i++) {
                    Integer itemI = itemset.get(i);
                    // get the map of occurences of that item
                    Map<Integer, PISOccurence> occurences = mapItemsOccur.get(itemI);
                    // if this map is null, create a new one
                    if (occurences == null) {
                        occurences = new HashMap<Integer, PISOccurence>();
                        // add occurences to the map
                        mapItemsOccur.put(itemI, occurences);
                    }

                    // get occurence of that item in sequence k
                    PISOccurence occurence = occurences.get(k);
                    // if occurence is null, create new occurence
                    // and add first occurence of an item in sequence k
                    if (occurence == null) {
                        List<Integer> frstlastocc = new ArrayList<Integer>();
                        frstlastocc.add(j);
                        occurence = new PISOccurence(frstlastocc, frstlastocc);
                        occurences.put(k, occurence);
                    } else {
                        // update the occurence by adding 
                        // j as occurence in sequence k
                        occurence.lastItemset.add(j);
                    }
                }
            }

            // Remove infrequent items in sequence k from the database 
            // and create list of infrequent items in sequence k  
            for (int item : mapItemsOccur.keySet()) {
                if (mapItemsOccur.get(item).get(k) != null
                        && ((((double) mapItemsOccur.get(item).get(k).firstItemset.size() / (double) sequence.getItemsets().size()) <= minSupp)
                        || (double) mapItemsOccur.get(item).get(k).firstItemset.size() < 2)) {

                    for (int tid : mapItemsOccur.get(item).get(k).firstItemset) {
                        //System.out.println(" remove item "+item+" k "+k+" tid "+tid);
                        List<Integer> itemset = sequence.get(tid);
                        itemset.remove((Integer) item);
                        itemSeqRemove.add(item);
                    }
                }
            }
            //System.out.println("NUMBER OF DIFFERENT ITEMS : " + mapItemsOccur.size());
            
            // Remove infrequent items in sequence k from map of occurences of items
            for (int l : itemSeqRemove) {
                mapItemsOccur.get(l).remove(k);
                if (mapItemsOccur.get(l).size() < 1) {
                    mapItemsOccur.remove(l);
                }
            }
        }

        System.out.println("NUMBER OF DIFFERENT ITEMS : " + mapItemsOccur.size());
        
        // create list of items with support below minRa threshold and remove 
        // them from map of occurences
        List<Integer> itemsToRemove = new ArrayList<Integer>();
        for (int item : mapItemsOccur.keySet()) {
            if (((double) mapItemsOccur.get(item).size() / database.size()) <= minRa) {
                itemsToRemove.add(item);
            }
        }
        for (int item : itemsToRemove) {
            mapItemsOccur.remove(item);
        }
        
        return mapItemsOccur;
    }

    /**
     * Print statistics about the last algorithm execution to System.out.
     */
    public void printStats() {
        System.out.println("=============  SEQRULEGROWTH - STATS ========");
        System.out.println("Sequential rules count: " + ruleCount);
        System.out.println("Total time: " + (timeEnd - timeStart) + " ms");
        System.out.println("Max memory: " + MemoryLogger.getInstance().getMaxMemory());
        System.out.println("==========================================");
    }

    /**
	 * This method returns occurences of items of a sequential rule 
         * with 1 item as antecedent and consequent in requested sequence.
	 * 
	 * @param occi          occurences of itemI in a sequence
	 * @param occj          occurences of itemJ in a sequence
	 * @param lenghtOfSeq   number of all tids in a sequence
         * @return              array of occurences of items I and J in rule IJ
         *                      in a sequence
	 */
    public PISOccurence[] getSingleItemRuleOccur(List<Integer> occi, List<Integer> occj, int lenghtOfSeq) {
        // initiallize last index of itemJ (marks ending of one rule)
        int previousIndex = -1;
        // create list of tids of items
        List<Integer> newocci = new ArrayList<Integer>();
        List<Integer> newoccj = new ArrayList<Integer>();
        // reduces number of iterations
        int counterJ = 0;
        // for each occurence of itemI
        for (int i = 0; i < occi.size(); i++) {
            // if tids of item I is after ending of previuos rule
            if (occi.get(i) > previousIndex) {
                for (int j = counterJ; j < occj.size(); j++) {
                    counterJ = j + 1;
                    // if item J is after item I
                    if (occj.get(j) > occi.get(i)) {
                        previousIndex = occj.get(j);
                        // add occurence to list of occurences
                        newocci.add(occi.get(i));
                        newoccj.add(occj.get(j));
                        break;
                    }
                }
            }
        }
        // if IJ appears only once or doesn't meet the minimum support threshold
        if ((newocci.size() < 2) || (((double) newocci.size() / lenghtOfSeq) <= minSupp)) {
            return null;
        }
        // create occurences of items I and J in rule IJ in a sequence
        PISOccurence IJocci = new PISOccurence(newocci, newocci);
        PISOccurence IJoccj = new PISOccurence(newoccj, newoccj);
        PISOccurence[] listOcc = new PISOccurence[2];
        listOcc[0] = IJocci;
        listOcc[1] = IJoccj;
        return listOcc;
    }

    /**
	 * This method checks if the rule is periodic and frequent in a sequence.
         * It checks if the rule meets parameter thresholds in a sequence.
	 * @param occInew       occurence of itemsetI in rule IJ in a sequence
         * @param occJnew       occurence of itemsetJ in rule IJ in a sequence
	 * @param lengthOfSeq   number of all the tids in a sequence
	 * @return  An array of values of the rule IJ for a sequence: 
         *          (support, confidence, standard deviation)
	 */
    Double[] isPisRule(PISOccurence occInew, PISOccurence occJnew, int occI, int lengthOfSeq) {
        double suppIJ = (double) occInew.firstItemset.size() / lengthOfSeq;
        if (suppIJ <= minSupp) {
            return null;
        }

        Double confIJ = (double) (occInew.firstItemset.size()) / (double) (occI);

        if (confIJ <= minConf) {
            confIJ=null;
        }
        // calculate standard deviation of the rule IJ in a sequence
        Double stdIJ = getStd(occInew, occJnew, lengthOfSeq);

        if (stdIJ >= maxStd) {
            stdIJ=null;
        }

        Double[] values = new Double[3];
        values[0] = suppIJ;
        values[1] = confIJ;
        values[2] = stdIJ;
        return values;
    }
    
/**
	 * This method calculates standard deviation of a rule in a sequence.
	 * @param oci          occurence of itemsetI in rule IJ in a sequence
         * @param ocj          occurence of itemsetJ in rule IJ in a sequence
	 * @param seqLen       number of all the tids in a sequence
	 * @return  value of the standard deviation of a rule in a sequence.
	 */
    public double getStd(PISOccurence oci, PISOccurence ocj, int seqLen) {
       
        int preTID = 0;
        double stanDev = 0;
        for (int tid : ocj.lastItemset) {

            tid += 1;

            int perI = tid - preTID;

            stanDev = stanDev + Math.pow(perI, 2);
            preTID = tid;
        }
        stanDev = stanDev + Math.pow(seqLen - preTID, 2);
        
        stanDev = stanDev / (double) (oci.firstItemset.size() + 1);
       
        stanDev = stanDev - Math.pow((double) seqLen / (double) (oci.firstItemset.size() + 1), 2);
        
        stanDev = Math.sqrt(stanDev);

        return stanDev;
    }

     /**
	 * This method returns occurences of items of a sequential rule IJ in 
         * requested sequences.
	 * 
	 * @param itemsetI   items contained in the left part of the rule
         * @param itemsetJ   items contained in the right part of the rule
         * @param seq        sids of sequences already containing itemsetIC
         * @return map of occurences of rule IJ, where key= sid,value = array of
         *         occurences of itemsetI and itemsetJ in the rule IJ in that sid.                
	 */
    public Map<Integer, PISOccurence[]> getRuleOccur(List<Integer> itemsetI, List<Integer> itemsetJ, List<Integer> seq) {
        // create map of occurences of itemsets I and J in the rule IJ in each sequence
        Map<Integer, PISOccurence[]> occIJAll = new HashMap<Integer, PISOccurence[]>();
        // for each sequence
        for (int s : seq) {
            // initialize default value of condition, if it changes to false, 
            // IJ has no more occurences in that sequence and we move to the next sequence 
            boolean cond = true;
            // initialize index of the current last occurence of the rule IJ in that sequence
            int previousPosition = -1;
            // create lists of frst and last occurences of itemsets I and J in the rule IJ in that sequence
            List<Integer> frstOcci = new ArrayList<Integer>();
            List<Integer> lastOcci = new ArrayList<Integer>();

            List<Integer> frstOccj = new ArrayList<Integer>();
            List<Integer> lastOccj = new ArrayList<Integer>();

            while (cond) {
                // initializes values of first and last occurence next of itemsetI
                int frstI = -1;
                int lastI = -1;
                // create array of frst and last next occurence of itemsetI
                int[] frstlastI = getNextIC(itemsetI, s, previousPosition);
                if (frstlastI == null) {
                    cond = false;
                    break;
                } else {
                    frstI = frstlastI[0];
                    lastI = frstlastI[1];
                }

                // initialize index of the current last occurence of the itemsetI 
                // in the rule IJ in that sequence
                int positionI = lastI;
                // initializes values of first and last occurence next of itemsetJ
                int frstJ = -1;
                int lastJ = -1;
                // create array of frst and last next occurence of itemsetJ
                int[] frstlastJ = getNextIC(itemsetJ, s, positionI);
                if (frstlastJ == null) {
                    cond = false;
                    break;
                } else {
                    frstJ = frstlastJ[0];
                    lastJ = frstlastJ[1];
                }

                frstOcci.add(frstI);
                lastOcci.add(lastI);
                frstOccj.add(frstJ);
                lastOccj.add(lastJ);
                previousPosition = lastJ;
            }
            // check if rule occurs minimum 2 times
            if (frstOcci.size() < 2) {
                continue;
            }
            // check support of the rule IJ in that sequence
            if (((double) frstOcci.size() / seqLenghts[s]) <= minSupp) {
                continue;
            }
            // create array of occurences of IJ rule for that sequence
            PISOccurence occurI = new PISOccurence(frstOcci, lastOcci);
            PISOccurence occurJ = new PISOccurence(frstOccj, lastOccj);
            PISOccurence[] occurences = new PISOccurence[2];
            occurences[0] = occurI;
            occurences[1] = occurJ;
            occIJAll.put(s, occurences);
        }
        return occIJAll;
    }

    /**
	 * This method calculates number of occurences of IC in each requested sequence.
	 * @param itemsetIC an itemset containing items from itemsetI and item c
	 * @param seq       sids of sequences that we want to calculate occurence for.
	 * @return a map, where key= sid, value = number of occurences of IC in that sequence
	 */
    public Map<Integer, Integer> countIC(List<Integer> itemsetIC, List<Integer> seq) {
        // creates map of number of occurences of IC in all requested sequences
        Map<Integer, Integer> mapCount = new HashMap<Integer, Integer>();
        
        // for each sequence
        for (int s : seq) {
            // initiallize number of occurence of IC in that sequence
            int count = 0;
            // initialize default value of condition, if it changes to false, 
            // IC has no more occurences in that sequence and we move to the next sequence 
            boolean cond = true;
            // initialize index of the current last occurence of itemsetIC in that sequence
            int previousPosition = -1;
            
            while (cond) {
                // create array of frst and last occurence of itemsetIC
                int[] frstlastI = getNextIC(itemsetIC, s, previousPosition);
                if (frstlastI == null) {
                    cond = false;
                    break;
                } else {

                    previousPosition = frstlastI[1];
                    count += 1;
                }

            }
            // check if itemsetIC has enough occurences in that sequence that 
            // meets the minimum support treshold, if it does, add occurence to the map
            if ((count > 1) && (((double) count / seqLenghts[s]) > minSupp)) {
                mapCount.put(s, count);
            }
        }
        return mapCount;
    }

     /**
	 * This method searches for next occurence of itemset IC, 
         * from assigned index in requested sequence.
	 * @param itemsetI          an itemset containing items from itemsetI and item c
	 * @param s                 sid of sequences that we want to calculate occurence for.
         * @param previousPosition  index that we want to search occurence from.
	 * @return array of the frst and last next occurence of itemsetIC.
	 */
    public int[] getNextIC(List<Integer> itemsetIC, int s, int previousPosition) {
        // create list of next indexes(positions) of each item contained in itemsetIC
        List<Integer> posList = new ArrayList<Integer>();

        // for each item in itemsetIC
        for (Integer item : itemsetIC) {
            // get occurences of an item
            List<Integer> occurItem = mapItemsOccur.get(item).get(s).firstItemset;
            // initialize counter of index of each item position in occurItemList
            int index = -1;
            // initialize variable that holds next occurence of the item
            int tempItem = -1;
            
            // for each occurence of the item
itemLoop:       for (int i : occurItem) {
                index++;
                // checks if occurence is next, 
                // e.g. if it has greater value than previous occurence
                if (i > previousPosition) {
                    tempItem = i;
                    break;
                }
            }

            // checks if next occurence of the item exists, 
            // if it does, adds them to the position list
            if (tempItem != -1) {
                posList.add(occurItem.get(index));
            } else {
                return null;
            }

        }
        
        // initializes values of first and last occurence of itemsetIC
        int frst = seqLenghts[s];
        int last = -1;

        // calculates first and last occurence of itemsetIC, 
        // among occurences of all contained items in that itemset
        for (int k = 0; k < (posList.size()); k++) {
            if (posList.get(k) < frst) {
                frst = posList.get(k);
            }
            if (posList.get(k) > last) {
                last = posList.get(k);
            }
        }
        // create array of the next first and last occurence of itemsetIC
        int[] frstlastIC = new int[2];
        frstlastIC[0] = frst;
        frstlastIC[1] = last;
        return frstlastIC;
    }
}
