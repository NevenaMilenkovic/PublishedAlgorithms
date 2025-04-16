package ca.pfv.spmf.algorithms.sequential_rules.pisrulegrowth;

import ca.pfv.spmf.input.sequence_database_list_integers.Sequence;
import ca.pfv.spmf.input.sequence_database_list_integers.SequenceDatabase;
import java.util.List;

/**
 * This class represent the list of first and last occurences of an itemset in a
 * sequence, as defined in the PISRuleGrowth algorithm.
 *
 * @see AlgoPISRuleGrowth
 * @see Sequence
 * @see SequenceDatabase
 * @author Nevena Milenkovic
 */
public class PISOccurence {

    /**
     * the list of first occurences of an itemset <br/>
     */
    public List<Integer> firstItemset;
    /**
     * the list of last occurences of an itemset
     */
    public List<Integer> lastItemset;

    /**
     * Constructor
     *
     * @param firstItemset the list of first occurences of an itemset
     * @param lastItemset the list of last occurences of an itemset
     */
    public PISOccurence(List<Integer> firstItemset, List<Integer> lastItemset) {
        this.firstItemset = firstItemset;
        this.lastItemset = lastItemset;
    }
}
