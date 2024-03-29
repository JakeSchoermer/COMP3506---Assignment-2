package a2;


import java.util.Arrays;

public class Alignment {

	private final PerfMeter perf; // monitor of performance
	private final Sequence[] dna; // all sequences that are used in the
  								  // alignment
	private int N = 0; 	// sequence length (each sequence must be the same)
	private int W; 		// width of alignment; N <= W <= N*2-1
	private boolean reverse = false;// whether the reverse strand should be
									// considered

    private int currentBest = 0;

	/**
	 * Construct an alignment from a list of DNA sequences.
	 *
	 * @param perf
	 *            performance meter
	 * @param dna
	 *            DNA sequences
	 */
	public Alignment(PerfMeter perf, Sequence[] dna) {
		this.perf = perf;
		this.dna = dna;
		for (int i = 0; i < dna.length; i++) {
			if (N == 0)
				N = dna[i].getLength();
			else if (N != dna[i].getLength())
				throw new RuntimeException("Different lengths of sequences");
		}
		this.W = N * 2 - 1;
		this.reverse = false;
	}

	/**
	 * Construct an alignment from a list of DNA sequences. Note that these
	 * sequences can come from either strand and that one or the other may align
	 * better to the other strings.
	 *
	 * @param perf
	 *            performance meter
	 * @param dna
	 *            DNA sequences
	 * @param reverse
	 *            use the reverse strand in addition to the one provided in DNA
	 *            sequence
	 * @param W
	 *            the width of the alignment (must be between N and N*2-1 where
	 *            N is the length of each sequence being aligned.
	 */
	public Alignment(PerfMeter perf, Sequence[] dna, boolean reverse, int W) {
		this(perf, dna);
		// make sure the width of the alignment is at least as great as the
		// length of each sequence
		// and not wider than two when they overlap
		if (W < N || W > N * 2 - 1)
			this.W = N * 2 - 1;
		else
			this.W = W;
		this.reverse = reverse;
	}

	/**
	 * Get number of symbols in each sequence (N)
	 */
	public int getN() {
		return N;
	}

	/**
	 * Determine the level of search (i.e. number of offset indices that have
	 * been assigned).
	 *
	 * @param s
	 *            the offset indices in alignment (-1 or end-of-array
	 *            terminates)
	 * @return number of sequences in current alignment
	 */
	public int getLevel(int[] s) {
		for (int i = 0; i < s.length; i++)
			if (s[i] == -1)
				return i;
		return s.length;
	}

	/**
	 * Find the symbol with the max count for each column of a profile
	 *
	 * @param profile
	 *            the profile for an alignment
	 * @return an array with symbol indices, identifying the consensus
	 */
	public int[] getConsensus(int[][] profile) {


        int[] sym = new int[W];
		// Problem 2: Your code here

        for (int j=0; j <profile[0].length; j++) {  //Column
            for (int a=0; a<profile.length; a++) {   //Letter
                if (a == 0) {
                    sym[j] = a+1;
                }
                else {
                    if (profile[a][j] > profile[a-1][j]) {
                        sym[j] = a+1;
                    }
                }
            }
        }
		return sym;
	}

	/**
	 * Determine the profile (M) of a specified alignment, i.e. the counts of
	 * each of the four symbols for all of W positions.
	 *
	 * @param s
	 *            the alignment represented by offset indices
	 * @return the profile
	 */
	public int[][] getProfile(int[] s) {
		int[][] profile = new int[Sequence.alpha.length][W]; //4, 2 * N - 1
		//Problem 1: Your code here

        char nucleotide;
        int nucleotideIdx;

        //If s[i] >= 0 AND s[i] > W-N, then do getSymbolIndex(j, false)
        //else, if s[i] >= 0, then do getSymbolIndex(j, true)
        //You also need to +=1 profile[a-1][j + s[i] - (W - N + 1) (this goes in the first if statement)
        //In the else if, update profile[a-1][j+s[i]] = profile[a-1][j+s[i]] + 1
         for (int i = 0; i<dna.length; i++) {
            for (int j=0; j<N;j++) {
                int offset = s[i];
                if (offset >= 0 && offset > W - this.getN()) {
                    nucleotideIdx = dna[i].getSymbolIndex(j, false);
                    profile[nucleotideIdx-1][j+offset-(W-N+1)] = profile[nucleotideIdx-1][j+offset - (W-N+1)] + 1;
                }
                else if (offset >=0) {
                    nucleotideIdx = dna[i].getSymbolIndex(j, true);
                    profile[nucleotideIdx-1][j+offset] +=1;
                }
            }
        }

		return profile;
	}

	/**
	 * Determine the score of a consensus given a profile
	 *
	 * @param profile
	 *            the profile of an alignment
	 * @param consensus
	 *            an array of symbol indices
	 * @return the score indicating the quality of the consensus
	 */
	public int getScore(int[][] profile, int[] consensus) {
		int score = 0;
		// Problem 3: Your code here
        for (int j = 0; j<profile[0].length;j++) {
            score += profile[consensus[j]-1][j];
        }
		return score;
	}

	/**
	 * Given a partial alignment, expand it to all possible alignments
	 * incorporating one more sequence.
	 * Note this version does NOT consider the reverse/complementary strand.
	 *
	 * @param s
	 *            the partial alignment (s.length is total number of sequences,
	 *            -1 marks end)
	 * @return alignments including one more sequence
	 */
	public int[][] expand(int[] s) {
		int level = getLevel(s);
        int mult;


        if (!this.reverse) {
            mult =1;
        }
        else {
            mult =2;
        }
        int[][] s_copies = new int[mult * (W - N + 1)][s.length];

		if (level == s.length) {
            return null;
        }

        for (int i = 0; i < (W - N + 1); i++) { // we can shift pattern this far
            for (int j = 0; j < s.length; j ++)
                s_copies[i][j] = s[j];
            s_copies[i][level] = i;
        }

        if (this.reverse) {
		    for (int i = (W-N+1); i < 2 * (W - N + 1); i++) { // we can shift pattern this far
                for (int j = 0; j < s.length; j ++)
                    s_copies[i][j] = s[j];
                s_copies[i][level] = i;
            }
        }

        return s_copies;

	}

    /**
     * Implement Branch-and-Bound for AlignmentScore
     *
     * @param current
     *      the current score at this point in the findAlignment process.
     *
     * @param length
     *      the the length of the next node in the tree.
     *
     *
     * @return true if we think we can get a better comparison, else return false
     */

    private boolean chanceOfImproving(AlignmentScore current, int length) {
        if (current.actual + length < this.currentBest) {
            return false;
        }
        return true;
    }

    /**
	 * Search for the optimal alignment.
	 *
	 * @param s
	 *            the offset indices for the proposed and possibly partial
	 *            alignment
	 * @return the best as far as we know or null
	 */

    public AlignmentScore findAlignment(int[] s) {

        perf.countFind(); // ********DO NOT REMOVE********//

		int level = getLevel(s); 	// Will be 0 first call when the s[0] == -1
									// before any offsets have been set for ANY
									// sequence.
									// Will be 1 for first sequence, s.length
									// for last sequence.

		// calculate the score for the current level
		int[][] profile = getProfile(s);
		int[] maxsym = getConsensus(profile);
		int score = getScore(profile, maxsym);
		AlignmentScore current = new AlignmentScore(score, s);

        if (score > this.currentBest) {
            this.currentBest = score;
        }

		if (level == s.length) { // At leaf node
            perf.countLeaf(); // ********DO NOT REMOVE********//
			return current;
		}

		// Make assessment of how good things can be, from here on,
		// because if there's no chance of it improving on an optimistic estimate...
		// then we give up.

        int length = (N*(s.length - level));

		if (!chanceOfImproving(current, length)) {
			perf.countBreak(); // ********DO NOT REMOVE********//
            return null;
		}

		// Generate all child nodes.
		// Consider ordering of the branches?
		int[][] extensions = expand(s);

		AlignmentScore bestScore = null;
		for (int i = 0; i < extensions.length; i++) {
            int[][] ext_prefix = extensions;
			AlignmentScore nextScore = findAlignment(extensions[i]);
			if (nextScore != null) {
				if (bestScore == null)
					bestScore = nextScore;
				else if (nextScore.actual > bestScore.actual)
					bestScore = nextScore;
			}
		}
		perf.countPropagate(); // ********DO NOT REMOVE********//

        return bestScore;
	}

	/**
	 * Holder of score and the path (offset indices defining the alignment)
	 * leading to those scores. You may modify the code for this but keep the
	 * original constructor signature.
	 */
	class AlignmentScore {

		final int actual; // the actual score for this config
		final int[] path; // the offset indices to which the distance applies

		/**
		 * Constructs an instance that combines a distance and the applicable
		 * offset indices.
		 *
		 * @param actual
		 *            the actual distance
		 * @param s
		 *            the offset indices that apply
		 */
		public AlignmentScore(int actual, int[] s) {
			this.actual = actual;
			this.path = s;
		}

		public String toString() {
			StringBuffer sbuf = new StringBuffer();
			for (int i = 0; i < path.length; i++) {
				if (path[i] == -1)
					break;
				sbuf.append(path[i] + ";");
			}
			sbuf.append(":" + actual);
			return sbuf.toString();
		}
	}
}
