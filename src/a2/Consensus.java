package a2;

import java.util.Arrays;

/**
 * Class for searching for a consensus in DNA sequence data.
 */
public class Consensus {

	private final PerfMeter perf; // a performance monitor class
	private final Sequence[] dna; // the sequence data that are searched
	private int N = 0; // number of symbols in each sequence
	private boolean reverse = false;// whether the consensus considers the
									// reverse strand
	private int W = -1; // the width of the sought consensus

	/**
	 * Constructs an instance of the class, prepares for search and checks for
	 * problems. Note if this constructor is used then the reverse strand is NOT
	 * used, and the width of the consensus is 2N - 1 (where N is the length of
	 * each sequence)
	 * 
	 * @param perf
	 *            the performance monitor
	 * @param seqs
	 *            the sequence data
	 */
	public Consensus(PerfMeter perf, Sequence[] dna) {
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
	 * Constructs an instance of the class, prepares for search and checks for
	 * problems.
	 * 
	 * @param perf
	 *            the performance monitor
	 * @param seqs
	 *            the sequence data
	 * @param reverse
	 *            use the reverse strand in addition to that in data set when
	 *            finding the consensus sequence
	 * @param W
	 *            the width of the consensus (set to 2N - 1 if W < N or W > 2N - 1)
	 */
	public Consensus(PerfMeter perf, Sequence[] dna, boolean reverse, int W) {
		this(perf, dna);
		if (W < N || W > N * 2 - 1)
			this.W = N * 2 - 1;
		else
			this.W = W;
		this.reverse = reverse;
	}

	/**
	 * Get number symbols in each sequence
	 */
	public int getN() {
		return N;
	}

	/**
	 * Determines the "level of search", i.e. the number of symbols that are
	 * included in the consensus. We refer to an incomplete consensus as a
	 * "prefix", an array as long as the consensus with all elements set to -1
	 * except those in the prefix, in which symbols are identified by their
	 * index (A=1, C=2, G=3 and T=4), e.g. int[] prefix = [1, 2, 3, -1, -1] for
	 * a W=5 consensus with the prefix 'ACG'.
	 * 
	 * @param prefix
	 *            an array containing assigned symbols
	 * @return the number of symbols in the prefix
	 */
	public int getLevel(int[] prefix) {
		for (int i = 0; i < prefix.length; i++) {
			if (prefix[i] == -1)
				return i;
		}
		return prefix.length;
	}

	/**
	 * Score each sequence against the specified consensus sequence (or prefix).
	 * The consensus must be at least as long as the sequence.
	 * 
	 * @param consensus
	 *            the sequence of symbols in the consensus
	 * @return the score (number of sequence positions that overlap with the
	 *         consensus sequence)
	 */
	public int getScore(int[] consensus) {
		int level = getLevel(consensus); // find how many symbols we are
											// checking
		int totscore = 0;
		for (int j = 0; j < dna.length; j++) { // we sum the score over all
												// sequences
			int[] seq_orig = dna[j].getSymbolIndices(true); // the original
															// strand
			int[] seq_reverse = dna[j].getSymbolIndices(false); // the reverse
																// strand
			int bestscore = 0;
			for (int i = 0; i < W - N + 1; i++) { // for each start position in
													// the consensus
				int score = 0; // note: we stop counting when the sequence does
								// not fit
				for (int w = 0; w < N; w++) { // walk through the whole sequence
												// (orig strand)
					if (i + w < level) // only count as far as symbols have been
										// assigned
						// add one if the sequence and consensus share symbol
						score += (seq_orig[w] == consensus[i + w]) ? 1 : 0;
				}
				if (score > bestscore) // if this is the best we've seen...
					bestscore = score; // remember it
				// if we are looking at the reverse strand sequence we do the
				// same again
				if (reverse) {
					score = 0;
					for (int w = 0; w < N; w++) { // walk through the reverse
													// strand
						if (i + w < level) // only count as far as symbols have
											// been assigned
							// add one if the reverse sequence and consensus
							// share the symbol
							score += (seq_reverse[w] == consensus[i + w]) ? 1 : 0;
					}
					if (score > bestscore)
						bestscore = score;
				}
			}
			totscore += bestscore;
		}
		return totscore;
	}

	/**
	 * Determine the alignment of sequences from a consensus (possibly partial).
	 * Note: this method is almost identical to @see {@link #getScore(int[])}.
	 * 
	 * @param consensus
	 *            the symbols making up the consensus
	 * @return the offset indices that define the alignment
	 */
	public int[] getAlignment(int[] consensus) {
		int level = getLevel(consensus);
		int totscore = 0;
		int[] s = new int[dna.length]; // offsets for alignment
		for (int j = 0; j < dna.length; j++) {
			int[] seq_orig = dna[j].getSymbolIndices(true);
			int[] seq_reverse = dna[j].getSymbolIndices(false);
			int bestscore = 0;
			for (int i = 0; i < W - N + 1; i++) { // for each start position in
													// the consensus
				int score = 0;
				for (int w = 0; w < N; w++) { // walk through the whole sequence
												// (orig strand)
					if (i + w < level)
						score += (seq_orig[w] == consensus[i + w]) ? 1 : 0;
				}
				if (score > bestscore) {
					bestscore = score;
					s[j] = i; // alignment
				}
				if (reverse) {
					score = 0;
					for (int w = 0; w < N; w++) { // walk through the reverse
													// strand too
						if (i + w < level)
							score += (seq_reverse[w] == consensus[i + w]) ? 1 : 0;
					}
					if (score > bestscore) {
						bestscore = score;
						s[j] = i + (W - N + 1); // alignment reverse strand
					}
				}
			}
			totscore += bestscore;
		}
		return s;
	}

	/**
	 * Expands the current prefix (aka partial consensus) into all prefixes that
	 * can be constructed by adding a single symbol from the alphabet. The order
	 * of prefixes follows the order defined for the alphabet.
	 * 
	 * @return an array of prefixes extending the current by one symbol
	 */
	public int[][] expand(int[] prefix) {
		int level = getLevel(prefix); // check length of prefix
		if (level == prefix.length) // if complete, fail
			return null;
		// we need to make separate copies for each new prefix to explore
		int[][] extensions = new int[Sequence.alpha.length][];
		for (int i = 0; i < Sequence.alpha.length; i++) {
			extensions[i] = Arrays.copyOf(prefix, prefix.length);
			extensions[i][level] = i + 1; // symbol index is i + 1 (A, C, G, T)
		}
		return extensions;
	}

	/**
	 * Search for consensus of a specified length, starting from the given
	 * prefix.
	 * 
	 * @param the
	 *            prefix from which the search is started, e.g. if prefix is
	 *            "AC", we explore "ACA", "ACC", "ACG" and "ACT" (and beyond
	 *            through recursive calls). Note: to understand the performance
	 *            of the search, we use the performance monitor supplied to the
	 *            constructor to count the number of "nodes" we search, and how
	 *            they are "exited."
	 * @param cutoff
	 *            is a "bound" that can be used to break branching
	 * @return the score of the best consensus (below this point in the tree)
	 */
	public ConsensusScore findConsensus(int[] prefix, int cutoff) {
		// We count the number of nodes examined (i.e. number of prefixes
		// explored)
		perf.countFind();

		int nSym = getLevel(prefix); // number of symbols in current string
		if (nSym == W) { // full string: we can determine score
			int actual = getScore(prefix);
			perf.countLeaf(); // we count this as a "leaf" exit
			if (actual < cutoff) // disregard the score if we're
				return null; // doing worse than we can do elsewhere
			else
				// if good, then return actual score
				return new ConsensusScore(actual, prefix);
		} else { // we are looking at an incomplete consensus (i.e. a prefix)
			// make an assessment on whether we can do better than "cutoff"...
			int actual = getScore(prefix);
			// we've just scored a shortened version... determine the best we
			// can do
			if (actual + dna.length * (W - nSym) < cutoff) {
				perf.countBreak(); // if there is no hope of improving, we give
									// up
				return null;
			}

			int[][] extend = expand(prefix); // find all prefixes that extend
												// the current
			ConsensusScore best = null; // remember the best because that is the
										// only one we return
			for (int i = 0; i < extend.length; i++) { // go through all extended
													// prefixes
				// Notice that the recursive call (below) will involve scoring
				// the same sequences again
				// but this time with ONE more symbol specified. Consider
				// passing this distance vector
				// to improve performance?
				ConsensusScore current = findConsensus(extend[i], cutoff);

				if (current != null) { // if we got a result (that was better than
										// the "cutoff")
					cutoff = current.actual; // update the cutoff for next
												// "branch"
					best = current; // and remember the score and path
				}
			}
			perf.countPropagate(); // exit by returning the value
									// "up the search tree"
			return best;
		}
	}

	/**
	 * Holder of score and the path (the symbols making up the consensus)
	 * leading to those scores.
	 */
	public class ConsensusScore {

		final int actual; // the actual distance of this k-mer
		final int[] path; // the k-mer to which the distance applies

		/**
		 * Constructs an instance that combines a distance and the applicable
		 * k-mer.
		 * 
		 * @param actual
		 *            the actual distance
		 * @param kmer
		 *            the k-mer
		 */
		public ConsensusScore(int actual, int[] kmer) {
			this.actual = actual;
			this.path = kmer;
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
