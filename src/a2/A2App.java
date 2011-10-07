package a2;

import java.io.IOException;

import a2.Alignment.AlignmentScore;
import a2.Consensus.ConsensusScore;

public class A2App {

	/**
	 * A command line application that accepts a number of parameters. 
	 * -f <filename> (load sequences from FASTA file) 
	 * -m alignment|consensus (the search method to use) 
	 * -t <#> (limit number of sequences used; only effective for alignment) 
	 * -w <#> (limit alignment or consensus width to this value) 
	 * -r (use reverse complement) 
	 * See usage message for more information.
	 * 
	 * @param args command line parameters
	 */
	public static void main(String[] args) {

		String file = null; // FASTA file
		String method = null; // search method
		Sequence[] seqs = null; // loaded sequences
		int N = -1; // length of each sequence (set after sequence data has be
					// read)
		int T = -1; // limit number of sequences (default is as long as the data
					// indicates)
		int W = -1; // limit width of consensus/alignment (default is N * 2 - 1
					// where N is the length of sequence)
		boolean reverse = false; // use reverse complement

		// parse the parameters
		for (int i = 0; i < args.length; i++) {
			if (args[i].charAt(0) == '-') // option
			{
				switch (args[i].charAt(1)) {
                case 'f':
					if (i + 1 < args.length)
						file = args[++i];
					break;
				case 'm':
					if (i + 1 < args.length)
						method = args[++i];
					break;
				case 't':
					if (i + 1 < args.length)
						T = Integer.parseInt(args[++i]);
					break;
				case 'w':
					if (i + 1 < args.length)
						W = Integer.parseInt(args[++i]);
					break;
				case 'r':
					reverse = true;
					break;
				default:
					System.err.println("Unknown option \"-" + args[i].charAt(1)
							+ "\"");
				}
			}
		}

		if (file == null || method == null) {
			System.err
					.println("Usage: A2App -f <sequence-filename> -m alignment|consensus [-t <limit-sequences>] [-w <limit-width>] [-r]");
			System.exit(1);
		}

		if (file != null) {
			try {
				seqs = Sequence.readFile(file); // read a FASTA file with
												// sequences
				if (seqs.length > 0) {
					N = seqs[0].getLength();
					if (T < 1 || T > seqs.length)
						T = seqs.length;
					else {
						Sequence[] tmp = seqs;
						seqs = new Sequence[T];
						for (int i = 0; i < T; i++)
							seqs[i] = tmp[i];
					}
					if (W < N || W > N * 2 - 1)
						W = N * 2 - 1; // the widest consensus we can search for
										// (any two sequences overlap by at
										// least one)
				}
			} catch (IOException e) {
				System.err.println(e.getMessage());
				System.exit(2);
			}
		}

		if (method != null) {
			PerfMeter perf = new PerfMeter();
			if (method.startsWith("c")) { // consensus
				Consensus problem = new Consensus(perf, seqs, reverse, W);
				int[] start = new int[W];
				for (int i = 0; i < start.length; i++)
					start[i] = -1; // marker for end of solution
				// Start searching
				ConsensusScore score = problem.findConsensus(start, 0); 
				printConsensus(score.path);
				printAlignment(problem.getAlignment(score.path), seqs, W);
				System.out
						.println(String.format("Score %d (%4.1f%%)",
								score.actual, score.actual * 100.0
										/ (seqs.length * N)));
			} else if (method.startsWith("a")) { // alignment
				Alignment problem = new Alignment(perf, seqs, reverse, W);
				int[] start = new int[T];
				for (int i = 0; i < start.length; i++)
					start[i] = -1; // marker for end of solution
				// Start searching
				AlignmentScore score = problem.findAlignment(start); 
				printAlignment(score.path, seqs, W);
				printConsensus(problem.getConsensus(problem.getProfile(score.path)));
				System.out.println(String.format("Score %d (%4.1f%%)",
						score.actual, score.actual * 100.0 / (T * N)));
			}
			perf.printReport(); // performance report
		}
	}

	/**
	 * Print the consensus sequence
	 * 
	 * @param consensus
	 */
	public static void printConsensus(int[] consensus) {
		System.out.print('[');
		for (int i = 0; i < consensus.length; i++)
			System.out.print(Sequence.toChar(consensus[i]));
		System.out.println(']');
	}

	/**
	 * Print the alignment
	 * 
	 * @param aln
	 *            the alignment specified by offset indices
	 * @param dna
	 *            the aligned sequences
	 * @param W
	 *            the width of the alignment
	 */
	public static void printAlignment(int[] aln, Sequence[] dna, int W) {
		for (int i = 0; i < aln.length; i++) {
			int N = dna[0].getLength();
			if (aln[i] == -1)
				break;
			boolean strand = aln[i] < W - N + 1;
			int offset = strand ? aln[i] : aln[i] - (W - N + 1);
			System.out.print('\'');
			for (int j = 0; j < offset; j++)
				System.out.print(' ');
			System.out.print(dna[i].getSymbolChars(strand));
			for (int j = offset + N; j < W; j++)
				System.out.print(' ');
			System.out.println('\'');
		}
	}

}
