package a2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

/**
 * A class for representing biological sequence data
 */
public class Sequence {

	public final static char[] alpha = {'A', 'C', 'G', 'T'}; // valid characters, mapping to index + 1
	public static int  toIndex(char ch)  { for (int i = 0; i < alpha.length; i ++) { if (alpha[i] == ch) return i + 1; } return 0; } 
	public static char toChar(int index) { if (index > 0 && index <= alpha.length) return alpha[index - 1]; return ' '; }
	
	private final String name; 			// name of sequence
	private final int[]  seq_orig; 		// original sequence of symbols represented by index
	private final int[]  seq_reverse; 	// sequence of symbols of reverse strand represented by index

	/**
	 * Constructs a DNA sequence instance.
	 * 
	 * @param name
	 *            the name of the sequence
	 * @param seq_orig
	 *            the character string representing the symbols of the sequence
	 *            (must be instances of alphabet).
	 * @see DNAAlphabet
	 */
	public Sequence(String name, char[] string) {
		this.name = name;
		// convert to index and check that the sequence is valid
		this.seq_orig = new int[string.length];
		this.seq_reverse = new int[string.length];
		for (int i = 0; i < seq_orig.length; i++) {
			this.seq_orig[i] = Sequence.toIndex(string[i]);
			this.seq_reverse[seq_orig.length - 1 - i] = alpha.length + 1 - Sequence.toIndex(string[i]);
			if (seq_orig[i] == 0 || seq_reverse[i] == 5)
				throw new SequenceRuntimeException("Invalid character in sequence: "+string[i]); 
		}
	}

	/**
	 * Retrieves the index of the symbol found at the specified position.
	 * 
	 * @param position
	 *            position of symbol 0..n-1 where n is the length of the sequence
	 * @param strand
	 * 			  the original (True) or reverse strand (False)
	 * @return the index of the symbol
	 * @throws SequenceRuntimeException
	 *             if an invalid position is given
	 */
	public int getSymbolIndex(int position, boolean strand) {
		if (position >= 0 && position < seq_orig.length)
			return strand ? seq_orig[position] : seq_reverse[position];
		else
			throw new SequenceRuntimeException(this,
					"Attempt to retrieve invalid index " + position + " in \"" + name + "\"");
	}

	/**
	 * Retrieves the indices of all the symbols in the sequence 0..n-1 where n
	 * is the length of the sequence
	 * 
	 * @return the indices
	 */
	public int[] getSymbolIndices(boolean strand) {
		return strand ? seq_orig : seq_reverse;
	}

	/**
	 * Retrieves the character representation of all the symbols in the sequence
	 * 0..n-1 where n is the length of the sequence
	 * 
	 * @return the printable characters
	 */
	public char[] getSymbolChars(boolean strand) {
		char[] str = new char[seq_orig.length];
		for (int i = 0; i < seq_orig.length; i ++) 
			str[i] = strand ? toChar(seq_orig[i]) : toChar(seq_reverse[i]);
		return str;
	}

	/**
	 * Retrieve the printable character of the symbol found at the specified
	 * position.
	 * 
	 * @param position
	 *            position of symbol 0..n-1 where n is the length of the
	 *            sequence
	 * @return the character of the symbol
	 */
	public char getSymbolChar(int position, boolean strand) {
		return toChar(getSymbolIndex(position, strand));
	}

	/**
	 * Retrieve the length of the DNA sequence.
	 * 
	 * @return the length (number of symbols)
	 */
	public int getLength() {
		return seq_orig.length;
	}

	/**
	 * Printable representation of sequence
	 */
	public String toString() {
		return name + " (" + seq_orig.length + ")";
	}

	/**
	 * Reads DNA sequences from a file on the FASTA standard format
	 * 
	 * @param filename
	 *            the name of the file
	 * @return an array of instance of {@link #DNASequence}
	 * @throws IOException
	 *             if the file operation fails
	 */
	public static Sequence[] readFile(String filename) throws IOException {
		List<Sequence> seqs = new ArrayList<Sequence>();
		BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
		// buffer variables to hold recently read data
		String name = null;
		StringBuffer buf = null;

		int row = 0;
		String line = br.readLine();
		while (line != null) {
			row++;
			line = line.trim(); // remove any spaces, tabs etc at the ends
			if (line.startsWith(">")) {
				if (buf != null) // there is data in the buffer, we need to
				// store it before processing the new entry
				{
					try {
						seqs.add(new Sequence(name, buf.toString().toCharArray()));
					} catch (SequenceRuntimeException e) {
						System.err.println("Ignored " + name + ": " + e.getMessage());
					}
					buf = null;
					name = null;
				}
				try {
					StringTokenizer stok = new StringTokenizer(line, " \t");
					name = stok.nextToken().substring(1);
				} catch (NoSuchElementException e) {
					throw new RuntimeException("Invalid format in file " + filename + " at row " + row);
				}
				buf = new StringBuffer();
			} else {
				if (buf != null) {
					buf.append(line);
				}
			}
			line = br.readLine();
		}
		if (buf != null) // there is data in the buffer, we need to store it
						 // before processing the new entry
		{
			try {
				seqs.add(new Sequence(name, buf.toString().toCharArray()));
			} catch (SequenceRuntimeException e) {
				System.err.println("Ignored " + name + ": " + e.getMessage());
			}
			buf = null;
			name = null;
		}
		Sequence[] all = new Sequence[seqs.size()];
		seqs.toArray(all);
		return all;
	}

	/**
	 * Example application that simply loads a FASTA file and prints out the
	 * sequences in it.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		for (int i = 0; i < args.length; i++) {
			try {
				Sequence[] seqs = Sequence.readFile(args[i]);
				System.out.println("Read " + seqs.length + " sequences from " + args[i]);
				for (Sequence s : seqs)
					System.out.println(s);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
	}
}

/**
 * A class that holds information about sequence-related runtime exceptions
 */
class SequenceRuntimeException extends RuntimeException {
	private static final long serialVersionUID = 1L;
	public Sequence s = null;

	public SequenceRuntimeException(Sequence s, String msg) {
		super(msg);
		this.s = s;
	}

	public SequenceRuntimeException(String msg) {
		super(msg);
	}
}
