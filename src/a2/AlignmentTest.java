package a2;

import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;

public class AlignmentTest {

	String[] cs7x5 = { "AGCTG", "AGCAG", "CAGCC", "CACAG", "GCAGC", "GATAA",
			"CAGGC" }; // orig
	String[] cs5x3 = { "AAG", "GCC", "CGC", "AGC", "GCT" }; //
	String[] cs5x10 = { "GCGCGCGAAG", "GCCGGGCGCG", "AACGCAGCGC", "ATTAGCAAAA", "TTTTGCTAAG" }; //
	Sequence[] dna1 = new Sequence[cs7x5.length];
	Sequence[] dna2 = new Sequence[cs5x3.length];
	Sequence[] dna3 = new Sequence[cs5x3.length];
	PerfMeter perf = new PerfMeter();
	Alignment testme1, testme2;
	Alignment testme3, testme4;
	Alignment testme5, testme6;

	public AlignmentTest() {
		for (int i = 0; i < cs7x5.length; i++)
			dna1[i] = new Sequence("S" + Integer.toString((i + 1)),
					cs7x5[i].toCharArray());
		testme1 = new Alignment(perf, dna1);
		testme2 = new Alignment(perf, dna1, true, 0);
		for (int i = 0; i < cs5x3.length; i++)
			dna2[i] = new Sequence("S" + Integer.toString((i + 1)),
					cs5x3[i].toCharArray());
		testme3 = new Alignment(perf, dna2);
		testme4 = new Alignment(perf, dna2, true, 0);
		for (int i = 0; i < cs5x10.length; i++)
			dna3[i] = new Sequence("S" + Integer.toString((i + 1)),
					cs5x10[i].toCharArray());
		testme5 = new Alignment(perf, dna3);
		testme6 = new Alignment(perf, dna3, true, 0);
	}

	/* Part A tests below */

	@Test
	public void testGetConsensus() {
		int[][] profile = {
		/* A */{ 2, 0, 0, 0, 9, 0, 2, 2, 0 },
		/* C */{ 1, 1, 1, 1, 0, 0, 2, 2, 0 },
		/* G */{ 1, 2, 0, 2, 0, 0, 2, 2, 0 },
		/* T */{ 1, 0, 0, 3, 0, 0, 2, 2, 0 }};
		int[] c = testme1.getConsensus(profile);
		assertEquals(1, c[0]);
		assertEquals(3, c[1]);
		assertEquals(2, c[2]);
		assertEquals(4, c[3]);
		assertEquals(1, c[4]);
	}

	@Test
	public void testGetProfile() {
		int[] aln = { 1, 1, 0, 0, 2, 2, 0 };
		// -AGCTG-
		// -AGCAG-
		// CAGCC--
		// CACAG--
		// --GCAGC
		// --GATAA
		// CAGGC--
		int[][] profile = testme1.getProfile(aln);
		assertEquals(5, profile[0][1]); // A pos 1
		assertEquals(4, profile[1][3]); // C pos 3
		assertEquals(3, profile[2][5]); // G pos 5
		assertEquals(0, profile[3][6]); // T pos 6
	}
//
//	@Test
//	public void testGetScore() {
//		int[][] profile = {
//		/* A */{ 2, 0, 0, 0, 9, 0, 0, 0, 0 },
//		/* C */{ 1, 1, 1, 1, 0, 4, 0, 0, 0 },
//		/* G */{ 1, 2, 0, 2, 0, 0, 0, 0, 0 },
//		/* T */{ 1, 0, 0, 3, 0, 0, 5, 0, 0 } };
//		int[] consens = { 1, 3, 2, 4, 1, 2, 4, 3, 4 };
//		int sc = testme1.getScore(profile, consens);
//		assertEquals(2 + 2 + 1 + 3 + 9 + 4 + 5, sc);
//	}

	/* Part B tests below */

//	@Test
//	public void testFindAlignment1() {
//		int[] s = new int[cs5x3.length];
//		for (int i = 0; i < s.length; i++)
//			s[i] = -1;
//		Alignment.AlignmentScore score = testme3.findAlignment(s);
//		// AAG-- 3
//		// --GCC 3
//		// -CGC- 2
//		// -AGC- 3
//		// --GCT 2
//		assertEquals(13, score.actual);
//		score = testme4.findAlignment(s);
//		// AAG-- 3
//		// --GCC 3
//		// -CGC- 2
//		// -AGC- 3
//		// -AGC- 3 [reverse strand]
//		assertEquals(14, score.actual);
//	}
//
//	@Test
//	public void testFindAlignment2() {
//		int[] s = new int[cs5x10.length];
//		for (int i = 0; i < s.length; i++)
//			s[i] = -1;
//		Alignment.AlignmentScore score = testme5.findAlignment(s);
//		//---------GCGCGCGAAG
//		//------GCCGGGCGCG---
//		//---AACGCAGCGC------
//		//ATTAGCAAAA---------
//		//TTTTGCTAAG---------
//		assertEquals(40, score.actual);
//		score = testme6.findAlignment(s);
//		//-GCGCGCGAAG--------
//		//CGCGCCCGGC---------
//		//---GCGCTGCGTT------
//		//---------ATTAGCAAAA
//		//---------CTTAGCAAAA
//		assertEquals(42, score.actual);
//	}
//
//	public static int[] allocateStartState(int n) {
//		int[] s = new int[n];
//		for (int i = 0; i < n; i ++)
//			s[i] = -1;
//		return s;
//	}
//
//	@Test
//	/**
//	 * Test on combinations of random sequences. Both original and reverse strand.
//	 * Compare scores with those of Consensus (consensus scores should be equivalent).
//	 * Note that this test takes a long time for larger values of n (>7) and t (>7).
//	 */
//	public void testFindAlignmentVsConsensus() {
//		char[] chars = {'A','C','G','T'};
//		Alignment aln_orig, aln_rev;
//		Consensus cons_orig, cons_rev;
//		Random rand = new Random();
//		for (int n = 2; n < 7; n ++) { // test different sequence widths
//			for (int t = 2; t < 7; t ++) { // test different sequence numbers
//				// System.out.println("====== Test n="+n+" t="+t+" ======");
//				Sequence[] dna4 = new Sequence[t];
//				for (int i = 0; i < t; i ++) {
//					StringBuffer sb = new StringBuffer();
//					for (int j = 0; j < n; j ++)
//						sb.append(chars[rand.nextInt(4)]);
//					dna4[i] = new Sequence("S"+i, sb.toString().toCharArray());
//				}
//				PerfMeter pm_aln_orig  = new PerfMeter();
//				PerfMeter pm_aln_rev   = new PerfMeter();
//				PerfMeter pm_cons_orig = new PerfMeter();
//				PerfMeter pm_cons_rev  = new PerfMeter();
//				int W = 2 * n - 1;
//				aln_orig  = new Alignment(pm_aln_orig,  dna4, false, 0);
//				aln_rev   = new Alignment(pm_aln_rev,   dna4, true,  0);
//				cons_orig = new Consensus(pm_cons_orig, dna4, false, 0);
//				cons_rev  = new Consensus(pm_cons_rev,  dna4, true,  0);
//
//				Alignment.AlignmentScore sc_aln_orig = aln_orig.findAlignment(allocateStartState(t));
//				// A2App.printAlignment(sc_aln_orig.path, dna4, W);
//				Consensus.ConsensusScore sc_cons_orig = cons_orig.findConsensus(allocateStartState(W), 0);
//				// A2App.printAlignment(cons_orig.getAlignment(sc_cons_orig.path), dna4, W);
//				assertEquals(sc_cons_orig.actual, sc_aln_orig.actual);
//
//				Alignment.AlignmentScore sc_aln_rev = aln_rev.findAlignment(allocateStartState(t));
//				// A2App.printAlignment(sc_aln_rev.path, dna4, W);
//				Consensus.ConsensusScore sc_cons_rev = cons_rev.findConsensus(allocateStartState(W), 0);
//				// A2App.printAlignment(cons_rev.getAlignment(sc_cons_rev.path), dna4, W);
//				assertEquals(sc_cons_rev.actual, sc_aln_rev.actual);
//			}
//		}
//	}
}
