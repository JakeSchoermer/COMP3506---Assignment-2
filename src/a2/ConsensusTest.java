package a2;

import static org.junit.Assert.*;

import org.junit.Test;

public class ConsensusTest {

	String[] cs7x5 = {"AGCTG","AGCAG","CAGCC","CACAG","GCAGC","GATAA","CAGGC"}; // orig
	String[] cs5x3 = {"AAG", "GCC", "CGC", "AGC", "GCT"}; // 
	Sequence[] dna1 = new Sequence[cs7x5.length];
	Sequence[] dna2 = new Sequence[cs5x3.length];
	PerfMeter perf = new PerfMeter();
	Consensus testme1, testme2;
	Consensus testme3, testme4;
	
	public ConsensusTest() {
		for (int i = 0; i < cs7x5.length; i ++)
			dna1[i] = new Sequence("S"+Integer.toString((i+1)), cs7x5[i].toCharArray());
		testme1 = new Consensus(perf, dna1, false, 0); // original strand only
		testme2 = new Consensus(perf, dna1, true,  0); // reverse
		for (int i = 0; i < cs5x3.length; i ++)
			dna2[i] = new Sequence("S"+Integer.toString((i+1)), cs5x3[i].toCharArray());
		testme3 = new Consensus(perf, dna2);
		testme4 = new Consensus(perf, dna2, true, 0);
	}

	@Test
	public void testGetLevel() {
		int[][] paths = {{-1}, {1, -1}, {1, 2}, {1, 2, 3, -1, -1}, {1, 2, 3, 4, -1}, {1, 2, 4, 1, 1}};
		for (int i = 0; i < paths.length; i ++) {
			int level1 = testme1.getLevel(paths[i]);
			assertEquals(i, level1);
		}
	}

	@Test
	public void testGetScore() {
		int[] mypat1 = {4,1,3,2,4,1,1,1,1}; 	
		// 'TAGCTAAAA' 
		int sc1o = testme1.getScore(mypat1); 	
		// ' AGCTG   '	(4)
		// ' AGCAG   '	(3)
		// 'CAGCC    ' 	(3)
		// '   CACAG ' 	(2)
		// '  GCAGC  '	(2)
		// '  GATAA  '	(4)
		// 'CAGGC    '	(2)  SUM = 20
		assertEquals(20, sc1o);
		// 						'TTTTAGCTA' (reverse)
		int sc1r = testme2.getScore(mypat1); 	
		// ' AGCTG   '	(4)	or	'    AGCTG'	(4) 
		// ' AGCAG   '	(3)		'    AGCAG'	(3)
		// 'CAGCC    ' 	(3)		'   CAGCC '	(3)
		// 'CACAG    ' 	(1)		'   CACAG '	(1)
		// '  GCAGC  '	(2)		'  GCAGC  '	(3)
		// '  GATAA  '	(4)		' GATAA   '	(2)
		// 'CAGGC    '	(2)  	'   CAGGC ' (2)	SUM = 21
		assertEquals(21, sc1r);
		int[] mypat3 = {3,3,3,3,3,1,4,1,1}; 	
		// 'GGGGGATAA' 		or 	'TTATCCCCC' 
		int sc3o = testme1.getScore(mypat3);	
		// '   AGCTG '	(2)	or	'  AGCTG  '	(2) 
		// 'AGCAG    '	(2)		'  AGCAG  '	(2)
		// 'CAGCC    ' 	(1)		'    CAGCC'	(3)
		// '    CACAG' 	(2)		;    CACAG'	(2)
		// 'GCAGC    '	(2)		'GCAGC    '	(2)
		// '    GATAA'	(5)		' GATAA   '	(2)
		// 'CAGGC    '	(2) 16 	' CAGGC   ' (2)	SUM = 18
		assertEquals(16, sc3o);
		int sc3r = testme2.getScore(mypat3);
		assertEquals(18, sc3r);
	}

	@Test
	public void testGetAlignment() {
		int[] mypat0 = {3,3,2,1,2,2,2,1,1}; 		
		// 'GGCACCCGA' or 'TCCGGTGCC' (reverse)
		int[] aln0o = testme1.getAlignment(mypat0);
		assertEquals(aln0o[2], aln0o[6]); 			
		//    CAGCC vs         CAGCC vs
		//    CAGGC		    CAGGC
		int[] mypat1 = {3,1,4,1,1,3,2,4,3}; 		
		// 'GATAAGCTG' or 	'CAGTTTATC' (reverse)
		int[] aln1o = testme1.getAlignment(mypat1);
		assertEquals(0, aln1o[5]);		 			
		//  CACAG (2x3) vs 	 CAGCC  (3) vs
		//  GATAA (5)		 GATAA(RC) (5)
		int[] aln1r = testme2.getAlignment(mypat0);
		assertEquals(true, aln1o[3] != aln1r[3]); 	
		assertEquals(aln1o[5], aln1r[5]); 			
	}

	@Test
	public void testFindConsensus() {
		int[] prefix = new int[dna2[0].getLength() * 2 - 1];
		for (int i = 0; i < prefix.length; i ++)
			prefix[i] = -1;
		Consensus.ConsensusScore score = testme3.findConsensus(prefix, 0);
		// AAG--	3
		// --GCC	3
		// -CGC-	2
		// -AGC-	3
		// --GCT	2
		assertEquals(13, score.actual);
		score = testme4.findConsensus(prefix, 0);
		// AAG--	3
		// --GCC	3
		// -CGC-	2
		// -AGC-	3
		// -AGC-*  	3	[reverse strand]
		assertEquals(14, score.actual);
	}

}
