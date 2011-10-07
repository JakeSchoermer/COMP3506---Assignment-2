package a2;

/**
 * Class for monitoring the performance of search algorithms
 * Do NOT change any of this code. 
 */
public final class PerfMeter {

	private final long start;
	private long end = -1;
	private static int CNT_FIND = 0;
	private static int CNT_EXIT_BREAK = 0;
	private static int CNT_EXIT_LEAF = 0;
	private static int CNT_EXIT_PROPAGATE = 0;
 
	public PerfMeter() {
		start = System.currentTimeMillis();
	}
	
	public int countFind() {
		if (end > 0)
			throw new RuntimeException("Already ended");
		return ++ CNT_FIND;
	}
	
	public int countBreak() {
		if (end > 0)
			throw new RuntimeException("Already ended");
		return ++ CNT_EXIT_BREAK;
	}
	
	public int countLeaf() {
		if (end > 0)
			throw new RuntimeException("Already ended");
		return ++ CNT_EXIT_LEAF;
	}
	
	public int countPropagate() {
		if (end > 0)
			throw new RuntimeException("Already ended");
		return ++ CNT_EXIT_PROPAGATE;
	}
	
	public void exit() {
		if (end > 0)
			throw new RuntimeException("Already ended");
		end = System.currentTimeMillis();
	}
	
	public void printReport() {
		if (end < 0) 
			end = System.currentTimeMillis();
		System.out.println("Started at "+new java.util.Date(start));
		System.out.println("Finished at "+new java.util.Date(end));
		System.out.println(String.format("Time elapsed: \t%9.2f secs", (end - start) / 1000.0));
		System.out.println("#ENTRY\t \t"+CNT_FIND);
		System.out.println("#EXIT by");
		System.out.println("  \tleaf \t"+CNT_EXIT_LEAF);
		System.out.println("  \tbreak\t"+CNT_EXIT_BREAK);
		System.out.println("  \tpropg\t"+CNT_EXIT_PROPAGATE);
	}
}

