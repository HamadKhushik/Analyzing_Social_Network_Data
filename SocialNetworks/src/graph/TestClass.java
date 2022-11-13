package graph;

import java.util.HashMap;
import java.util.HashSet;

import graph.GirvanNewman.GirvanNode;

// test HashMap with reference variables 
// test updating reference variables in HashMap while looping through it
// test updating variables in a list while looping through it
public class TestClass {
	
	private HashMap<TestNode, HashSet<TestNode>> testMap;
	private static int count = 0;
	private int arrayId;
	
	public TestClass() {
		testMap = new HashMap<TestNode, HashSet<TestNode>>();
		arrayId = count++;
	}
	
	public int getArrayId() {
		return arrayId;
	}
	
	public static void main(String[] args) {
//		TestClass test = new TestClass();
//		TestNode hamad = test.new TestNode("Hamad","0");
//		TestNode khan = test.new TestNode("khan","1");
//		TestNode khushk = test.new TestNode("khushk","2");
//		
//		HashSet<TestNode> neighbours = new HashSet<TestNode>();
//		neighbours.add(khan);
//		neighbours.add(khushk);
//		test.testMap.put(hamad, neighbours);
//		System.out.println(test.testMap);
//		
//		HashSet<TestNode> neighbour = new HashSet<TestNode>();
//		neighbour.add(hamad);
//		neighbour.add(khan);
//		
//
//		test.testMap.put(khushk, neighbour);
//		
//		System.out.println(test.testMap);
//		hamad.name = "HAMMAD";
//		System.out.println(test.testMap);
		
		// test auto increment
		for (int i = 0; i < 10; i++) {
			TestClass curr = new TestClass();
			System.out.println(" loop count: " + i + " arrayId: " + curr.getArrayId());
		}
	}
	
	
	//************************************************************************
	
	class TestNode{
		private String name;
		private String id;
		
		public TestNode (String name, String id) {
			this.name = name;
			this.id = id;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			TestNode other = (TestNode) obj;
			if (!getEnclosingInstance().equals(other.getEnclosingInstance()))
				return false;
			if (name == null) {
				if (other.name != null)
					return false;
			} else if (name.equals(other.name))
				return false;
			return true;
		}
		
		private TestClass getEnclosingInstance() {
			return TestClass.this;
		}


		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getEnclosingInstance().hashCode();
			result = prime * result + ((name == null) ? 0 : name.hashCode());
			return result;
		}

		@Override
		public String toString() {
			// TODO Auto-generated method stub
			return this.name;
		}
		
		
		
		
	}

}
