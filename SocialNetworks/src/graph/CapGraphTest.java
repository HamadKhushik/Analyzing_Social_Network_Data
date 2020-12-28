/**
 * 
 */
package graph;

import static org.junit.Assert.*;
import java.util.*;

import org.junit.Before;
import org.junit.Test;

/**
 * @author Hammad
 *
 */
public class CapGraphTest {

	/**
	 * @throws java.lang.Exception
	 */
	
	CapGraph empty;
	CapGraph test;
	
	@Before
	public void setUp() throws Exception {
		
		empty = new CapGraph();
		test = new CapGraph();
		test.addVertex(32);
		test.addVertex(50);
		test.addVertex(44);
		test.addVertex(32);
		
		test.addEdge(32, 50);
		test.addEdge(50, 32);
		test.addEdge(32, 44);
		test.addEdge(44, 32);
		test.addEdge(50, 44);
		test.addEdge(44, 50);
		test.addEdge(32, 50);
	}

	/**
	 * Test method for {@link graph.CapGraph#addVertex(int)}.
	 */
	@Test
	public void testAddVertex() {
		
		assertEquals("check total number of vertices", 3, test.getNumOfVertices());
		assertEquals("check total number of vertices in an empty graph", 0, empty.getNumOfVertices());
		HashSet<Integer> expected = new HashSet<Integer>(Arrays.asList(32,50,44));
		assertEquals("check the vertices", expected, test.getVertices());
	}

	/**
	 * Test method for {@link graph.CapGraph#addEdge(int, int)}.
	 */
	@Test
	public void testAddEdge() {
		
		assertEquals("check total number of edges", 6, test.getNumOfEdges());
		assertEquals("check total number of edges in an empty graph", 0, empty.getNumOfEdges());
		HashSet<Integer> expected = new HashSet<Integer>(Arrays.asList(50, 44));
		assertEquals("check the edges for vertex 32", expected, test.getNeighbours(32));
	}

	/**
	 * Test method for {@link graph.CapGraph#getEgonet(int)}.
	 */
	@Test
	public void testGetEgonet() {
		fail("Not yet implemented");
	}

	/**
	 * Test method for {@link graph.CapGraph#getSCCs()}.
	 */
	@Test
	public void testGetSCCs() {
		fail("Not yet implemented");
	}

}
