/**
 * 
 */
package graph;

//import java.util.ArrayList;
import java.util.HashMap;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import util.GraphLoader;

/**
 * @author Your name here.
 * 
 * For the warm up assignment, you must implement your Graph in a class
 * named CapGraph.  Here is the stub file.
 *
 */
public class CapGraph implements Graph {

	/* (non-Javadoc)
	 * @see graph.Graph#addVertex(int)
	 */
	
	private HashMap<Integer, HashSet<Integer>> map; 
	private int numVertices;  	// total number of vertices in graph
	private int numEdges;		// total number of edges in graph
	
	public CapGraph() {
		
		map = new HashMap<Integer, HashSet<Integer>>();
		numVertices = 0;
		numEdges = 0;
	}
	
	@Override
	public void addVertex(int num) {
		// TODO Auto-generated method stub
		
		if (map.containsKey(num)) {
			System.out.println("The map already contains the vertex: " + num);
			return;
		}
		
		HashSet<Integer> neighbours = new HashSet<Integer>();
		map.put(num, neighbours);
		numVertices++;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#addEdge(int, int)
	 */
	@Override
	public void addEdge(int from, int to) {
		// TODO Auto-generated method stub
		
		HashSet<Integer> neighbours = map.get(from);
		if (neighbours.contains(to)) {
			System.out.println("This edge: " + from +" - " + to + " already exists");
			return;
		}
		if (!map.containsKey(to)) {
			System.out.println("Vertex: " + to + " does not exist");
			return;
		}
		
		neighbours.add(to);
		map.put(from, neighbours);
		numEdges++;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#getEgonet(int)
	 */
	@Override
	public Graph getEgonet(int center) {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#getSCCs()
	 */
	@Override
	public List<Graph> getSCCs() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#exportGraph()
	 */
	@Override
	public HashMap<Integer, HashSet<Integer>> exportGraph() {
		// TODO Auto-generated method stub
		return map;
	}
	
	/**
	 * @return
	 */
	public int getNumOfVertices() {
		return numVertices;
	}
	
	/**
	 * @return
	 */
	public int getNumOfEdges() {
		return numEdges;
	}
	
	/**
	 * @return the vertices in a graph
	 */
	public Set<Integer> getVertices() {
		return map.keySet();
	}
	
	/**
	 * @param vertex of the graph
	 * @return set of neighbouring vertices 
	 */
	public HashSet<Integer> getNeighbours(int vertex){
		return map.get(vertex);
	}
	
	public static void main(String[] args) {
		CapGraph graph = new CapGraph();
		GraphLoader.loadGraph(graph, "data/small_test_graph.txt");
		HashMap<Integer, HashSet<Integer>> result = graph.exportGraph();
		System.out.println(result);
		
	}

}
