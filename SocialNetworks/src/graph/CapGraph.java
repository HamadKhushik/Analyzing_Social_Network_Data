/**
 * 
 */
package graph;

import java.util.ArrayList;
import java.util.Arrays;
//import java.util.ArrayList;
import java.util.HashMap;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

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
		
		if (center < 0) {
			throw new IllegalArgumentException("negative numbers are not allowed") ;
		}
		
		Graph egoGraph = initEgoGraph(center, map.get(center));
		HashSet<Integer> neighbours = map.get(center);		
		for (int i : neighbours) {
			HashSet<Integer> currNeighbours = map.get(i);
			for (int j : currNeighbours) {
				if (neighbours.contains(j)) {
					egoGraph.addEdge(i, j);
				}
			}
		}
		return egoGraph;
	}
	
	/**
	 * Helper function for getEgonet to initialize the graph with center and its neighbours
	 * @param neighbours - neighbours for vertex 'center'
	 * @param center - vertex 'Center'
	 * @return returns a graph object initialized with center and its neighbours
	 */
	private Graph initEgoGraph(int center, HashSet<Integer> neighbours) {
		
		Graph result = new CapGraph();
		result.addVertex(center);
		
		for (int i : neighbours) {
			result.addVertex(i);
			result.addEdge(center, i);
			// add edges in both directions
			result.addEdge(i, center);
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#getSCCs()
	 */
	@Override
	public List<Graph> getSCCs() {
		// TODO Auto-generated method stub
		
		List<Graph> scc = new ArrayList<Graph>();
		
		// first iteration of dfs
		Stack<Integer> vertices = new Stack<Integer>();
		vertices.addAll(map.keySet());
		Stack<Integer> dfs1 = dfs(null, vertices);
		System.out.println("Stack after dfs iteration " + dfs1);
		
		// transpose of graph
		Graph tr = transpose();
		
		// second iteration of dfs
		vertices.clear();
		vertices.addAll(((CapGraph) tr).getVertices());
		dfs(scc, vertices);
		
		
		
		return scc;
	}
	
	/**
	 * @param g graph
	 * @return list of graphs which are strongly connected components
	 */
	private Stack<Integer> dfs (List<Graph> g, Stack<Integer> vertices){
		
		HashSet<Integer> visited = new HashSet<Integer>();
		Stack<Integer> finished = new Stack<Integer>();
		
		while (!vertices.empty()) {
			int v = vertices.pop();
			System.out.println("dfs - pop " + v);
			if (!visited.contains(v)) {
				if (g != null) {
					Graph subGraph = new CapGraph();
					dfsVisit(subGraph, v, visited, finished);
				}
				else {
					dfsVisit(null, v, visited, finished);
				}
			}
		}

		return finished;
	}
	
	private Graph dfsVisit(Graph g, int v, HashSet<Integer> visited, Stack<Integer> finished){
		
		visited.add(v);
		for (int i : getNeighbours(v)){
			if (!visited.contains(i)) {
				dfsVisit(g, i, visited, finished);
			}
		}
		finished.push(v);
		System.out.println("dfs - visit push " + v);
		if (g != null) {
			g.addVertex(v);
		}
		return g;
	}
	
	/*
	 * private HashMap<Integer, HashSet<Integer>> transpose(HashMap<Integer,
	 * HashSet<Integer>> graph){
	 * 
	 * HashMap<Integer, HashSet<Integer>> tr = new HashMap<Integer,
	 * HashSet<Integer>>(); //Graph tr = new CapGraph();
	 * 
	 * for (int i : graph.keySet()) { for (int j : map.get(i)) { if
	 * (tr.containsKey(j)) { HashSet<Integer> currNeighbours = tr.get(j);
	 * currNeighbours.add(i); tr.put(j, currNeighbours); } else { tr.put(j, new
	 * HashSet<Integer>(Arrays.asList(i))); } } }
	 * System.out.println("Transpose of Graph: " + tr); return tr; }
	 */
	
	/**
	 * @param graph overloaded method to return the transpose of graph
	 * @return
	 */
	private Graph transpose () {
		
		Graph tr = new CapGraph();
		for (int i : exportGraph().keySet()) {
			for (int j : getNeighbours(i)) {
				if (tr.exportGraph().containsKey(j) && tr.exportGraph().containsKey(i)) {
					tr.addEdge(j, i);
				}
				else {
					tr.addVertex(j);
					tr.addVertex(i);
					tr.addEdge(j, i);
				}
			}
		}
		System.out.println("Overloaded transpose: " + tr.exportGraph());
		return tr;
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
	
	public String toString() {
		return map.toString();
	}
	
	
	public static void main(String[] args) {
		CapGraph graph = new CapGraph();
		GraphLoader.loadGraph(graph, "data/egoNetTest.txt");
		HashMap<Integer, HashSet<Integer>> result = graph.exportGraph();
		System.out.println(result);
		Graph egoNet = graph.getEgonet(50);
		System.out.println(egoNet.exportGraph());
		
		//test Transpose of graph
		CapGraph testTr = new CapGraph();
		GraphLoader.loadGraph(testTr, "data/sccTestFile.txt");
		System.out.println("\r\n" + testTr.exportGraph());
		//System.out.println(testTr.transpose(testTr.exportGraph()));
		testTr.transpose();
		
		// test dfs
		Stack<Integer> vertices = new Stack<Integer>();
		vertices.addAll(testTr.exportGraph().keySet());
		System.out.println("\r\n dfs first iteration: " + testTr.dfs(null, vertices));
		
		// test SCC
		List<Graph> scc= new ArrayList<Graph>();
		Graph testScc = new CapGraph();
		GraphLoader.loadGraph(testScc, "data/sccTestFile.txt");
		scc = testScc.getSCCs();
		for (int i = 0; i < scc.size(); i++){
			System.out.println(scc.get(i));
		}
		
		
		// test stack
		HashSet<Integer> test1 = new HashSet<Integer>(Arrays.asList(1,2,3,4));
		Stack<Integer> test = new Stack<Integer>();
		test.addAll(test1);
		System.out.println(test);
		for (int i = 0; i < 4; i++) {
			System.out.println(test.pop());
		}
	}

}
