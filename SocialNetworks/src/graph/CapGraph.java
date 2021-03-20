/**
 * 
 */
package graph;

import java.util.ArrayList;
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
		
		// if vertex already exists
		if (map.containsKey(num)) {
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
			return;
		}
		if (!map.containsKey(to)) {
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
	 * Helper function for getEgonet to initialize the graph with given vertex and its neighbours
	 * @param neighbours - neighbours for given vertex
	 * @param center - given vertex 'Center'
	 * @return returns a graph object initialized with center and its neighbours
	 */
	private Graph initEgoGraph(int center, HashSet<Integer> neighbours) {
		
		Graph result = new CapGraph();
		result.addVertex(center);
		
		for (int i : neighbours) {
			result.addVertex(i);
			
			// add edges in both directions
			result.addEdge(center, i);
			result.addEdge(i, center);
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#getSCCs() returns strongly connected components
	 */
	@Override
	public List<Graph> getSCCs() {
		// TODO Auto-generated method stub
		
		if (this.exportGraph().isEmpty()) {
			return null;
		}
		
		List<Graph> scc = new ArrayList<Graph>();
		
		// first iteration of dfs
		Stack<Integer> vertices = new Stack<Integer>();
		vertices.addAll(map.keySet());
		// perform first iteration of dfs, list not required- so passed as null
		Stack<Integer> finished1 = dfs(this, vertices, null);
		
		// transpose of graph
		Graph tr = transpose();
		
		// second iteration of dfs, with list passed in to be populated
		dfs(tr, finished1, scc);
		
		return scc;
	}
	
	
	/**
	 * @param g -> graph object to be explored
	 * @param vertices -> stack of vertices to be explored in the graph
	 * @param list -> list of scc(Strongly connected components) to be populated by the method
	 * @return stack of vertices in the order algorithm finished visiting them, only required in the first iteration of dfs
	 */
	private Stack<Integer> dfs (Graph g, Stack<Integer> vertices, List<Graph> list){
		
		// initialization
		HashSet<Integer> visited = new HashSet<Integer>();
		Stack<Integer> finished = new Stack<Integer>();
		
		while (!vertices.empty()) {
			int v = vertices.pop();
			if (!visited.contains(v)) {
				
				// if second iteration of dfs, pass in subGraph to be populated
				if (list != null) {
					Graph subGraph = new CapGraph();
					dfsVisit(g, v, visited, finished, subGraph);
					list.add(subGraph);
				}
				// if first iteration of dfs, subGraph not required
				else {
					dfsVisit(g, v, visited, finished, null);
				}
			}
		}
		return finished;
	}
	
	/**
	 * called from dfs() method. it explores all the neighbours of current vertex in the graph
	 * and populates a graph(subGraph) with the strongly connected components 
	 * @param g -> graph object to be explored
	 * @param currVertex -> current vertex to be explored
	 * @param visited -> set of visited nodes/vertices
	 * @param finished -> stack of completely visited nodes in the order they were finished, populated for the dfs() method. used in first iteration of dfs
	 * @param subGraph -> graph of strongly connected components, populated for the dfs() method. used in second iteration of dfs
	 */
	private void dfsVisit(Graph g, int currVertex, HashSet<Integer> visited, Stack<Integer> finished, Graph subGraph){
		
		visited.add(currVertex);
		for (int i : ((CapGraph) g).getNeighbours(currVertex)){
			if (!visited.contains(i)) {
				dfsVisit(g, i, visited, finished, subGraph);
			}
		}
		finished.push(currVertex);
		
		// if second iteration of dfs, populate the graph with strongly connected components
		if (subGraph != null) {
			subGraph.addVertex(currVertex);
		}
	}
	
	/**
	 * @param method to transpose a graph, called from getScc() method
	 * @return
	 */
	protected Graph transpose () {
		
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
		//System.out.println("\r\n dfs first iteration: " + testTr.dfs(null, vertices));
		
		// test SCC
		List<Graph> scc= new ArrayList<Graph>();
		Graph testScc = new CapGraph();
		GraphLoader.loadGraph(testScc, "data/sccTestFile.txt");
		scc = testScc.getSCCs();
		for (int i = 0; i < scc.size(); i++){
			System.out.println("Scc" + scc.get(i));
		}
	}

}
