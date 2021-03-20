package graph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

public class GirvanNewman2 extends CapGraph{
	
	
	private HashMap<Integer, HashSet<Integer>> girvanMap; 
	private int girvanNumVertices;  	// total number of vertices in graph
	private int girvanNumEdges;		// total number of edges in graph
	private HashMap<String, GirvanEdge> edges;
	
	private static final int INFINITY = (int) Double.POSITIVE_INFINITY;
	
	
	public GirvanNewman2() {
		
		girvanMap = new HashMap<Integer, HashSet<Integer>>();
		edges = new HashMap<String, GirvanEdge>();
	}
	
	@Override
	public void addEdge(int from, int to) {
		// TODO Auto-generated method stub
		super.addEdge(from, to);
		GirvanEdge edge = new GirvanEdge(from, to);
		edges.put(edge.getId(), edge);
	}
	
	public HashMap getEdges() {
		return edges;
	}

	/**
	 * recursive bfs method to check connectivity
	 * performs a breadth first search on the graph
	 * @return visited set
	 */
	public HashSet<Integer> bfRecursive(Graph g, int vertex){
		
		Queue<Integer> queue = new LinkedList<Integer>();
		HashSet<Integer> visited = new HashSet<Integer>();
		queue.add(vertex);
		while (!queue.isEmpty()) {
			int curr = queue.remove();
			if (!visited.contains(curr)) {
				bfsVisit(g, curr, visited, queue);
			}
		}
		
		return visited;
	}
	
	private void bfsVisit(Graph g, int vertex, HashSet<Integer> visited, Queue<Integer> queue ){
		
		visited.add(vertex);
		for (int i : ((CapGraph) g).getNeighbours(vertex)) {
			if (!visited.contains(i)) {
				bfsVisit(g, i, visited, queue);
			}
		}
	}
	
	/**
	 * performs bfs on a graph and returns a HashSet of visited vertices
	 * -> used in checkConnectivityBfs() method to check the connectivity of graph
	 * 
	 * @param g -> graph to perform bfs on
	 * @param vertex -> the vertex to start the bfs from
	 * @return -> HashSet of all the visited vertices
	 */
	public HashSet<Integer> bfs(Graph g, int vertex){

		Queue<Integer> queue = new LinkedList<Integer>();
		HashSet<Integer> visited = new HashSet<Integer>();

		queue.add(vertex);
		visited.add(vertex);

		while (!queue.isEmpty()) {
			int curr = queue.remove();
			for (int i : ((CapGraph) g).getNeighbours(curr)) {
				if (!visited.contains(i)) {
					visited.add(i);
					queue.add(i);
				}
			}
		}

		return visited;

	}

	/**
	 * finds the shortest paths in a graph based on one vertex
	 * @param g -> graph of the object
	 * @param vertex -> vertex to find the shortest paths from
	 * @param paths -> an array to update the number of paths from the vertex to all the other vertices
	 * @param dist -> distance of the shortest path from vertex
	 * @return -> not important at the minute
	 */
	public HashSet<Integer> bfs(int vertex, int[] paths, int[] dist){

		paths[vertex] = 1;
		dist[vertex] = 0;
		
		
		Queue<Integer> queue = new LinkedList<Integer>();
		HashSet<Integer> visited = new HashSet<Integer>();

		queue.add(vertex);
		visited.add(vertex);

		while (!queue.isEmpty()) {
			int curr = queue.remove();

			for (int i : this.getNeighbours(curr)) {

				if (!visited.contains(i)) {
					visited.add(i);
					queue.add(i);
				}
				// if this is the first shortest path found, update dist[] and path[]
				if (dist[i] > dist[curr] + 1) {
					dist[i] = dist[curr] + 1;
					paths[i] = paths[curr];
				}

				// if more shortest paths are found
				else if (dist[i] == dist[curr] + 1) {
					paths[i] += paths[curr];
				}
			}
		}

		return visited;

	}
	
	
	public int[] shortestPaths(int vertex) {
		
		int numOfVertices = this.getNumOfVertices();
		
		int[] paths = new int[numOfVertices];
		int[] dist = new int[numOfVertices];
		
		Arrays.fill(paths, 0);
		Arrays.fill(dist, INFINITY);
		
		bfs(vertex, paths, dist);
		
		return paths;
	}
	
	/**
	 * checks if the graph is connected using the recursive bfs (and bfsVisit) methods
	 * @param vertex -> vertex to check the connectivity from
	 * @return true or false depending on if the graph is connected
	 * 
	 */
	public boolean checkConnectivity(int vertex) {
		
		Graph tr =  (Graph) this.transpose();
		
		Set visited = bfRecursive(this, vertex);
		Set visited2 = bfRecursive(tr, vertex);
		Set vertices = this.getVertices();
		System.out.println("Check Connectivity GirvanNewman");
		System.out.println("Visited " + visited);
		System.out.println("Vertices " + vertices);
		System.out.println("Visited2 " + visited2);
		
		if (visited.containsAll(vertices)) {
			return true;
		}
		return false;
	}
	
	/**
	 * checks connectivity of graph using simple bfs. Two iterations of bfs, one of which is on the inverse of graph
	 * @param vertex -> to check connectivity from
	 * @return
	 */
	public boolean checkConnectivityBfs(int vertex) {
		
		Graph tr = (CapGraph) this.transpose();
		
		Set visited1 = bfs(this, vertex);
		Set visited2 = bfs(tr, vertex);
		Set vertices = this.getVertices();
		System.out.println("Check Connectivity GirvanNewman checkConnectivityBfs()");
		System.out.println("Visited " + visited1);
		System.out.println("Vertices " + vertices);
		System.out.println("Visited2 " + visited2);
		
		if (visited1.containsAll(vertices)) {
			return true;
		} else if (visited1.equals(visited2)) {
			return false; 
		}
		
		return false;
	}
	
	public double[] brandes(int vertex1) {
		
		// Initialization
		Queue<Integer> queue = new LinkedList<Integer>();
		Stack<Integer> stack = new Stack<Integer>();
		Stack<Integer> stack2 = new Stack<Integer>();  // for edge betweenness run
		double[] dist = new double[this.getNumOfVertices()];
		List<ArrayList<Integer>> pred = new ArrayList<ArrayList<Integer>>();
		int[] numShortestPaths = new int[this.getNumOfVertices()];
		double[] sourceDependency = new double[this.getNumOfVertices()];

		double[]  vertexBetweenness = new double[this.getNumOfVertices()];
		HashMap<GirvanEdge, Double> edgeBetweenness = new HashMap<GirvanEdge, Double>();
		double[] totalSourceDependency = new double[this.getNumOfVertices()];
		
		for (int i = 0; i < totalSourceDependency.length; i++) {
			totalSourceDependency[i] = 0;
		}

		// setting initial values for variables
		for (int w=0; w < vertexBetweenness.length; w++) {
			vertexBetweenness[w] = 0;
		}
		
		for (String id : edges.keySet()) {
			edgeBetweenness.put(edges.get(id), 0.0);
		}
		
		// for all vertices
		for (int vertex : this.getVertices()) {
			
			// initialize for each bfs run
			for(int i=0; i < this.getNumOfVertices(); i++) {
				pred.add(i, new ArrayList<Integer>());
				dist[i] = INFINITY;
				numShortestPaths[i] = 0;
			}

			dist[vertex] = 0;
			numShortestPaths[vertex] = 1;
			queue.add(vertex);
			
			// bfs
			while (!queue.isEmpty()) {
				int currVertex = queue.remove();
				stack.push(currVertex);
				stack2.push(currVertex);

				for (int i : this.getNeighbours(currVertex)) {
					// path discovery
					if (dist[i] == INFINITY) {
						dist[i] = dist[currVertex] + 1;
						queue.add(i);
					}
					// path counting
					if (dist[i] == dist[currVertex] + 1) {
						numShortestPaths[i] += numShortestPaths[currVertex];
						//pred.get(i).add(currVertex);
						ArrayList<Integer> temp = pred.get(i);
						temp.add(currVertex);
						pred.set(i, temp);
					}
				}
			}

			// accumulation - back propagation
			for (int i = 0; i <sourceDependency.length; i++ ) {
				sourceDependency[i] = 0;
			}
			
			// vertex betweenness
			while (!stack.isEmpty()) {
				int currVertex = stack.pop();

				for(int j : pred.get(currVertex)) {
					double temp =  (double)numShortestPaths[j] / numShortestPaths[currVertex] * (1 + sourceDependency[currVertex]);
					sourceDependency[j] = sourceDependency[j] + (double)numShortestPaths[j] / numShortestPaths[currVertex] * (1 + sourceDependency[currVertex]);
					double tempC = (double)numShortestPaths[j]/numShortestPaths[currVertex] * (1 + sourceDependency[currVertex]);
					System.out.println("Edge betweenness temp = " + tempC);
					//GirvanEdge currEdge = new GirvanEdge(j, currVertex);
					String edgeId = j + "-" + currVertex; 
					double betweennessValue = edgeBetweenness.get(edges.get(edgeId)) + tempC;
					edgeBetweenness.put(edges.get(edgeId), betweennessValue);
				}
				if (currVertex != vertex) {
					vertexBetweenness[currVertex] += sourceDependency[currVertex];
				}
			}
			
			for (int i = 0; i < sourceDependency.length; i++) {
				totalSourceDependency[i] += sourceDependency[i];
			}
		}

//		// edge betweenness
//		while (!stack2.isEmpty()) {
//			int currVertex = stack2.pop();
//			
//			for (int j : pred.get(currVertex)) {
//				double temp = (double)numShortestPaths[j]/numShortestPaths[currVertex] * (1 + sourceDependency[currVertex]);
//				System.out.println("Edge betweenness temp = " + temp);
//				GirvanEdge currEdge = new GirvanEdge(j, currVertex);
//				double betweennessValue = edgeBetweenness.get(currEdge) + temp;
//				edgeBetweenness.put(currEdge, betweennessValue);
//				
//			}
//		}
		
		
		// checking values/testing
		for (int i = 0; i < sourceDependency.length; i++) {
			System.out.println("source dependency : " + vertexBetweenness[i]);
		}

		System.out.println("***************************************************************************");
		for (GirvanEdge edge : edgeBetweenness.keySet()) {
			System.out.println("Edge Betweenness : " + edge + " Value = " + edgeBetweenness.get(edge));
		}
		System.out.println(edgeBetweenness.size());
		
		return totalSourceDependency;
	}
	
	
	public class GirvanEdge {
		
		private int source;
		private int destination;
		private String id;  // id is source - destination
		
		public GirvanEdge(int source, int destination) {
			this.source = source;
			this.destination = destination;
			id = source + "-" + destination;
		}
		
		public String getId() {
			return this.id;
		}
		
		public String getId(int source, int destination) {
			
			return source + "-" + destination;
		}
		
		public String toString() {
			return id;
		}
		
		@Override
		public boolean equals(Object o) {
			// TODO Auto-generated method stub
			
			if (this == o) {
				return true;
			} 
			
			if (!(o instanceof GirvanEdge)) {
				return false;
			}
			
			GirvanEdge edge = (GirvanEdge) o;
			if (this.source == edge.source && this.destination == edge.destination && this.id.equals(edge.id)) {
				return true;
			}
			return false;
		}
		
		

	}
	
	
	public static void main(String[] args) {
		GirvanNewman2 girvan = new GirvanNewman2();
		girvan.addVertex(0);
		girvan.addVertex(1);
		girvan.addVertex(2);
		girvan.addVertex(3);
		girvan.addVertex(4);
		girvan.addVertex(5);
		girvan.addVertex(6);
		
		girvan.addEdge(0, 1);
		girvan.addEdge(0, 2);
		girvan.addEdge(1, 2);
		girvan.addEdge(1, 3);
		girvan.addEdge(2, 3);
		girvan.addEdge(3, 4);
		girvan.addEdge(3, 5);
		girvan.addEdge(4, 6);
		girvan.addEdge(5, 6);
		
		int vertex = 0;
		
		// check connectivity of a graph
//		if (girvan.checkConnectivityBfs(vertex)) {
//			System.out.println("Graph is connected");
//		} else {
//			System.out.println("Graph is not connected");
//		}
		
		//int[] path = girvan.shortestPaths(0);
		//System.out.println("Shortest paths from vertex 0: " + Arrays.toString(path));
		
		System.out.println("*************************************************************");
		
		System.out.println("Brandes betweenness: ");
		//double[] between = girvan.brandes(0);
		
//		for (int i = 0; i < between.length; i++) {
//			System.out.println(between[i]);
//		}	
		
		GirvanNewman2 girvan2 = new GirvanNewman2();
		girvan2.addVertex(0);
		girvan2.addVertex(1);
		girvan2.addVertex(2);
		girvan2.addVertex(3);
		girvan2.addVertex(4);
		girvan2.addVertex(5);
		girvan2.addVertex(6);
		
		girvan2.addEdge(0, 1);
		girvan2.addEdge(1, 0);
		girvan2.addEdge(0, 3);
		girvan2.addEdge(3, 0);
		girvan2.addEdge(3, 4);
		girvan2.addEdge(4, 3);
		girvan2.addEdge(1, 2);
		girvan2.addEdge(2, 1);
		girvan2.addEdge(1, 4);
		girvan2.addEdge(4, 1);
		girvan2.addEdge(4, 5);
		girvan2.addEdge(5, 4);
		girvan2.addEdge(2, 5);
		girvan2.addEdge(5, 2);
		
		double[] betweenness = girvan2.brandes(0);
	}
}
