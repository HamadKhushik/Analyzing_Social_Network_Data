package graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

import util.GraphLoader;


public class GirvanNewman extends CapGraph{

	//private HashMap<Integer, HashSet<Integer>> girvanMap; 
	private HashMap<GirvanNode, HashSet<GirvanNode>> girvanMap; 
	private int girvanNumVertices;  	// total number of vertices in graph
	private int girvanNumEdges;		// total number of edges in graph
	private HashMap<String, GirvanEdge> edges;

	private static final int INFINITY = (int) Double.POSITIVE_INFINITY;


	public GirvanNewman() {

		girvanMap = new HashMap<GirvanNode, HashSet<GirvanNode>>();
		edges = new HashMap<String, GirvanEdge>();
	}

	@Override
	public void addEdge(int from, int to) {
		// TODO Auto-generated method stub
		super.addEdge(from, to);
		GirvanEdge edge = new GirvanEdge(from, to);
		edges.put(edge.getId(), edge);
		
		GirvanNode source = this.getGirvanVertex(Integer.toString(from));
		GirvanNode dest = this.getGirvanVertex(Integer.toString(to));
		HashSet<GirvanNode> neighbours = girvanMap.get(source);
		if (neighbours.contains(dest)) {
			return;
		}
		if (!girvanMap.containsKey(dest)) {
			return;
		}
		
		neighbours.add(dest);
		girvanMap.put(source, neighbours);
		girvanNumEdges++;
	}
	
	@Override
	public void addVertex(int num) {
		// TODO Auto-generated method stub
		super.addVertex(num);
		GirvanNode vertex = new GirvanNode(num);
		if (girvanMap.containsKey(vertex)) {
			return;
		}
		HashSet<GirvanNode> neighbours = new HashSet<GirvanNode>();
		girvanMap.put(vertex, neighbours);
		girvanNumVertices++;
	}
	
	public GirvanNode getGirvanVertex(String id) {
		for (GirvanNode vertex : this.girvanMap.keySet()) {
			if (vertex.vertexId.equals(id)) {
				return vertex;
			}
		}
		return null;
	}
	
	public Set<GirvanNode> getGirvanVertices(){
		return girvanMap.keySet();
	}

	public HashMap getEdges() {
		return edges;
	}
	
	public int getGirvanNumVertices() {
		return girvanNumVertices;
	}

	public int getGirvanNumEdges() {
		return girvanNumEdges;
	}

	public void setGirvanNumEdges(int girvanNumEdges) {
		this.girvanNumEdges = girvanNumEdges;
	}
	
	public Set<GirvanNode> getGirvanNeighbours(GirvanNode vertex){
		return girvanMap.get(vertex);
	}
	
	public HashMap<GirvanNode, HashSet<GirvanNode>> exportGirvanGraph(){
		return girvanMap;
	}
		
	/**
	 * initialize the graph for each bfs run
	 */
	public void bfsInit() {
		
		for (GirvanNode currVertex : girvanMap.keySet()) {
			
			HashSet<GirvanNode> neighbours = girvanMap.get(currVertex);
			HashSet<GirvanNode> temp = new HashSet<GirvanNode>();
			for (GirvanNode neighbour : neighbours) {
				neighbour.pred.clear();
				neighbour.distance = INFINITY;
				neighbour.numShortestPaths = 0;
				temp.add(neighbour);
			}
			currVertex.pred.clear();
			currVertex.distance = INFINITY;
			currVertex.numShortestPaths = 0;
			girvanMap.put(currVertex, temp);
			
		}
	}
	
	/**
	 *  initialize source dependency of graph to 0.0
	 */
	private void sourceDependencyInit() {
		
		for (GirvanNode i : girvanMap.keySet()) {
			i.sourceDependency = 0.0;
			HashSet<GirvanNode> temp = new HashSet<GirvanNode>();
			for (GirvanNode w : girvanMap.get(i)) {
				w.sourceDependency = 0.0;
				temp.add(w);
			}
			girvanMap.put(i, temp);
		}
		
	}

	public HashMap<GirvanEdge, Double> brandes() {

		// Initialization
		Queue<GirvanNode> queue = new LinkedList<GirvanNode>();
		Stack<GirvanNode> stack = new Stack<GirvanNode>();
//		double[] dist = new double[this.getNumOfVertices()];
//		List<ArrayList<Integer>> pred = new ArrayList<ArrayList<Integer>>();
//		int[] numShortestPaths = new int[this.getNumOfVertices()];
//		double[] sourceDependency = new double[this.getNumOfVertices()];
//
//		double[]  vertexBetweenness = new double[this.getNumOfVertices()];
		HashMap<GirvanEdge, Double> edgeBetweenness = new HashMap<GirvanEdge, Double>();

		// setting initial values for variables
//		for (int w=0; w < vertexBetweenness.length; w++) {
//			vertexBetweenness[w] = 0;
//		}

		for (String id : edges.keySet()) {
			edgeBetweenness.put(edges.get(id), 0.0);
		}

		// for all vertices
		for (GirvanNode vertex : girvanMap.keySet()) {
			//GirvanNode vertex = this.getGirvanVertex("0");

			// initialize for each bfs run
			bfsInit();
			vertex.setDistance(0.0);
			vertex.setNumShortestPaths(1);
			queue.add(vertex);

			// bfs
			while (!queue.isEmpty()) {
				GirvanNode currVertex = queue.remove();
				stack.push(currVertex);

				for (GirvanNode i : girvanMap.get(currVertex)) {

					// path discovery
					if (i.distance == INFINITY) {
						i.distance = currVertex.distance + 1;
						queue.add(i);
					}

					// path counting
					if (i.distance == currVertex.distance + 1) {
						i.numShortestPaths += currVertex.numShortestPaths;
						i.pred.add(currVertex);
					}

				}
			}

			// vertex betweenness / back propagation
			sourceDependencyInit();
			
			while (!stack.isEmpty()) {
				GirvanNode currVertex = stack.pop();

				for(GirvanNode j : currVertex.pred) {
					double temp = (double) j.numShortestPaths / currVertex.numShortestPaths * (1 + currVertex.sourceDependency);
					j.sourceDependency = j.sourceDependency + temp;
					String edgeId = j + "->" + currVertex; 
					double betweennessValue = edgeBetweenness.get(edges.get(edgeId)) + temp;
					edgeBetweenness.put(edges.get(edgeId), betweennessValue);
				}
				
				if (currVertex != vertex) {
					currVertex.vertexBetweenness += currVertex.sourceDependency;
				}
			}
		}

		// checking values/testing
		System.out.println("***************************************************************************");
		for (GirvanEdge edge : edgeBetweenness.keySet()) {
			System.out.println("Edge Betweenness : " + edge + " Value = " + edgeBetweenness.get(edge));
		}
		System.out.println(edgeBetweenness.size());

		return edgeBetweenness;
	}
		
	@Override
	public List<Graph> getSCCs() {
		
		if (this.exportGirvanGraph().isEmpty()) {
			return null;
		}
		
		List<Graph> scc = new ArrayList<Graph>();
		
		// first iteration of dfs
		Stack<GirvanNode> vertices = new Stack<GirvanNode>();
		vertices.addAll(this.girvanMap.keySet());
		// perform first iteration of dfs, list not required- so passed as null
		Stack<GirvanNode> finished1 = dfs(this, vertices, null);
		
		// transpose of graph
		Graph tr = this.transpose();
		
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
	private Stack<GirvanNode> dfs (Graph g, Stack<GirvanNode> vertices, List<Graph> list){
		
		// initialization
		HashSet<GirvanNode> visited = new HashSet<GirvanNode>();
		Stack<GirvanNode> finished = new Stack<GirvanNode>();
		
		while (!vertices.empty()) {
			GirvanNode v = vertices.pop();
			if (!visited.contains(v)) {
				
				// if second iteration of dfs, pass in subGraph to be populated
				if (list != null) {
					Graph subGraph = new GirvanNewman();
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
	private void dfsVisit(Graph g, GirvanNode currVertex, HashSet<GirvanNode> visited, Stack<GirvanNode> finished, Graph subGraph){
		
		visited.add(currVertex);
		//g = (GirvanNewman) g;

		for (GirvanNode i : ((GirvanNewman) g).getGirvanNeighbours(currVertex)){
			if (!visited.contains(i)) {
				dfsVisit(g, i, visited, finished, subGraph);
			}
		}
		finished.push(currVertex);
		
		// if second iteration of dfs, populate the graph with strongly connected components
		if (subGraph != null) {
			subGraph.addVertex(Integer.parseInt(currVertex.getId()));
		}
	}
	
	/**
	 * @param method to transpose a graph, called from getScc() method
	 * @return  GirvanNewman graph object
	 */
	@Override
	protected Graph transpose() {
				
			GirvanNewman tr = new GirvanNewman();
			for (GirvanNode i : this.girvanMap.keySet()) {
				int iId = Integer.valueOf(i.getId());
				for (GirvanNode j : this.girvanMap.get(i)) {
					int jId = Integer.valueOf(j.getId());
					if (tr.exportGirvanGraph().containsKey(j) && tr.exportGirvanGraph().containsKey(i)) {
						tr.addEdge(jId, iId);
					}
					else {
						tr.addVertex(jId);
						tr.addVertex(iId);
						tr.addEdge(jId, iId); 
					}
				}
			}
			return tr;
		}
	
	//****************************************************************
	
	public class GirvanNode{
		
		private String vertexId;
		private double distance;
		private int numShortestPaths;
		private double sourceDependency;
		private double vertexBetweenness;
		private List<GirvanNode> pred;
		
		
		public GirvanNode(int vertex) {
			
			vertexId = String.valueOf(vertex);
			pred = new ArrayList<GirvanNode>();
		}
		
		public List<GirvanNode> getPred() {
			return pred;
		}

		public void setPred(List<GirvanNode> pred) {
			this.pred = pred;
		}

		public String getId() {
			return vertexId;
		}
		
		public double getDistance() {
			return distance;
		}

		public int getNumShortestPaths() {
			return numShortestPaths;
		}

		public void setNumShortestPaths(int numShortestPaths) {
			this.numShortestPaths = numShortestPaths;
		}

		public double getSourceDependency() {
			return sourceDependency;
		}

		public void setSourceDependency(double sourceDependency) {
			this.sourceDependency = sourceDependency;
		}

		public double getVertexBetweenness() {
			return vertexBetweenness;
		}

		public void setVertexBetweenness(double vertexBetweenness) {
			this.vertexBetweenness = vertexBetweenness;
		}

		public void setDistance(double distance) {
			this.distance = distance;
		}
		
		public String toString() {
			return vertexId;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result; // + getEnclosingInstance().hashCode();
			result = prime * result + ((vertexId == null) ? 0 : vertexId.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			GirvanNode other = (GirvanNode) obj;
//			if (!getEnclosingInstance().equals(other.getEnclosingInstance()))
//				return false;
			if (vertexId == null) {
				if (other.vertexId != null)
					return false;
			} else if (!vertexId.equals(other.vertexId))
				return false;
			return true;
		}

		private GirvanNewman getEnclosingInstance() {
			return GirvanNewman.this;
		}
		
		
		
	}
		
	
	// ***************************************************************

	public class GirvanEdge {

		private int source;
		private int destination;
		private String id;  // id is source - destination
		private double edgeBetweenness;

		public GirvanEdge(int source, int destination) {
			this.source = source;
			this.destination = destination;
			id = source + "->" + destination;
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
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getEnclosingInstance().hashCode();
			result = prime * result + ((id == null) ? 0 : id.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			GirvanEdge other = (GirvanEdge) obj;
			if (!getEnclosingInstance().equals(other.getEnclosingInstance()))
				return false;
			if (id == null) {
				if (other.id != null)
					return false;
			} else if (!id.equals(other.id))
				return false;
			return true;
		}

		private GirvanNewman getEnclosingInstance() {
			return GirvanNewman.this;
		}
	}
	
	// ***************************************************************
	
	public static void main(String[] args) {
		GirvanNewman girvan = new GirvanNewman();
		girvan.addVertex(0);
		girvan.addVertex(1);
		girvan.addVertex(2);
		girvan.addVertex(3);
		girvan.addVertex(4);
		girvan.addVertex(5);
		girvan.addVertex(6);
		
		girvan.addEdge(0, 1);
		girvan.addEdge(1, 0);
		girvan.addEdge(0, 3);
		girvan.addEdge(3, 0);
		girvan.addEdge(3, 4);
		girvan.addEdge(4, 3);
		girvan.addEdge(1, 2);
		girvan.addEdge(2, 1);
		girvan.addEdge(1, 4);
		girvan.addEdge(4, 1);
		girvan.addEdge(4, 5);
		girvan.addEdge(5, 4);
		girvan.addEdge(2, 5);
		girvan.addEdge(5, 2);
		
		GirvanNewman girvan2 = new GirvanNewman();
		GraphLoader.loadGraph(girvan2, "data/smallTest.txt");
		
		
		HashMap<GirvanEdge, Double> edgeBetweenness = girvan2.brandes();
		HashMap<GirvanNode, HashSet<GirvanNode>> graph = girvan2.exportGirvanGraph();
		System.out.println("Graph = " + graph);
		
		// test Scc()
//		girvan2 = new GirvanNewman();
//		GraphLoader.loadGraph(girvan2, "data/sccTestFile.txt");
//		List<Graph> scc = girvan2.getSCCs();
//		for (int i = 0; i < scc.size(); i++) {
//			System.out.println(scc.get(i));
//		}
		
		// transpose check
//		GirvanNewman tr = (GirvanNewman) girvan2.transpose();
//		System.out.println(tr instanceof GirvanNewman);
//		System.out.println("Girvan Graph: " + girvan2.exportGirvanGraph());
//		System.out.println("Transpose of Girvan Graph: " + ((GirvanNewman)tr).girvanMap);
//		
//		GirvanNode testNode = girvan2.new GirvanNode(65);
//		System.out.println(tr.exportGirvanGraph().containsKey(testNode));
//		System.out.println(girvan2.getGirvanNeighbours(testNode).hashCode());
//		
//		System.out.println("***********************************************************");
//		GirvanNode testNode2 = tr.getGirvanVertex("65");
//		System.out.println(testNode2.hashCode());
//		System.out.println(girvan2.exportGirvanGraph().containsKey(testNode2));
//		System.out.println(testNode2 instanceof GirvanNode);
		
//		HashSet<GirvanNode> test = new HashSet<GirvanNode>();
//		HashMap<GirvanNode, HashSet<GirvanNode>> map = new HashMap<GirvanNode, HashSet<GirvanNode>>();
//		GirvanNode testNode = girvan2.new GirvanNode(100);
//		GirvanNode testNode2 = girvan2.new GirvanNode(100);
//		map.put(testNode, test);
//		System.out.println(map.get(testNode2));
//		for (GirvanNode i : map.get(testNode)) {
//			System.out.println("test loop for null set");
//		}
		
//		GirvanNode vertex = girvan2.getGirvanVertex("65");
//		System.out.println("getNeighboursTest : " + girvan2.getGirvanNeighbours(vertex));
		
//		GirvanNewman tr = new GirvanNewman();
//		//for (GirvanNode i : this.girvanMap.keySet()) {
//		GirvanNode i = girvan2.getGirvanVertex("32");
//			int iId = Integer.valueOf(i.getId());
//			System.out.println("Transpose : " + iId + " String: " + i.getId());
//			for (GirvanNode j : girvan2.girvanMap.get(i)) {
//				int jId = Integer.valueOf(j.getId());
//				System.out.println("Transpose : " + jId + " String: " + j.getId());
//				if (tr.exportGirvanGraph().containsKey(j) && tr.exportGirvanGraph().containsKey(i)) {
//					tr.addEdge(jId, iId);
//					System.out.println("Hello");
//				}
//				else {
//					tr.addVertex(jId);
//					tr.addVertex(iId);
//					tr.addEdge(jId, iId);
//				}
//			}
//		//}
//			
//			System.out.println(tr.getGirvanVertex("44").hashCode());
//			System.out.println(girvan2.getGirvanVertex("44").hashCode());
//			
//			Graph tr2 = tr;
//			System.out.println(tr2);
//			System.out.println(((GirvanNewman) tr2).getGirvanVertex("44").hashCode());
		 
	}
}



