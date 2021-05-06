package graph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
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
	private HashMap<GirvanNode, Integer> arrayIdMap;  // parallel map to keep track of nodeIndex in adjacency matrix for modularization 
	private int girvanNumVertices;  	// total number of vertices in graph
	private int girvanNumEdges;		// total number of edges in graph
	private HashMap<String, GirvanEdge> edges;	// edges in the graph. id is 'source -> destination'

	private static final int INFINITY = (int) Double.POSITIVE_INFINITY;


	public GirvanNewman() {

		girvanMap = new HashMap<GirvanNode, HashSet<GirvanNode>>();
		edges = new HashMap<String, GirvanEdge>();
		arrayIdMap = new HashMap<GirvanNode, Integer>();
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
		source.incrementDegree();
		//dest.incrementDegree();
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

	/** prints any two dimensional double array
	 * 
	 */
	private void printDoubleMatrix(double[][] matrix) {
		
		if (matrix == null) {
			return;
		}

		for (int row = 0; row < matrix.length; row++) {
			for (int column = 0; column < matrix[row].length; column++) {
				System.out.printf("%.3f  ", matrix[row][column]);
			}
			System.out.println();
		}
	}
	
	/** prints the community structure/structures found
	 * @param optimalMap
	 */
	private void printCommunities(HashMap<List<Graph>, Double> optimalMap) {
		
		System.out.println("Community structures found with maximum Modularity are: " + optimalMap.size());
		int count = 1;
		
		for (List<Graph> community : optimalMap.keySet()) {
			System.out.println("************************************************");
			System.out.println("Community Structure #  " + count++ + " : ");
			System.out.printf("Modularity Score for this structure is: " + "%.4f",optimalMap.get(community));
			System.out.println("\r\nCommunities detected within this structure are: " + community.size());
			for (int i = 0; i < community.size(); i++) {
				System.out.println("Community No: " + (i+1) + " -> " + community.get(i));
			}
		}
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
		HashMap<GirvanEdge, Double> edgeBetweenness = new HashMap<GirvanEdge, Double>();

		// setting initial values for variables
		for (String id : edges.keySet()) {
			edgeBetweenness.put(edges.get(id), 0.0);
		}

		// for all vertices
		for (GirvanNode vertex : girvanMap.keySet()) {

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
//		System.out.println("***************************************************************************");
//		for (GirvanEdge edge : edgeBetweenness.keySet()) {
//			System.out.println("Edge Betweenness : " + edge + " Value = " + edgeBetweenness.get(edge));
//		}
//		System.out.println("brandes(): total number of edges in original graph: " + edgeBetweenness.size());

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
			if (this.girvanMap.get(i).size() == 0) {
				tr.addVertex(iId);
			}
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

	/** removes the edges with highest betweenness, if tie - removes all
	 * @param eBetweenenss : betweenness value of al the edges
	 */
	private void removeEdges (HashMap<GirvanEdge, Double> eBetweenness) {

		HashMap<GirvanEdge, Double> max = maxBetweenness(eBetweenness);

		for (GirvanEdge curr : max.keySet()) {
			int source = curr.getSource();
			int dest = curr.getDestination();
			GirvanNode sourceNode = getGirvanVertex(Integer.toString(source));
			GirvanNode destNode = getGirvanVertex(Integer.toString(dest));
			girvanMap.get(sourceNode).remove(destNode);
			girvanNumEdges--;
		}
//		System.out.println("Graph after removing edges: " + girvanMap);
//		System.out.println("Remove Edges -> edges remaining = " + girvanNumEdges);
	}


	/** returns edges with max betweenness
	 * @param eBetweeness edge betweenness for the graph
	 * @return
	 */
	private HashMap<GirvanEdge, Double> maxBetweenness(HashMap<GirvanEdge, Double> eBetweenness) {

		HashMap<GirvanEdge, Double> maxBetweenness = new HashMap<GirvanEdge, Double>();
//		double temp = 0.0;
		double THRESHOLD = 1e-6;
//
//		for (GirvanEdge i : eBetweenness.keySet()) {
//
//			int comp = Double.compare(eBetweenness.get(i), temp);
//
//			if (Math.abs(eBetweenness.get(i) - temp) < THRESHOLD) {
//				maxBetweenness.put(i, temp);
//			}
//			if (Math.abs(eBetweenness.get(i) - temp) > THRESHOLD) {
//				if (comp > 0) {
//					maxBetweenness.clear();
//					temp = eBetweenness.get(i);
//					maxBetweenness.put(i, temp);
//				}
//			} 
//		}
		
		double maxValueInMap = Collections.max(eBetweenness.values());
		
		for (GirvanEdge i : eBetweenness.keySet()) {
			
			if (Math.abs(eBetweenness.get(i) - maxValueInMap) <  THRESHOLD) {
				maxBetweenness.put(i, eBetweenness.get(i));
			}
		}
		
		return maxBetweenness;
	}

	/** find communities until no edges are left
	 * 
	 */
	private HashMap<List<Graph>, Double> getCommunities() {

		List<Graph> communities = new ArrayList<Graph>();
		this.assignArrayIds();
		
		// to map the modularity score of the community structure
		HashMap<List<Graph>, Double> communityModularityMap = new HashMap<List<Graph>, Double>();
//		double [][] probabilityMat = this.getProbabilityMatrix();
		double[][] modularityMatrix = this.getModularityMatrix();
		double modularity = -1.0;
		// get total number of edges before removing edges -> required for calculating modularity
		int numOfEdges = this.getGirvanNumEdges();
		// to save top three modularity values while removing edges
		Double [] topThree = new Double[3];
		Arrays.fill(topThree, 0.0);
		
		// initial community structure and its modularity score
		System.out.println("getCommunities() -> getSCCs() = " + getSCCs() );
		System.out.println(this.getModularity(getSCCs(), modularityMatrix, numOfEdges));
		communityModularityMap.put(getSCCs(), this.getModularity(getSCCs(), modularityMatrix, numOfEdges));
		
		while (this.getGirvanNumEdges() > 0) {
			HashMap<GirvanEdge, Double> edgeBetweenness = this.brandes();
			removeEdges(edgeBetweenness);
//			System.out.println("getCommunities(): no of edges = " + this.getGirvanNumEdges());
			communities = getSCCs();
			modularityMatrix = this.getModularityMatrix();
			modularity = this.getModularity(communities, modularityMatrix, numOfEdges);
			System.out.println("Modularity = " + modularity);

			if (communityModularityMap.size() >= 3) {
//				System.out.println("Find communities >= 3 ");
				if (modularity > Collections.min( communityModularityMap.values())) {
					topThreeModularities(communityModularityMap, Collections.min(communityModularityMap.values()));
					System.out.println("Hello===");
//					System.out.println(modularity);
//					System.out.println(communities);
//					System.out.println("=================================================");
					communityModularityMap.put(communities, modularity);
				}
			}
			else {
				communityModularityMap.put(communities, modularity);
				System.out.println("Find Communities: size less than 3");				
			}

		}
		
		// get the community structure/structures with maximum modularity
//		communityModularityMap = this.getOptimalCommunities(communityModularityMap);
//		this.printCommunities(communityModularityMap);
		
		return communityModularityMap;
	}
	
	/** this function maintains communities with top three modularity scores, the communities with top 3 modularity scores will then be returned
	 * @param map -> map of community structure and modularity score
	 * @param minValue -> minimum modularity score in the map (the one that needs to be removed)
	 */
	public void topThreeModularities(HashMap<List<Graph>, Double> map, double minValue) { 

		List<Graph> toRemove = new ArrayList<Graph>();
		for (List<Graph> k : map.keySet()) {
			if (map.get(k) == minValue) {
				toRemove = k;
			}
		}
		map.remove(toRemove, minValue);
	}
	
	/** finds the community structures/structures with maximum modularity value and returns its map
	 * @param communityModularityMap -> mapping of community structures found and their corresponding modularity value
	 * @return optimalCommunityMap -> mapping of community structure/structures with maximum modularity value
	 */
	public HashMap<List<Graph>, Double> getOptimalCommunities(){

		HashMap<List<Graph>, Double> optimalCommunitiesMap = new HashMap<List<Graph>, Double>();
		double max = -1.0;
		double THRESHOLD = 1e-6;  // tolerance for double calculation
		
		HashMap<List<Graph>, Double> communityModularityMap = this.getCommunities();
		
		for (List<Graph> communities : communityModularityMap.keySet()) {

			int comp = Double.compare(communityModularityMap.get(communities), max);

			if (Math.abs(communityModularityMap.get(communities) - max) < THRESHOLD) {
				optimalCommunitiesMap.put(communities, communityModularityMap.get(communities));
			}

			if (Math.abs(communityModularityMap.get(communities) - max) > THRESHOLD) {
				if (comp > 0) {
					optimalCommunitiesMap.clear();
					optimalCommunitiesMap.put(communities, communityModularityMap.get(communities));
					max = communityModularityMap.get(communities);
				}
			}
		}
		return optimalCommunitiesMap;
	}

	/** METHOD NOT USED 
	 * when we remove all the edges, we return list of graph with individual vertices as individual communities
	 * as there are no edges. every vertex is a community
	 * @return list<Graph> with every vertex as individual community
	 */
	private List<Graph> getIndividualCommunities(){

		List<Graph> iCommunity = new ArrayList<Graph>();
		for (GirvanNode node : getGirvanVertices()) {
			Graph curr = new GirvanNewman();
			curr.addVertex(Integer.parseInt(node.getId()));
			iCommunity.add(curr);
		}
		return iCommunity;
	}

	/**
	 * assign arrayIds for every node in the graph for the Adjacency matrix
	 * dataset dont always start from '0', so we have to assign ids for the graph nodes to be compatible
	 */
	private void assignArrayIds() {

		int id = 0;
		for (GirvanNode curr : this.getGirvanVertices()) {
			arrayIdMap.put(curr, id++);
		}
	}

	/**calculates the Modularity matrix i-e (AdjacencyMatrix - ProbabilityMatrix) and returns it
	 * 
	 * @return returns the modularity matrix
	 */
	private double[][] getModularityMatrix() {

		double[][] modularityMatrix = new double[this.getGirvanNumVertices()][this.getGirvanNumVertices()];
		double[][] adjacencyMatrix = this.getAdjacencyMatrix();
		double[][] probabilityMatrix = this.getProbabilityMatrix();
		
		if (adjacencyMatrix == null || probabilityMatrix == null) {
			return null;
		}
		
		for (int row = 0; row < modularityMatrix.length; row++) {
			for (int column = 0; column < modularityMatrix[row].length; column++) {
				modularityMatrix[row][column] = adjacencyMatrix[row][column] - probabilityMatrix[row][column];
			}
		}
		return modularityMatrix;
	}

	/** generates the adjacency matrix for the graph
	 * @return
	 */
	private double[][] getAdjacencyMatrix(){

		// create a n x n adjacency matrix
		double[][] adjacencyMatrix = new double[this.getGirvanNumVertices()][this.getGirvanNumVertices()];

		for (GirvanNode i : this.getGirvanVertices()) {

			for (GirvanNode w : this.getGirvanNeighbours(i)) {

				adjacencyMatrix[arrayIdMap.get(i)][arrayIdMap.get(w)] = 1.0;
			}
		}
		return adjacencyMatrix;
	}

	/** probability of edge (1->2) = (degree of '1' * degree of '2') / 2* total number of edges
	 * @return return the probability matrix of the graph
	 */
	private double[][] getProbabilityMatrix() {

		double[][] probabilityMatrix = new double [this.getGirvanNumVertices()][this.getGirvanNumVertices()];
		double numEdges = this.getGirvanNumEdges();
//		System.out.println("getProbability");
		
		for (GirvanNode i : this.getGirvanVertices()) {
			for (GirvanNode w : this.getGirvanVertices()){
//				System.out.println("getProbabilityMatrix() -> " + i.getId() + " " + i.getDegree());
//				System.out.println("getProbabilityMatrix() -> " + w.getId() + " " + w.getDegree());
				probabilityMatrix[arrayIdMap.get(i)][arrayIdMap.get(w)] = (double)(i.getDegree() * w.getDegree()) / numEdges; // (2*this.getGirvanNumEdges());
			}
		}
		return probabilityMatrix;
	}
	
	/** calculates the total modularity of the graph by computing the modularity of individual communities and adding them up
	 * @param communities -> communities returned by eliminating the edges
	 * @param modularityMatrix -> modularity matrix of the graph
	 * @return total modularity of the graph
	 */
	private double getModularity(List<Graph> communities, double[][] modularityMatrix, int numOfEdges) {
		
		double totalModularity = 0.0;
		
		for (Graph community : communities) {
			
			Set<GirvanNode> vertices = ((GirvanNewman) community).getGirvanVertices();
//			System.out.println("getModularity -> : no of communities" + communities);
			
			for (GirvanNode node : vertices) {
//				System.out.println("getModularity -> no of vertices");

				for (GirvanNode neighbour : vertices) {
				//for (GirvanNode neighbour : ((GirvanNewman) this).getGirvanNeighbours(node)) {
					
					// if 'node' and 'neighbour' are neighbours in given community
					// then calculate the modularity
					//if (vertices.contains(neighbour)) {
//						System.out.println(arrayIdMap.get(node) + "      " + arrayIdMap.get(neighbour));
//						System.out.println("getModularity -> " + (double) modularityMatrix[arrayIdMap.get(node)][arrayIdMap.get(neighbour)]);
						totalModularity += (double) modularityMatrix[arrayIdMap.get(node)][arrayIdMap.get(neighbour)];
					//}
				}
			}
		}
		
		totalModularity = (double) totalModularity/(2*numOfEdges);
		
		return totalModularity;
	}
	
	

	//****************************************************************

	public class GirvanNode{

		private String vertexId;
		private double distance;
		private int numShortestPaths;
		private double sourceDependency;
		private double vertexBetweenness;
		private List<GirvanNode> pred;
		private int degree;


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

		public int getDegree() {
			return degree;
		}

		//		public void setDegree(int degree) {
		//			this.degree = degree;
		//		}

		public void incrementDegree() {
			this.degree++;
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

		public int getSource() {
			return source;
		}

		public void setSource(int source) {
			this.source = source;
		}

		public int getDestination() {
			return destination;
		}

		public void setDestination(int destination) {
			this.destination = destination;
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
		GraphLoader.loadGraph(girvan2, "data/small_test_graph.txt");


		//HashMap<GirvanEdge, Double> edgeBetweenness = girvan2.brandes();
		HashMap<GirvanNode, HashSet<GirvanNode>> graph = girvan2.exportGirvanGraph();
		System.out.println("Graph = " + graph);
		System.out.println("-------------------------------------------");
		girvan2.assignArrayIds();
		System.out.println(Arrays.deepToString(girvan2.getModularityMatrix()));
		System.out.println(Arrays.deepToString(girvan2.getProbabilityMatrix()));
		System.out.println(Arrays.deepToString(girvan2.getAdjacencyMatrix()));
		System.out.println(girvan2.getSCCs());
		System.out.printf("%.4f" , girvan2.getModularity(girvan2.getSCCs(), girvan2.getModularityMatrix(), girvan2.getGirvanNumEdges()));

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


		// test maxBetweenness()

		//HashMap<GirvanEdge, Double> between = girvan2.(edgeBetweenness);
		//girvan2.removeEdges(edgeBetweenness);
		//System.out.println("==========================================");
		//System.out.println(between);

		// test find communities
//		List<Graph> communities = girvan2.findCommunities(4);
//		for (int i = 0; i < communities.size(); i++) {
//			System.out.println(((GirvanNewman) communities.get(i)).getGirvanVertices());
//		}

		// test arrayId, AdjacencyMatrix and probability matrix
//		girvan2.assignArrayIds();
//		double[][] adjacencyMatrix = girvan2.getAdjacencyMatrix();
//		System.out.println("-----------------------------------------------------------");
//
//		double [][] probabilityMatrix = girvan2.getProbabilityMatrix();
//		double [][] modularityMatrix = girvan2.getModularityMartix();
//		System.out.println("Adjacency Matrix");
//		girvan2.printDoubleMatrix(adjacencyMatrix);
//		System.out.println();
//		System.out.println("Probability Matrix");
//		girvan2.printDoubleMatrix(probabilityMatrix);
//		System.out.println();
//		System.out.println("Modularity Matrix");
//		girvan2.printDoubleMatrix(modularityMatrix);
//
//		System.out.println("Array Ids: ");
//		System.out.println(girvan2.arrayIdMap);

		// test degree
		//		for (GirvanNode curr : girvan2.getGirvanVertices()) {
		//			System.out.println("Degree of vertex: " + curr + " is = " + curr.getDegree());
		//		}
		
		
		// test getOptimalCommunities
//System.out.println("main method");
//System.out.println("communities getSccs: " + girvan2.getSCCs());
//HashMap<List<Graph>, Double> communityMap = girvan2.getCommunities();
//girvan2.printCommunities(communityMap);
//		for (List<Graph> curr : communityMap.keySet()) {
//			System.out.println(curr + " Modularity " + communityMap.get(curr));
//		}
//		HashMap<List<Graph>, Double> currStructure = new HashMap<List<Graph>, Double>();
//		double value = 1.000001;
//		for (List<Graph> community : communityMap.keySet()) {
//			currStructure.put(community, value);
//			System.out.println("Value: " + value);
//			value += 0.000001;
//		}
//		System.out.println("******************** optimal Structure********************************");
//		System.out.println("communityMap size: " + communityMap.size());
//		HashMap<List<Graph>, Double> optimalStructure = girvan2.getOptimalCommunities(currStructure);
//		for (List<Graph> curr : optimalStructure.keySet()) {
//			System.out.println( curr);
//			System.out.println(optimalStructure.get(curr));
//		}
		
//		// test modularity for individual communities
//		girvan2.assignArrayIds();
//		List<Graph> individualCommunities = girvan2.getIndividualCommunities();
//		double modularity = girvan2.getModularity(individualCommunities, girvan2.getModularityMartix());
//		System.out.println(modularity);
//		
		
	}
}



