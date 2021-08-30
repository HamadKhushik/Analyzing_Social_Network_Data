package graph;

import java.io.File;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import edu.uci.ics.jung.algorithms.cluster.EdgeBetweennessClusterer;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Hypergraph;
import edu.uci.ics.jung.graph.SparseGraph;
import util.GraphLoader;

public class GirvanNewmanJung {
	
	
	public static void main(String[] args) {
		
		Graph<Integer, String> jungGraph = new DirectedSparseGraph<Integer, String>();
		JungGraphLoader.loadGraph(jungGraph, "data/facebook_1000.txt");
		System.out.println("Total Number of edges in the graph are: " + jungGraph.getEdgeCount());
		EdgeBetweennessClusterer<Integer, String> clusterer = new EdgeBetweennessClusterer<Integer, String>(2000);
		System.out.println(clusterer.transform(jungGraph));

	}

}
