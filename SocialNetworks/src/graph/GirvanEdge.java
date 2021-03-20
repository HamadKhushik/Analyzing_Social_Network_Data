package graph;

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

}
