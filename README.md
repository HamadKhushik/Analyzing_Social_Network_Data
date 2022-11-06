

  
    
      
        
          
          


-------------------------------------------------------------------------------------------------------------------

Overview: 
==========

This project analyzes data to find sub-communities/clusters within the data based on the connections in data. A detailed description of community structure can be found [on this wikipedia page](https://en.wikipedia.org/wiki/Community_structure)

The picture below shows that there are two sub-communities within the data, pointed out by blue and purple colours. In picture 2 a more complex structure of sub-communities is shown.

![Alt text](images/community%20finding%201.PNG?raw=true "Picture 1")
Picture 1

![](images/community%20finding%202.PNG?raw=true "Picture 2")
Picture 2


This program detects sub-communities/clusters within a community using Girvan-Newman algorithm.    
Girvan-Newman Algorithm is implemented in two ways

1. without using external libraries(from scratch)
2. using external library i-e JUNG library (GirvanNewmanJung class implements this)

----------------------------------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------------------------------------

Data Structure: 
---------------
The main Data Structure used is graph using the adjacency list.   
However, adjacency matrix is also implemented (which is used only for modularity calculation). Modularity is calculated using the Modularity Matrix - which is evaluated using the Adjacency Matrix.


Algorithm: 
-----------
The main algorithm used is Girvan-Newman Edge Betweenness Algorithm for community finding. In order to find the edge betweenness, Brande’s edge betweenness algorithm is used.   And Girvan-Newman method is used to find the modularity of the graph.


Algorithm analysis: 
--------------------
Total running time of the algorithm is O(n2m2), where n is the number of vertices and m is the total number of edges. 
This time complexity is calculated from the getCommunities() method.  

In getCommunities() method, I have a while loop which loops until no edges are left i-e O(m). within this loop, methods called   are as follows with their complexities  
brandes() -> 			O(n2 * m)  
removeEdges() -> 		O(m)  
getSCCs()->			O(n * m)  
getModularity() -> 		O(n * m)  
  
combining all the above complexities gives us the following  
O(m) { O(n2 *m) + O(m) + O(n*m) + O(n*m)} -> ignoring the smaller terms  
= O(m) {O(n2*m)}  
= O(m) * O(n2*m)  
= O(n2*m2)  
  
--------------------------------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------------------------------------
Overall Design: 
=================
The design is simple and focused towards readability.  
GirvanNewman class is used to implement the Girvan Newman Algorithm without any libraries  
GirvanNewmanJung class implements the algorithm using the JUNG library.  
GirvanNewman class has the major chunk of code. It extends the CapGraph class to make it compatible with provided Graph Interface and also to re-use the code written in CapGraph class.  
GirvanNewman class has two inner classes to represent the vertex and edge in the graph. Method brandes() calculates the edge-betweenness in the graph using the Brande’s algorithm, based on which, the edges with highest betweenness are then removed. Method getModularity() calculates the modularity of the community structure found by removing the edges.   

Method getCommunities() performs all the work in a single call, finds edge-betweenness, removes edges, calculates modularity for every community structure found until no more edges are left, records communities with top three modularity scores and returns them.

--------------------------------------------------------------------------------------------------------------------------
Class Design
=============

GirvanNewman.java
--------------------------
This class extends CapGraph class to make it compatible with the provided code. This is the main class to find communities in a graph using Girvan-Newman Edge Betweenness algorithm. This class performs all the functions to find the communities in the graph. It has two inner classes for Vertex and Edge in the graph.

GirvanNode.java  
--------------------------
This is a inner class of Girvan Newman class and serves as the vertex for the graph to be analysed by Girvan Newman class.

GirvanEdge.java  
------------------------
This is a inner class of Girvan Newman class and serves as the edge for the graph to be analysed by Girvan Newman class.

GirvanNewmanJung.java  
--------------------------------
This class do not have lot of code and uses JUNG library to find communities in the graph. This class depends on JungGraphLoader class to load the graph for GirvanNewmanJung class to work on.

JungGraphLoader.java  
--------------------------------
Purpose and description of class: The only purpose of this class is to load the graph for GiranNewmanJung class. It has only one static method which loads the graph.

-------------------------------------------------------------------------------------------------------------------------------
 
 
