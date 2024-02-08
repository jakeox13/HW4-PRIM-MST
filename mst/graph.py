import numpy as np
import heapq
import random
from typing import Union

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """
        
        
        # Intilize mst to a numpy array of zeros
        self.mst = np.zeros(np.shape(self.adj_mat))
        
        # Intiize to a random node
        start=random.randint(0, np.shape(self.adj_mat)[0]-1)
        # Create s and T
        s=[start]
        t=dict()
        dist=dict()
        prev=dict()
        
        # Set all starting distances to inf
        for x in range(0,np.shape(self.adj_mat)[0]):
            dist[x]=float("inf")
            prev[x]=None
        dist[start]=0
        # Find intial distances to node
        for x in range(0,np.shape(self.adj_mat)[0]):
            if self.adj_mat[start,x] !=0:
                dist[x]=self.adj_mat[start,x]
                prev[x]=start
        # Set up prority queue
        pq=[]
        for x in range(0,np.shape(self.adj_mat)[0]):
            if x !=start:
                heapq.heappush(pq, (dist[x],x))
        
        # While there are things in the queue
        while pq:
            # Pull the node of shortest distance
            dist_val,u=heapq.heappop(pq)
            # Check if value is in s
            if u not in s:
                # Add node to searched
                s.append(u)
                
                # Add node to end graph
                self.mst[prev[u],u]=dist_val
                self.mst[u,prev[u]]=dist_val
                
                # Check if the distancace from new node to all nodes not in s are shorter than the current distance from s
                for node in range(0,np.shape(self.adj_mat)[0]):
                    if (node not in s) & (self.adj_mat[u,node] !=0) & (self.adj_mat[u,node] < dist[node]) :
                        # Update the shortest distance value in the matrix
                        dist[node]=self.adj_mat[u,node]
                        # Note that edge orginates fromo the new node
                        prev[node]=u
                        heapq.heappush(pq, (dist[node],node))
        
        # If inf remains then there is a disconnect
        if float("inf") in self.mst:
            print("Not all nodes visited graph may be disconnected")
            #Replace with zero(unconnected)
            self.mst[self.mst == float("inf")] = 0