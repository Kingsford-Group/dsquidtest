# Python program to print all paths from a source to destination. 

from collections import defaultdict 
import time
import sys
import random
import tqdm
import networkx as nx

#This class represents a directed graph 
# using adjacency list representation 
class Graph: 

    def __init__(self):         
        # default dictionary to store graph 
        self.graph = defaultdict(set) 
        self.vertices = set()

    # function to add an edge to graph 
    def add_edge(self,u,v): 
        self.graph[u].add(v) 
        self.graph[v].add(u)
        self.vertices.add(u)
        self.vertices.add(v)

    def removeEdge(self, u,v):
        self.graph[u].remove(v)
        self.graph[v].remove(u)

    '''A recursive function to print all paths from 'u' to 'd'. 
    visited[] keeps track of vertices in current path. 
    path[] stores actual vertices and path_index is current 
    index in path[]'''
    def printAllPathsUtil(self, u, d, visited, path, s, start_time):
        if time.time() - start_time > 0.5:
            return -1, None
        # print(u)
        # Mark the current node as visited and store in path 
        visited[u]= True
        path.append(u) 

        result = 0

        # If current vertex is same as destination, then print 
        # current path[] 
        if u == d :
            path.append(s)
            adj = 0
            for i in range(0, len(path)-2):
                if path[i][:-1] != path[i+1][:-1] and path[i+1][:-1] != path[i+2][:-1]:
                    adj+=1
            if path[-2][:-1] != path[0][:-1] and path[0][:-1] != path[1][:-1]:
                adj+=1
            #if d[:-1] == '893':
            #    print(path)
            #    print(adj)
            if adj != 2:
                return 1, path
        else: 
            # If current vertex is not destination 
            #Recur for all the vertices adjacent to this vertex 
            neighbors = list(self.graph[u])
            random.shuffle(neighbors)
            for i in neighbors: 
                if visited[i]==False: 
                    result, path1 = self.printAllPathsUtil(i, d, visited, path, s, start_time)
                    if result == 1:
                        return 1, path1  
                    if result == -1:
                        return -1, path1
        # Remove current vertex from path[] and mark it as unvisited 
        # print("==")
        path.pop()
        visited[u]= False
        return 0, None


    # Prints all paths from 's' to 'd' 
    def printAllPaths(self,s, d): 
        self.removeEdge(s,d)
        # print(self.graph)
        # Mark all the vertices as not visited 
        visited = {}
        for v in self.vertices:
            visited[v] = False

        # Create an array to store paths 
        path = [] 

        # Call the recursive helper function to print all paths 
        ret = self.printAllPathsUtil(s, d,visited, path, s, time.time()) 
        counter = 0
        while ret[0] == -1:
            if counter > 100:
                print(s +"," + d+" -- Time exceeded more than 20 times, skipping...")
                return 0, []
            visited = {}
            for v in self.vertices:
                visited[v] = False
            path=[]
            ret = self.printAllPathsUtil(s, d,visited, path, s, time.time())
            counter+=1
        return ret

def isDiscordant(n1, n2):
    if (n1[-1] == 'H' and n2[-1] == 'T') and (int(n1[:-1]) < int(n2[:-1])):
        return True
    elif n1[-1] == n2[-1]:
        return True
    else:
        return False

def read_graph(fname):
    g = Graph()
    disc_edges = set()
    for l in open(fname):
        if "edge" in l:
            line = l.rstrip("\n").split("\t")
            u1 = line[2]+line[3]
            u2 = line[4]+line[5]
            g.add_edge(u1,u2)
            g.add_edge(u1[:-1]+"H", u1[:-1]+"T")
            g.add_edge(u2[:-1]+"H", u2[:-1]+"T")
            if isDiscordant(u1, u2):
                if int(u1[:-1]) <= int(u2[:-1]):
                    disc_edges.add(u1+","+u2)
                else:
                    disc_edges.add(u2+","+u1)
    return g, disc_edges

if __name__ == "__main__":
    g, disc_edges = read_graph(sys.argv[1])
    disc_conf = defaultdict(int)
    #with tqdm.tqdm(total=len(disc_edges)) as pbar:
    for e in disc_edges:
        if disc_conf[e] == 0:
            r, path = g.printAllPaths(e.split(",")[0], e.split(",")[1])
            if r == 1:
                disc_conf[e] = 1
                for i in range(len(path)-2):
                    if isDiscordant(path[i], path[i+1]):
                        disc_conf[e] = 1
            g.addEdge(e.split(",")[0], e.split(",")[1])
           # pbar.update(1)

    dis_conf_true = [i for i in disc_conf.keys() if disc_conf[i] == 1]
    num_disc_conf = len(dis_conf_true)
    out = open(sys.argv[2], 'w')
    out.write("TOTAL\t"+str(len(disc_edges))+"\n")
    for e in dis_conf_true:
        u1 = e.split(",")[0]
        u2 = e.split(",")[1]
        out.write(u1[:-1]+"\t"+str(u1[-1]=='H')+"\t"+u2[:-1]+"\t"+str(u2[-1]=='H')+"\n")
