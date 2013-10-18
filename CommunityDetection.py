from collections import defaultdict
from decimal import Decimal
import sys, time
import copy
from operator import itemgetter

#reading input file  
def init(graph):
 data=[(line.strip().split(' ')[0],line.strip().split(' ')[1]) for line in file(graph)]
 data.remove(data[0])
 return data

print "enter graph file path"
input_graph=raw_input()
start_time = time.time()
data=init(input_graph)


#adjacency list creation
def getAdjacencylist(data):
 alist = defaultdict(list) 
 for e in data:
   alist[int(e[0])].append(int(e[1]))  
   alist[int(e[1])].append(int(e[0]))
 return alist
 
alist = getAdjacencylist(data)
size = len(alist)
p = [Decimal(1)/size] * size   #initialize page rank vector 

#calculate the page rank vector sum in the reducer
def sumRank(index,rankList,p3):
 rank = reduce(lambda x,y: x+y,rankList)
 p3[index]=rank
 
#emit neighbor and its rank 
def mapper(adjList,page_rank_vector):
 reducerList = defaultdict(list) 
 for k,v in adjList.items():
  res=map(lambda neighbor,rank:reducerList[neighbor].append(rank),adjList[k],[page_rank_vector[k]/len(adjList[k])]* len(adjList[k]))
 return reducerList     
    
def reducer(reducerList):
 p2 = [Decimal(0)] * size
 for k,v in reducerList.items():
  sumRank(k,v,p2)
 return p2

#page rank algorithm implementation
def sequentialPageRank(p):
 for i in range(1,30):
  iteration_start_time=time.time()
  reducerInputList = mapper(alist,p)
  p_new = reducer(reducerInputList)
  iteration_end_time=time.time()
  if abs(p[0]-p_new[0]) < sys.float_info.epsilon : break #if page rank is not changing significantly then break 
  p=p_new
 return p 
 

p=sequentialPageRank(p)
#pickle the page rank vector output
#pageRankVector=open('page_rank_vector.pkl','wb')
#pickle.dump(p,pageRankVector)
#pageRankVector.close()

#f=open('page_rank_vector.pkl','rb')
#p=pickle.load(f)

#calculate the weight metric - average degree 
def density(cluster,adjList, new_vertex=-1):
    if new_vertex ==-1:
     m_s= Decimal(cluster['num_edges'])
     #number of edges inside the community
     n_s= Decimal(len(cluster['vertices']))#vertices inside the community 
     if n_s == 0: return 0
     cluster_density= 2 *m_s /n_s
     
     return cluster_density
    else:
     #calculating density after adding new vertex
     new_links=[n for n in adjList[new_vertex] if n in cluster['vertices']]
     m_s= Decimal(cluster['num_edges'] + len(new_links)) #number of edges inside the community
     n_s= Decimal(len(cluster['vertices']) +1) #vertices inside the community 
     #According to the paper, we should use the weight metric cluster_weight= 2 *m_s /n_s
     #return cluster_density
     cluster_density = 2 * m_s/ n_s
     return cluster_density

#Implementation of Link Aggregation algorithm to generate seed communities 
def LinkAggregation(adjList, page_rank_vector):
 community=[] #list of communities
 dic={}
 counter=0 
 for k,v in enumerate(page_rank_vector): dic[k]=v #vertex id is key, page rank is value 
 for vertex,rank in sorted(dic.iteritems(), key=itemgetter(1), reverse=True): #iteration over all vertices, non increasing page rank order
  added=False
  counter+=1
  #add the vertex to each community and compare the new weight of the community with the previous weight.
  for C_i in community:
   newDensity=density(C_i,adjList,vertex); 
   if newDensity >0 and newDensity >= density(C_i,adjList):
    #append the vertex to the current community if the weight increases on adding the vertex.
    new_links=[n for n in adjList[vertex] if n in C_i['vertices']]
    C_i['num_edges']=C_i['num_edges']+len(new_links)
    C_i['vertices'].append(vertex)
    added=True
  if added == False:
   #making new cluster if the weight of any community doesn't increase on adding the given vertex to it.
   new_cluster={}
   new_cluster['vertices']=[]
   new_cluster['vertices'].append(vertex)
   new_cluster['num_edges']=0
   community.append(new_cluster)
 return community
 
#Implementation of IS2 to iteratively improve the quality of detected communities
def IS2 (seedCommunity,adjList):
    C=copy.deepcopy(seedCommunity)
    w = density(C,adjList);
    increased = True;
    while increased:
        N = copy.deepcopy(C);
        for v in C['vertices']:
            neighbours = adjList[v];
            #generate the neighbourhood of the current community.
            for vertex in neighbours:
              if ((vertex not in N['vertices'])):
                  new_links=[n for n in adjList[vertex] if n in N['vertices']]
                  N['num_edges']=N['num_edges']+len(new_links)
                  N['vertices'].append(vertex)
                  N['vertices'] = list(set(N['vertices']))
        for v in N['vertices']:
            C_prime= copy.deepcopy(C)
            if v in C['vertices']:
                #removing v from C_prime if v exists in the community.
                neighbours=[u for u in adjList[v] if u in C_prime['vertices']]
                C_prime['num_edges'] = C_prime['num_edges']-len(neighbours)
                C_prime['vertices'].remove(v)
            else:
                #add vertex to the C_prime community.
                new_links=[u for u in adjList[v] if u in C_prime['vertices']]
                C_prime['num_edges']=C_prime['num_edges']+len(new_links)
                C_prime['vertices'].append(v)
            if density(C_prime,adjList) > density(C,adjList):
                #if the weight increases, select the new community else discard it.
                C = C_prime;
        if density(C,adjList) == w:
            increased = False
        else:
            w = density(C,adjList);
    return C

#Write the output to a file
def writeResults(results):
    fp = open("output.txt","w");
    for result in results:
        for vertex in result:
            fp.write(str(vertex)+" ");
        fp.write("\n");

communities_seed=LinkAggregation(alist,p)
#laOutput=open('laoutput.pkl','wb')
#pickle.dump(communities_seed,laOutput)
#laOutput.close()

finalCommunities = []

for community in communities_seed:
 C=IS2(community,alist)
 finalCommunities.append(C['vertices'])
sortedFinalCommunities=[]
for i in range(len(finalCommunities)):
    sortedFinalCommunities.append(sorted(finalCommunities[i]))
sortedFinalCommunities = sorted(sortedFinalCommunities)
#remove duplicate communitites.
result = [sortedFinalCommunities[i] for i in range(len(sortedFinalCommunities)) if i == 0 or sortedFinalCommunities[i] != sortedFinalCommunities[i-1]]
print "Final communities = ",result
writeResults(result)
