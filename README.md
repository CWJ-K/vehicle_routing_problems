# Vehicle Routing Problems

Vehicle Routing Problems (VRP) are commonly used in the transportation industry. The goal of this problem is to find the optimal tour from the depot to the customers in order to meet customers’ demands. Capacitated VRP (CVRP) is one of the VRP problems in which the customers will be served by the identical or homogenous vehicles and each vehicle is limited in capacity. 

Geetha, et.al (2009) proposed one method which could be approached as the initial step while
solving CVRP with multiple vehicles. The method basically uses the concept of k-means clustering algorithm and improves that algorithm by adding priority term to assign customer in the cluster rather than using the distance variable alone. Therefore, this method ensures the vehicle capacity is optimally used since the clustering algorithm considers not only the distance but also the demand of the customers. The customers with larger demand and lower distances get the higher priority.

In this project, we deal with 63 problems of CVRP derived from several sources of VRP
problems. The CVRP problems used are the symmetrical problems with single depot and known location. The initial step used to solve the problem is clustering the customer using improved kmeans clustering algorithm. This algorithm is successfully optimizing the usage of the vehicles (the clusters formed are equal to the minimum number of the vehicles). 

Seven construction heuristics methods are applied to build the complete tours for each cluster.
Farthest Insertion method is chosen among the other methods based on the lowest percentage improvement over the optimal solution. Furthermore, six Improvement Heuristics methods are used to improve the initial TSP tour from construction heuristics. While implementing Improvement Heuristics methods, the intra-move (2-opt and 3-opt) and inter-move (node exchange and node relocation) are used to perform nodes movements. Three criteria are used to compare each of the Improvement Heuristic methods. Firstly, the percentage improvement over the initial solution. Secondly, the CPU time used to perform the algorithm. Thirdly, the percentage
improvement over the optimal solution.
