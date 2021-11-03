function [TSP_Cost ]=Cluster_NN(k,Iter_Param,Sample_Location,Depot)
%=============================NEAREST_NEIGHBOUR_PROCEDURES================
TSP_NN = cell(1,k);                 % Create Array to Store TSP
TSP_Distance = cell(1,k);           % Create Array to Store Total Distance
Cluster_Nodes = Iter_Param;         % Create Array Cluster Nodes from each cluster
Depot_Node = 1;                     % Initiate the Depot Node

% Update Nodes Number in Cluster
for it_x = 1:k
    for it_y = 1:size(Cluster_Nodes{it_x},2)
        Cluster_Nodes{it_x}(it_y) = Cluster_Nodes{it_x}(it_y)+1;
    end
end
Unvisited_Nodes = Cluster_Nodes;        % Initiate Unvisited Nodes
Location = [Depot;Sample_Location];     % Update Matrix by adding Location of Depot
DM = pdist2(Location,Location);         % Create Distance Matrix
% Solve the TSP using NN Procedure
for it = 1 : k
    go_sign = numel (Unvisited_Nodes{it});
    total_route = 0;
    Visited_Nodes = Depot_Node;
    while go_sign ~=0
        [min_distance,min_index] = min(DM(Visited_Nodes,Unvisited_Nodes{it}));
        Visited_Nodes = Unvisited_Nodes{it}(min_index);
        total_route = total_route+min_distance;
        TSP_NN{it}=[TSP_NN{it} Visited_Nodes];
        Unvisited_Nodes{it}(min_index)=[];
        go_sign=numel (Unvisited_Nodes{it});
    end
    TSP_NN{it}=[Depot_Node TSP_NN{it} Depot_Node];
    TSP_Distance{it}=total_route+DM(Depot_Node,TSP_NN{it}(end-1));
end
TSP_Cost = 0;
for it = 1 : k
    Cost = TSP_Distance{it};
    TSP_Cost = TSP_Cost + Cost;
end
TSP_Cost;
end



