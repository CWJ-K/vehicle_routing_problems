function [TSP_Cost]=Cheapest_Insertion(k,Iter_Param,Sample_Location,Depot)
%=============================CHEAPEST_INSERTION_PROCEDURES================
TSP_CI = cell(1,k);                 % Create Array to Store TSP
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
% Solve the TSP using Cheapest Insertion Procedure
for it = 1 : k
    go_sign = numel (Unvisited_Nodes{it});
    total_route = 0;
    Visited_Nodes = Depot_Node;
    % find the closest node to the depot
     [min_distance,min_index] = min(DM(Visited_Nodes,Unvisited_Nodes{it}));
      Visited_Nodes = Unvisited_Nodes{it}(min_index);
    % create close tour from depot to visited_nodes
      TSP_CI{it} = [TSP_CI{it} Depot_Node Visited_Nodes Depot_Node];
      TSP_Distance{it} = min_distance + DM(TSP_CI{it}(end-1),TSP_CI{it}(end));
    % update list Unvisited_Nodes
    if numel(Unvisited_Nodes{it})>1
        Unvisited_Nodes{it}(min_index)=[];
    end
      
    
    while go_sign ~=0
       MinInsertionCost = Inf;
       for y_axis = 1 : size(TSP_CI{it},2)-1 
            for x_axis = 1 : numel(Unvisited_Nodes{it})
                insertioncost = DM(TSP_CI{it}(y_axis),Unvisited_Nodes{it}(x_axis))+DM(Unvisited_Nodes{it}(x_axis),TSP_CI{it}(y_axis+1))- DM(TSP_CI{it}(y_axis),TSP_CI{it}(y_axis+1));
                if (insertioncost<MinInsertionCost)
                  inspos = y_axis+1;
                  insnode = Unvisited_Nodes{it}(x_axis);
                  axis_ref = x_axis;
                  MinInsertionCost = insertioncost;
                end
            end
       end
       
       Unvisited_Nodes{it}(axis_ref) = [];
       go_sign=numel(Unvisited_Nodes{it});
       TSP_CI{it}= [TSP_CI{it}(1:inspos-1) insnode TSP_CI{it}(inspos:end)];
       TSP_Distance{it} = TSP_Distance{it} + MinInsertionCost; 
    end

end
TSP_Cost = 0;
for it = 1 : k
    Cost = TSP_Distance{it};
    TSP_Cost = TSP_Cost + Cost;
end
TSP_Cost;

