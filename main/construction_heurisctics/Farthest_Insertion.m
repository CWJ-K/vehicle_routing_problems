function [TSP_Cost ]=Farthest_Insertion(k,Iter_Param,Sample_Location,Depot)
%=============================FARTHEST_INSERTION_PROCEDURES================
TSP_FI = cell(1,k);                 % Create Array to Store TSP
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

% Solve the TSP using FI Procedure
for it = 1 : k
    total_route = 0;
    TSP = [Depot_Node];
    [min_distance,min_index] = min(DM(TSP,Unvisited_Nodes{it}));
    TSP = [TSP Unvisited_Nodes{it}(min_index) TSP];
    Tour_size = 3;
    Unvisited_Nodes{it}(min_index)=[];
    total_route = total_route + min_distance + min_distance;
    go_sign = numel(Unvisited_Nodes{it});
    while go_sign ~=0
    % select the farthest node to insert based on Insertion Criteria =
    % minimum insertion cost
    if(go_sign>1)
        Dis_to_TSP=DM(Unvisited_Nodes{it},TSP(1,1:Tour_size-1));
        [Dis_to_TSP,i_uninserted]=max(Dis_to_TSP);
        [M,I]=max(Dis_to_TSP);
        nodeindex=i_uninserted(I);
        node=Unvisited_Nodes{it}(nodeindex);
    else
        node=Unvisited_Nodes{it}(1);
        nodeindex=1;
    end
    
    %calculate insertion cost of the chosen node in the TSP tour
    MinInsertionCost=inf;
    for j=1:size(TSP,2)-1
        insertionCost=DM(TSP(j),node)+DM(TSP(j+1),node)-DM(TSP(j),TSP(j+1));
        if(MinInsertionCost>insertionCost)
            inspos=j;
            MinInsertionCost=insertionCost;
        end  
    end
    
    % insert the node in the best position in the tour
    TSP=[TSP(1,1:inspos) node TSP(1,inspos+1:end)];
    total_route = total_route + MinInsertionCost;        
    Unvisited_Nodes{it}(nodeindex)=[];
    go_sign=numel (Unvisited_Nodes{it});
    Tour_size=Tour_size+1;
    end
    %update the TSP and Total Distance
    TSP_FI{it}= TSP;
    TSP_Distance{it}=total_route;
end
TSP_Cost = 0;
for it = 1 : k
    Cost = TSP_Distance{it};
    TSP_Cost = TSP_Cost + Cost;
end
TSP_Cost;
end



