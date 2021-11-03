function [TSP_Cost] = CWS(k,Iter_Param,Sample_Location,Depot)

%  Initialize the sub-tours to be merged
TSP_NM = cell(1,k);                 % Create Array to Store TSP
TSP_Distance = cell(1,k);           % Create Array to Store Total Distance
Cluster_Nodes = Iter_Param;         % Create Array Cluster Nodes from each cluster
Depot_Node = 1;                     % Initiate the Depot Node

%Update Nodes Number in Cluster
for it_x = 1:k
    for it_y = 1:size(Cluster_Nodes{it_x},2)
        Cluster_Nodes{it_x}(it_y) = Cluster_Nodes{it_x}(it_y)+1;
    end
end


Unvisited_Nodes = Cluster_Nodes;        % Initiate Unvisited Nodes

% A list can be empty:


% Create an back and forth subtour to the depot
for it=1:k
    Subtours={};
    %create a new matrix for each nodes in each k groups
    nodes_k=Iter_Param{it};   %use array without depot to locate each node in Sample_Location(without depot)
    nodes_location=Sample_Location(nodes_k,:); % pick out the (x,y) of each node in group k
    Location = [Depot;nodes_location];  % Update Matrix by adding Location of Depot
    Dis = pdist2(Location,Location);         % Create Distance Matrix
    Nb_Nodes = numel(Unvisited_Nodes{it})+1;
    for i = 1:Nb_Nodes
        if(i~=Depot_Node)
            Subtour=[Depot_Node i Depot_Node];
            % add new entries to the end of a list, use the append command or the . (dot) operator:
            Subtours(end+1,:) =  {Subtour};
            
        end
    end
    
    % Access Data in Cell Array
    % https://uk.mathworks.com/help/matlab/matlab_prog/access-data-in-a-cell-array.html
    Nb_Subtours=size(Subtours,1);
    
    
    % calculate the savings obtained from merging two subtours
    Savings_matrix=zeros(Nb_Subtours,Nb_Subtours);
    for S1 = 1:Nb_Subtours
        for S2 = 1:Nb_Subtours
            if(S1~=S2)
                % i & i+1 in S1 AND j  & j+1 in S2 should be broken
                Savings_matrix(S1,S2)=Dis(Subtours{S1}(2),Subtours{S1}(3))+Dis(Subtours{S2}(1),Subtours{S2}(2))-Dis(Subtours{S1}(2),Subtours{S2}(2));
            end
        end
        
    end
    
    
    while Nb_Subtours(1)~=1
        % Select the two sub-tours with maximum savings, where two sub-tours, say S and S?, qualify as the closest when the  distance between any pair of nodes,
        % one belonging to S and the other belonging to S?, is minimal amongst all pairs of sub-tours (S, S?).
        
        % find two Subtours with Maximum saving (S1,S2)
        maxSaving=0;
        for s1 = 1:Nb_Subtours
            
            [MAX,s2] = max(Savings_matrix(s1,:));
            if(maxSaving<MAX)
                
                S1=s1;
                S2=s2;
                maxSaving=MAX;
            end
        end
        
        % Merge S and S'
        NewSubtour=[Subtours{S1}(1:end-1) Subtours{S2}(2:end)];
        % Update list of subtours and savings
        if(S1>S2)
            % Update list of subtours
            Subtours(S1,:) = [];
            Subtours(S2,:) = [];
            
            % Update list of  savings
            Savings_matrix(S1, :) = [];
            Savings_matrix(S2, :) = [];
            
            Savings_matrix(:,S1) = [];
            Savings_matrix(:,S2) = [];
            
        else
            % Update list of subtours
            Subtours(S2,:) = [];
            Subtours(S1,:) = [];
            
            % Update list of  savings
            Savings_matrix(S2, :) = [];
            Savings_matrix(S1, :) = [];
            
            Savings_matrix(:,S2) = [];
            Savings_matrix(:,S1) = [];
        end
        
        
        Subtours(end+1,:) = {NewSubtour};
        Nb_Subtours=Nb_Subtours-1;
        r=zeros(1,Nb_Subtours-1);
        Savings_matrix=[Savings_matrix;r];% add one row
        
        c=zeros(Nb_Subtours,1);
        Savings_matrix=[Savings_matrix,c];% add one column
        
        % calculate savings for new subtour and the rest
        
        S1=Nb_Subtours;
        for S2 = 1:Nb_Subtours
            if(S1~=S2)
                % i & i+1 in S1 AND j  & j+1 in S2 should be broken
                Savings_matrix(S1,S2)=Dis(Subtours{S1}(end-1),Subtours{S1}(end))+Dis(Subtours{S2}(1),Subtours{S2}(2))-Dis(Subtours{S1}(end-1),Subtours{S2}(2));
                Savings_matrix(S2,S1)=Dis(Subtours{S2}(end-1),Subtours{S2}(end))+Dis(Subtours{S1}(1),Subtours{S1}(2))-Dis(Subtours{S2}(end-1),Subtours{S1}(2));
            end
            
        end
    end
    
    
    
    % Output : calculate the tour cost
    x=Location( Subtours{1},1);
    y=Location( Subtours{1},2);
    
    L = sum(sqrt(diff(x).^2 + diff(y).^2));
    Total_Distance{it} = L;
    Tour{it}=Subtours{1};
end
TSP_Cost = 0;
for it = 1 : k
    Cost = Total_Distance{it};
    TSP_Cost = TSP_Cost + Cost;
end
end
