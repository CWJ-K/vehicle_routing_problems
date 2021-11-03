function [NM_Cost ]= Nearst_Merger(k,Iter_Param,Sample_Location,Depot)
%=============================NEARST_MERGER_PROCEDURES================
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

% Solve the TSP using NM Procedure
isdepot=1;  % for the first node is depot

for it=1:k
    %create a new matrix for each nodes in each k groups
    nodes_k=Iter_Param{it};   %use array without depot to locate each node in Sample_Location(without depot)
    nodes_location=Sample_Location(nodes_k,:); % pick out the (x,y) of each node in group k
    Location = [Depot;nodes_location];  % Update Matrix by adding Location of Depot
    DM = pdist2(Location,Location);         % Create Distance Matrix
    
    Nb_Nodes = numel(Unvisited_Nodes{it})+1;  % add depot
    
    %Initialize the sub-tours to be merged
    Subtours=cell(Nb_Nodes,1);
    
    for i = 1:Nb_Nodes
        Subtours(i)={i};
    end
    
    Nb_Subtours=Nb_Nodes;
    
    Subtours_Dist=DM;
    for S1 = 1:Nb_Subtours
        Subtours_Dist(S1,S1)=inf;
    end
    
    % Repeate until all subtours are merged
    while Nb_Subtours~=1
        [MIN_Distances,S] = min(Subtours_Dist);
        [MinMergingCost,S1] = min(MIN_Distances);
        S2=S(S1);
        %%% Merge S1 and S2
        
        s1_Size=size(Subtours{S1},2);
        s2_Size=size(Subtours{S2},2);
        %
        if(s1_Size==1 && s2_Size==1 ) % if both subtours have length of one
            NewSubtour=[Subtours{S1}(1) Subtours{S2}(1) Subtours{S1}(1)];
        else % if one subtour has length of one and the other has length of more than one
            if((s1_Size==1 && s2_Size~=1 )||(s1_Size~=1 && s2_Size==1 ))
                if(s2_Size~=1)
                    S_Temp=S2;
                    S2=S1; % subtour 2 is with lenght 1
                    S1=S_Temp;
                    
                    S_Temp_Size= s1_Size;
                    s1_Size=s2_Size;
                    s2_Size=S_Temp_Size;
                end
                
                Cheapest_Merging_Cost=Inf;
                
                %%%
                for i=1:s1_Size-1
                    % find the cheapeast insertion cost of S1 into S2 and Merge them
                    Merging_Cost=DM(Subtours{S1}(i),Subtours{S2}(1))+DM(Subtours{S2}(1),Subtours{S1}(i+1))-DM(Subtours{S1}(i),Subtours{S1}(i+1));
                    
                    if(Merging_Cost< Cheapest_Merging_Cost)
                        Cheapest_Merging_Cost=Merging_Cost;
                        Position_i=i;
                    end
                    
                    
                end
                
                NewSubtour = [Subtours{S1}(1:Position_i) Subtours{S2}(1) Subtours{S1}(Position_i+1:end)];
                
            else % if both subtours have length of more than one
                %%% Merge S1 and S2
                
                Cheapest_Merging_Cost=Inf;
                % find cheapest merging cost of S1 and S2 by calculating all combinations of breaking an edge in each subtour and merging them
                for i=1:s1_Size-1
                    for j=1:s2_Size-1
                        % calculate the breaking cost of edges
                        % i & i+1 in S1 AND j  & j+1 in S2 should be broken
                        BreakingCost=DM(Subtours{S1}(i),Subtours{S1}(i+1))+DM(Subtours{S2}(j),Subtours{S2}(j+1));
                        % No Reversing
                        Merging_Cost1=DM(Subtours{S1}(i),Subtours{S2}(j+1))+DM(Subtours{S2}(j),Subtours{S1}(i+1))-BreakingCost;
                        
                        % Reversing
                        % i & j and i+1 & j+1 are the nodes to join
                        Merging_Cost2=DM(Subtours{S1}(i),Subtours{S2}(j))+DM(Subtours{S2}(j+1),Subtours{S1}(i+1))-BreakingCost;
                        
                        if(Merging_Cost1<= Merging_Cost2 && Merging_Cost1<Cheapest_Merging_Cost)
                            Reverse=false ;
                            Cheapest_Merging_Cost=Merging_Cost1;
                            Position_i=i;
                            Position_j=j;
                        else
                            if(Merging_Cost1> Merging_Cost2 && Merging_Cost2<Cheapest_Merging_Cost)
                                Reverse=true ;
                                Cheapest_Merging_Cost=Merging_Cost2;
                                Position_i=i;
                                Position_j=j;
                            end
                        end
                    end
                end
                
                
                
                %%% merge two subtours given best merging option
                if(Reverse)
                    NewSubtour = [Subtours{S1}(1:Position_i) fliplr( Subtours{S2}(1:Position_j)) fliplr(Subtours{S2}(Position_j+1:end-1)) Subtours{S1}(Position_i+1:end)];
                    
                else
                    NewSubtour = [Subtours{S1}(1:Position_i) Subtours{S2}(Position_j+1:end-1) Subtours{S2}(1:Position_j) Subtours{S1}(Position_i+1:end)];
                end
                
            end
            
            
        end
        
        %Update list of subtours
        if(S1>S2)
            % Update list of subtours by removing S1 and S2 from the list
            Subtours(S1,:) = [];
            Subtours(S2,:) = [];
            
            % Update the subtours distance matrix
            % by deleting the row ans colomn corresponding to S1 & S2
            Subtours_Dist(S1, :) = [];
            Subtours_Dist(S2, :) = [];
            
            Subtours_Dist(:,S1) = [];
            Subtours_Dist(:,S2) = [];
            
        else
            % Update list of subtours by removing S1 and S2 from the list
            Subtours(S2,:) = [];
            Subtours(S1,:) = [];
            
            % Update the subtours distance matrix
            % by deleting the row ans colomn corresponding to S1 & S2
            Subtours_Dist(S2, :) = [];
            Subtours_Dist(S1, :) = [];
            
            Subtours_Dist(:,S2) = [];
            Subtours_Dist(:,S1) = [];
        end
        
        
        % add the new subtour to list of subtours
        Subtours(end+1,:) = {NewSubtour};
        
        % update Nb_Subtours
        Nb_Subtours=Nb_Subtours-1;
        
        % Update the subtours distance matrix by adding a row and column
        % corresponding to the new subtour at the end
        
        % add one row
        r=Inf(1,Nb_Subtours-1);
        Subtours_Dist=[Subtours_Dist;r];
        % add one column
        c=Inf(Nb_Subtours,1);
        Subtours_Dist=[Subtours_Dist,c];
        
        % calculate distance between new subtour and the rest
        S1=Nb_Subtours;
        
        for S2 = 1:Nb_Subtours-1
            
            Min=min(DM(Subtours{S1},Subtours{S2}));
            
            if(size(Min,2)>1 ||size(Min,1)>1 )
                Min=min(Min);
            end
            
            Subtours_Dist(S1,S2)=Min;
            Subtours_Dist(S2,S1)=Subtours_Dist(S1,S2);
        end
    end

TSP_Tour= Subtours{1};
% Calculate the total cost of the TSP tour
x=Location( Subtours{1},1);
y=Location( Subtours{1},2);

KCost = sum(sqrt(diff(x).^2 + diff(y).^2));
TSP_NM {it}=KCost;

end

NM_Cost = 0;
for it = 1 : k
    Cost = TSP_NM {it};
    NM_Cost = NM_Cost + Cost;
end

NM_Cost;

end
