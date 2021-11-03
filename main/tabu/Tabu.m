clear
%======================================= Retrieve Data %==============================
Filelist = dir('*.mat');
FileName = {Filelist.name};
StoreTable = {'No', 'Problem','TSP Cost','Time (in second)'};
for length = 1 : length(Filelist)
    newData = load(FileName{length}, 'Demand','Sample_Location','Capacity','Depot');
    Sample_Location = newData.Sample_Location;
    Demand =newData.Demand;
    Capacity = newData.Capacity;
    Depot = newData.Depot;
    % ============================= Clustering ===============================%
    k = ceil(sum(Demand)/Capacity) ;                                % NB of Cluster
    [Sort_Demand,Index_Sort_Demand] = sort(Demand,'descend');       % sort_Demand
    ListofNodes = [1:size(Sample_Location,1)];                      % List of Nodes
    NbofNodes=size(ListofNodes,2);                                  % Nb of Nodes
    max_capacity = Capacity;                                        % maximum capacity
    available_nodes = ListofNodes;                                  % available_nodes to be clustered
    available_capacity = cell(1,k);                                 % capacity parameter
    Iter_Param = cell(1,k);                                         % final cluster member
    T_Demand = cell(1,k);                                           % final total_demand for each cluster
    NP = cell(1,k);                                                 % final cluster size
    TD = cell(1,k);                                                 % Total Demand for each cluster
    CL = cell(1,k);                                                 % list of cluster
    
    %========================= ARRAY INIT VALUE ============================
    for it = 1 : k
        Iter_Param{it} = Inf;
        NP{it} = size(Iter_Param{it},2);
        TD{it} = 0;
        available_capacity{it} = max_capacity;
    end
    tic;
    %========================= K_MEANS_CLUSTERING =============================
    iteration = 0       ; % set iteration criteria
    k_means_stop = 1    ; % 1 = go loop ; 0 = stop loop
    while (k_means_stop ~= 0)
        % centroids Selection
        if iteration == 0
            centroids = [];
            for it = 1 : k
                selected_nodes = Index_Sort_Demand(it);
                centroids = [centroids selected_nodes];
            end
        else
            centroids = cell(1,k);
            for c = 1 : k
                cl = Iter_Param{c};
                x_init_cl = [];
                y_init_cl = [];
                for nbcl = 1 : size(cl,2)
                    x_cl = Sample_Location(cl(nbcl),1);
                    y_cl = Sample_Location(cl(nbcl),2);
                    x_init_cl = [x_init_cl x_cl];
                    y_init_cl = [y_init_cl y_cl];
                end
                x_mean = mean(x_init_cl);
                y_mean = mean(y_init_cl);
                centroids{c} = [x_mean y_mean];
            end
        end
        
        % create distance matrix from centroids to all nodes
        k_lists = [];
        for it = 1 : k
            % calculate distance for each centroid to all nodes
            if iteration == 0
                distance = pdist2(Sample_Location((centroids(it)),:),Sample_Location);
            else
                distance = pdist2(centroids{it},Sample_Location);
            end
            % put distance value to matrix k_lists
            if it == 1
                k_lists = [distance];
            else
                k_lists = [k_lists;distance];
            end
        end
        
        % create priority matrix from centroids to all nodes
        priority_lists = [];
        col_list_array = cell(1,1);       % helper array to put list value
        for row_i = 1:k
            % calculate priority value for each centroid to all nodes
            for col_j = 1:size(k_lists,2)
                col_list = k_lists(row_i,col_j)/Demand(col_j);
                col_list_array{1}(:,col_j) = col_list;
            end
            % put priority value to matrix priority_lists
            if row_i==1
                priority_lists=[col_list_array{1}(1:end)];
            else
                priority_lists=[priority_lists;col_list_array{1}(1:end)];
            end
        end
        
        % Adjust available nodes & available capacity
        if iteration == 0
            available_nodes = setdiff(ListofNodes,centroids);
            for it = 1 : k
                available_capacity{it} = available_capacity{it}-Demand(centroids(it));
            end
        else
            available_nodes = ListofNodes;
            for it = 1 : k
                available_capacity{it} = max_capacity;
            end
        end
        
        % start nodes assignments procedures
        j = 1;
        while size(available_nodes,2)~=0
            [k_value_sort,k_row_index] = sort(k_lists(:,available_nodes(j)),'ascend');
            m = k_row_index(1);         % set m as cluster referrence
            % group unassigned nodes which have the lowest priority value refer to m
            nodes_lists = [];           % empty cell to store node
            candidate_lists = [];       % empty cell to store priority value
            for x = 1 : size(available_nodes,2)
                % Find if the node has the lowest priority value to m
                [candidate_value,candidate_row] = min(priority_lists(:,available_nodes(x)));
                if candidate_row == m
                    nodes_lists = [nodes_lists available_nodes(x)]; % put the node into node list
                    candidate_lists = [candidate_lists candidate_value]; % put the priority value into priority list
                end
            end
            % sort the priority for node lists in ascending order (node which we want to put into m)
            [prior_value,prior_ref] = sort(candidate_lists,'ascend') ;
            node_choose = nodes_lists(prior_ref);
            % sort the cluster priority for the each node in the shortlist in ascending
            [value_prio_cn,ind_prio_cn] = sort(priority_lists(:,node_choose),'ascend');
            % assign chosen node to the cluster
            for loop_x = 1 : size(ind_prio_cn,2)
                i = 1;
                stop = 0;
                while (stop ~= 1)
                    % if demand meets the capacity then include the nodes to cluster i
                    if (Demand(node_choose(:,loop_x)) <= available_capacity{ind_prio_cn(i,loop_x)})
                        id_cluster = ind_prio_cn(i,loop_x);
                        CL_choose = node_choose(:,loop_x);
                        available_nodes = setdiff(available_nodes, CL_choose);
                        stop = 1;
                    else
                        % if demand over the capacity then select the next nearest cluster
                        if i < k
                            i = i+1;
                            stop = 0;
                        else
                            id_cluster = ind_prio_cn(i,loop_x);
                            CL_choose = node_choose(:,loop_x);
                            available_nodes = setdiff(available_nodes, CL_choose);
                            stop = 1;
                        end
                    end
                end
                % Update total demand & available for each cluster as well as the cluster member
                TD{id_cluster} = TD{id_cluster} + Demand(CL_choose);
                CL{id_cluster} = [CL{id_cluster}(1:end) CL_choose];
                available_capacity{id_cluster} = available_capacity{id_cluster}-Demand(CL_choose);
            end
        end
        % insert the initial centroids to the cluster (just for iteration 0)
        if iteration ==0
            for cl = 1 : k
                TD{cl} = TD{cl} + Demand(centroids(cl));
                CL{cl} = [CL{cl}(1:end) centroids(cl)];
            end
        end
        
        % store node for each cluster to iterparam
        % if iteration = 0 then automatically store the cluster value
        if iteration == 0
            for it = 1 : k
                Iter_Param{it}=CL{it};
                T_Demand{it} = TD{it};
                CL{it} = [];
                TD{it} = 0;
                NP{it} = size(Iter_Param{it},2);
            end
        else
            % if iteration > 0 do each step
            for it = 1:k
                if (size(CL{it},2)== size(Iter_Param{it},2))
                    if (sort(CL{it}) == sort(Iter_Param{it}))
                        Iter_Param{it}=CL{it};
                        T_Demand{it} = TD{it};
                        CL{it}=[];
                        TD{it}=0;
                        NP{it} = 0;
                    else
                        Iter_Param{it}=CL{it};
                        T_Demand{it} = TD{it};
                        CL{it}=[];
                        TD{it}=0;
                        NP{it}=NP{it};
                    end
                else
                    Iter_Param{it}=CL{it};
                    T_Demand{it} = TD{it};
                    CL{it}=[];
                    TD{it}=0;
                    NP{it}=size(Iter_Param{it},2);
                end
            end
        end
        % calculate stopping condition
        pattern_checker = [];
        for it = 1:k
            checker = NP{it};
            pattern_checker = [pattern_checker checker];
        end
        k_means_stop = sum(pattern_checker);
        % update iteration
        if k_means_stop > 0
            iteration = iteration +1;
        else
            iteration = iteration;
        end
        if iteration >=3000
            k_means_stop =0;
        end
    end
    
    %=============================FARTHEST_INSERTION_PROCEDURES================
    TSP_FI = cell(1,k);                 % Create Array to Store TSP
    TSP_Distance = cell(1,k);           % Create Array to Store Total Distance
    Cluster_Nodes = Iter_Param;         % Create Array Cluster Nodes from each cluster
    Depot_Node = 1;                     % Initiate the Depot Node
    Depot = newData.Depot;              % Retrieve Depot Location
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
    
    %===============================Tabu====================================================
    %BestImp=TSP_Cost;
    Tour = cell(1,k);                 % Create Array to Store TSP
    Distance = cell(1,k);
    
    %calculate the improvements acheived by all combination of 2-opt and pick the best amongst them
    for cl=1:k
        BestImp=-1;
        TSP = TSP_FI{cl};
        Nb_Nodes= size(TSP,2);
        TSPSize=size(TSP,2);
        
        Cost_all = TSP_Distance{cl};
        Cost_Input = Cost_all;
        
        if (Nb_Nodes > 4)
            [TBCost, TSP_TB] = Two_opt_TS(TSP,Cost_Input,TSPSize,DM);
    
        NewTB_Tour{cl} = TSP_TB;
        NewTB_Distance{cl} = TBCost;
        else
            NewTB_Tour{cl} = TSP;
            NewTB_Distance{cl} = Cost_all;
    end
end
NewTB_Cost=0;
for it = 1 : k
    Cost_TB = NewTB_Distance{it};
    NewTB_Cost = NewTB_Cost + Cost_TB;
end




%======================================================================================================
timer = toc;
StoreTable{length+1,1} = length;
StoreTable{length+1,2} = FileName{length};
StoreTable{length+1,3} = NewTB_Cost;
StoreTable{length+1,4} = timer;
sprintf ("Finished %d", length)
end


xlsFile = 'output_TB.xlsx';
sheetName='TB';
data = StoreTable;
[status, message] = xlswrite(xlsFile, data, sheetName);
dos(['start ' xlsFile]);