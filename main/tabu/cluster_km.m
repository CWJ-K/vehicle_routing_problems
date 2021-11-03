function [k,Iter_Param,Sample_Location,Depot] = cluster_km(FileName)
load dat;
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
end

