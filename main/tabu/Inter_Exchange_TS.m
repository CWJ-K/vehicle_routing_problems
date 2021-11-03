function [TSP_IR_Cost ]=Inter_Exchange_TS(Demand,TBCost,TSP_TB)
try
    Demand_Ref = [0 Demand];
    cost_order = [];
    Tot_Demand_Cl = T_Demand;
    TSP_IR_Route = TSP_New_Route;
    TSP_IR_Cost = TSP_New_Cost;
    for it = 1:k
        cost_insert = TSP_IR_Cost{it};
        cost_order = [cost_order cost_insert];
    end
    [cost_val_ord,cl_idx_order] = sort(cost_order,'descend');
    for cl_idx = 1 : k-1
        x = 0;
        % Create array to store the value of reducement
        C1_Idx=[];
        C2_Idx=[];
        P1_Idx=[];
        P2_Idx=[];
        Cost_Reduce=[];
        vertex = 0;
        P1 = TSP_IR_Route{cl_idx_order(cl_idx)};
        D1 = Tot_Demand_Cl{cl_idx_order(cl_idx)};
        C1 = TSP_IR_Cost(cl_idx);
        step = cl_idx+1;
        for cl_idy = step : k
            D2 = Tot_Demand_Cl{cl_idx_order(cl_idy)};
            P2 = TSP_IR_Route{cl_idx_order(cl_idy)};
            C2 = TSP_IR_Cost(cl_idy);
            % Perform inter movement
            % node exchange
            cost_difference = 0;
            for p1_idx = 2 : size(P1,2)-1
                node1=P1(p1_idx);
                node1_demand = Demand_Ref(node1);
                for p2_idx = 2 : size(P2,2)-1
                    node2=P2(p2_idx);
                    node2_demand = Demand_Ref(node2);
                    demand1 = D1 - node1_demand + node2_demand ;
                    demand2 = D2-node2_demand+node1_demand;
                    go_sign=0;
                    go_sign_1=0;
                    if (demand1<=Capacity && demand2<=Capacity)
                        go_sign = 1;
                        breakcost_P1 = DM(P1(p1_idx-1),P1(p1_idx))+DM(P1(p1_idx),P1(p1_idx+1));
                        addcost_P1 = DM(P1(p1_idx-1),node2)+DM(node2,P1(p1_idx+1));
                        deltacost_P1 = addcost_P1 - breakcost_P1;
                        
                        breakcost_P2 = DM(P2(p2_idx-1),P2(p2_idx))+DM(P2(p2_idx),P2(p2_idx+1));
                        addcost_P2 = DM(P2(p2_idx-1),node1)+DM(node1,P2(p2_idx+1));
                        deltacost_P2 = addcost_P2 - breakcost_P2;
                        deltacost = deltacost_P1 + deltacost_P2;
                        initdelta_cost = 0;
                        if (deltacost< 0)
                            if deltacost < initdelta_cost
                                initdelta_cost = deltacost;
                                cl1_index = p1_idx;
                                cl2_index = p2_idx;
                                go_sign_1 = 1;
                            end
                        end
                    end
                end
            end
            if (go_sign==1)
                if (go_sign_1==1)
                    C1_Idx=cl_idx_order(cl_idx);
                    C2_Idx=cl_idx_order(cl_idy);
                    P1_Idx(x+1)=cl1_index;
                    P2_Idx(x+1)=cl2_index;
                    Cost_Reduce(x+1)=deltacost;
                    go_sign_1=0;
                    x = x+1;
                end
            end
        end
        
        if (size(Cost_Reduce,2)>0)
            % perform moves k = 1
            [value_reducement,idx_moves] = min(Cost_Reduce);
            nodes1 = TSP_IR_Route{C1_Idx}(P1_Idx(idx_moves));
            nodes2 = TSP_IR_Route{C2_Idx}(P2_Idx(idx_moves));
            TSP_IR_Route{C1_Idx}=[TSP_IR_Route{C1_Idx}(1:P1_Idx(idx_moves)-1) nodes2 TSP_IR_Route{C1_Idx}(P1_Idx(idx_moves)+1:end)] ;
            TSP_IR_Route{C2_Idx}=[TSP_IR_Route{C2_Idx}(1:P2_Idx(idx_moves)-1) nodes1 TSP_IR_Route{C2_Idx}(P2_Idx(idx_moves)+1:end)] ;
            TSP_IR = 0;
            TSP_IR_2=0;
            for iter = 1 : size(TSP_IR_Route{C1_Idx},2)-1
                IR_Cost = DM(TSP_IR_Route{C1_Idx}(iter),TSP_IR_Route{C1_Idx}(iter+1));
                TSP_IR = TSP_IR+IR_Cost;
            end
            for iter = 1 : size(TSP_IR_Route{C2_Idx},2)-1
                IR_Cost_2 = DM(TSP_IR_Route{C2_Idx}(iter),TSP_IR_Route{C2_Idx}(iter+1));
                TSP_IR_2 = TSP_IR_2+IR_Cost_2;
            end
            TSP_IR_Cost{C1_Idx}=TSP_IR;
            TSP_IR_Cost{C2_Idx}=TSP_IR_2;
            Tot_Demand_Cl{C1_Idx}=Tot_Demand_Cl{C1_Idx}-Demand_Ref(nodes1)+Demand_Ref(nodes2);
            Tot_Demand_Cl{C2_Idx}=Tot_Demand_Cl{C2_Idx}-Demand_Ref(nodes2)+Demand_Ref(nodes1);
        end
    end
end
end