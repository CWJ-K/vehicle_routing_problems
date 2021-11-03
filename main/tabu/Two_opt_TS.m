function [BestCost, TSP ] = Two_opt_TS(TSP,TSPCost,TSPSize,Dis)
% Note that the initial feasible solution is TSP and the corresponding objective function
% value is TSPCost
try
    % the neighborhood structure to use is
    % whole neighborhood of TSP except the tabu moves, for its best neighbor using Two_Opt move
    
    
    
    
    %Initialize iteration counter iter
    iter=1;
    %Initialize  node1 & node2 to -1
    node1=-1;
    node2=-1;
    CurrentCost=TSPCost;
    %Initialize BestCost to TSPCost
    BestCost=TSPCost;
    BestTSP=TSP;
    
    %initialise the Tabu list
    TSize=10;           % change from 10 to 7
    TabuList=zeros(TSize,2);
    % set TCounter to one
    TCounter=1;
    % set the stopping criterion
    Stop=0;
    No_Imp_iter=0;
    Max_NoImp_iter=1000;  % change from 5 to 1000
    
    %REPEAT until stopping condition = true //
    
    while(Stop==0)
        
        Tabuidx=0;% initialise the id of the tabu edges to 0
        
        % search the neighborhood of TSP given the neighborhood structure
        [node1,node2,Improvement_in_Cost,Tabuidx]=Two_Opt_with_TL(TSP,CurrentCost,BestCost,TSPSize,Dis,TabuList);
        
        BestImp=round(Improvement_in_Cost,3);
        % update the tabu list
        
        if(Tabuidx>0)
            Tabuidx= sort(Tabuidx,'descend');
            % delete the moves corresponding to Tabuidx
            TabuList(Tabuidx,:)=[];
            % if only one edge is deleted from the tabu list
            if(size(Tabuidx,2)==1)
                % add an empty row at the end of the tabu list
                TabuList(end+1,:)=[0 0];
                % decrease the TCounter by one
                TCounter=TCounter-1;
                
            else% if two edges are deleted from the tabu list
                
                % add two empty rows at the end of the tabu list
                TabuList(end+1,:)=[0 0];
                TabuList(end+1,:)=[0 0];
                % decrease the TCounter by two
                TCounter=TCounter-2;
            end
        end
        
        % Add the new edges to the tabu list
        TabuList(TCounter:(TCounter+1),:)=[TSP(node1) TSP(node2) ;TSP(node1+1) TSP(node2+1)];
        TCounter=TCounter+2;
        
        % if the tabu list is full
        if(TCounter-1>TSize)
            % delete the oldest ones (first ones)
            TabuList(1:(TCounter-TSize-1),:)=[];
            
            %set the Tcounter to TSize+1
            TCounter=TSize+1;
        end
        
        %update the current tour by performing the move given the best
        %combination
        if BestImp < 0
         NewTSP=[TSP(1:node1) fliplr(TSP(node1+1:node2)) TSP(node2+1:end) ];
        
        %update the current cost
        CurrentCost=CurrentCost+Improvement_in_Cost;
        end
        % if the BestImp is less than zero (i.e. a better neighbor is found
        % that improves the tour by Improvement_in_Cost, update the BestCost
        % and current neighbor TSP
        
        if(  BestCost>CurrentCost)
            BestCost=CurrentCost;

            BestTSP=NewTSP;
            
            
        else % in no improvement found increase No_Imp_iter by one
            No_Imp_iter=No_Imp_iter+1;
            % if No_Imp_iter is greater than Max_NoImp_iter THEN stop
            if(No_Imp_iter>Max_NoImp_iter)
                break;
            end
        end
        
        
        TSP=NewTSP;
        NewTSP=[];
        
    end
end
end

