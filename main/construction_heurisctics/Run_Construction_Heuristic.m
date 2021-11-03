clear
Files = dir('*.mat'); 
FileName = {Files.name};
StoreTable = {'No','Problem','Optimal TSP Cost','NN Cost','Over Optimal(%)','AI Cost','Over Optimal(%)','FI Cost','Over Optimal(%)','NI Cost',...
    'Over Optimal(%)','CI Cost','Over Optimal(%)','NM Cost','Over Optimal(%)','CWS Cost','Over Optimal(%)'};
StoreTime = {'No','Problem','NN time','AI time', 'FI time','NI time','CI time','NM time','CWS time','Cluster time'};
for r = 1:length(Files)
    newData = load(FileName{r}, 'Demand','Sample_Location','Capacity','Depot','Best_Value');
    Sample_Location = newData.Sample_Location;
    Demand =newData.Demand;
    Capacity = newData.Capacity;
    Depot = newData.Depot;  
    Best_Value=newData.Best_Value;
    %save dat;
    tic;
    save dat;
%=======================Cluster
    [k,Iter_Param,Sample_Location,Depot] = cluster_km(FileName);
    cluster_t=toc;
%=========================================================================%    
%=======================Construction Heuristics
    %Cluster_NN
    tic;
    [NN_Cost]=Cluster_NN(k,Iter_Param,Sample_Location,Depot);
    NN_t=toc;
    %Arbitrary_Insertion
    tic;
    [AI_Cost]=Arbitrary_Insertion(k,Iter_Param,Sample_Location,Depot);
    AI_t=toc;
    %Farthest_Insertion
    tic;
    [FI_Cost]=Farthest_Insertion(k,Iter_Param,Sample_Location,Depot);
    FI_t=toc;
    %Nearest Insertion
    tic;
    [NI_Cost]= Nearest_Insertion(k,Iter_Param,Sample_Location,Depot);
    NI_t=toc;
    %Cheapest Insertion
    tic;
    [CI_Cost]=Cheapest_Insertion(k,Iter_Param,Sample_Location,Depot);
    CI_t=toc;
    %Nearst_Merger
    tic;
    [NM_Cost]= Nearst_Merger(k,Iter_Param,Sample_Location,Depot);
    NM_t=toc;
    %Clarck & Wright Saving
    tic;
    [CWS_Cost] = CWS(k,Iter_Param,Sample_Location,Depot);
    CWS_t=toc;
%=========================================================================%   
%=======================Classical Local Search $ Iterative Local Search 
    %Steepest Descent    
     
    %Random Descent 
    
    %Percentage Random Descent 
    
%=========================================================================%   
%=======================Metaheuristics 
    %Simulated Annealing
    
    %Tabu Search 
    
    %Variable Neighborhood Search

%=========================================================================%   
%=======================Store Results
%=====Cost
 Optimal_Cost{r}=Best_Value;
%=====Cost
%Construction Heuristics
NNCost{r}=NN_Cost;
Improvement_NN{r}=(NNCost{r}(1)/Optimal_Cost{r}(1))-1;
AICost{r}=AI_Cost;
Improvement_AI{r}=(AICost{r}(1)/Optimal_Cost{r}(1))-1;
FICost{r}=FI_Cost;
Improvement_FI{r}=(FICost{r}(1)/Optimal_Cost{r}(1))-1;
NICost{r}=NI_Cost;
Improvement_NI{r}=(NICost{r}(1)/Optimal_Cost{r}(1))-1;
CICost{r}=CI_Cost;
Improvement_CI{r}=(CICost{r}(1)/Optimal_Cost{r}(1))-1;
NMCost{r}=NM_Cost;
Improvement_NM{r}=(NMCost{r}(1)/Optimal_Cost{r}(1))-1;
CWSCost{r}=CWS_Cost;
Improvement_CWS{r}=(CWSCost{r}(1)/Optimal_Cost{r}(1))-1;

% %Classical Local Search $ Iterative Local Search 
% 
% %Metaheuristics 
% 
% 
%%=====Time
 cluster_time{r}=cluster_t;
%%=====Time
% % %Construction Heuristics
NN_time{r}=NN_t;
AI_time{r}=AI_t;
FI_time{r}=FI_t;
NI_time{r}=NI_t;
CI_time{r}=CI_t;
NM_time{r}=NM_t;
CWS_time{r}=CWS_t;
%%=====Excel Part: write result into Excel
%Cost
StoreTable{r+1,1} = r;
StoreTable{r+1,2} = FileName{r};
StoreTable{r+1,3} = Best_Value;
StoreTable{r+1,4} = NN_Cost;
StoreTable{r+1,5} = Improvement_NN{r};
StoreTable{r+1,6} = AI_Cost;
StoreTable{r+1,7} = Improvement_AI{r};
StoreTable{r+1,8} = FI_Cost;
StoreTable{r+1,9} = Improvement_FI{r};
StoreTable{r+1,10} = NI_Cost;
StoreTable{r+1,11} = Improvement_NI{r};
StoreTable{r+1,12} = CI_Cost;
StoreTable{r+1,13} = Improvement_CI{r};
StoreTable{r+1,14} = NM_Cost;
StoreTable{r+1,15} = Improvement_NM{r};
StoreTable{r+1,16} = CWS_Cost;
StoreTable{r+1,17} = Improvement_CWS{r};
%Time
StoreTime{r+1,1} = r;
StoreTime{r+1,2} = FileName{r};
StoreTime{r+1,3} = NN_t;
StoreTime{r+1,4} = AI_t;
StoreTime{r+1,5} = FI_t;
StoreTime{r+1,6} = NI_t;
StoreTime{r+1,7} = CI_t;
StoreTime{r+1,8} = NM_t;
StoreTime{r+1,9} = CWS_t;
StoreTime{r+1,10} = cluster_t;
%
sprintf ("Finished %d", r)
end


xlsFile = 'Construction_Heuristics_old_cluster.xlsx';
sheetName1='Cost';
sheetName2='Time';
data1 = StoreTable;
data2=StoreTime;
[status, message] = xlswrite(xlsFile, data1, sheetName1);
[status, message] = xlswrite(xlsFile, data2, sheetName2);


dos(['start ' xlsFile]);

    
    

