clear
Files = dir('*.txt'); 
numfiles = length(Files);
mydata = cell(1, numfiles);
FileName = {Files.name}

for k = 1:numfiles 

[c1 c2 c3] = textread(FileName{k},'%s %s %s','headerlines',3);


n=str2num(char(c3(1)));
x=str2num(char(c2(6:5+n-1)));
y=str2num(char(c3(6:5+n-1)));
Sample_Location=[x y];
dx=str2num(char(c2(7+n:5+2*n)));
Demand=[dx'];
xd=str2num(char(c2(5)));
yd=str2num(char(c3(5)));
Depot=[xd yd];
Number_Nodes=n-1;
Capacity=str2num(char(c3(3)));

%'Best_Value','Number_Vehicle'
[b] = textread(FileName{k},'%s','delimiter',':');
Best_Value=str2num(char(erase(b(6),")")));



FileMat = [extractBefore(FileName{k},'.') +'.mat'];
save(FileMat,'Sample_Location','Demand','Depot','Number_Nodes','Best_Value','Capacity');

end


    
    
    
  
