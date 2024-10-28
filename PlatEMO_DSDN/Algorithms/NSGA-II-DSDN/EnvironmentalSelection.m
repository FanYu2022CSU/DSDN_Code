function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,R)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort_Mod(Population.objs,Population.cons,N,Population.decs,R);
    [FrontNo] = Move(Population,FrontNo);
    for f = 1:MaxFNo+1
        if sum(FrontNo<=f)>=N
            MaxFNo = f;
            break;
        end
    end
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end



function [FrontNo] = Move(Population,FrontNo)
       value=Population.objs;
    [T,M]=size(value);
    First = FrontNo <=1;%记录第一层个体的序号
    %在record当中，第一个表示个体的序号，第二行表示对应的目标值
    Record=zeros(4,T);%记录第一层个体的序号
    
    
    %下面表示记录相关的数据
    num=0;
    for i=1:T
        if First(1,i)==1
            num=num+1;
            Record(1,num)=i;
        end
    end
    Record=Record(:,1:num);
    
    %计算第一层中所有个体的目标函数值之和，放到Record当中的第二行
    for i=1:num
        s=0;
        for k=1:M
            s=s+value(Record(1,i),k);
        end
        Record(2,i)=s;
    end
    small=min(Record(2,:));%计算最小的目标函数值之和
    
    for i=1:num
        Record(3,i)=Record(2,i)/small;
        if Record(3,i)>4% 大于5倍的个体放到第二行去
            Record(4,i)=2;
        end
    end
    
    % 修改原来的矩阵
    for i=1:num
        if Record(4,i)==2
            n=Record(1,i);
            FrontNo(1,n)=2;
        end
    end
end
