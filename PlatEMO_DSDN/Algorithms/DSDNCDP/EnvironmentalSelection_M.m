function [Population,FrontNo] = EnvironmentalSelection_M(Global,Population,N,Z,Zmin,R)
% The environmental selection of NSGA-III

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort_Mod(Population.objs,Population.cons,N,Population.decs,R);
    [FrontNo] = Move(Population,FrontNo);
    % 更新临界层
    for f = 1:MaxFNo+1
        if sum(FrontNo<=f)>=N
            MaxFNo = f;
            break;
        end
    end
    %% Select the solutions in the last front
    if Global.M == 2
        Next = FrontNo <= MaxFNo;
        Del  = Truncation(Population(Next).objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    else
        Next = FrontNo < MaxFNo;
        Last   = find(FrontNo==MaxFNo);
        Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Z,Zmin);
        Next(Last(Choose)) = true;
    end
    % Population for next generation
    Population = Population(Next);
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end
function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
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

