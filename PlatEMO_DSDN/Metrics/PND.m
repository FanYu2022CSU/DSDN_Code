function Score = PND(PopObj,PF)
% <metric> <min>
% Generational distance

%------------------------------- Reference --------------------------------
% D. A. Van Veldhuizen, Multiobjective evolutionary algorithms:
% Classifications, analyses, and new innovations, Ph.D. thesis, Department
% of Electrical and Computer Engineering, Graduate School of Engineering,
% Air Force Institute of Technology, Wright Patterson Air Force Base, 1999.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Score    = Cal_PND(PopObj);
end

function PND =  Cal_PND(PopObj)
    [FrontNo,~]=NDSort(PopObj,inf);    %快速非支配排序
    ND=size(find(FrontNo==1),2);    %非支配层个体数量
    PND=ND/length(PopObj)  ;       %非支配个体在种群中所占比率
end