function CMOEA_SDE(Global)
% <algorithm> <A>
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
%J. Zhou, Y. Zhang, J. Zheng and M. Li, "Domination-based Selection and Shift-based Density Estimation for Constrained Multiobjective Optimization," 
%in IEEE Transactions on Evolutionary Computation, 2022, doi: 10.1109/TEVC.2022.3190401.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
     %% Parameter setting
    type = Global.ParameterSet(1);
    %% Generate random population
    Population = Global.Initialization();
    CV=sum(max(Population.cons,0),2);
    r=Global.evaluated/Global.evaluation;
    Fitness = CalFitness(Population.objs,CV,r);
    arch = ArchiveUpdate(Population,Global.N);
    %保留数据
    I1 = [];
    S1 = GHD(Population.objs,Global.PF);
    I1 = [I1,S1];
    I2 = [];
    S2 = IGDp(Population.objs,Global.PF);
    I2 = [I2,S2];
    %% Optimization
    while Global.NotTermination(Population)
%         Draw(Global.PF,'ok','Markeredgecolor',[.7 .0 .0],'Markerfacecolor',[.0 .0 .0],'MarkerSize',1);
%         Draw(Population.objs,'sk','Markeredgecolor',[.1 .1 .1],'Markerfacecolor',[.0 .6 .6]);
        if mod(Global.evaluated,Global.N)==0
            S1 = GHD(Population.objs,Global.PF);
            I1 = [I1,S1];
            S2 = IGDp(Population.objs,Global.PF);
            I2 = [I2,S2];
        end
         r=Global.evaluated/Global.evaluation;
         if type == 1
           MatingPool = TournamentSelection(2,Global.N,Fitness);  
           Offspring  = GA(Population(MatingPool));
         elseif type == 2
           Mat1 = TournamentSelection(2,Global.N,Fitness);
           Mat2 = TournamentSelection(2,Global.N,Fitness);
           Offspring = DE(Population,Population(Mat1),Population(Mat2));
         end
        [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Global.N,r);
		 % Output the non-dominated and feasible solutions.
         arch = ArchiveUpdate([arch,Population],Global.N);
        if Global.evaluated >= Global.evaluation
            Population = arch;
            problem_name = class(Global.problem); % 获取类名
            fileName1 = sprintf('D:\\experiment\\DSDN\\Curve\\CMOEASDE_%s_GHD.mat', problem_name);
            fileName2 = sprintf('D:\\experiment\\DSDN\\Curve\\CMOEASDE_%s_IGDp.mat', problem_name);
            save(fileName1, 'I1', '-V7.3'); % 保存第一个文件
            save(fileName2, 'I2', '-V7.3'); % 保存第二个文件
        end
    end
end