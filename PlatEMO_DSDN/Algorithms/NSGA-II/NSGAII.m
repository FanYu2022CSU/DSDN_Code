function NSGAII(Global)
% <algorithm> <N>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
%     Arch = ArchiveUpdate(Population,Global.N);
    %保留数据
    I1 = [];
    S1 = GHD(Population.objs,Global.PF);
    I1 = [I1,S1];
    I2 = [];
    S2 = IGDp(Population.objs,Global.PF);
    I2 = [I2,S2];
%     Arch = ArchiveUpdate(Population,Global.N);
    %% Optimization
    while Global.NotTermination(Population)
        if mod(Global.evaluated,Global.N)==0
            S1 = GHD(Population.objs,Global.PF);
            I1 = [I1,S1];
            S2 = IGDp(Population.objs,Global.PF);
            I2 = [I2,S2];
        end 
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = GA(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
%         Arch = ArchiveUpdate([Arch,Population],Global.N);
        if Global.evaluated >= Global.evaluation
            problem_name = class(Global.problem); % 获取类名
            fileName1 = sprintf('D:\\experiment\\DSDN\\Curve1\\NSGAII_%s_GHD.mat', problem_name);
            fileName2 = sprintf('D:\\experiment\\DSDN\\Curve1\\NSGAII_%s_IGDp.mat', problem_name);
            save(fileName1, 'I1', '-V7.3'); % 保存第一个文件
            save(fileName2, 'I2', '-V7.3'); % 保存第二个文件
        end
    end
end