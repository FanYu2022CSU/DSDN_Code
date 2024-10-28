function NSGAIII(Global)
% <algorithm> <N>
% Nondominated sorting genetic algorithm III

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    [Z,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    Zmin         = min(Population(all(Population.cons<=0,2)).objs,[],1);
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
%         Draw(Global.PF,'ok','Markeredgecolor',[.7 .0 .0],'Markerfacecolor',[.0 .0 .0],'MarkerSize',1);
%         Draw(Population.objs,'sk','Markeredgecolor',[.1 .1 .1],'Markerfacecolor',[.0 .6 .6]);
        MatingPool = TournamentSelection(2,Global.N,sum(max(0,Population.cons),2));
        Offspring  = GA(Population(MatingPool));
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        Population = EnvironmentalSelection([Population,Offspring],Global.N,Z,Zmin);
%         Arch = ArchiveUpdate([Arch,Population],Global.N);
        if Global.evaluated >= Global.evaluation
            problem_name = class(Global.problem); % 获取类名
            fileName1 = sprintf('D:\\experiment\\DSDN\\Curve\\NSGAIII_%s_GHD.mat', problem_name);
            fileName2 = sprintf('D:\\experiment\\DSDN\\Curve\\NSGAIII_%s_IGDp.mat', problem_name);
            save(fileName1, 'I1', '-V7.3'); % 保存第一个文件
            save(fileName2, 'I2', '-V7.3'); % 保存第二个文件
        end
    end
end