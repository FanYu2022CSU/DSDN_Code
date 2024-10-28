function NSGAIIDSDN(Global)
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
    [P_pre,ROC_P] =  Global.ParameterSet(1,1);
    isCov = false;
    mu = Global.maxgen/2;
    sigma = 0;
    Population = Global.Initialization();
    Archive = ArchiveUpdate(Population,Global.N);
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N,100);
    I1 = [];
    S1 = GHD(Archive.objs,Global.PF);
    I1 = [I1,S1];
    I2 = [];
    S2 = IGDp(Archive.objs,Global.PF);
    I2 = [I2,S2];
    I3 = [];
    S3 = GHD(Population.objs,Global.PF);
    I3 = [I3,S3];
    I4 = [];
    S4 = IGDp(Population.objs,Global.PF);
    I4 = [I4,S4];
    %% Optimization
    while Global.NotTermination(Archive)
%         Draw(Population.objs,'sk','Markeredgecolor',[.5 .5 .5],'Markerfacecolor',[.8 .5 .5]);
        if mod(Global.evaluated,Global.N)==0
            S1 = GHD(Archive.objs,Global.PF);
            I1 = [I1,S1];
            S2 = IGDp(Archive.objs,Global.PF);
            I2 = [I2,S2];
            S3 = GHD(Population.objs,Global.PF);
            I3 = [I3,S3];
            S4 = IGDp(Population.objs,Global.PF);
            I4 = [I4,S4];
        end 
        R = (Global.D^(1/3)) * (1 - exp(-((Global.gen-mu)^2) / (2 * sigma^2)));
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        if rand() > 0.5
            Offspring  = DE(Population,Archive,Population(MatingPool),{0.5,0.5,0.5,0.75});
        else
            Offspring  = GA(Population(MatingPool));
        end
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,R);
        if Global.evaluated >= Global.evaluation
            problem_name = class(Global.problem); 
            fileName1 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNNSGAII_%s_GHD.mat', problem_name);
            fileName2 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNNSGAII_%s_IGDp.mat', problem_name);
            fileName3 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNNSGAIIP_%s_GHD.mat', problem_name);
            fileName4 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNNSGAIIP_%s_IGDp.mat', problem_name);
            save(fileName1, 'I1', '-V7.3'); 
            save(fileName2, 'I2', '-V7.3'); 
            save(fileName3, 'I3', '-V7.3'); 
            save(fileName4, 'I4', '-V7.3'); 
        end
        Archive = ArchiveUpdate([Archive,Population],Global.N);
        mean_P = sum(sum(Population.objs,1))/Global.N;
        if mod(Global.gen,10)==0
            [ROC_P,P_pre] = Detect(P_pre,mean_P);
        end
        PND = Cal_PND(Population); 
        if PND > 0.99 && ROC_P<1e-3 && ~isCov
            sigma = floor(abs(mu-Global.gen)/3);
            isCov = true;
        end
    end
end
function [ROC,pre] = Detect(pre,now)
    ROC = abs((now-pre)/pre);
    pre = now; 
end
function PND =  Cal_PND(Population)
    [FrontNo,~]=NDSort(Population.objs,inf);  
    ND=size(find(FrontNo==1),2);    
    PND=ND/length(Population)  ;       
end
