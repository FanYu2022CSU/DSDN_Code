function CCMODSDN(Global)
% <algorithm> <C>
% Coevolutionary constrained multi-objective optimization framework
%------------------------------- Reference --------------------------------
% Y. Tian, T. Zhang, J. Xiao, X. Zhang, and Y. Jin, A coevolutionary
% framework for constrained multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [P_pre,ROC_P] =  Global.ParameterSet(1,1);
    isCov = false;
    mu = Global.maxgen/2;
    sigma = 0;
    %% Generate random population
    Population1 = Global.Initialization();
    Population2 = Global.Initialization();
    Archive = ArchiveUpdate([Population1,Population2],Global.N);
    Fitness1    = CalFitness1(Population1.objs,100,Population1.decs,Population1.cons);
    Fitness2    = CalFitness(Population2.objs);
    I1 = [];
    S1 = GHD(Archive.objs,Global.PF);
    I1 = [I1,S1];
    I2 = [];
    S2 = IGDp(Archive.objs,Global.PF);
    I2 = [I2,S2];
    I3 = [];
    S3 = GHD(Population1.objs,Global.PF);
    I3 = [I3,S3];
    I4 = [];
    S4 = IGDp(Population1.objs,Global.PF);
    I4 = [I4,S4];
    %% Optimization
    while Global.NotTermination(Archive)
%         Draw(Population1.objs,'sk','Markeredgecolor',[.5 .5 .5],'Markerfacecolor',[.8 .5 .5]);
        if mod(Global.evaluated,Global.N)==0
            S1 = GHD(Archive.objs,Global.PF);
            I1 = [I1,S1];
            S2 = IGDp(Archive.objs,Global.PF);
            I2 = [I2,S2];
            I2 = [I2,S2];
            S3 = GHD(Population1.objs,Global.PF);
            I3 = [I3,S3];
            S4 = IGDp(Population1.objs,Global.PF);
            I4 = [I4,S4];
        end
        if rand() > 0.5
            MatingPool1 = TournamentSelection(2,Global.N,Fitness1);
            MatingPool2 = TournamentSelection(2,Global.N,Fitness2);
            Offspring1  = DE(Archive(1:end/2),Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)),{0.5,0.5,0.5,0.75});
            Offspring2  = DE(Archive(1:end/2),Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)),{0.5,0.5,0.5,0.75});
        else
            MatingPool1 = TournamentSelection(2,Global.N,Fitness1);
            MatingPool2 = TournamentSelection(2,Global.N,Fitness2);
            Offspring1  = GAhalf(Population1(MatingPool1));
            Offspring2  = GAhalf(Population2(MatingPool2));
        end
        R = (Global.D^(1/3)) * (1 - exp(-((Global.gen-mu)^2) / (2 * sigma^2)));
        [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Global.N,true,R);
        [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Global.N,false,R);
        if Global.evaluated >= Global.evaluation
            I1 = [I1,S1];
            I2 = [I2,S2];
            I3 = [I3,S3];
            problem_name = class(Global.problem); 
            fileName1 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNCCMO_%s_GHD.mat', problem_name);
            fileName2 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNCCMO_%s_IGDp.mat', problem_name);
            fileName3 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNCCMOP_%s_GHD.mat', problem_name);
            fileName4 = sprintf('D:\\experiment\\DSDN\\Curve\\DSDNCCMOP_%s_IGDp.mat', problem_name);
            save(fileName1, 'I1', '-V7.3'); 
            save(fileName2, 'I2', '-V7.3'); 
            save(fileName3, 'I3', '-V7.3'); 
            save(fileName4, 'I4', '-V7.3'); 
        end
        Archive = ArchiveUpdate([Archive,Population1,Population2],Global.N);
        Archive(length(Archive)+1:Global.N) = Population1(randsample(Global.N,Global.N-length(Archive)));
        mean_P = sum(sum(Population1.objs,1))/Global.N;
        if mod(Global.gen,10)==0
            [ROC_P,P_pre] = Detect(P_pre,mean_P);
        end
        PND = Cal_PND(Population1); 
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