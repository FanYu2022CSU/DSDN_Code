function CTAEA(Global)
% <algorithm> <C>
% Two-archive evolutionary algorithm for constrained MOPs

%------------------------------- Reference --------------------------------
% K. Li, R. Chen, G. Fu, and X. Yao, Two-archive evolutionary algorithm for
% constrained multi-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2018, 23(2): 303-315.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
   
    %% Generate random population
    Population = Global.Initialization();
    CA=UpdateCA([],Population,W);            % initial CA
    DA=UpdateDA(CA,[],Population,W);         % initial DA
      %保留数据
    I1 = [];
    S1 = GHD(CA.objs,Global.PF);
    I1 = [I1,S1];
    I2 = [];
    S2 = IGDp(CA.objs,Global.PF);
    I2 = [I2,S2];
    %% Optimization
    while Global.NotTermination(CA)
        Draw(DA.objs,'sk','Markeredgecolor',[.5 .5 .5],'Markerfacecolor',[.8 .5 .5]);
%         if mod(Global.evaluated,Global.N)==0
%             S1 = GHD(CA.objs,Global.PF);
%             I1 = [I1,S1];
%             S2 = IGDp(CA.objs,Global.PF);
%             I2 = [I2,S2];
%         end
        %% mating pool choosing
        % calculate the ratio of non-dominated solutions of CA and DA in Hm
        Hm=[CA,DA];                         
        [FrontNo,~]=NDSort(Hm.objs,inf);
        FrontNo_C=FrontNo(1:ceil(length(Hm)/2));
        Nc=size(find(FrontNo_C==1),2);      
        Pc=Nc/length(Hm);
        FrontNo_D=FrontNo(ceil(length(Hm)/2)+1:length(Hm));
        Nd=size(find(FrontNo_D==1),2);      
        Pd=Nd/length(Hm);

        % calculate the proportion of non-dominated solutions in CA
        [FrontNo,~]=NDSort(CA.objs,inf);
        NC=size(find(FrontNo==1),2);         
        PC=NC/length(CA);                     % PC denotes the proportion of non-dominated solutions in CA,it is different from Pc

        %reproduction
        Q=[];
        for i=1:size(W,1)
            if Pc>Pd
                P1=MatingSelection(CA); 
            else
                P1=MatingSelection(DA);
            end
            pf=rand();
            if pf<PC
                P2=MatingSelection(CA);
            else
                P2=MatingSelection(DA);
            end
            MatingPool=[P1,P2];
            Offspring=GAhalf(MatingPool);
            Q=[Q,Offspring];
        end
        if Global.evaluated >= Global.evaluation
            problem_name = class(Global.problem); % 获取类名
            fileName1 = sprintf('D:\\experiment\\DSDN\\Curve\\CTAEA_%s_GHD.mat', problem_name);
            fileName2 = sprintf('D:\\experiment\\DSDN\\Curve\\CTAEA_%s_IGDp.mat', problem_name);
%             save(fileName1, 'I1', '-V7.3'); % 保存第一个文件
%             save(fileName2, 'I2', '-V7.3'); % 保存第二个文件
        end
       %% update CA and DA
        CA=UpdateCA(CA,Q,W);
        DA=UpdateDA(CA,DA,Q,W);
    end
end