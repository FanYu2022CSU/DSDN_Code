classdef RWMOP15 < PROBLEM
% <problem> <RWMOPs>

%------------------------------- Reference --------------------------------
% Abhishek Kumar, Guohua Wu, Mostafa Ali, Qizhang Luo, Rammohan Mallipeddi,
% Ponnuthurai Suganthan, and Swagatam Das, A Benchmark-Suite of Real-World 
% Constrained Multi-Objective Optimization Problems and some Baseline Resu-
% -lts, submitted to Swarm and Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    methods
        %% Initialization
        function obj =RWMOP15()
            par = Cal_par(15);
            obj.Global.M = par.fn;
            obj.Global.D = par.n;
            obj.Global.lower    = par.xmin;
            obj.Global.upper    = par.xmax;
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [PopObj,~,~] = CEC2021_func(X,15);    
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            [~,g,h] = CEC2021_func(X,15);
            PopCon = [g,abs(h)-1e-4];
        end 
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = load('RefPoints/nadir_15.txt');
        end
    end
end