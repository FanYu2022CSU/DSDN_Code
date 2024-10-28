function [par] = Cal_par(prob_k)
% CEC2021 Real world multi-objective Constrained Optimization Test Suite 
% Abhishek Kumar (email: abhishek.kumar.eee13@iitbhu.ac.in, Indian Institute of Technology (BHU), Varanasi) 

% prob_k -> Index of problem.
% par.n  -> Dimension of the problem.
% par.fn -> Number of objective.
% par.g  -> Number of inequility constraints.
% par.h  -> Number of equality constraints.
% par.xmin -> lower bound of decision variables.
% par.xmax -> upper bound of decision variables.


D        = [4,5,3,4,4,7,4,7,4,2,3,4,7,5,3,2,6,3,10,4,6,9,6,9,2,3,3,7,7,25,25,25,30,30,30,28,28,28,28,34,34,34,34,34,34,34,18,18,18,6];
par.n    = D(prob_k); 
O        = [2,2,2,2,2,2,2,3,2,2,5,2,3,2,2,2,3,2,3,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,3,2,2,3,3,4,2,2,3,2];
par.fn   = O(prob_k);
gn       = [2,5,3,4,4,11,1,9,1,2,7,1,11,8,8,2,9,3,10,7,4,2,1,0,2,1,3,4,9,24,24,24,29,29,29,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
hn       = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,6,0,1,0,4,0,0,0,0,0,0,0,24,24,24,24,26,26,26,26,26,26,26,12,12,12,1];
par.gn    = gn;
par.hn    = hn;
par.g     = gn(prob_k);
par.h     = hn(prob_k);
%% range
% bound constraint definitions for all 18 test functions
xmin1    = [0.51,0.51,10,10];
xmax1    = [99.49,99.49,200,200];
xmin2    = [0.05,0.2,0.2,0.35,3];
xmax2    = [0.5,0.5,0.6,0.5,6];
xmin3    = [1e-5,1e-5,1];
xmax3    = [100, 100, 3];
xmin4    = [0.125,0.1,0.1,0.125];
xmax4    = [5,10,10,5];
xmin5    = [55,75,1000,11];
xmax5    = [80,110,3000,20];
xmin6    = [2.6,0.7,16.51,7.3,7.3,2.9,5];
xmax6    = [3.6,0.8,28.49,8.3,8.3,3.9,5.5];
xmin7    = [11.51,11.51,11.51,11.51];
xmax7    = [60.49,60.49,60.49,60.49];
xmin8    = [0.5,0.45,0.5,0.5,0.875,0.4,0.4];
xmax8    = [1.5,1.35,1.5,1.5,2.625,1.2,1.2];
xmin9    = [1,sqrt(2),sqrt(2),1];
xmax9    = [3,3,3,3];
xmin10   = [0.1,0.5];
xmax10   = [2,2.5];
xmin11   = [0.01 0.01 0.01];
xmax11   = [0.45 0.1 0.1];
xmin12   = [10 10 0.9 0.9];
xmax12   = [80 50 5 5];  
xmin13    = [2.6,0.7,16.51,7.3,7.3,2.9,5];
xmax13    = [3.6,0.8,28.49,8.3,8.3,3.9,5.5];
xmin14    =  [60 90 1 0 2];
xmax14    = [80 110 3 1000 9];
xmin15    = [0.51,0.6,0.51];
xmax15    = [70.49,3,42.49];
xmin16    = [0.01 0.20];
xmax16    = [0.05 1];
xmin17    = [150.0 20.0 13.0 10.0 14.0 0.63];
xmax17    = [274.32 32.31 25.0 11.71 18.0 0.75];
xmin18    = [136,56,1.4];
xmax18    = [146,68,2.2];
xmin19   = [ 0.51,0.51,0.51,250,250,250,6,4,40,10];
xmax19   = [3.49,3.49,3.49,2500,2500,2500,20,16,700,450];
xmin20   = [ 1, 1,  1e-6,1];
xmax20   = [16, 16, 16*1e-6,16];
xmin21 = [1.3, 2.5, 1.3, 1.3, 1.3, 1.3];
xmax21 = [1.7, 3.5, 1.7, 1.7, 1.7, 1.7];
xmin22 = zeros(1,9);
xmax22 = [100,200,100,100,100,100,200,100,200];
xmin23 = [0,0,0,0,0.00001,0.00001];
xmax23 = [ 1,1,1,1,16,16];
xmin24 = [0,0,0,0,1000,0,100,100,100];
xmax24 = [10,200,100,200,2000000,600,600,600,900];
xmin25 = [0,-0.49];
xmax25 = [1.6,1.49];
xmin26 = [0.5,0.5,-0.49];
xmax26 = [1.4,1.4,1.49];
xmin27 = [0.2,-2.22554,-0.49];
xmax27 = [1,-1,1.49];
xmin28 = [0,0,0,0,-0.49,-0.49,0];
xmax28 = [20,20,10,10,1.49,1.49,40];
xmin29 = [0,0,0,-0.49,-0.49,-0.49,-0.49];
xmax29 = [100,100,100,1.49,1.49,1.49,1.49];
xmin30   = -0*ones(1,par.n);
xmax30   = +90*ones(1,par.n);
xmin31   = -0*ones(1,par.n);
xmax31   = +90*ones(1,par.n);
xmin32   = -0*ones(1,par.n);
xmax32   = +90*ones(1,par.n);
xmin33   = -0*ones(1,par.n);
xmax33   = +90*ones(1,par.n);
xmin34   = -0*ones(1,par.n);
xmax34   = +90*ones(1,par.n);
xmin35   = -0*ones(1,par.n);
xmax35   = +90*ones(1,par.n);
xmin36   = -1*ones(1,par.n);xmin36(25:28) = 0;
xmax36   = +1*ones(1,par.n);
xmin37   = -1*ones(1,par.n);xmin37(25:28) = 0;
xmax37   = +1*ones(1,par.n);
xmin38   = -1*ones(1,par.n);xmin38(25:28) = 0;
xmax38   = +1*ones(1,par.n);
xmin39   = -1*ones(1,par.n);xmin39(25:28) = 0;
xmax39   = +1*ones(1,par.n);
xmin40   = -1*ones(1,par.n);xmin40(27:34) = 0;
xmax40   = +1*ones(1,par.n);
xmin41   = -1*ones(1,par.n);xmin41(27:34) = 0;
xmax41   = +1*ones(1,par.n);
xmin42   = -1*ones(1,par.n);xmin42(27:34) = 0;
xmax42   = +1*ones(1,par.n);
xmin43   = -1*ones(1,par.n);xmin43(27:34) = 0;
xmax43   = +1*ones(1,par.n);
xmin44   = -1*ones(1,par.n);xmin44(27:34) = 0;
xmax44   = +1*ones(1,par.n);
xmin45   = -1*ones(1,par.n);xmin45(27:34) = 0;
xmax45   = +1*ones(1,par.n);
xmin46   = -1*ones(1,par.n);xmin46(27:34) = 0;
xmax46   = +1*ones(1,par.n);
xmin47   = -1*ones(1,par.n);xmin47(11:12) = 0;xmin47(13:18) = 0;
xmax47   = +1*ones(1,par.n);xmax47(11:12) = 2;xmax47(13:18) = 500;
xmin48   = -1*ones(1,par.n);xmin48(11:12) = 0;xmin48(13:18) = 0;
xmax48   = +1*ones(1,par.n);xmax48(11:12) = 2;xmax48(13:18) = 500;
xmin49   = -1*ones(1,par.n);xmin49(11:12) = 0;xmin49(13:18) = 0;
xmax49   = +1*ones(1,par.n);xmax49(11:12) = 2;xmax49(13:18) = 500;
xmin50    = [10 10 35 35 125 130];
xmax50    = [125 150 210 225 315 325];

eval(['par.xmin=xmin' int2str(prob_k) ';']);
eval(['par.xmax=xmax' int2str(prob_k) ';' ]);
end