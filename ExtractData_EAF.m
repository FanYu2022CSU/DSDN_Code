
problem_name = 'MW13';  
dim = 15;              
obj_num = 2;          
num_runs = 30;         
num_solutions = 200;  
algorithms = {'NSGAII','NSGAIIDSDN','NSGAIII','NSGAIIIDSDN','SPEA2','SPEA2DSDN','CCMO','CCMODSDN'};  

save_dir = 'D:\experiment\DSDN_EAF\'; 

for alg_idx = 1:length(algorithms)
    algorithm_name = algorithms{alg_idx};
    
    all_results = zeros(num_runs, num_solutions, obj_num);
    
    for run = 1:num_runs
        mat_file = sprintf('D:\\experiment\\DSDN\\PlatEMO\\Data\\%s\\%s_%s_M%d_D%d_%d.mat', algorithm_name, algorithm_name, problem_name, obj_num, dim, run);
        loaded_data = load(mat_file);
        all_results(run, :, :) = loaded_data.result{1, 2}.objs; 
    end
    
    save_file = sprintf('%s%s_%s_results.mat', save_dir, algorithm_name, problem_name);
    save(save_file, 'all_results');
    
end

disp('Complicated');
