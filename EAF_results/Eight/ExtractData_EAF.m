% MATLAB 脚本：提取算法运行结果并保存为合适的格式

% 设置相关路径和参数
problem_name = 'MW';  % 问题名称
dim = 15;              % 问题的维度
obj_num = 2;           % 目标函数的数量
num_runs = 30;         % 运行次数
num_solutions = 200;   % 每次运行的解的数量
algorithms = {'NSGAII','NSGAIIDSDN','NSGAIII','NSGAIIIDSDN','SPEA2','SPEA2DSDN','CCMODSDN'};  % 算法名称列表

% 创建保存结果的目录
save_dir = 'D:\experiment\DSDN_EAF\';  % 保存目录

% 循环提取每个算法的运行结果
for alg_idx = 1:length(algorithms)
    algorithm_name = algorithms{alg_idx};
    
    % 创建一个存储当前算法所有运行结果的矩阵，尺寸为 (num_runs, num_solutions, obj_num)
    all_results = zeros(num_runs, num_solutions, obj_num);
    
    % 遍历每次运行，提取结果
    for run = 1:num_runs
        % 生成对应的文件名（假设文件命名格式：算法_问题_M目标数_D维度_运行次数.mat）
        mat_file = sprintf('D:\\experiment\\PlatEMO-4.8\\PlatEMO-master\\PlatEMO\\Data\\%s\\%s_%s_M%d_D%d_%d.mat', algorithm_name, algorithm_name, problem_name, obj_num, dim, run);
        
        % 加载 .mat 文件
        loaded_data = load(mat_file);
        
        % 提取结果：假设结果存储在 result{1,2}.objs 中
        all_results(run, :, :) = loaded_data.result{1, 2}.objs;  % 200个解，2个目标函数值
    end
    
    % 将当前算法的运行结果保存为 .mat 文件，方便Python读取
    save_file = sprintf('%s%s_%s_results.mat', save_dir, algorithm_name, problem_name);
    save(save_file, 'all_results');
    
    % 也可以导出为 Python 友好的 .csv 或 .npz 格式
    % csvwrite() 或者使用 mat2npz 工具来保存为 Python 可直接读取的格式
end

disp('所有算法的运行结果提取并保存完成。');
