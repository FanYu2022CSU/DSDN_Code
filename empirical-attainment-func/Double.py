import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from eaf import get_empirical_attainment_surface, EmpiricalAttainmentFuncPlot
import os


algorithms = ['CCMO', 'CCMODSDN']
algorithms1 = ['CCMO', 'DSDN-CCMO']
problem_name = 'MW13'  
data_dir = 'D:\\experiment\\DSDN_EAF'  
dim, n_samples, n_runs = 2, 200, 30  

def load_algorithm_data(algorithm: str, problem: str) -> np.ndarray:
    file_path = os.path.join(data_dir, f'{algorithm}_{problem}_results.mat')
    try:
        mat_data = scipy.io.loadmat(file_path)
        return mat_data['all_results']  
    except Exception as e:
        print(f"Error loading data for {algorithm}: {e}")
        return None

def add_random_noise(costs: np.ndarray, noise_level=0.01) -> np.ndarray:
    noise = np.random.normal(loc=0.0, scale=noise_level, size=costs.shape)
    return costs + noise

if __name__ == "__main__":
    n_algorithms = len(algorithms)

    all_costs = []
    for alg in algorithms:
        data = load_algorithm_data(alg, problem_name)
        if data is not None:
            all_costs.append(data)
        else:
            print(f"No data for {alg}, skipping.")

    assert len(all_costs) == n_algorithms, "Not all algorithms' data were loaded correctly."

    transformed_costs = []
    for costs in all_costs:
        std_per_obj = np.std(costs, axis=(0, 1))  

        scaled_costs = np.log(1 + (costs / std_per_obj))

        transformed_costs.append(scaled_costs)

    min_value = np.min([np.min(costs) for costs in transformed_costs])
    adjusted_costs = [costs - min_value for costs in transformed_costs]

    combined_costs = np.concatenate(adjusted_costs, axis=0)  
    print("Combined costs shape:", combined_costs.shape)

    combined_costs_with_noise = add_random_noise(combined_costs, noise_level=0.01)


    shared_levels = [1, n_algorithms * n_runs]  
    shared_surfs = get_empirical_attainment_surface(costs=combined_costs_with_noise, levels=shared_levels)


    algorithm_colors = ["red", "green"]  


    _, ax = plt.subplots()

    eaf_plot = EmpiricalAttainmentFuncPlot()
    shared_labels = ["shared best attainment", "shared worst attainment"]  
    shared_colors = ["black", "gray"]  
    eaf_plot.plot_multiple_surface(ax, colors=shared_colors, labels=shared_labels, surfs=shared_surfs)

    for i in range(n_algorithms):
        levels = [n_runs // 4, n_runs // 2, 3 * n_runs // 4] 
        surfs = get_empirical_attainment_surface(costs=adjusted_costs[i], levels=levels)

        eaf_plot.plot_surface_with_band(ax, color=algorithm_colors[i], label=f"{algorithms1[i]} median attainment", surfs=surfs)

    # ax.set_xlim([2.6, 2.7])
    # ax.set_ylim([1.4, 1.7])
    ax.legend()  
    ax.grid() 
    plt.show()  
