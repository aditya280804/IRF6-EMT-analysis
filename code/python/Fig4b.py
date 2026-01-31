import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import os
import sys
import matplotlib.pyplot as plt
from scipy import stats

# JCO color palette
jco_colors = [
    '#0073C2FF',  # blue
    '#EFC000FF',  # yellow
    '#868686FF',  # gray
    '#CD534CFF',  # red
    '#7AA6DC',  # light blue
    '#003C67',  # dark blue
    '#8F7700',  # dark yellow
    '#3B3B3B',  # dark gray
    '#A73030',  # dark red
]
plt.rcParams.update({
    "font.size": 10,              # base font
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 10
})


def process_file(file_path, cfg_file_path):
    # Load the data
    df = pd.read_csv(file_path, delimiter='\t', header=None)
    exclude_cols = [0, 1, 2]
    cols_to_normalize = df.drop(columns=exclude_cols)
    cols_to_exclude = df[exclude_cols]

    # Normalize the data
    scaler = StandardScaler()
    normalized_data = scaler.fit_transform(cols_to_normalize)
    normalized_df = pd.DataFrame(normalized_data, columns=cols_to_normalize.columns)

    # Recombine the normalized columns with the excluded columns
    final_df = pd.concat([cols_to_exclude, normalized_df], axis=1)

    # Read gene names from the config file
    gene_names = []
    with open(cfg_file_path, 'r') as file:
        lines = file.readlines()
        gene_section = False
        for line in lines:
            if line.startswith("NumberOfGenes"):
                gene_section = True
            elif gene_section:
                if line.strip():  # Avoid empty lines
                    parts = line.split()
                    if len(parts) == 2 and parts[0].isdigit():
                        gene_names.append(parts[1])
                    else:
                        break  # Stop if the format is not as expected

    original_columns = ['Model', 'Parameter Set', 'Stable State']
    final_df.columns = original_columns + gene_names

    # Create an empty DataFrame for the states
    data = np.zeros((10000, 3))
    df_2 = pd.DataFrame(data, columns=['E', 'M', 'H'])

    # Populate df_2 based on conditions
    for i, row in final_df.iterrows():
        model_number = int(row['Model']) - 1  # Adjusting for 0-based index in df
        if row['ZEB'] < 0 and row['mir200'] > 0:
            df_2.at[model_number, 'E'] = 1
        elif row['ZEB'] > 0 and row['mir200'] > 0:
            df_2.at[model_number, 'H'] = 1
        elif row['ZEB'] > 0 and row['mir200'] < 0:
            df_2.at[model_number, 'M'] = 1
        if i in (1,2):
            print(i,row)

    # Count the occurrences of each state
    counts = {
        'E': 0,
        'M': 0,
        'H': 0,
        'EM': 0,
        'EH': 0,
        'HM': 0,
        'EHM': 0
    }
    states = {
        'E': 0,
        'M': 0,
        'H': 0
    }

    for i, row in df_2.iterrows():
        if i in (1,2):
            print(i,row)
        e = row['E']
        m = row['M']
        h = row['H']

        if e == 1 and m == 0 and h == 0:
            counts['E'] += 1
            states['E'] +=1
        elif e == 0 and m == 1 and h == 0:
            counts['M'] += 1
            states['M'] +=1
        elif e == 0 and m == 0 and h == 1:
            counts['H'] += 1
            states['H'] +=1
        elif e == 1 and m == 1 and h == 0:
            counts['EM'] += 1
            states['E'] +=1
            states['M'] +=1
        elif e == 1 and m == 0 and h == 1:
            counts['EH'] += 1
            states['E'] +=1
            states['H'] +=1
        elif e == 0 and m == 1 and h == 1:
            counts['HM'] += 1
            states['M'] +=1
            states['H'] +=1
        elif e == 1 and m == 1 and h == 1:
            counts['EHM'] += 1
            states['E'] +=1
            states['M'] +=1
            states['H'] +=1
    
    total_counts = sum(states.values())
    counts_ratio = {key: value / total_counts for key, value in counts.items()}
    states_ratio = {key: value / total_counts for key, value in states.items()}
    
    return counts_ratio, states_ratio


# Function to calculate statistical significance
def calculate_significance(ts_data_1, ts_data_2):
    # Perform t-tests or other statistical tests
    ts_data_1 = ts_data_1.apply(pd.to_numeric, errors='coerce')
    ts_data_2 = ts_data_2.apply(pd.to_numeric, errors='coerce')
    p_values = {}
    for col in ts_data_1.columns:
        t_stat, p_value = stats.ttest_ind(ts_data_1[col], ts_data_2[col], equal_var=False)
        p_values[col] = p_value

    # Determine significance based on p-value (e.g., p < 0.05)
    significance = {}
    for col, p_value in p_values.items():
        significance[col] = p_value < 0.05

    return significance

def plot_mean_values_with_significance(counts_df, groups):
    # Calculate number of groups
    num_groups = len(counts_df) // 3

    # Initialize the figure
    plt.figure(figsize=(5.3, 3))

    # Bar width
    bar_width = 0.15
    bar_containers = []
    # Iterate through groups
    for i in range(num_groups):
        # Separate data for current group
        ts_data = counts_df.iloc[i*3:(i+1)*3, 1:]

        # Calculate mean and standard deviation
        ts_mean_values = ts_data.mean()
        ts_std_values = ts_data.std()

        # Positions of the bars
        bar_positions = [pos + i * bar_width for pos in range(len(ts_mean_values))]

        # Plotting bars for current group
        bars = plt.bar(bar_positions, ts_mean_values, bar_width, yerr=ts_std_values, label=f'{groups[i*3]}', capsize=5, alpha=0.8, color=jco_colors[i % len(jco_colors)] )
        # Store the bar containers
        bar_containers.append(bars)
        # Calculate statistical significance between the first and current group
        if i > 0 and len(groups) > 1:
            significance = calculate_significance(ts_data, counts_df.iloc[0:3, 1:])

            # Add significance markers (asterisks)
            for j, bar in enumerate(bars):
                if significance[ts_data.columns[j]]:
                    plt.text(bar.get_x() + bar.get_width()/2., 1.05 * bar.get_height(), '*', ha='center', va='bottom', fontsize=12)

    # Adding labels, title, and ticks
    plt.ylim(0,0.6)
    plt.xlabel('Category')
    plt.ylabel('Probability')
    plt.xticks([pos + (num_groups - 1) * bar_width / 2 for pos in range(len(ts_mean_values))], ts_data.columns)
    #plt.legend()
    #custom_legend_labels = ['WT ', 'IRF6↑', 'ELF3↑', 'KLF4↑']
    custom_legend_labels = ['WT ', 'IRF6↓', 'ELF3↓', 'KLF4↓']
    plt.legend(
        bar_containers,
        custom_legend_labels,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
        handlelength=1.2,
        handletextpad=0.4,
        borderaxespad=0.2
    )

    plt.tight_layout()
    plt.savefig(
        "./figures/Fig4/FigB_Down.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


if __name__ == "__main__":
    directory_path = "./data/RACIPE_practice"
    groups = sys.argv[1:]
    
    # Initialize empty DataFrame for counts
    counts_df = pd.DataFrame(columns=['File', 'E', 'M', 'H', 'EM', 'EH', 'HM', 'EHM'])
    states_df = pd.DataFrame(columns=['File', 'E', 'M', 'H'])
    # Process each group argument
    for topo_name in groups:
        print("Processing:", topo_name)
        if "over_IRF6" in topo_name:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_OE_7_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}_OE_7.cfg")
        elif "under_IRF6" in topo_name:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_DE_7_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}_DE_7.cfg")
        elif "over_KLF4" in topo_name:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_OE_3_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}_OE_3.cfg")
        elif "under_KLF4" in topo_name:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_DE_3_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}_DE_3.cfg")
        elif "over_ELF3" in topo_name:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_OE_5_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}_OE_5.cfg")
        elif "under_ELF3" in topo_name:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_DE_5_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}_DE_5.cfg")
        else:
            file_path = os.path.join(directory_path, topo_name, f"{topo_name}_solution.dat")
            cfg_file_path = os.path.join(directory_path, f"{topo_name}.cfg")

        
        counts, states = process_file(file_path, cfg_file_path)
        
        # Create DataFrame for current counts
        counts_df_current = pd.DataFrame([{'File': topo_name, **counts}])
        
        # Append current counts to main DataFrame
        counts_df = pd.concat([counts_df, counts_df_current], ignore_index=True)
        # Create DataFrame for current counts
        states_df_current = pd.DataFrame([{'File': topo_name, **states}])
        
        # Append current counts to main DataFrame
        states_df = pd.concat([states_df, states_df_current], ignore_index=True)
    
    # Plot mean values with statistical significance indicators
    # plot_mean_values_with_significance(counts_df, groups)
    plot_mean_values_with_significance(states_df, groups)
