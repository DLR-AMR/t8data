import sys
import matplotlib.pyplot as plt

def extract_data(file_path, columns):
# The file consists of multiple exectuings with differen number of procs. 
# Each execution starts with a line containing:
# -------------  Running: <command> <args> with <procs> procs -------------
# In args we find the argument -n which states how often the command is executed.
# Each executin of the command starts with a line containing:
# [t8] #################### Run i of n ####################
# At the end of each run there is a line containing:
# [t8] Summary = [ time time time .... time ];
# To extract the data we need to find the lines containing:
# -------------  Running: <command> <args> with <procs> procs -------------
# until the end of the file.
    data = []
    with open(file_path, 'r') as file:
        lines = file.readlines()  # Read all lines into a list
        for line in lines:
            if "-------------  Running:" in line and "with" in line and "procs" in line:
                parts = line.split("-------------  Running:")[1].strip().split("with")
                args = parts[0].strip()
                procs = parts[1].strip().split()[0]
                data.append({"procs": int(procs)})
                # Extract the number of runs from the args
                # Assuming args is in the form of '-n <number>'
                n_value = None
                if any(arg.startswith('-n') and arg[2:].isdigit() for arg in args.split()):
                    n_value = next(arg[2:] for arg in args.split() if arg.startswith('-n') and arg[2:].isdigit())
                    args_list = args.split()  # Split args into a list of arguments
                else:
                    print("No -n argument found in args.")
                # for each run we need to find the line containing:
                # [t8] Summary = [ time time time .... time ];
                current_run_summaries = []
                for run_line in lines[lines.index(line) + 1:]:
                    if "[t8] #################### Run" in run_line:
                        # find the number of the current run
                        run_number = run_line.split("[t8] #################### Run")[1].split("of")[0].strip()
                        continue  # Skip the run header lines
                    if "[t8] Summary = [" in run_line:
                        # Extract the summary data
                        summary_data = run_line.split("[t8] Summary = [")[1].split("]")[0].strip().split()
                        # Filter the times based on the specified columns
                        summary_data = [summary_data[col] for col in columns if col < len(summary_data)]
                        # Convert the summary data to floats and filter based on columns
                        current_run_summaries.append([float(value) for value in summary_data])
                    if "-------------  Running:" in run_line:
                        break  # Exit the inner loop to reprocess the line in the outer loop
                # Compute the average of the times, given by the columns
                if current_run_summaries:
                    avg_summary = [sum(x) / len(x) for x in zip(*current_run_summaries)]
                    # Append the average summary to the data
                    data[-1]["summaries"] = avg_summary
                else:
                    print("No summaries found for this run.")
    return data

def create_graphics(num_files, names, graph_name_base, graph_name_compare, data_base, data_compare, pyra_flag=True):
    for i in range(len(names)):
        plt.figure(figsize=(10, 6))
        plt.xlabel('Number of Processes')
        plt.ylabel('Average Time (s)')
        plt.xscale('log', base=2)
        plt.yscale('log', base=2)
        plt.title('Performance Comparison')
        for ifile in range(int(num_files)):
            procs = [entry["procs"] for entry in data_base[ifile]]
            base_times = [entry["summaries"][i] for entry in data_base[ifile]]
            compare_times = [entry["summaries"][i] for entry in data_compare[ifile]]
            ideal_scaling = [base_times[0] / 2**iproc for iproc in range(len(procs))]
            if ifile == 0:
                plt.plot(procs, base_times, label=f"{names[i]} {graph_name_base} ", color='orange')
                plt.plot(procs, compare_times, label=f"{names[i]} {graph_name_compare} ", color='blue')
                plt.plot(procs, ideal_scaling, label=f"{names[i]} Ideal Scaling", color='black', linestyle='dashed')
                for iproc in range(len(procs)):
                    val = data_base[ifile][iproc]["summaries"][i]  # Use the first value in summaries as the reference
                    ideal_weak_scaling= []
                    if pyra_flag:
                        ideal_weak_scaling = [val * ((2 * (8**i) - 6**i)/(8**i)) for i in range((int(num_files)))]
                    else:
                        ideal_weak_scaling = [val for _ in range(len(procs))]
                    shifted_procs = [procs[iproc] * 8**i for i in range(int(num_files))]
                    if iproc == 0:
                        plt.plot(shifted_procs, ideal_weak_scaling, label=f"{names[i]} Ideal Weak Scaling", color='black', linestyle='dotted')
                    else:
                        plt.plot(shifted_procs, ideal_weak_scaling, color='black', linestyle='dotted')
            else:   
                plt.plot(procs, base_times, color='orange')
                plt.plot(procs, compare_times, color='blue')
                plt.plot(procs, ideal_scaling, color='black', linestyle='dashed')
            
        plt.grid()
        plt.legend()
        plt_name = f"graph_{names[i]}.png"
        plt.savefig(plt_name)
        print(f"Graph saved as {plt_name}")


def main():  
    num_files = sys.argv[1]

    base = sys.argv[2:2 + int(num_files)]
    compare = sys.argv[2 + int(num_files):2 + 2 * int(num_files)]

    additional_args_index = 2 + 2 * int(num_files)
    names = sys.argv[additional_args_index].split(',')
    indices = list(map(int, sys.argv[additional_args_index + 1].split(',')))

    graph_name_base = sys.argv[additional_args_index + 2]
    graph_name_compare = sys.argv[additional_args_index + 3]

    if len(names) != len(indices):
        print("Error: The number of names and indices must be the same.")
        sys.exit(1)
    

    data_base = []
    data_compare = []


    for file in range(int(num_files)):
        data_base.append(extract_data(base[file], indices))
        data_compare.append(extract_data(compare[file], indices))

    create_graphics(num_files, names, graph_name_base, graph_name_compare, data_base, data_compare)


if __name__ == "__main__":
    main()