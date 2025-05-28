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
                # Ensure the last element_type is associated with the current procs
                element_type = data[-1]["element_type"] if data and "element_type" in data[-1] else "Unknown"
                data.append({"element_type": element_type, "procs": int(procs)})
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

def extract_data_elems(file_path, columns):
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
        element_type = None
        for line in lines:
            if "TEST " in line:
                element_type = line.split("TEST")[1].strip()
            elif "-------------  Running:" in line and "with" in line and "procs" in line:
                parts = line.split("-------------  Running:")[1].strip().split("with")
                args = parts[0].strip()
                procs = parts[1].strip().split()[0]
                data.append({"element_type": element_type, "procs": int(procs)})
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
                        ideal_weak_scaling = [val for _ in range(int(num_files))]
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

def create_graphics_elem(names, num_files, data):
    for i in range(len(names)):
        plt.figure(figsize=(10, 6))
        plt.xlabel('Number of Processes')
        plt.ylabel('Average Time (s)')
        plt.xscale('log', base=2)
        plt.yscale('log', base=2)
        plt.title(f'Performance Comparison for {names[i]}')
        color_map = {
            "PYRAMID": "orange",
            "HEXAHEDRON": "blue",
            "TETRAHEDRON": "green",
            "PRISM": "red"
        }
        name_map = {
            "PYRAMID": "Pyramid",
            "HEXAHEDRON": "Hexahedron",
            "TETRAHEDRON": "Tetrahedron",
            "PRISM": "Prism"
        }
        for ifile in range(int(num_files)):
            element_types = set(entry["element_type"] for entry in data[ifile])
            for element_type in element_types:
                element_data = [entry for entry in data[ifile] if entry["element_type"] == element_type]
                procs = [entry["procs"] for entry in element_data]
                times = [entry["summaries"][i] for entry in element_data]
                color = color_map.get(element_type, "gray")
                plt.plot(procs, times, label=name_map[element_type] if ifile == 0 else None, marker='o', color=color)

                if element_type == "PYRAMID":
                    ideal_scaling = [times[0] / 2**iproc for iproc in range(len(procs))]
                    plt.plot(procs, ideal_scaling, label="Ideal Strong Scaling" if ifile == 0 else None, color='black', linestyle='dashed')
                    if ifile != int(num_files) - 1:
                        for iproc, proc in enumerate(procs):
                            val = times[iproc]
                            ideal_weak_scaling = [val * ((2 * (8**i) - 6**i) / (8**i)) for i in range(int(num_files))]
                            shifted_procs = [proc * 8**i for i in range(int(num_files))]
                            plt.plot(shifted_procs, ideal_weak_scaling, color='black', linestyle='dotted', label="Ideal Weak Scaling" if iproc == 0 and ifile == 0 else None)

        plt.grid()
        plt.legend()
        plt_name = f"graph_{names[i]}.png"
        plt.savefig(plt_name)
        print(f"Graph saved as {plt_name}")

def compare_versions():
    num_files = sys.argv[2]

    base = sys.argv[3:3 + int(num_files)]
    compare = sys.argv[3 + int(num_files):3 + 2 * int(num_files)]

    additional_args_index = 3 + 2 * int(num_files)
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

def compare_elements():
    if (len(sys.argv) < 6):
        print("Usage: python create_graphic.py elements <num_files> <base> <names> <indices> <graph_name_base>")
        sys.exit(1)

    num_files = sys.argv[2]
    base = sys.argv[3:3 + int(num_files)]

    additional_args_index = 3 + int(num_files)

    names = sys.argv[additional_args_index].split(',')
    indices = list(map(int, sys.argv[additional_args_index + 1].split(',')))

    graph_name_base = sys.argv[additional_args_index + 2]

    data_base = []
    for file in range(int(num_files)):
        data_base.append(extract_data_elems(base[file], indices))

    print(f"Data extracted from {base}: {data_base}")

    create_graphics_elem(names, num_files, data_base)


def main():
    if len(sys.argv) < 2:
        print("Usage: python create_graphic.py <mode> [args]")
        sys.exit(1)
    
    mode = sys.argv[1]

    if mode == "compare":
        compare_versions()
    elif mode == "elements":
        compare_elements()
    else:
        print("Error: Invalid mode. Use 'compare' or 'elements'.")
        sys.exit(1)



if __name__ == "__main__":
    main()