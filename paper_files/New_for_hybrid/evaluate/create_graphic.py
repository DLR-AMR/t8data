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
                        print(f"Extracted summary data: {summary_data}")
                        print("\n")
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
    print(f"Extracting data from {file_path} for columns {columns}")
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
        num_elements = 0
        for line in lines:
            if "TEST " in line:
                element_type = line.split("TEST")[1].strip()
                num_elements = 0  # Reset num_elements for each new element type
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
                    if "Done t8_forest_balance with" in run_line and num_elements == 0:
                        # Extract the number of elements from the line
                        num_elements = int(run_line.split("with")[1].split("global elements")[0].strip())
                        data[-1]["num_elements"] = num_elements
                    if "[t8] Summary = [" in run_line:
                        # Extract the summary data
                        summary_data = run_line.split("[t8] Summary = [")[1].split("]")[0].strip().split()
                        # Filter the times based on the specified columns
                        summary_data = [summary_data[col] for col in columns if col < len(summary_data)]
                        print(f"Extracted summary data for type {element_type}: {summary_data}")
                        # Convert the summary data to floats and filter based on columns
                        current_run_summaries.append([float(value) for value in summary_data])
                        print(f"Current run summaries: {current_run_summaries}")
                    if "-------------  Running:" in run_line:
                        break  # Exit the inner loop to reprocess the line in the outer loop
                # Compute the average of the times, given by the columns
                avg_summary = [sum(x) / len(x) for x in zip(*current_run_summaries)]
                # Append the average summary to the data
                data[-1]["summaries"] = avg_summary
    return data

def create_graphics(num_files, names, graph_name_base, data_base):
    plt.figure(figsize=(10, 6))
    plt.xlabel('Number of Processes', fontsize=16)
    plt.ylabel('Average Time (s)', fontsize=16)
    plt.xscale('log', base=2)
    plt.yscale('log', base=10)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.title('Performance Comparison', fontsize=16)
    print(data_base)
    for i in range(int(num_files)):
        entry = data_base[i][0]
        procs = entry["procs"]
        summaries = entry["summaries"]
        total = sum(summaries)
        color = plt.cm.tab10(i % 10)
        # Use a line plot (single-point plotted as a marker; use sequences to connect points if desired)
        if i == 0:
            procs_list = []
            totals_list = []
        procs_list.append(procs)
        totals_list.append(total)
        # draw/update a line that connects the accumulated points
        plt.plot(procs_list, totals_list, marker='o', linestyle='-', color='black', label=graph_name_base if i == 0 else None)
        #bottom = 0
        #for j, value in enumerate(summaries):
        #    color = plt.cm.tab10(j % 10)  # Use a consistent color for each section
        #    plt.bar(procs, value, bottom=bottom, label=f"{names[j]}" if i == 0 else "",  width=procs / 10, color=color)
        #    bottom += value

    
    plt_name = f"graph_{graph_name_base}.png"
    plt.legend(fontsize=14)
    plt.savefig(plt_name)
    print(f"Graph saved as {plt_name}")

def create_graphics_elem(names, num_files, data):
    compute_ghosts = False
    for i in range(len(names)):
        if "Ghost" in names[i]:
            compute_ghosts = True
        if "Ghost_sent" in names[i]:
            continue
        
        plt.figure(figsize=(10, 6))
        plt.xlabel('Number of Processes', fontsize=16)
        plt.ylabel('Average Time (s)', fontsize=16)
        plt.xscale('log', base=2)
        plt.yscale('log', base=10)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title(f'Performance Comparison for {names[i]}', fontsize=16)
        color_map = {
            "TETRAHEDRON": "green",
            "HEXAHEDRON": "blue",
            "PRISM": "red",
            "PYRAMID": "orange"
        }
        name_map = {
            "TETRAHEDRON": "Tetrahedron",
            "HEXAHEDRON": "Hexahedron",
            "PRISM": "Prism",
            "PYRAMID": "Pyramid"
        }
        for ifile in range(int(num_files)):
            element_types = set(entry["element_type"] for entry in data[ifile])
            for element_type in element_types:
                element_data = [entry for entry in data[ifile] if entry["element_type"] == element_type]
                procs = [entry["procs"] for entry in element_data]
                times = [entry["summaries"][i] for entry in element_data]
                if compute_ghosts:
                    ghost_sent = [entry["summaries"][i+1] for entry in element_data]
                    parallel_efficiency = [times[0] * ghost_sent[j] / (times[j] * ghost_sent[0]) for j in range(len(times))]
                    print(f"Parallel Efficiency for {element_type}: {parallel_efficiency}")
                    if element_type == "PYRAMID":
                        print(f"Times for PYRAMID: {times}")
                        print(f"Ghost Sent for PYRAMID: {ghost_sent}")
                color = color_map.get(element_type, "gray")
                marker_map = {
                    "TETRAHEDRON": "^",  # triangle_up
                    "HEXAHEDRON": "s",   # square
                    "PRISM": "8",        # octagon
                    "PYRAMID": "D"       # diamond
                }
                marker = marker_map.get(element_type, "o")
                plt.plot(procs, times, label=name_map.get(element_type, element_type) if ifile == 0 else None,
                         marker=marker, color=color, linestyle='-')

                if element_type == "PYRAMID":
                    if compute_ghosts:
                        ghost_sent = [entry["summaries"][i+1] for entry in element_data]
                        ideal_scaling = [times[0]] * len(procs)
                        for j in range(1, len(procs)):
                            ideal_scaling[j] = ideal_scaling[j-1] * (ghost_sent[j] / ghost_sent[j-1])
                        #print(f"Times for PYRAMID: {times}")
                        #print(f"Ghost Sent for PYRAMID: {ghost_sent}")
                        #print(f"Ideal Scaling for PYRAMID: {ideal_scaling}")
                        plt.plot(procs, ideal_scaling, label="Ideal Strong Scaling" if ifile == 0 else None, color='black', linestyle='dashed')
                    else:
                        ideal_scaling = [times[0] / 2**iproc for iproc in range(len(procs))]
                        plt.plot(procs, ideal_scaling, label="Ideal Strong Scaling" if ifile == 0 else None, color='black', linestyle='dashed')
                    if ifile == 0:
                        #if compute_ghosts:
                        #    elem_data = [[entry for entry in data[j] if entry["element_type"] == "PYRAMID"] for j in range(int(num_files))]
                        #    ghost_sent = [[entry["summaries"][i+1] for entry in elem_data[j]] for j in range(int(num_files))]
                        #    print(f"Ghost Sent for PYRAMID: {ghost_sent}")
                        #    for iproc, proc in enumerate(procs):
                        #        val = times[iproc]
                        #        ideal_weak_scaling = [val ] * int(num_files)
                        #        for j in range(1, int(num_files)):
                        #            ideal_weak_scaling[j] = ideal_weak_scaling[j-1] * (ghost_sent[j][iproc] / ghost_sent[j-1][iproc])
                        #        shifted_procs = [proc * 8**i for i in range(int(num_files))]
                        #        plt.plot(shifted_procs, ideal_weak_scaling, color='black', linestyle='dotted', label="Ideal Weak Scaling" if iproc == 0 and ifile == 0 else None)
                        #else:
                        num_elements = element_data[ifile].get("num_elements", 1)
                        num_elements_next_file = next(
                            (entry.get("num_elements", 1) for entry in data[ifile + 1] if entry.get("element_type") == "PYRAMID"),
                            1
                        )
                        element_scaling = num_elements_next_file / num_elements
                        for iproc, proc in enumerate(procs):
                            val = times[iproc]
                            ideal_weak_scaling = [val * (element_scaling**i) / (8**i) for i in range(int(num_files))]
                            shifted_procs = [proc * 8**i for i in range(int(num_files))]
                            plt.plot(shifted_procs, ideal_weak_scaling, color='black', linestyle='dotted', label="Ideal Weak Scaling" if iproc == 0 and ifile == 0 else None)

        plt.grid()
        plt.legend()
        plt.legend(loc='lower left', fontsize=14)
        handles, labels = plt.gca().get_legend_handles_labels()
        order = []
        label_order = ["Pyramid", "Prism", "Tetrahedron", "Hexahedron", "Ideal Strong Scaling", "Ideal Weak Scaling"]
        for lbl in label_order:
            if lbl in labels:
                order.append(labels.index(lbl))
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower left')
        plt_name = f"graph_{names[i]}.png"
        plt.savefig(plt_name)
        print(f"Graph saved as {plt_name}")
        if "Ghost" in names[i]:
            compute_ghosts = False

def compare_mesh():
    num_files = sys.argv[2]

    base = sys.argv[3:3 + int(num_files)]

    additional_args_index = 3 + int(num_files)
    names = [name.strip() for name in sys.argv[additional_args_index].split(',')]

    indices = list(map(int, sys.argv[additional_args_index + 1].split(',')))

    if len(names) != len(indices):
        print("Error: The number of names and indices must be the same.")
        sys.exit(1)
    
    graph_name_base = sys.argv[additional_args_index + 2]

    data_base = []

    for file in range(int(num_files)):
        data_base.append(extract_data(base[file], indices))

    create_graphics(num_files, names, graph_name_base, data_base)

def compare_elements():
    if (len(sys.argv) < 6):
        print("Usage: python create_graphic.py elements <num_files> <base> <names> <indices> <graph_name_base>")
        sys.exit(1)

    num_files = sys.argv[2]
    print (num_files)
    base = sys.argv[3:3 + int(num_files)]
    print(base)
    additional_args_index = 3 + int(num_files)

    names = [name.strip() for name in sys.argv[additional_args_index].split(',')]
    print (names)
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

    if mode == "mesh":
        compare_mesh()
    elif mode == "elements":
        compare_elements()
    else:
        print("Error: Invalid mode. Use 'mesh' or 'elements'.")
        sys.exit(1)



if __name__ == "__main__":
    main()