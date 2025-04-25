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
        for line in file:
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
                for run_line in file:
                    if "[t8] #################### Run" in run_line:
                        continue  # Skip the run header lines
                    if "[t8] Summary =" in run_line:
                        # Extract the summary data
                        summary_data = run_line.split("[t8] Summary = [")[1].split("]")[0].strip().split()
                        # Filter the times based on the specified columns
                        summary_data = [summary_data[col] for col in columns if col < len(summary_data)]
                        # Convert the summary data to floats and filter based on columns
                        current_run_summaries.append([float(value) for value in summary_data])
                    if "-------------  Running: " in run_line:
                        break  # Stop when the next run starts
                # Compute the average of the times, given by the columns
                if current_run_summaries:
                    avg_summary = [sum(x) / len(x) for x in zip(*current_run_summaries)]
                    # Append the average summary to the data
                    data[-1]["summaries"] = avg_summary
                else:
                    print("No summaries found for this run.")
    return data

def create_graph(base_data, compare_data, names):
# Create a graph comparing the base and compare data
    plt.figure(figsize=(10, 6))
    for i in range(len(names)):
        procs = [entry["procs"] for entry in base_data]
        base_times = [entry["summaries"][i] for entry in base_data]
        compare_times = [entry["summaries"][i] for entry in compare_data]
        plt.plot(procs, base_times, label=f"{names[i]} Base")
        plt.plot(procs, compare_times, label=f"{names[i]} Compare")
        plt.xlabel('Number of Processes')
        plt.ylabel('Average Time (s)')
        plt.title('Performance Comparison')
        plt.legend()
        plt.grid()
        plt_name = f"graph_{names[i]}.png"
        plt.savefig(plt_name)
        plt.figure(figsize=(10, 6))  # Start a new plot for the next comparison
        print(f"Graph saved as {plt_name}")


def main():  
    if len(sys.argv) != 5:
        print("Usage: python3 create_graph.py <base_file> <compare_file> <names> <indices>")
        print("Example: python3 create_graph.py base.txt compare.txt 'New,Adapt' '0,1'")
        sys.exit(1)

    base = sys.argv[1]
    compare = sys.argv[2]
    names = sys.argv[3].split(',')
    indices = list(map(int, sys.argv[4].split(',')))

    if len(names) != len(indices):
        print("Error: The number of names and indices must be the same.")
        sys.exit(1)

    print(f"Base file: {base}")
    print(f"Compare file: {compare}")

    data_base = extract_data(base, indices)
    data_compare = extract_data(compare, indices)

    create_graph(data_base, data_compare, names)

if __name__ == "__main__":
    main()