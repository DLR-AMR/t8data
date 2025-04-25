import sys

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
    print(f"Reading file: {file_path}")
    with open(file_path, 'r') as file:
        for line in file:
            print(f"Processing line: {line.strip()}")
            if "-------------  Running:" in line and "with" in line and "procs" in line:
                print(f"Found line: {line.strip()}")
                parts = line.split("-------------  Running:")[1].strip().split("with")
                args = parts[0].strip()
                procs = parts[1].strip().split()[0]
                data.append({"args": args, "procs": int(procs)})
                print(f"Extracted args: {args}, procs: {procs}")

    return data


def main():  
    if len(sys.argv) != 3:
        print("Usage: python3 create_graph.py <base_file> <compare_file>")
        sys.exit(1)

    base = sys.argv[1]
    compare = sys.argv[2]

    print(f"Base file: {base}")
    print(f"Compare file: {compare}")

    extract_data(base, [0, 1])

if __name__ == "__main__":
    main()