import sys

def read_assigned_count_from_summary(filename):
    with open(filename, 'r') as file:
        for line in file:
            if "Assigned" in line:
                return int(line.split('\t')[1])
    return 0

def copy_file(source_file, destination_file):
    with open(source_file, 'r') as source:
        with open(destination_file, 'w') as destination:
            for line in source:
                destination.write(line)

def main():
    summary_strand1 = snakemake.input.counts1_summary + '.summary'
    summary_strand2 = snakemake.input.counts2_summary + '.summary'
    output_file = snakemake.output.best_strand_file
    output_file_2 = snakemake.output.best_strand_number_file

    # read assigned counts from summary files
    counts_strand1 = read_assigned_count_from_summary(summary_strand1)
    counts_strand2 = read_assigned_count_from_summary(summary_strand2)

    # determine which strand specificity yields higher count
    chosen_strand = 'strand1' if counts_strand1 > counts_strand2 else 'strand2'

    # output the chosen strand to a file
    with open(output_file_2, 'w') as file:
        file.write(chosen_strand + "\n")

    chosen_summary = summary_strand1 if counts_strand1 > counts_strand2 else summary_strand2
    chosen_file = chosen_summary.rsplit('.summary', 1)[0]
    copy_file(chosen_file, output_file)


if __name__ == "__main__":
    main()
