"""
This is the test file of the project.
It loads the example files (from data folder) and runs the algorithm on them,
storing the results in the test_results folder.
"""
import networkx as nx
import matplotlib.pyplot as plt
from alleleconsolidator import *


def main():
    """
    This is the main function of the test file.
    It loads the example files (from data folder) and runs the algorithm on them,
    storing the results in the test_results folder.
    """
    file_mapping_1 = "blast_1.txt"
    file_mapping_2 = "blast_2.txt"
    file_expression_in = "expression.csv"
    file_expression_out = "expression_2.csv"
    folder_in = "data/"
    folder_out = "test_results/"
    prefix_restricted = "restricted"
    g = allele_transform(folder_in + file_mapping_1, folder_in + file_mapping_2,
                         folder_in + file_expression_in, folder_out,
                         file_expression_out, prefix_restricted, verbose=True)
    print("Results from the allele transform algorithm (first 20 lines):")
    with open(f"{folder_out}{prefix_restricted}.txt", "r") as f:
        for line in f.readlines()[:20]:
            print(line.strip())

    # Plotting the resulting graph
    draw_multigraph(g, ["SS", "qseqid", "sseqid", "ST"], x_offset=10,
                    savefig=True, filename="test_results/graph.pdf")


if __name__ == "__main__":
    main()
