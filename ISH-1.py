from Bio import SeqIO
from fuzzywuzzy import fuzz
from Bio.Seq import reverse_complement
from collections import Counter
import matplotlib.pyplot as plt
import gzip
from tqdm import tqdm
import random
import tkinter as tk
from tkinter import filedialog, messagebox, Text, Label, Button, Entry, Toplevel
import os

def fuzzy_match(target, read_sequence, threshold, use_reverse_complement=True):
    target_to_match = target if not use_reverse_complement else reverse_complement(target)
    partial_ratio = fuzz.partial_ratio(target_to_match, read_sequence)
    return partial_ratio >= threshold, partial_ratio

def display_result_and_progress(output_file, file_label):
    with open(output_file, 'r') as txt_file:
        result_window = Toplevel(root)
        result_window.title("Analysis Result")
        result_text = Text(result_window, height=20, width=60)
        result_text.pack(pady=10)
        result_text.insert(tk.END, "Analysis Result:\n\n")
        result_text.insert(tk.END, txt_file.read())

    file_label.config(text="Selected FASTQ File: None")

def plot_results(target_counts):
    names = [f"{name} ({'Reverse Complement' if is_reverse_complement else 'Original'})" for (name, is_reverse_complement) in target_counts.keys()]
    counts = list(target_counts.values())

    plt.bar(names, counts, color='skyblue')
    plt.xlabel('Sequence Names')
    plt.ylabel('Read Counts')
    plt.title('Target Read Quantification')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

def detect_and_quantify_target_reads(input_file, targets, threshold, sampling_percentage, output_file_path, file_label):
    target_counts = Counter({(name, False): 0 for name, _ in targets})
    total_matched_reads = 0
    total_reads = 0

    # Update the label to display the selected FASTQ file
    file_label.config(text=f"Selected FASTQ File: {input_file}")

    with gzip.open(input_file, "rt") as handle:
        all_records = list(tqdm(SeqIO.parse(handle, "fastq"), desc="Reading Records", unit="reads", dynamic_ncols=True, position=0, leave=True))
        num_reads_to_sample = int(len(all_records) * (sampling_percentage / 100))
        sampled_records = random.sample(list(tqdm(all_records, desc="Sampling", unit="reads", dynamic_ncols=True, position=0, leave=True)), min(num_reads_to_sample, len(all_records)))

        for record in tqdm(sampled_records, desc="Processing", unit="reads", dynamic_ncols=True, position=0, leave=True):
            read_sequence = str(record.seq)
            total_reads += 1
            for name, target in targets:
                # Match with the original target
                match_original, match_percentage_original = fuzzy_match(target, read_sequence, threshold, use_reverse_complement=False)
                if match_original:
                    target_counts[(name, False)] += 1
                    total_matched_reads += 1

                # Match with the reverse complement
                match_reverse_complement, match_percentage_reverse_complement = fuzzy_match(target, read_sequence, threshold, use_reverse_complement=True)
                if match_reverse_complement:
                    target_counts[(name, True)] += 1
                    total_matched_reads += 1

    # Create text file with results
    with open(output_file_path, 'w') as txt_file:
        txt_file.write("Metric\tValue\n")
        txt_file.write(f"Total Selected Reads\t{len(sampled_records)}\n")
        txt_file.write(f"Total Matched Reads\t{total_matched_reads}\n")
        txt_file.write(f"Average Match Percentage\t{(total_matched_reads / total_reads) * 100:.2f}%\n")
        txt_file.write("\n")
        txt_file.write("Sequence Name\tRead Count\tTotal Match Count\tAverage Match Percentage\n")
        for (name, is_reverse_complement), count in target_counts.items():
            total_match_count = target_counts[(name, is_reverse_complement)]
            average_match_percentage = (total_match_count / total_matched_reads) * 100 if total_matched_reads > 0 else 0
            match_type = "Reverse Complement" if is_reverse_complement else "Original"
            txt_file.write(f"{name} ({match_type})\t{count}\t{total_match_count}\t{average_match_percentage:.2f}%\n")

    # Display the result and progress in separate windows
    display_result_and_progress(output_file_path, file_label)

    # Plot the results
    plot_results(target_counts)

def browse_file(file_var, file_label):
    file_path = filedialog.askopenfilename(filetypes=[("FASTQ Files", "*.fastq.gz")])
    file_var.set(file_path)

    # Update the label to display the selected FASTQ file
    file_label.config(text=f"Selected FASTQ File: {file_path}")

def start_analysis(file_var, targets_text_widget, threshold_entry, percentage_entry, output_file_entry, file_label):
    input_file = file_var.get()
    if not input_file:
        messagebox.showwarning("Warning", "Please select a FASTQ file.")
        return

    targets_text = targets_text_widget.get("1.0", tk.END)
    if not targets_text.strip():
        messagebox.showwarning("Warning", "Please enter target sequences in FASTA format.")
        return

    targets = []
    fasta_sequences = targets_text.split('>')
    for fasta_sequence in fasta_sequences:
        lines = fasta_sequence.strip().split('\n')
        if not lines:
            continue
        name = lines[0].strip()
        sequence = ''.join(line.strip() for line in lines[1:])
        if name and sequence:
            targets.append((name, sequence))

    threshold = threshold_entry.get()
    sampling_percentage = percentage_entry.get()

    # Extract the directory and file name from the user input
    output_file_path = output_file_entry.get()
    output_directory, output_file_name = os.path.split(output_file_path)

    # Use the extracted directory and file name in the analysis function
    try:
        detect_and_quantify_target_reads(input_file, targets, int(threshold), int(sampling_percentage), output_file_path, file_label)
        messagebox.showinfo("Analysis Complete", "Analysis completed successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred during analysis:\n{str(e)}")

# GUI setup
root = tk.Tk()
root.title("Target Read Detection")

file_var = tk.StringVar()

# Browse file button
browse_button = tk.Button(root, text="1. Select a FASTQ file", command=lambda: browse_file(file_var, file_label))
browse_button.pack(pady=5)

# Label to display the selected FASTQ file
file_label = Label(root, text="Selected FASTQ File: None")
file_label.pack(pady=5)

# Target sequences input text window
targets_label = Label(root, text="2. Enter target sequences in FASTA format:")
targets_label.pack(pady=5)

targets_text_widget = Text(root, height=10, width=50)
targets_text_widget.pack(pady=5)

# Fuzzy match threshold input
threshold_label = Label(root, text="3. Enter the threshold for fuzzy match percentage (e.g., 80):")
threshold_label.pack(pady=5)

threshold_entry = Entry(root)
threshold_entry.pack(pady=5)

# Percentage of reads to sample input
percentage_label = Label(root, text="4. Enter the percentage of reads to sample (1-10 recommended, e.g., 5):")
percentage_label.pack(pady=5)

percentage_entry = Entry(root)
percentage_entry.pack(pady=5)

# Output text file input with directory and file name
output_file_label = Label(root, text="5. Enter the output text file path (e.g., /path/to/directory/output.txt):")
output_file_label.pack(pady=5)

output_file_entry = Entry(root)
output_file_entry.pack(pady=5)

# Start Analysis button
start_button = Button(root, text="Start Analysis", command=lambda: start_analysis(file_var, targets_text_widget, threshold_entry, percentage_entry, output_file_entry, file_label))
start_button.pack(pady=10)

root.mainloop()
