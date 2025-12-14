#!/usr/bin/env python3
import sys

def add_column(input_file, output_file, value="3"):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            line = line.rstrip('\n')
            if line:
                f_out.write(f"{line}\t{value}\n")
            else:
                f_out.write("\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python add_column.py <input.tsv> <output.tsv> [value]")
        print("Example: python add_column.py input.tsv output.tsv 3")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    value = sys.argv[3] if len(sys.argv) > 3 else "3"
    
    add_column(input_file, output_file, value)
    print(f"Done: {output_file}")