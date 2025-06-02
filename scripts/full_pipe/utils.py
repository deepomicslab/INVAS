import os

def check_graphg_file(graphg_file):
    if not os.path.exists(graphg_file):
        raise FileNotFoundError(f"GraphG file {graphg_file} does not exist")
    segs_lines=[]
    with open(graphg_file, 'r') as f:
        for line in f:
            if line.startswith("SEG"):
                segs_lines.append(line.strip())
    if len(segs_lines[0].split()) == 7:
        return False
    elif len(segs_lines[-1].split()) == 7:
        return False
    else:
        return True
    