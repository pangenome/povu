import sys

def parse_tree_file(input_path):
    node_info = {}  # maps node_id to label
    edges = []
    with open(input_path, 'r') as f:
        lines = f.readlines()
    for line in lines[1:]:  # skip header
        parts = line.strip().split('\t')
        if len(parts) < 4:
            continue
        label = parts[0].strip()
        node_id = parts[1].strip()
        children_field = parts[3].strip()
        if node_id not in node_info:
            node_info[node_id] = label

        if children_field != '.':
            children_list = [c.strip() for c in children_field.split(',') if c.strip()]
            for child_id in children_list:
                edges.append((node_id, child_id))

    return node_info, edges

def write_gml(output_path, node_info, edges):
    with open(output_path, 'w') as g:
        g.write("graph [\n  directed 1\n")
        for node_id in sorted(node_info.keys(), key=lambda x: int(x) if x.isdigit() else x):
            label = node_id
            annotation = node_info[node_id]  # escape quotes for GML
            g.write(f"  node [\n    id {label}\n    label \"{label}\"\n    annotation \"{annotation}\"\n  ]\n")
        for src, dst in edges:
            g.write(f"  edge [ source {src} target {dst} ]\n")
        g.write("]\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pvst_to_gml.py /path/to/input.pvst /path/to/output.gml")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    node_info, edges = parse_tree_file(input_file)
    write_gml(output_file, node_info, edges)
    print(f"GML file written to {output_file}")