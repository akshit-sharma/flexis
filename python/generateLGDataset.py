import argparse
import random
from pathlib import Path

def commandParser():
    parser = argparse.ArgumentParser(description='Generate LG dataset')
    parser.add_argument('--input-dir', type=Path, default='./data/extracted/', help='Input directory')
    parser.add_argument('--verbose', action='store_true', help='Verbose mode')
    parser.add_argument('--vertex-labels', type=int, default=5, help='Number of vertex labels')
    parser.add_argument('--edge-labels', type=int, default=5, help='Number of edge labels')
    parser.add_argument('--seed', type=int, default=0, help='Random seed')
    return parser.parse_args()

def readInputFile(input_file):
    vertices = set()
    edges = set()
    with open(input_file, 'r') as read:
        for line in read:
            if line[0] == '#':
                continue
            src, dst = line.split()
            src = int(src)
            dst = int(dst)
            vertices.add(src)
            vertices.add(dst)
            edges.add((src, dst))
    return vertices, edges


def generateLGFile(input_file, args):
    assert input_file.exists(), f'Input file {args.input_file} does not exist'
    # input_file will have ./data/temp/som-filename.txt
    # output_file will have ./data/lg/filename.lg
    vertices, edges = readInputFile(input_file)
    
    def vertexLabel():
        return random.randint(0, args.vertex_labels - 1)
    def edgeLabel():
        return random.randint(0, args.edge_labels - 1)

    random.seed(args.seed)
    vertices = [(vertex, vertexLabel()) for vertex in vertices]
    random.seed(args.seed)
    edges = [(src, dst, edgeLabel()) for src, dst in edges]

    # changes the mapping of vertices from 0 to n-1 if not already there
    # and remap edges accordingly
    vertices = sorted(vertices)
    mapping = {vertex: i for i, (vertex, _) in enumerate(vertices)}

    vertices = [(mapping[vertex], label) for vertex, label in vertices]
    edges = [(mapping[src], mapping[dst], label) for src, dst, label in edges]

    assert len(vertices) == max(vertices)[0] + 1, f"len(vertices) = {len(vertices)}, max(vertices) = {max(vertices)}"

    return vertices, edges

def outputLGFile(output_file, vertices, edges):
    with open(output_file, 'w') as write:
        write.write(f'# t 1\n')
        for vertex, label in vertices:
            write.write(f'v {vertex} {label}\n')
        for src, dst, label in edges:
            write.write(f'e {src} {dst} {label}\n')

def main(args):
    for input_file in args.input_dir.glob('**/*.txt'):
        filename = input_file.stem
        output_file = input_file.parent.parent / 'lg' / (filename + '.lg')
        if output_file.exists():
            continue
        if args.verbose:
            print(f"{input_file} to {output_file}")
        vertices, edges = generateLGFile(input_file, args)
        outputLGFile(output_file, vertices, edges)

if __name__ == '__main__':
    args = commandParser()
    main(args)
