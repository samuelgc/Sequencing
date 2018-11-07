"""

Sequencing

1. build deBruijn graph
2. get contigs
3. maybe do something fancy to deal with errors in the reads

"""

# All the data segments are divisible by 5 so probably use kmer of size 5

def debruijn(dna):
    result = {}
    for seq in dna:
        suffix = seq[:-1]
        result[suffix] = []
    for seq in dna:
        suffix = seq[:-1]
        prefix = seq[1:]
        if prefix not in result[suffix]:
            result[suffix].append(prefix)
    return result


def get_hub_node(graph):
    keys = []
    for key in graph:
        if check_deg(key, graph):
            keys.append(key)
    return keys


def check_deg(node, graph):
    in_deg = 0
    for _, value in graph.items():
        for val in value:
            if val == node:
                in_deg += 1
    if in_deg > 1 or len(graph[node]) > 1 or in_deg != len(graph[node]):
        return True
    return False


def euler_path(node, graph, cycle):
    cycle += [node]
    if node not in graph:
        return cycle
    if len(graph[node]) == 0:
        return cycle
    while len(graph[node]) > 0:
        tmp_node = graph[node][0]
        graph[node].remove(tmp_node)
        sub_cycle = euler_path(tmp_node, graph, [])
        cycle = cycle[:1] + sub_cycle + cycle[1:]
    return cycle


def get_start_node(graph):
    for key in graph:
        in_deg = 0
        for _, value in graph.items():
            for val in value:
                if val == key:
                    in_deg += 1
        if in_deg < len(graph[key]):
            return key


def contigs(hubs, graph):
    paths = []
    for node in hubs:
        for p in graph[node]:
            path = node
            path += p[-1]
            q = p
            while q not in hubs and q in graph.keys():
                q = graph[q][0]
                path += q[-1]
            paths.append(path)
    return paths


def main():
    k = 25
    with open('./fasta/example.data.fasta') as f:
        reads = []
        for i, line in enumerate(f.readlines()):
            if i % 2 == 1:
                reads.append(line.strip())
        dna = []
        for read in reads:
            dna += [read[i:i+k] for i in range(len(read) - k + 1)]
        graph = debruijn(dna)
        hubs = get_hub_node(graph)
        contig = contigs(hubs, graph)
        contig.sort()

        start = get_start_node(graph)
        path = euler_path(start, graph, [])

    with open('./result/result.txt', 'w+') as f:
        f.write(' '.join(contig))
        f.write('\n')
        for i, seq in enumerate(path):
            if i == 0:
                f.write(seq)
            else:
                f.write(seq[-1])


if __name__ == '__main__':
    main()
