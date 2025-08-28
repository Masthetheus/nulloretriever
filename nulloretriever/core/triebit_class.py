"""Defines the TrieBit class for nullomer retrieving pipeline"""

from bitarray import bitarray

class TrieBitNode:
    def __init__(self, m):
        self.children = [None] * 4  # 0:A, 1:T, 2:C, 3:G
        self.v2_set = bitarray(m)
        self.v2_set.setall(0)

    def iterate(self, path=None):
        if path is None:
            path = []
        yield path, self
        for i, child in enumerate(self.children):
            if child is not None:
                yield from child.iterate(path + [i])

class TrieBit:
    def __init__(self, m, l):
        self.root = TrieBitNode(m)
        self.m = m
        self.l = l

        def build(node, depth):
            if depth == l:
                return
            for i in range(4):
                if node.children[i] is None:
                    node.children[i] = TrieBitNode(m)
                build(node.children[i], depth + 1)
        build(self.root, 0)

    def insert(self, v1, v2):
        node = self.root
        for value in v1:
            if node.children[value] is None:
                node.children[value] = TrieBitNode(self.m)
            node = node.children[value]
        node.v2_set[v2] = 1

    def iterate(self):
        yield from self.root.iterate([])

    def count_nullomers(self, target_length=None):
        l = target_length or self.l
        def dfs(node, depth):
            if depth == l:
                return sum(1 for genome_id in range(len(node.v2_set)) 
                        if not node.v2_set[genome_id])
            
            return sum(dfs(child, depth + 1) 
                    for child in node.children 
                    if child is not None)
        return dfs(self.root, 0)

