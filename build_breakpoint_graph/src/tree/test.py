from ete3 import Tree


tree_file = '/Users/alexey/PycharmProjects/true-dist-infer/src/tree/hg38.100way.tree'
# infercars_file = '/Users/alexey/PycharmProjects/true-dist-infer/real_data/ucsc_11_diff_res/540k/Conserved.Segments'

# ers = count_tree_errors(tree_file, infercars_file,
#                         error_func=lambda y, y_tree: (y - y_tree)
#                         )
# print(ers)

t = Tree(tree_file)
t.prune(['hg38', 'mm10', 'rn6', 'panTro6', 'galGal6'])

# t.set_outgroup('galGal6')
t.unroot()
print(t)