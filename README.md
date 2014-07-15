PolytomySolver
==============

Given a non-binary gene tree (rooted or unrooted), a species tree and a gene distance matrix, outputs a binarization of the gene tree
that minimizes the duplication+loss score.  Since the number of solutions can easily get exponential,
the algorithm uses a NJ distance criterion to join the 'nearest' subtrees first.
If the tree is treated as unrooted (using the -r argument), the program tries every possible root (this takes a while)
in order to find the one that yields the lowest DL-score after correction, or to output every rooted correction.

Example calls :
Correct the 15th binary gene tree in star_trees.trees using the specified species tree (-s) and distances file  (-d), verbose mode.
```
polytomysolver -s "speciestree.newick" -g "star_trees.trees" -gl 15 -d "distances/famille_15.dist" -o "famille_15.corrected" -v
```

Try every rooting of the 1st tree in star_trees.trees, and output every possible correction along with the rooted-binary version
of each tree and their DL-score.  The -n flags makes the program replace negative distances by 99999.
```
polytomysolver -gl 1 -s "speciestree.newick" -g "star_trees.trees" -d "famille_1.dist" -r outputallroots -n
```

Try every rooting in the same manner, but only output the gene tree correction yielding the smallest DL-score
```
polytomysolver -gl 1 -s "speciestree.newick" -g "star_trees.trees" -d "famille_1.dist" -r findbestroot -n
```

Solve the gene tree found at line 15 in star_trees.trees, using the species tree in speciestree.newick and distances found in distances/famille_15.dist.
Outputs the final tree in famille_15.corrected. -v is for verbose mode, to print out info in the std output.
If -v is replaced by -v2, then ultra-verbose mode rushes in and outputs all the gory details of the NJ calculation.

Command line arguments
----------------------

-s [filename] : name of the file containing the species tree newick
-sn [species newick]: newick string of the species tree
[One of -s or -sn must be specified.  If both are, -s takes precedence]

-g [filename] : name of the file containing the gene tree
-gl [line #] : index of the line in the file specified by -g that contains the gene tree newick.  If unspecified, the while file is considered to contain the newick.
-gn [gene newick] : newick string of the gene tree
[One of -g or -gn must be specified.  If both are, -g takes precedence]

-d [filename] : name of the file containing the distances between each pair of genes (leaves in the gene tree)

-o [filename] : file in which to output the corrected tree.  If omitted, the gene tree newick is printed on the std out.

-v : verbose mode.  Optional.

-v2 : ultra-verbose mode.  Outputs all the details about the NJ distance calculations, for each cell of the dynamic programming table.

-r : reroot mode.  Must be one of {none,outputallroots,findbestroot}.  Defaults to none.
     The 3 modes act as such :
     none : reads the input gene tree newick as a rooted tree, corrects its polytomies and outputs the result.
     outputallroots : reads the input gene tree as unrooted, tries every possible root on the tree, and outputs all
                      polytomy corrected version of each rooted tree, along with two comment lines before the corrected tree :
                      the non-binary rooted tree before correction, and the resulting dup-loss score after correction.
     findbestroot : reads the input gene tree as unrooted, tries every possible root on the tree, binarizes them all and outputs the
                    tree with the smallest DL-score (the first one encountered if not unique).

-n : non-negative flag.  If present, will treat negative distances as large distances (99999).

-e : Edge flag.  If -r is outputallroots or findbestroot, the default behavior is to test every node as a root, as opposed to
     also testing every edge as a root (meaning adding a vertex on an edge and setting it as the root).
     However, this would be useless, as argued by Eric Tannier.
-cno : No cache flag.  After resolving a given clade C (a clade is a set of leaf labels), the program caches the resolution
       for C.  If C is met again, the cached resolution is used instead of recomputing the resolution.
       This should not be a problem, but in event that it is, -cno disables the cache feature.

File formats
------------

Each leaf of the species tree must have a distinct label in the newick string.
Not only each leaf of the gene tree must also have a distinct label, it must also specify a mapping with the species tree.
Both fields must be part of the label, separated by ";;" (without quotes).
Say the gene labeled G0000001 is mapped to the species pelodiscus_sinensis.  Then the leaf corresponding to this gene
must have this label :
G0000001;;pelodiscus_sinensis

As for the distances file, it must specify a distance between each pair of leaves, in the form of a matrix with fields separated by a space.
The first line of the dist file is the # of genes it contains.
Say the gene tree has n leaves.  Then the rest of the dist file specifies a n * (n + 1) matrix, denoted D.
The first column must contain the gene names from the gene tree (e.g. G0000001), and
D[i,j + 1] is the distance between the i-th gene and the j-th gene.

```
polytomysolver -gl 1 -s "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/Compara.73.species_tree" -g "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/famille_2.start_tree" -d "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/famille_1.dist" -r findbestroot -n
```

Change log
----------
1.2.3
- Corrected a dumb bug that would wrongly parse a Newick string, sometimes including the last ';' to the label of the root of the tree.

1.2.2
- Added cache feature to speed up multiple root testing (see -cno argument).
- The program only tests rooting at nodes of the tree, and does not test rooting on edges (see -e argument).
