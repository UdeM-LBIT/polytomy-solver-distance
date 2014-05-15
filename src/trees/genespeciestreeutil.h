#ifndef GENESPECIESTREEUTIL_H
#define GENESPECIESTREEUTIL_H

#include "trees/node.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <iostream>

using namespace std;

class GeneSpeciesTreeUtil
{
private:
    GeneSpeciesTreeUtil();

public:
    static GeneSpeciesTreeUtil* Instance()
    {
        static GeneSpeciesTreeUtil instance;

        return &instance;
    }

    unordered_map<Node*, Node*> GetLCAMapping(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> &geneLeavesSpeciesMapping);

    unordered_map<Node*, Node*> GetLCAMapping(Node *geneTree, Node *speciesTree, string geneLabelSeparator, int speciesIndex);

    unordered_set<Node*> GetGeneTreeSpecies(Node *geneTree, unordered_map<Node*, Node*> &lcaMapping);

    vector<Node*> GetGenesSpecies(vector<Node*> genes, unordered_map<Node*, Node*> &lcaMapping);

    unordered_map<Node*, Node*> GetGeneSpeciesMappingByLabel(Node* geneTree, Node* speciesTree, string separator = "_", int speciesIndex = 0);

    Node* CopyTreeWithNodeMapping(Node* tree, unordered_map<Node*, Node*> &yourMapping, unordered_map<Node*,Node*> &mappingToFill);

    bool HaveCommonSpecies(Node* tree1, Node* tree2, unordered_map<Node*, Node*> &mapping);

    Node* GetSingleNodeLCAMapping(Node* n, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);

    void PrintMapping(Node* g, unordered_map<Node*, Node*> &mapping);

    vector<Node*> GetNADNodes(Node* g, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);


    void RelabelGenes(Node* geneTree, string search = ";;", string replace = "__");

    void RelabelGenesByIndex(Node* geneTree, string separator, int indexToKeep);

    int GetDLScore(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> lcaMapping);
};

#endif // GENESPECIESTREEUTIL_H
