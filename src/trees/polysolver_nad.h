#ifndef NADCORRECTOR_H
#define NADCORRECTOR_H

#include "trees/node.h"
#include "trees/genespeciestreeutil.h"

#include "trees/newicklex.h"

#include <map>
#include <unordered_map>

#include <iostream>

#include <utility>

using namespace std;

class SADNADNode;
/**
  Note : this is a naive, slow implementation
  **/
class SADNADGraph
{
public:
    SADNADGraph();
    ~SADNADGraph();
    void AddNode(Node* n);
    void AddEdge(Node* n1, Node* n2, string edgeType);
    //set<Node*> GetADComponentOf(Node* n);
    Node* GetADComponentRepresentantOf(Node* n);

    bool HaveSameADComponent(Node* n1, Node* n2);

    string GetEdgeType(Node* n1, Node* n2, unordered_map<Node*, Node*> &lcaMapping);
    void BuildGraph(vector<Node*> nodes, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);
    void MergeNodes(Node* n1, Node* n2, Node* newNode, unordered_map<Node*, Node*> &lcaMapping);



    const map<Node*, SADNADNode*> GetNodes(){return nodes;}

    int GetNbADComponents();
    int GetNbADComponents(set<pair<Node*, Node*> > &additionalSEdges);

    vector<pair<Node *, Node *> > GetUsefulSpeciationEdges();

    void PrintGraph();

    bool HasSEdge(Node* n1, Node* n2);

private:
    map<Node*, SADNADNode*> nodes;

    //map<Node*, set<Node*> > ADComponents;   //the key is a representant of the component

};


class SADNADNode
{
public:
    set<Node*> S_Neighbors;
    set<Node*> AD_Neighbors;
    set<Node*> NAD_Neighbors;

    Node* ADComponent;

    string GetEdgeTypeWith(Node* n)
    {
        if (S_Neighbors.find(n) != S_Neighbors.end())
            return "S";
        if (AD_Neighbors.find(n) != AD_Neighbors.end())
            return "AD";
        if (NAD_Neighbors.find(n) != NAD_Neighbors.end())
            return "NAD";
        throw "No edge defined !";
    }

    void RemoveNeighbor(Node* n)
    {
        S_Neighbors.erase(n);
        AD_Neighbors.erase(n);
        NAD_Neighbors.erase(n);
    }
};



class PolySolverCorrectionInfo
{
public:
    Node* correction;
    int firstPolySize;
    int secondPolySize;
    vector<string> nadCladeGenes;
};


class PolySolverNAD
{
public:
    PolySolverNAD();


    Node* SolvePolytomies(Node *geneTree, Node *speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping);

    Node* SolvePolytomy(vector<Node*> leaves, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);
    Node* SolvePolytomy(Node* polytomyRoot, Node* speciesTree, unordered_map<Node*, Node*> &lcaMapping);

    //TODO : does this really belong here ?  Asking the question is answering it...
    /**
      Returns a copy of geneTree in which the highest NAD is corrected, or NULL if no correction is applied.
      The returned value must be deleted by caller, unless NULL.
      **/
    PolySolverCorrectionInfo CorrectHighestNAD(Node* geneTree, Node* speciesTree, unordered_map<Node*, Node*> geneLeavesSpeciesMapping);

    void PerformRandomTest(int k, int verbose = 0);

    int FindBestSolution(SADNADGraph &graph, vector<pair<Node*, Node*> > availableSpecs, set<pair<Node*, Node*> > chosenSpecs, int verbose = 0);

    int bestSolSoFar;
    set<pair<Node*, Node*> > bestChoiceSoFar;

     int last_nb_ad_components;  //use with care after calling SolvePolytomy
private:
    pair<Node*, unordered_map<Node*, Node*> > PolytomizeNAD(Node* nadNode, Node* speciesTree, unordered_map<Node*, Node*> lcaMapping);

    void SpeciateCleverly(SADNADGraph &graph, set<Node*> &gLeftGuys, set<Node*> &gRightGuys, Node* curSpecies, unordered_map<Node*, Node*> &lcaMapping);

    vector<vector<Node*> > GetSortedLocalADComponents(SADNADGraph &graph, set<Node*> guys, set<Node*> otherGuys);





};

#endif // NADCORRECTOR_H
