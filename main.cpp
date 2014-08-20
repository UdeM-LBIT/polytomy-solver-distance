#include "trees/node.h"
#include "trees/polysolver.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"
#include "trees/polysolver_nad.h"
#include "trees/polysolver_distance.h"

#include "trees/paralogycorrector.h"

#include "div/tinydir.h"
//#include "trees/autester.h"

#include <Python.h>

#include <fstream>

#include <string>
#include <unordered_map>
#include <iostream>


using namespace std;







int verbose = 0;
string rerootMode = "none"; //{none,findbestroot,outputallroots}
map<string, unordered_map<string, unordered_map<string, double> > > cachedDistances;
bool hasNonnegativeDistanceFlag = false;
bool testEdgeRoots = false;
bool useCache = true;
bool outputOnlyGenes = false;
int runCounter = 0;

/**
 What's going on here :
 Every tree we want to correct is in a .parc file (here, the hard-coded test.parc),
 whcih is a pseudo-xml file.
 This file includes, for every tree, the list of gene names,
 their corresponding species, the gene tree, the species tree,
 and the orthology constraints.  Most of the main() is dedicated to parsing this file,
 and the magic happens on these two lines :
    ParalogyCorrector pc;
    ctree = pc.CorrectGeneTree(genetree, speciestree, geneSpeciesMapping, orthologs);
 The .parc file is not necessary, but here it is used to fill up
 the genetree, speciestree, geneSpeciesMapping, orthologs parameters.
 */
void DoParalogyCorrection()
{

    //leave this set empty to parse all trees
    set<string> restrictToTrees;

    std::ifstream infile("C:\\Users\\Manuel\\Desktop\\test.parc");

    //Parse the parc file, extract everythng
    //allModeStrings is the list of instances to correct (demarcated by <INSTANCE>...</INSTANCE>)
    //The map<string, string> for a given instance in allModeStrings,
    //contains specific values of the instance, the key being the xml tag and the value, well, the value
    //eg. string newickk = my_map["<GENETREE>"]
    vector<map<string, string> > allModeStrings;
    map<string, string> curModeStrings;

    bool inInstance = false;
    std::string line;
    string mode = "";
    string modestr = "";
    while (std::getline(infile, line))
    {
        cout<<line<<endl;
        if (mode == "")
        {
            if (line == "<INSTANCE>")
                inInstance = true;
            else if (line == "</INSTANCE>" && inInstance)
            {

                if (restrictToTrees.size() == 0 || restrictToTrees.find(curModeStrings["<TREEID>"]) != restrictToTrees.end())
                {
                    allModeStrings.push_back(curModeStrings);
                }

                curModeStrings = map<string, string>();
                mode = "";
                modestr = "";
                inInstance = false;
            }
            else if (line == "<CONSTRAINTS>" || line == "<GENETREE>" || line == "<SPECIESTREE>" || line == "<GENESPECIESMAPPING>" || line == "<TREEID>")
                mode = line;
        }
        else
        {

            string ender = Util::ReplaceAll(mode, "<", "</");
            if (line == ender)
            {
                curModeStrings[mode] = modestr;
                mode = "";
                modestr = "";
            }
            else
            {
                if (line != "")
                {
                    if (modestr != "")
                        modestr += "\n";
                    modestr += line;
                }
            }

        }
    }

    infile.close();



    //For each <INSTANCE> in the .parc file, build up the parameters and correct the tree
    for (int i = 0; i < allModeStrings.size(); i++)
    {
        map<string, string> modeStrings = allModeStrings[i];

        string treeID = modeStrings["<TREEID>"];
        Node* genetree = NewickLex::ParseNewickString(modeStrings["<GENETREE>"], false);
        Node* speciestree = NewickLex::ParseNewickString(modeStrings["<SPECIESTREE>"], true);
        Node* ctree = NULL;


        //Here we parse the <GENESPECIESMAPPING> to build, you guessed it, geneSpeciesMapping
        map<Node*, Node*> geneSpeciesMapping;
        vector<string> geneSpeciesLines = Util::Split(modeStrings["<GENESPECIESMAPPING>"], "\n");

        //TODO : this lazy loop really sucks - we re-search the tree everytime we need a specific node
        //       Consider using treeInfo, which keeps a map of label -> node
        for (int j = 0; j < geneSpeciesLines.size(); j++)
        {
            vector<string> gs = Util::Split(geneSpeciesLines[j], ":");
            Node* g = genetree->GetNodeWithLabel(gs[0], true);
            Node* s = speciestree->GetNodeWithLabel(gs[1], true);

            geneSpeciesMapping[g] = s;
        }

        set<Node*> leaves = genetree->GetLeafSet();
        set<Node*> unmapped;
        for (set<Node*>::iterator it = leaves.begin(); it != leaves.end(); it++)
        {
            if (geneSpeciesMapping.find((*it)) == geneSpeciesMapping.end())
            {
                unmapped.insert(*it);
            }
        }

        cout<<"GENELEAVES="<<genetree->GetLeafSet().size()<<" NBMAPPINGS="<<geneSpeciesMapping.size()<<endl;

        //And here we parse the "<CONSTRAINTS>", the pairs of gene names required to be orthologs
        vector<pair<string, string> > orthologs;
        vector<string> orthoLines = Util::Split(modeStrings["<CONSTRAINTS>"], "\n");

        for (int j = 0; j < orthoLines.size(); j++)
        {
            vector<string> xz = Util::Split(orthoLines[j], ":");
            pair<string, string> p = make_pair(xz[0], xz[1]);
            Node* g1 = genetree->GetNodeWithLabel(p.first);
            Node* g2 = genetree->GetNodeWithLabel(p.second);

            if (unmapped.find(g1) == unmapped.end() && unmapped.find(g2) == unmapped.end())
                orthologs.push_back(p);
        }

        cout<<"NBORTHOLOGS="<<orthologs.size()<<endl;

        for (set<Node*>::iterator it = unmapped.begin(); it != unmapped.end(); it++)
        {
            geneSpeciesMapping[(*it)] = speciestree->GetNodeWithLabel("Gorilla_gorilla");
        }

        ParalogyCorrector pc;
        ctree = pc.CorrectGeneTree(genetree, speciestree, geneSpeciesMapping, orthologs);



        //ctree might be NULL if no correction was applied
        if (ctree)
        {
            cout<<"TREEID : "<<treeID;
            cout<<"BEFORE : "<<NewickLex::ToNewickString(genetree);
            cout<<"AFTER : "<<NewickLex::ToNewickString(ctree);

            delete ctree;
        }
        else
        {
            cout<<"TREEID : "<<treeID<<" WAS OK";
        }


        geneSpeciesMapping.clear();
        geneSpeciesLines.clear();
        orthologs.clear();
        orthoLines.clear();


        delete genetree;
        delete speciestree;



    }
}



unordered_map<string, unordered_map<string, double> > loadDistances(string strDistances, string cacheID = "")
{
    if (cachedDistances.find(cacheID) != cachedDistances.end())
    {
        return cachedDistances[cacheID];
    }

    unordered_map<string, unordered_map<string, double> > distances;


    //std::ifstream ifs( filename );
    //std::string content( (std::istreambuf_iterator<char>(ifs) ),
    //                     (std::istreambuf_iterator<char>()    ));


    vector<vector<string> > splitLines;
    vector<string> lines = Util::Split(strDistances, "\n");

    //first line is irrelevant - it has no genes
    //we make a first pass to split every thing and store the splits
    //If we try to pare the file as it's read, we don't know the upcoming gene-gene correspondance
    for (int i = 1; i < lines.size(); i++)
    {
        vector<string> split = Util::Split(Util::Trim(lines[i]), " ", false);
        splitLines.push_back(split);
    }


    for (int i = 1; i < lines.size(); i++)
    {
        if (verbose > 0 && i % 100 == 0)
        {
            cout<<"distline lines parsed="<<i<<endl;
        }
        if (lines[i] != "")
        {
            vector<string> split = splitLines[i - 1];
            string lineGene = split[0];

            unordered_map<string, double> geneDists;
            if (distances.find(lineGene) != distances.end())
            {
                geneDists = distances[lineGene];
            }

            for (int j = 1; j < split.size(); j++)
            {
                vector<string> otherSplit = splitLines[j - 1];
                string otherGene = otherSplit[0];

                double dist = Util::ToDouble(split[j]);
                if (dist < 0 && hasNonnegativeDistanceFlag)
                    dist = 99999.9;
                geneDists[otherGene] = dist;
            }

            distances[lineGene] = geneDists;
        }

    }

    if (useCache)
        cachedDistances[cacheID] = distances;

    return distances;
}




Node* CorrectPolytomyByDistances(Node* speciesTree, Node* geneTree, string strDistances, string runID = "")
{
    Node* solved = NULL;

    unordered_map<string, unordered_map<string, double> > rawDists = loadDistances(strDistances, runID);


    if (verbose > 0)
        cout<<"DISTANCES FILE LOADED"<<endl;


    unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree, ";;", 1);
    unordered_map<Node*, Node*> lcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, geneSpeciesMapping);

    if (verbose > 0)
        cout<<"LCA MAPPING DONE"<<endl;

    //FROM HERE WE PUT THE DISTANCES INTO A MAP
    map<Node*, string> geneIDs;
    TreeIterator* it = geneTree->GetPostOrderIterator(true);
    while (Node* g = it->next())
    {
        vector<string> gstr = Util::Split(g->GetLabel(), ";;");
        geneIDs[g] = gstr[0];
    }
    geneTree->CloseIterator(it);


    map<Node*, map<Node*, double> > distances;
    vector<Node*> leaves = geneTree->GetLeafVector();
    for (int i = 0; i < leaves.size(); i++)
    {
        for (int j = i + 1; j < leaves.size(); j++)
        {
            string name1 = geneIDs[leaves[i]];
            string name2 = geneIDs[leaves[j]];

            if (rawDists.find(name1) != rawDists.end() && rawDists.find(name2) != rawDists.end())
            {
                distances[leaves[i]][leaves[j]] = rawDists[name1][name2];
                distances[leaves[j]][leaves[i]] = rawDists[name2][name1];
            }
            else
            {
                cout<<"DISTANCE PROBLEM !  RESULTS MIGHT NOT BE TRUSTWORTHY !"<<endl;
            }

        }
    }

    if (verbose > 0)
        cout<<"DISTANCE PARSING DONE"<<endl;


    PolySolverDistance solver;
    solver.verbose = verbose;
    solver.useCache = useCache;
    //Node* solved = solver.SolvePolytomy(leaves, speciesTree, lcaMapping, distances);
    solved = solver.SolvePolytomies(geneTree, speciesTree, lcaMapping, distances);


    rawDists.clear();
    geneSpeciesMapping.clear();
    lcaMapping.clear();
    distances.clear();
    leaves.clear();



    return solved;

}





/**
  Find one or many corrections for the given gene tree (in NewickFormat), and returns the string corresponding to the correction (format depends on _rerootMode)

  string speciesTreeString : species tree Newick
  string geneTreeString : gene tree Newick
  string distfilename : pairwise gene distances file
  string _rerootMode :
    one of {"none","findbestroot","outputallroots"}
    none : does standard correction - returns a single gene tree newick
    findbestroot : tries every possible rerooting of the gene tree and returns one with the minimum DL score (return value is a single newick line)
    outputallroots : returns the result for every rerooting, along with the dlscore in this format
        #before=[newick string 1]
        #dlscore=[newdlscore 1]
        [correctedtree 2]
        #before=[newick string 2]
        #dlscore=[newdlscore 2]
        [correctedtree 2]
        ...

        where #before gives the rooted tree before correction, #dlscore gives the dlscore after correction and the next line is the corrected tree

  bool _testEdgeRoots : true to try rooting on the gene tree edges, as opposed to rooting on nodes only.  Rooting on edges should not be necessary with respect to dlscores
  bool _hasNonnegativeDistanceFlag : if true, all negative distances will be replaced by a high distance (99999)
  bool _useCache : in the cases where we test all roots, cache can be applied to remember subtrees previously resolved.  Recommended, set to false only if you don't trust this feature and are patient enough.

  **/
extern "C" string CorrectPolytomyByDistancesStrings(string speciesTreeString, string geneTreeString, string strDistances, string _rerootMode, bool _testEdgeRoots, bool _hasNonnegativeDistanceFlag, bool _useCache)
{
    //TODO : not clean
    runCounter++;
    string runID = "CorrectPolytomyByDistancesStrings" + Util::ToString(runCounter);
    rerootMode = _rerootMode;
    useCache = _useCache;
    testEdgeRoots = _testEdgeRoots;
    hasNonnegativeDistanceFlag = _hasNonnegativeDistanceFlag;

    Node* geneTree = NewickLex::ParseNewickString(geneTreeString);
    Node* speciesTree = NewickLex::ParseNewickString(speciesTreeString, true);

    if (verbose)
    {
        cout<<"SpeciesTree="<<NewickLex::ToNewickString(speciesTree)<<endl<<endl;
        cout<<"GeneTree="<<NewickLex::ToNewickString(geneTree)<<endl<<endl;
    }


    vector<Node*> allSolved;
    vector<string> allBefore;
    vector<int> allDLScore;
    int bestDLScore = 99999;
    Node* bestSolved = NULL;
    if (rerootMode == "none")
    {
        Node* solved = CorrectPolytomyByDistances(speciesTree, geneTree, strDistances, runID);
        if (solved)
        {
            allSolved.push_back(solved);
            bestSolved = solved;
        }
    }
    else if (rerootMode == "findbestroot" || rerootMode == "outputallroots")
    {
        TreeIterator* it = geneTree->GetPostOrderIterator();
        vector<Node*> allTreesToTry;
        while (Node* n = it->next())
        {
            if (!n->IsLeaf())
            {
                Node* r = n->SetAsRootInCopy();

                allTreesToTry.push_back(r);
            }


            if (testEdgeRoots)
            {
                if (!n->IsRoot())
                {
                    Node* r = n->SetRootOnParentEdgeInCopy();

                    allTreesToTry.push_back(r);
                }
            }

        }
        geneTree->CloseIterator(it);


        for (int i = 0; i < allTreesToTry.size(); i++)
        {
            Node* r = allTreesToTry[i];
            Node* solved = CorrectPolytomyByDistances(speciesTree, r, strDistances, runID);
            if (solved)
            {
                allSolved.push_back(solved);
                allBefore.push_back(NewickLex::ToNewickString(r));

                //losing time here
                unordered_map<Node*, Node*> solvedmap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(solved, speciesTree, ";;", 1);
                int dlScore = GeneSpeciesTreeUtil::Instance()->GetDLScore(solved, speciesTree, solvedmap);
                allDLScore.push_back(dlScore);

                if (!bestSolved || dlScore < bestDLScore)
                {
                    bestDLScore = dlScore;
                    bestSolved = solved;
                }
            }


            delete r;
        }

    }

    string strout = "";
    if (allSolved.size() > 0)
    {
        for (int i = 0; i < allSolved.size(); i++)
        {
            Node* solved = allSolved[i];

            if (rerootMode == "outputallroots" || solved == bestSolved)
            {
                GeneSpeciesTreeUtil::Instance()->RelabelGenes(solved, ";;", "__");
                //GeneSpeciesTreeUtil::Instance()->RelabelGenesByIndex(solved, ";;", 0);

                if (strout != "")
                    strout += "\n";

                if (rerootMode == "outputallroots")
                {
                    strout += "#before=" + allBefore[i] + "\n" + "#dlscore=" + Util::ToString(allDLScore[i]) + "\n";
                }

                strout += NewickLex::ToNewickString(solved);

                if (verbose)
                {
                    cout<<NewickLex::ToNewickString(solved)<<endl;
                }
                /*if (outfilename == "")
                {
                    cout<<strout<<endl;
                    strout = "";
                }*/
            }

            delete solved;
        }
        allSolved.clear();

        /*if (outfilename != "")
        {#before=
            ofstream outfile;
            outfile.open (outfilename);

            outfile<<strout;
            outfile.close();
        }*/
    }



    delete geneTree;
    delete speciesTree;

    return strout;
}


extern "C" string CorrectParalogy(string geneTreeString, string speciesTreeString, boost::python::list geneSpeciesMapping, boost::python::list orthologs){

    ParalogyCorrector pc;
    string ctree = pc.CorrectGeneTree(geneTreeString, speciesTreeString, geneSpeciesMapping, orthologs);
    return ctree;
}



int main(int argc, char *argv[])
{
    //-gl 1 -s "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/Compara.73.species_tree" -g "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/famille_2.start_tree" -d "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/famille_2.dist" -r findbestroot -n -o "/u/lafonman/Projects/PolytomySolverDistance_1.2.2/example_files/bad/test_2.newick"
    //-gl 1 -s "/u/lafonman/Projects/PolytomySolverDistance_1.2.3/example_files/bad/Compara.73.species_tree" -g "/u/lafonman/Projects/PolytomySolverDistance_1.2.3/example_files/bad/famille_2.start_tree" -d "/u/lafonman/Projects/PolytomySolverDistance_1.2.3/example_files/bad/famille_2.dist" -r findbestroot -n -v
    //-gl 1 -s "/u/lafonman/Projects/PolytomySolverDistance_1.2.3/example_files/joseg/jospecies.txt" -g "/u/lafonman/Projects/PolytomySolverDistance_1.2.3/example_files/joseg/jogene.txt" -d "/u/lafonman/Projects/PolytomySolverDistance_1.2.3/example_files/joseg/jodist.txt" -r findbestroot -n -v


    string geneTreeString = "((A1,B1),C1);";
    string speciesTreeString = "((A,C)x,B)y;";

    boost::python::list geneSpeciesMapping;
    geneSpeciesMapping.append(boost::python::make_tuple("A1","A"));
    geneSpeciesMapping.append(boost::python::make_tuple("B1","B"));
    geneSpeciesMapping.append(boost::python::make_tuple("C1","C"));

    boost::python::list orthologs;
    orthologs.append(boost::python::make_tuple("A1","C1"));

    string corrected = CorrectParalogy(geneTreeString, speciesTreeString, geneSpeciesMapping, orthologs);
    cout << "corrected : " << corrected <<endl;


    map<string, string> args;

    string prevArg = "";
    for (int i = 0; i < argc; i++)
    {
        if (string(argv[i]) == "-v")
        {
            verbose = 1;
            prevArg = "";
        }
        else if (string(argv[i]) == "-v2")
        {
            verbose = 2;
            prevArg = "";
        }
        else if (string(argv[i]) == "-n")
        {
            hasNonnegativeDistanceFlag = true;
        }
        else if (string(argv[i]) == "-e")
        {
            testEdgeRoots = true;
        }
        else if (string(argv[i]) == "-cno")
        {
            useCache = false;
        }
        else
        {
            if (prevArg != "" && prevArg[0] == '-')
            {
                args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
            }

            prevArg = string(argv[i]);
        }
    }

    if (args.find("r") != args.end())
    {
        rerootMode = args["r"];
    }

    if (args.find("m") == args.end() || args["m"] == "polydistance")
    {
        bool isok = true;
        if (args.find("s") == args.end() && args.find("sn") == args.end())
        {
            cout<<"Specify a species tree file with -s or a species tree newick with -sn"<<endl;
            isok = false;
        }
        if (args.find("g") == args.end() && args.find("gn") == args.end())
        {
            cout<<"Specify a gene tree file with -g or a gene tree newick with -gn"<<endl;
            isok = false;
        }
        if (args.find("d") == args.end())
        {
            cout<<"Specify a distances matrix file with -d"<<endl;
            isok = false;
        }
        else
        {
            if (!Util::FileExists(args["d"]))
            {
                cout<<"Distances file does not exist"<<endl;
                isok = false;
            }
        }

        if (isok)
        {
            string scontent = "";
            string gcontent = "";

            if (args.find("s") != args.end())
            {
                scontent = Util::GetFileContent(args["s"]);
                //speciesTree = NewickLex::ParseNewickString(scontent, true);
            }
            else
            {
                scontent   = args["sn"];
                //speciesTree = NewickLex::ParseNewickString(args["sn"], true);
            }



            int geneTreeLine = -1;
            if (args.find("gl") != args.end())
            {
                geneTreeLine = Util::ToInt(args["gl"]);
            }

            if (args.find("g") != args.end())
            {
                gcontent = "";

                if (geneTreeLine > 0)
                {
                    gcontent = Util::GetFileLine(args["g"], geneTreeLine - 1);
                }
                else
                {
                    gcontent = Util::GetFileContent(args["g"]);
                }

                //geneTree = NewickLex::ParseNewickString(gcontent, false);
            }
            else
            {
                gcontent = args["gn"];
                //geneTree = NewickLex::ParseNewickString(args["gn"], false);
            }




            string distfile = args["d"];

            string strDistances = "";

            strDistances = Util::GetFileContent(distfile);



            if (verbose)
            {
                cout<<"DistFile="<<distfile<<endl<<endl;
            }

            string outfilename = "";
            if (args.find("o") != args.end())
            {
                outfilename = args["o"];
            }


            if (verbose)
            {
                cout<<"OutFile="<<outfilename<<endl<<endl;
            }


            string result = CorrectPolytomyByDistancesStrings(scontent, gcontent, strDistances, rerootMode, testEdgeRoots, hasNonnegativeDistanceFlag, useCache);

            if (outfilename != "")
            {
                ofstream outfile;
                outfile.open (outfilename);

                //outfile<<"dlscore="<<((PSDNodeInfo*)solved->nodeInfo)->dlScore<<endl;

                outfile<<result;
                outfile.close();
            }
            else
            {
                cout<<result<<endl;
            }


            /*vector<Node*> allSolved;
            vector<string> allBefore;
            vector<int> allDLScore;
            int bestDLScore = 9999;
            Node* bestSolved = NULL;
            if (rerootMode == "none")
            {
                Node* solved = CorrectPolytomyByDistances(speciesTree, geneTree, distfile, outfilename);
                if (solved)
                {
                    allSolved.push_back(solved);
                    bestSolved = solved;
                }
            }
            else if (rerootMode == "findbestroot" || rerootMode == "outputallroots")
            {
                TreeIterator* it = geneTree->GetPostOrderIterator();
                vector<Node*> allTreesToTry;
                while (Node* n = it->next())
                {
                    if (!n->IsLeaf())
                    {
                        Node* r = n->SetAsRootInCopy();

                        allTreesToTry.push_back(r);
                    }


                    if (testEdgeRoots)
                    {
                        if (!n->IsRoot())
                        {
                            Node* r = n->SetRootOnParentEdgeInCopy();

                            allTreesToTry.push_back(r);
                        }
                    }

                }
                geneTree->CloseIterator(it);


                for (int i = 0; i < allTreesToTry.size(); i++)
                {
                    Node* r = allTreesToTry[i];
                    Node* solved = CorrectPolytomyByDistances(speciesTree, r, distfile, outfilename);
                    if (solved)
                    {
                        allSolved.push_back(solved);
                        allBefore.push_back(NewickLex::ToNewickString(r));

                        //losing time here
                        unordered_map<Node*, Node*> solvedmap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(solved, speciesTree, ";;", 1);
                        int dlScore = GeneSpeciesTreeUtil::Instance()->GetDLScore(solved, speciesTree, solvedmap);
                        allDLScore.push_back(dlScore);

                        if (!bestSolved || dlScore < bestDLScore)
                        {
                            bestDLScore = dlScore;
                            bestSolved = solved;
                        }
                    }


                    delete r;
                }

            }

            string strout = "";
            if (allSolved.size() > 0)
            {
                for (int i = 0; i < allSolved.size(); i++)
                {
                    Node* solved = allSolved[i];

                    if (rerootMode == "outputallroots" || solved == bestSolved)
                    {
                        GeneSpeciesTreeUtil::Instance()->RelabelGenes(solved, ";;", "__");
                        //GeneSpeciesTreeUtil::Instance()->RelabelGenesByIndex(solved, ";;", 0);

                        if (strout != "")
                            strout += "\n";

                        if (rerootMode == "outputallroots")
                        {
                            strout += "#before=" + allBefore[i] + "\n" + "#dlscore=" + Util::ToString(allDLScore[i]) + "\n";
                        }

                        strout += NewickLex::ToNewickString(solved);

                        if (outfilename == "")
                        {
                            cout<<strout<<endl;
                            strout = "";
                        }
                    }

                    delete solved;
                }
                allSolved.clear();

                if (outfilename != "")
                {
                    ofstream outfile;
                    outfile.open (outfilename);

                    //outfile<<"dlscore="<<((PSDNodeInfo*)solved->nodeInfo)->dlScore<<endl;

                    outfile<<strout;
                    outfile.close();
                }
            }



            delete geneTree;
            delete speciesTree;*/
        }


    }



    return -1;
}

//Python wrapper
#include <boost/python.hpp>
using namespace boost::python;
BOOST_PYTHON_MODULE(libpolytomysolver)
{
    def("PolytomySolver", CorrectPolytomyByDistancesStrings);
    def("ParalogyCorrector", CorrectParalogy);
}
