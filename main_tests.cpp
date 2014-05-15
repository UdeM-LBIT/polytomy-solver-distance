#ifdef USE_MAIN_TEST

#include "trees/node.h"
#include "trees/polysolver.h"
#include "trees/newicklex.h"
#include "trees/genespeciestreeutil.h"
#include "trees/polysolver_nad.h"
#include "trees/polysolver_distance.h"


#include "div/tinydir.h"
#include "trees/autester.h"

#include <fstream>

#include <string>
#include <unordered_map>
#include <iostream>


using namespace std;




////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// FOR NAD CORRECTION
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

//IN MAIN :
//-gl 598 -s "C:\Users\Manuel\Desktop\School\PolytomySolver_home\data_eric\speciestree.newick" -g "C:\Users\Manuel\Desktop\School\PolytomySolver_home\data_eric\star_trees.trees" -d "C:\Users\Manuel\Desktop\School\PolytomySolver_home\data_eric\famille_598.dist" -v2

//-gl 1 -s "C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\speciestree_test.newick" -g "C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\test_tree.trees" -d "C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\famille_test.dist" -v -r outputallroots -o C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\test_out.txt
//-gl 1 -s "C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\speciestree.newick" -g "C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\star_trees.trees" -d "C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\famille_1.dist" -v -r outputallroots -o C:\Users\Manuel\Desktop\PolytomySolverDistance_1.2\example_files\test_out.txt
//HERE I SHAMELESSLY LEAVE MY JUNK CODE !

//EvaluateRFDistFile();
//    CorrectFolderNADs("/u/lafonman/Projects/all_protein_trees_70/testtrees/", "/u/lafonman/Projects/all_protein_trees_70/testtrees_corrected/");
//    FindTreesInWhichEnsemblCorrectedNADs();
//CorrectFolderNADs("C:/Users/Manuel/Desktop/School/PolytomySolver_home/data/fishtrees/", "C:/Users/Manuel/Desktop/School/PolytomySolver_home/data/fishtrees_corrected_clever/");

//CompareNADCounts("C:/Users/Manuel/Desktop/School/PolytomySolver_home/data/fishtrees_corrected_clever/", "C:/Users/Manuel/Desktop/School/PolytomySolver_home/data/fishtrees_corrected_v2/");
//    PolySolverNAD solver;
//    solver.PerformRandomTest(4);
//DoParalogyCorrection();
//return 0;


/*for (int k = 13; k <= 15; k++)
{
    for (int j = 0; j < 1000; j++)
    {
        PolySolverNAD* solver = new PolySolverNAD();
        cout<<"Test k="<<k<<", iter="<<j<<" : ";
        solver->PerformRandomTest(k); //, (k == 7) ? 1 : 0);
        delete solver;
    }
}

return 0;
//string sstr = "(((XL0, (XL3, XL1))S11, (XR0, XR1)S6)_0_1, (((XL4, ((AD_G4_G8, (AD_G3_G8, XL8)), XR3))S15, ((XL2, ((AD_G2_G4, (XL7, XL5)), (XR4, AD_G2_G3)))S4, (XR5, (XR7, XR2))S13)_2_5_7)_4, (((((XL6, XR8)S10, ((S7, S9), (S0, S16))), (S3, XR6))_6, ((S2, (S8, S5)), S1)), (S14, S17)))_8)_3;";
//string gstr = "((XL0, XR0)G0, (XL1, XR1)G1, (AD_G2_G4, (AD_G2_G3, (XL2, XR2)))G2, (AD_G2_G3, (AD_G3_G8, (XL3, XR3)))G3, (XR4, ((XL4, AD_G2_G4), AD_G4_G8))G4, (XL5, XR5)G5, (XL6, XR6)G6, (XL7, XR7)G7, (XL8, ((XR8, AD_G3_G8), AD_G4_G8))G8);";


string sstr = "(((XL2, (XL5, (XL0, AD_G2_G5)))S2, ((XL6, (XR2, XR0))S4, XR6)_6)_0_2, ((((AD_G1_G7, XL4), ((XL7, AD_G4_G5), (XL1, XR5)))S1, (XR7, (XR4, XR1))S5)_1_4_7, (XL3, XR3)_3))_5;";
string gstr = "((XL0, XR0)G0, (XR1, (XL1, AD_G1_G7))G1, (XR2, (XL2, AD_G2_G5))G2, (XL3, XR3)G3, (XL4, (AD_G4_G5, XR4))G4, (AD_G4_G5, (XL5, (AD_G2_G5, XR5)))G5, (XL6, XR6)G6, (AD_G1_G7, (XL7, XR7))G7);";




Node* s = NewickLex::ParseNewickString(sstr, true);
Node* g = NewickLex::ParseNewickString(gstr);

unordered_map<Node*,Node*> lcaMap = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(g, s, "*", 0);

cout<<"NADS B4="<<GeneSpeciesTreeUtil::Instance()->GetNADNodes(g,s,lcaMap).size()<<endl;

PolySolverNAD solverX;
Node* solved = solverX.SolvePolytomy(g, s, lcaMap);

cout<<"SOLVED="<<endl<<NewickLex::ToNewickString(solved)<<endl;
unordered_map<Node*,Node*> lcaMapSolved = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(solved, s, "*", 0);
vector<Node*> nads = GeneSpeciesTreeUtil::Instance()->GetNADNodes(solved,s,lcaMapSolved);
cout<<"NADS="<<nads.size()<<endl;

cout<<"COMPS="<<solverX.last_nb_ad_components;
return 0;



vector<Node*> forest;
for (int i = 0; i < g->GetNbChildren(); i++)
{
    forest.push_back(g->GetChild(i));
}


SADNADGraph graph;
graph.BuildGraph(forest, s, lcaMap);

vector<pair<Node*, Node*> > useful = graph.GetUsefulSpeciationEdges();
/*useful.clear();

useful.push_back( make_pair(g->GetNodeWithLabel("G7"), g->GetNodeWithLabel("G0")));
useful.push_back( make_pair(g->GetNodeWithLabel("G6"), g->GetNodeWithLabel("G2")));
useful.push_back( make_pair(g->GetNodeWithLabel("G6"), g->GetNodeWithLabel("G0")));
useful.push_back( make_pair(g->GetNodeWithLabel("G4"), g->GetNodeWithLabel("G1")));
//


cout<<"ALL USEFUL EDGES"<<endl;
for (vector<pair<Node*, Node*> >::iterator it = useful.begin(); it != useful.end(); it++)
{
    cout<<(*it).first->GetLabel()<<" "<<(*it).second->GetLabel()<<endl;
}
cout<<endl<<endl;

set<pair<Node*, Node*> > ch;

PolySolverNAD solver;
solver.bestSolSoFar = 9999;
solver.bestChoiceSoFar.clear();
solver.FindBestSolution(graph, useful, ch, 0);



cout<<"BEST SOL="<<solver.bestSolSoFar<<" NADS"<<endl;
for (set<pair<Node*, Node*> >::iterator it = solver.bestChoiceSoFar.begin(); it != solver.bestChoiceSoFar.end(); it++)
{
    cout<<(*it).first->GetLabel()<<" "<<(*it).second->GetLabel()<<endl;
}

cout<<"------------------------"<<endl;

return -1;
*/














map<string, string> FindTreesInWhichEnsemblCorrectedNADs()
{
    //CORRESPONDANCE FORMAT : V74FILE:V70FILE


    map<string, string> treePairsOfInterest;

    string speciesNewick = "((((((Xiphophorus maculatus, Oryzias latipes), Gasterosteus aculeatus),Oreochromis niloticus), (Takifugu rubripes,Tetraodon nigroviridis)), Gadus morhua), Danio rerio)";
    Node* speciesTree = NewickLex::ParseNewickString(speciesNewick, true);
    string corr = Util::GetFileContent("/u/lafonman/Projects/protein_trees_correspondance.txt");

    vector<string> lines = Util::Split(corr, "\n", false);

    int cptMoreNadsIn74 = 0;
    int cptLessNadsIn74 = 0;
    int cptSameNads = 0;
    for (int i = 0; i < lines.size(); i++)
    {
        vector<string> parts = Util::Split(lines[i],":", false);

        if (parts.size() == 2 && parts[1] != "NONE")
        {
            string t1str = Util::GetFileContent(parts[0]);
            Node* pre_tree1 = NewickLex::ParseNewickString(t1str);
            Node* tree1 = pre_tree1;

            //HACK ALERT !  The following is due to a bug in which the tree might have a single child root.
            if (pre_tree1->GetNbChildren() == 1)
            {
                tree1 = pre_tree1->GetChild(0);
            }

            unordered_map<Node*, Node*> geneSpeciesMapping1 = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(tree1, speciesTree, "___", 1);
            unordered_map<Node*, Node*> lcaMapping1 = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(tree1, speciesTree, geneSpeciesMapping1);

            string t2str = Util::GetFileContent(parts[1]);
            Node* pre_tree2 = NewickLex::ParseNewickString(t2str);
            Node* tree2 = pre_tree2;
            if (pre_tree2->GetNbChildren() == 1)
            {
                tree2 = pre_tree2->GetChild(0);
            }

            unordered_map<Node*, Node*> geneSpeciesMapping2 = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(tree2, speciesTree, "___", 1);
            unordered_map<Node*, Node*> lcaMapping2 = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(tree2, speciesTree, geneSpeciesMapping2);

            vector<Node*> nads1 = GeneSpeciesTreeUtil::Instance()->GetNADNodes(tree1, speciesTree, lcaMapping1);

            vector<Node*> nads2 = GeneSpeciesTreeUtil::Instance()->GetNADNodes(tree2, speciesTree, lcaMapping2);




            //same # of leaves, but less nads = interesting
            if (geneSpeciesMapping1.size() == geneSpeciesMapping2.size())
            {
                if (nads1.size() < nads2.size())
                    cptLessNadsIn74++;
                else if (nads1.size() > nads2.size())
                    cptMoreNadsIn74++;
                else
                    cptSameNads++;
                //cout<<nads1.size()<<" "<<nads2.size()<<endl;


                if (false && nads1.size() < nads2.size())
                {
                    cout<<parts[0]<<":"<<parts[1]<<" --> "<<nads2.size()<<" - "<<nads1.size()<<" = "<<(nads2.size() - nads1.size())<<endl;
                    cout<<geneSpeciesMapping1.size()<<" "<<geneSpeciesMapping2.size()<<endl;
                    treePairsOfInterest[parts[0]] = parts[1];

                    string v70corrected_file = Util::ReplaceAll(parts[1], "fishtrees", "fishtrees_corrected_v2");
                    /*string v70corrected_newick = Util::GetFileContent(v70corrected_file);

                    Node* pre_v70corrected = NewickLex::ParseNewickString(v70corrected_newick);
                    Node* v70_corrected = pre_v70corrected;
                    if (v70_corrected->GetNbChildren() == 1)
                        v70_corrected = v70_corrected->GetChild(0);
                    unordered_map<Node*, Node*> geneSpeciesMapping70_corrected = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(v70_corrected, speciesTree, "___", 1);
                    unordered_map<Node*, Node*> lcaMapping70_corrected = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(v70_corrected, speciesTree, geneSpeciesMapping70_corrected);
                    vector<Node*> nads70Corrected = GeneSpeciesTreeUtil::Instance()->GetNADNodes(v70_corrected, speciesTree, lcaMapping70_corrected);

                    v70corrected_newick = NewickLex::ToNewickString(v70_corrected);
                    delete pre_v70corrected;*/

                    string tree1newick = NewickLex::ToNewickString(tree1);
                    string tree2newick = NewickLex::ToNewickString(tree2);

                    string f74 = parts[0];
                    string f70 = parts[1];

                    string pairFile = "/u/lafonman/Projects/nad_data/v70_vs_v74/" + Util::GetPathFilename(f70) + "__vs__" + Util::GetPathFilename(f74);

                    string pairContent = "V70FILE=" + f70 + "\n" + "V70NADS=" + Util::ToString((int)nads2.size()) + "\n" + "V70TREE=" + tree2newick + "\n" +
                                          "V70CORRFILE=" + v70corrected_file + "\n" +  //"V70CORRNADS=" + Util::ToString((int)nads70Corrected.size()) + "\n" + "V70CORRTREE=" + v70corrected_newick + "\n" +
                                         "V74FILE=" + f74 + "\n" + "V74NADS=" + Util::ToString((int)nads1.size()) + "\n" + "V74TREE=" + tree1newick + "\n";

                    cout<<"WRITING "<<pairFile<<endl;
                    Util::WriteFileContent(pairFile, pairContent);
                }

            }



            delete pre_tree1;
            delete pre_tree2;
        }
        cout<<"MORE-NADS-IN-74="<<cptMoreNadsIn74<<endl<<"LESS-NADS-IN-74="<<cptLessNadsIn74<<endl<<"SAME-NADS="<<cptSameNads<<endl<<"TOTAL="<<(cptMoreNadsIn74 + cptLessNadsIn74 + cptSameNads)<<endl;
    }

    delete speciesTree;

    return treePairsOfInterest;
}







void EvaluateRFDistFile()
{
    /*tinydir_dir dir;
    tinydir_open(&dir, "/u/lafonman/Projects/all_protein_trees_70/fishtrees_corrected_v2/");
    int cnt = 0;
    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);


        if (!file.is_dir)
        {
            string fullname = "/u/lafonman/Projects/all_protein_trees_70/fishtrees_corrected_v2/" + string(file.name);


            string poly1 = "";
            string poly2 = "";
            string corrclade = "";
            string nadsBefore = "";
            string nadsAfter = "";
            string v70corrcontent = Util::GetFileContent(fullname);
            vector<string> v70lines = Util::Split(v70corrcontent, "\n", false);

            for (int v70l = 0; v70l < v70lines.size(); v70l++)
            {
                vector<string> kv = Util::Split(v70lines[v70l], "=", false);

                string key = kv[0];
                string val = kv[1];

                if (key == "POLYSIZE1")
                    poly1 = val;
                if (key == "POLYSIZE2")
                    poly2 = val;
                if (key == "CORRECTEDCLADE")
                    corrclade = val;
                if (key == "NADSBEFORE")
                    nadsBefore = val;
                if (key == "NADSAFTER")
                    nadsAfter = val;


            }

            if (Util::ToInt(poly1) > 3 || Util::ToInt(poly2) > 2)
            {
                cout<<"P1="<<poly1<<" P2="<<poly2<<" NB="<<nadsBefore<<" NA="<<nadsAfter<<" F="<<fullname<<endl;
                cnt++;
            }

        }
        tinydir_next(&dir);
    }
    tinydir_close(&dir);

    cout<<cnt<<endl;*/

    vector<string> fieldsToOutput;
    fieldsToOutput.push_back("V70C_POLYSIZE1");
    fieldsToOutput.push_back("V70C_POLYSIZE2");
    fieldsToOutput.push_back("V70C_NADSBEFORE");
    fieldsToOutput.push_back("V70C_NADSAFTER");
    fieldsToOutput.push_back("V74_NADS");
    fieldsToOutput.push_back("V70C_DLSCORE");
    fieldsToOutput.push_back("V74_DLSCORE");
    fieldsToOutput.push_back("V74_NADS");
    fieldsToOutput.push_back("V74_CLADETYPE");
    fieldsToOutput.push_back("V74_ISGONE");
    fieldsToOutput.push_back("NBLEAVES");
    fieldsToOutput.push_back("CORRECTED_CLADE_SIZE");
    fieldsToOutput.push_back("RF70_70CORR_V");
    fieldsToOutput.push_back("RF70_74_V");
    fieldsToOutput.push_back("RF70CORR_74_V");
    fieldsToOutput.push_back("RF70_70CORR");
    fieldsToOutput.push_back("RF70_74");
    fieldsToOutput.push_back("RF70CORR_74");
    fieldsToOutput.push_back("V70CORRFILE");
    fieldsToOutput.push_back("V74FILE");


    for (int k = 0; k < fieldsToOutput.size(); k++)
    {
        if (k != 0)
            cout<<" ";
        cout<<fieldsToOutput[k];
    }
    cout<<endl;

    string speciesNewick = "((((((Xiphophorus maculatus, Oryzias latipes), Gasterosteus aculeatus),Oreochromis niloticus), (Takifugu rubripes,Tetraodon nigroviridis)), Gadus morhua), (Astyanax mexicanus,Danio rerio))";
    Node* speciesTree = NewickLex::ParseNewickString(speciesNewick, true);


    int cptSpec = 0, cptNAD = 0, cptAD = 0, cptGone = 0;
    string content = Util::GetFileContent("/u/lafonman/Projects/nad_data/rfdists_v5.txt");

    vector<string> lines = Util::Split(content, "\n", false);

    for (int i = 0; i < lines.size(); i++)
    {
        map<string, string> theLineData;

        string rfdist = "";
        string v70corrfile = "";
        string nadsBefore = "";
        string nadsAfter = "";
        string v74nads = "";

        Node* v70Tree = NULL;
        Node* v74Tree = NULL;
        Node* v70correctedTree = NULL;

        vector<string> correctedClade;

        vector<string> params = Util::Split(lines[i], ";;", false);

        for (int p = 0; p < params.size(); p++)
        {
            vector<string> kv = Util::Split(params[p], "=");
            string key = kv[0];
            string val = kv[1];

            if (key == "RF")
                rfdist = val;
            if (key == "V70CORRFILE")
                v70corrfile = val;
            if (key == "V74NADS")
                v74nads = val;
            if (key == "V74FILE")
            {
                string ttttt = Util::GetFileContent(val);
                v74Tree = NewickLex::ParseNewickString(ttttt);

                //TODO : MEMORY LEAK HERE
                if (v74Tree->GetNbChildren() == 1)
                    v74Tree = v74Tree->GetChild(0);
            }
            if (key == "V70FILE")
            {
                string uuuuu = Util::GetFileContent(val);
                v70Tree = NewickLex::ParseNewickString(uuuuu);

                //TODO : MEMORY LEAK HERE
                if (v70Tree->GetNbChildren() == 1)
                    v70Tree = v70Tree->GetChild(0);
            }

            theLineData[key] = val;

        }


        string poly1 = "";
        string poly2 = "";

        if (v70corrfile != "")
        {
            string v70corrcontent = Util::GetFileContent(v70corrfile);
            vector<string> v70lines = Util::Split(v70corrcontent, "\n", false);

            for (int v70l = 0; v70l < v70lines.size(); v70l++)
            {
                vector<string> kv = Util::Split(v70lines[v70l], "=", false);

                string key = kv[0];
                string val = kv[1];

                if (key == "POLYSIZE1")
                    poly1 = val;
                if (key == "POLYSIZE2")
                    poly2 = val;
                if (key == "NADSBEFORE")
                    nadsBefore = val;
                if (key == "NADSAFTER")
                    nadsAfter = val;
                if (key == "CORRECTEDTREE")
                {
                    v70correctedTree = NewickLex::ParseNewickString(val);
                    //TODO : MEMORY LEAK HERE
                    if (v70correctedTree->GetNbChildren() == 1)
                        v70correctedTree = v70correctedTree->GetChild(0);
                }
                if (key == "CORRECTEDCLADE")
                {
                    correctedClade = Util::Split(val, ";");
                }

            }
        }

        theLineData["V70C_POLYSIZE1"] = poly1;
        theLineData["V70C_POLYSIZE2"] = poly2;
        theLineData["V70C_NADSBEFORE"] = nadsBefore;
        theLineData["V70C_NADSAFTER"] = nadsAfter;

        theLineData["RF70_70CORR_V"] = Util::GetSubstringBefore(theLineData["RF70_70CORR"], "/");
        theLineData["RF70_74_V"] = Util::GetSubstringBefore(theLineData["RF70_74"], "/");
        theLineData["RF70CORR_74_V"] = Util::GetSubstringBefore(theLineData["RF70CORR_74"], "/");


        theLineData["V74_NADS"] = v74nads;

        //theLineData["RFDIST"] = rfdist;

        theLineData["V70CORRFILE"] = v70corrfile;

        if (v70Tree && v74Tree && v70correctedTree)
        {
            set<Node*> v70leaves = v70correctedTree->GetLeafSet();
            theLineData["NBLEAVES"] = Util::ToString((int)v70leaves.size());


            unordered_map<Node*, Node*> v70correctedMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(v70correctedTree, speciesTree, "___", 1);
            int v70corrected_DL = GeneSpeciesTreeUtil::Instance()->GetDLScore(v70correctedTree, speciesTree, v70correctedMapping);
            theLineData["V70C_DLSCORE"] = Util::ToString(v70corrected_DL);

            unordered_map<Node*, Node*> v74mapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(v74Tree, speciesTree, "___", 1);
            int v74_DL = GeneSpeciesTreeUtil::Instance()->GetDLScore(v74Tree, speciesTree, v74mapping);
            theLineData["V74_DLSCORE"] = Util::ToString(v74_DL);

            vector<Node*> v74CladeLeaves;
            for (int k = 0; k < correctedClade.size(); k++)
            {
                Node* g = v74Tree->GetNodeWithLabel(correctedClade[k]);

                v74CladeLeaves.push_back(g);
            }

            theLineData["CORRECTED_CLADE_SIZE"] = Util::ToString((int)correctedClade.size());

            Node* lca = v74Tree->FindLCA(v74CladeLeaves);

            set<Node*> lcaLeaves = lca->GetLeafSet();


            string lcaType = "";
            if (v74mapping[lca->GetChild(0)] != v74mapping[lca] && v74mapping[lca->GetChild(1)] != v74mapping[lca])
            {
                lcaType = "SPEC";
                //cout<<"CLADE IS A SPEC"<<endl;

                theLineData["V74_CLADETYPE"] = "S";


                //-------------------------------------------------------------
                //SPECIAL CODE HERE : run au test on ensembl made 'S' and polysize is interesting
                //if (Util::ToInt(poly1) > 2 || Util::ToInt(poly2) > 2)
                //{
                //    RunAUTestOnTrees(Util::GetSubstringAfter(theLineData["V74FILE"], "/"), v74Tree, v70correctedTree);
                //}
                //-------------------------------------------------------------



                cptSpec++;
            }
            else
            {
                if (!GeneSpeciesTreeUtil::Instance()->HaveCommonSpecies(lca->GetChild(0), lca->GetChild(1), v74mapping))
                {
                    lcaType = "NAD";
                    theLineData["V74_CLADETYPE"] = "NAD";
                    //cout<<"CLADE IS A NAD"<<endl;
                    cptNAD++;
                }
                else
                {
                    lcaType = "AD";
                    theLineData["V74_CLADETYPE"] = "AD";
                    //cout<<"CLADE IS AN AD"<<endl;
                    cptAD++;
                }
            }

            if (lcaLeaves.size() != correctedClade.size())
            {
                //cout<<"CLADE IS GONE (type="<<lcaType<<")"<<endl;
                theLineData["V74_ISGONE"] = "1";
                cptGone++;
            }
            else
            {
                theLineData["V74_ISGONE"] = "0";
            }

            //cout<<"V70C_DL="<<v70corrected_DL<<" V74DL="<<v74_DL<<endl;

            for (int k = 0; k < fieldsToOutput.size(); k++)
            {
                string fval = "-";
                if (theLineData.find(fieldsToOutput[k]) != theLineData.end())
                    fval = Util::ReplaceAll(theLineData[fieldsToOutput[k]], " ", "_");
                if (k != 0)
                    cout<<" ";
                cout<<fval;
            }
            cout<<endl;

        }




 //       if (Util::ToInt(poly1) > 2 || Util::ToInt(poly2) > 2)
        {

            /*if (Util::ToInt(nadsBefore) == 1)
            {
                //cout<<i<<" : RF="<<rfdist<<" P1="<<poly1<<" P2="<<poly2<<" NB="<<nadsBefore<<" NA="<<nadsAfter<<" NA(Ensembl)="<<v74nads<<" V70="<<v70corrfile<<endl;
                cpt++;

                if (rfdist.find("0/") != string::npos)
                    cptrf0++;
            }*/
        }

        correctedClade.clear();
        delete v70Tree;
        delete v74Tree;
        delete v70correctedTree;



    }

    //cout<<"GONE="<<cptGone<<" SPEC="<<cptSpec<<" NAD="<<cptNAD<<" AD="<<cptAD<<endl;

    //cout<<cpt<<endl;
    //cout<<cptrf0<<endl;
}








void CorrectFolderNADs(string folderIn, string folderOut)
{
    string speciesNewick = "((((((Xiphophorus maculatus, Oryzias latipes), Gasterosteus aculeatus),Oreochromis niloticus), (Takifugu rubripes,Tetraodon nigroviridis)), Gadus morhua), (Astyanax mexicanus,Danio rerio))";
    Node* speciesTree = NewickLex::ParseNewickString(speciesNewick, true);

    tinydir_dir dir;
    tinydir_open(&dir, folderIn.c_str());

    int cpt = 0, cptCorrected = 0;
    int cptDidNothing = 0, cptKilledOne = 0, cptKilledMore = 0, cptCreatedMore = 0;

    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);


        if (file.is_dir)
        {
            ;
        }
        else
        {
            string filename = string(file.name);    //for some reason, calling file.name later causes crashes
            string fullname = folderIn + filename;


            if (fullname.find(".newick") != string::npos)
            {

                std::ifstream ifs( fullname );
                std::string content( (std::istreambuf_iterator<char>(ifs) ),
                                     (std::istreambuf_iterator<char>()    ));
                //cout<<fullname<<endl<<content<<endl;


                Node* geneTree = NewickLex::ParseNewickString(content);
                Node* realGeneTree = geneTree;

                //HACK ALERT !  The following is due to a bug in which the tree might have a single child root.
                //PolySolver expects its input gene trees to be binary
                if (realGeneTree->GetNbChildren() == 1)
                {
                    realGeneTree = realGeneTree->GetChild(0);
                }

                unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(realGeneTree, speciesTree, "___", 1);

                PolySolverNAD solver;
                PolySolverCorrectionInfo info = solver.CorrectHighestNAD(realGeneTree, speciesTree, geneSpeciesMapping);

                Node* solved = info.correction;
                if (solved)
                {
                    Node* realSolved = solved;
                    if (solved->GetNbChildren() == 1)
                        realSolved = solved->GetChild(0);

                    //TODO : lca mappings were already computed by solver
                    //cout<<NewickLex::ToNewickString(realSolved);
                    //return;

                    unordered_map<Node*, Node*> beforeLCAMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(realGeneTree, speciesTree, geneSpeciesMapping);
                    vector<Node*> nadsBefore = GeneSpeciesTreeUtil::Instance()->GetNADNodes(realGeneTree, speciesTree, beforeLCAMapping);


                    unordered_map<Node*, Node*> afterGeneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(realSolved, speciesTree, "___", 1);
                    unordered_map<Node*, Node*> afterLCAMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(realSolved, speciesTree, afterGeneSpeciesMapping);
                    vector<Node*> nadsAfter = GeneSpeciesTreeUtil::Instance()->GetNADNodes(realSolved, speciesTree, afterLCAMapping);

                    cout<<"BEFORE : "<<nadsBefore.size()<<", AFTER : "<<nadsAfter.size()<<endl;

                    if (nadsBefore.size() == nadsAfter.size())
                        cptDidNothing++;
                    if (nadsBefore.size() == nadsAfter.size() + 1)
                        cptKilledOne++;
                    if (nadsBefore.size() > nadsAfter.size() + 1)
                        cptKilledMore++;
                    if (nadsBefore.size() < nadsAfter.size())
                        cptCreatedMore++;

                    cptCorrected++;


                    string baseDir = "/u/lafonman/Projects/PolytomySolver/data/";
                    string treePairsDir = baseDir + "treepairs/";
                    string fastaDir = baseDir + "fasta/";
                    string alignementsDir = baseDir + "alignments/";
                    string conseloutDir = baseDir + "conselout/";
                    string treeID = Util::GetSubstringAfter(Util::ReplaceAll(fullname, ".newick", ""), "/");

                    //AUTester tester;
                    //tester.RunAUTest(realGeneTree, solved, treeID, "___", 0, 1, treePairsDir, fastaDir, alignementsDir, conseloutDir);


                    cout<<"WRITING "<<folderOut<<filename<<endl;
                    string outtreefile = folderOut + filename;
                    string newick = NewickLex::ToNewickString(solved);
                    string beforeNewick = NewickLex::ToNewickString(realGeneTree);


                    string txt = "ORIGINALFILE=" + fullname + "\n" +
                            "ORIGINALTREE=" + beforeNewick + "\n" +
                            "NADSBEFORE=" + Util::ToString((int)nadsBefore.size()) + "\n" +
                            "NADSAFTER=" + Util::ToString((int)nadsAfter.size()) + "\n" +
                            "POLYSIZE1=" + Util::ToString(info.firstPolySize) + "\n" +
                            "POLYSIZE2=" + Util::ToString(info.secondPolySize) + "\n" +
                            "CORRECTEDCLADE=";
                    for (int i = 0; i < info.nadCladeGenes.size();i++)
                    {
                        if (i != 0)
                            txt += ";";
                        txt += info.nadCladeGenes[i];
                    }
                    txt = txt + "\n" + "CORRECTEDTREE=" + newick;


                    Util::WriteFileContent(outtreefile, txt);

                    delete solved;
                }

                delete geneTree;

                cpt++;
            }
        }


        tinydir_next(&dir);
    }

    tinydir_close(&dir);

    cout<<"FILES="<<cpt<<endl<<"CORRECTED="<<cptCorrected<<endl<<
          "KILLED ONE="<<cptKilledOne<<endl<<"KILLED MORE="<<cptKilledMore<<endl<<
          "DID NOTHING="<<cptDidNothing<<endl<<"HARMFUL="<<cptCreatedMore<<endl;

    delete speciesTree;
}


void CompareNADCounts(string folder1, string folder2)
{
    string speciesNewick = "((((((Xiphophorus maculatus, Oryzias latipes), Gasterosteus aculeatus),Oreochromis niloticus), (Takifugu rubripes,Tetraodon nigroviridis)), Gadus morhua), (Astyanax mexicanus,Danio rerio))";
    Node* speciesTree = NewickLex::ParseNewickString(speciesNewick, true);

    tinydir_dir dir;
    tinydir_open(&dir, folder1.c_str());
    int cntDiff = 0;
    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);


        if (!file.is_dir)
        {
            string fullname1 = folder1 + string(file.name);
            string fullname2 = folder2 + string(file.name);

            if (!Util::FileExists(fullname1) || !Util::FileExists(fullname2))
                continue;

            string c1 = Util::GetFileContent(fullname1);
            Node* t1 = NewickLex::ParseNewickString(c1);

            //TODO : leak here, root not deleted
            if (t1->GetNbChildren() == 1)
                t1 = t1->GetChild(0);

            string c2 = Util::GetFileContent(fullname2);
            Node* t2 = NewickLex::ParseNewickString(c2);
            if (t2->GetNbChildren() == 1)
                t2 = t2->GetChild(0);

            unordered_map<Node*, Node*> geneSpeciesMapping1 = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(t1, speciesTree, "___", 1);
            unordered_map<Node*, Node*> lcaMapping1 = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(t1, speciesTree, geneSpeciesMapping1);
            vector<Node*> nads1 = GeneSpeciesTreeUtil::Instance()->GetNADNodes(t1, speciesTree, lcaMapping1);

            unordered_map<Node*, Node*> geneSpeciesMapping2 = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(t2, speciesTree, "___", 1);
            unordered_map<Node*, Node*> lcaMapping2 = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(t2, speciesTree, geneSpeciesMapping2);
            vector<Node*> nads2 = GeneSpeciesTreeUtil::Instance()->GetNADNodes(t2, speciesTree, lcaMapping2);

            cout<<nads1.size()<<" VS "<<nads2.size()<<endl;

            cntDiff += (int)(nads1.size() - nads2.size());

            delete t1;
            delete t2;
        }
        tinydir_next(&dir);
    }
    tinydir_close(&dir);

    cout<<"CNTDIFF="<<cntDiff<<endl;

    delete speciesTree;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// END NAD CORRECTION
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////























































map<string, string> FindCorrespondingTrees(string folder1, string folder2)
{
    map<string,string> correspondance;
    map<Node*, string> treesToFile1;
    map<Node*, string> treesToFile2;

    tinydir_dir dir;
    tinydir_open(&dir, folder1.c_str());
    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);


        if (!file.is_dir)
        {
            string fullname = folder1 + string(file.name);



            string content1 = Util::GetFileContent(fullname);
            Node* tree1 = NewickLex::ParseNewickString(content1);

            if (tree1->GetLeafVector().size() > 3)
            {
                treesToFile1[tree1] = fullname;
            }
        }
        tinydir_next(&dir);
    }
    tinydir_close(&dir);

    cout<<"LOADED DIR1 WITH "<<treesToFile1.size()<<" CONSIDERED"<<endl;


    tinydir_open(&dir, folder2.c_str());
    while (dir.has_next)
    {
        tinydir_file file2;
        tinydir_readfile(&dir, &file2);


        if (!file2.is_dir)
        {
            string fullname2 = folder2 + string(file2.name);

            string content2 = Util::GetFileContent(fullname2);
            Node* tree2 = NewickLex::ParseNewickString(content2);

            if (tree2->GetLeafVector().size() > 3)
            {
                treesToFile2[tree2] = fullname2;
            }
        }
        tinydir_next(&dir);
    }
    tinydir_close(&dir);



    cout<<"LOADED DIR2 WITH "<<treesToFile2.size()<<" CONSIDERED"<<endl;


    for (map<Node*, string>::iterator t1it = treesToFile1.begin(); t1it != treesToFile1.end(); t1it++)
    {
        Node* tree1 = (*t1it).first;
        string file1 = (*t1it).second;
        set<Node*> leaves1 = tree1->GetLeafSet();

        for (map<Node*, string>::iterator t2it = treesToFile2.begin(); t2it != treesToFile2.end(); t2it++)
        {
            Node* tree2 = (*t2it).first;
            string file2 = (*t2it).second;

            set<Node*> leaves2 = tree2->GetLeafSet();

            bool foundEverything = true;    //true until proven otherwise
            for (set<Node*>::iterator it1 = leaves1.begin(); it1 != leaves1.end() && foundEverything; it1++)
            {
                bool foundN1 = false;
                Node* n1 = (*it1);
                for (set<Node*>::iterator it2 = leaves2.begin(); it2 != leaves2.end() && !foundN1; it2++)
                {
                    Node* n2 = (*it2);
                    if (n1->GetLabel() == n2->GetLabel())
                    {
                        foundN1 = true;
                    }
                }

                if (!foundN1)
                    foundEverything = false;
            }

            if (foundEverything)
            {
                correspondance[file1] = file2;
                cout<<file1<<":"<<file2<<endl;
            }

        }

        if (correspondance.find(file1) == correspondance.end())
        {
            cout<<file1<<":NONE"<<endl;
        }
    }

    //TODO : delete trees

    return correspondance;
}




map<string, string> FindTreesInWhichEnsemblCorrectedNADs()
{
    //CORRESPONDANCE FORMAT : V74FILE:V70FILE


    map<string, string> treePairsOfInterest;

    string speciesNewick = "((((((Xiphophorus maculatus, Oryzias latipes), Gasterosteus aculeatus),Oreochromis niloticus), (Takifugu rubripes,Tetraodon nigroviridis)), Gadus morhua), Danio rerio)";
    Node* speciesTree = NewickLex::ParseNewickString(speciesNewick, true);
    string corr = Util::GetFileContent("/u/lafonman/Projects/protein_trees_correspondance.txt");

    vector<string> lines = Util::Split(corr, "\n", false);

    for (int i = 0; i < lines.size(); i++)
    {
        vector<string> parts = Util::Split(lines[i],":", false);

        if (parts.size() == 2 && parts[1] != "NONE")
        {
            string t1str = Util::GetFileContent(parts[0]);
            Node* pre_tree1 = NewickLex::ParseNewickString(t1str);
            Node* tree1 = pre_tree1;

            //HACK ALERT !  The following is due to a bug in which the tree might have a single child root.
            if (pre_tree1->GetNbChildren() == 1)
            {
                tree1 = pre_tree1->GetChild(0);
            }

            unordered_map<Node*, Node*> geneSpeciesMapping1 = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(tree1, speciesTree, "___", 1);
            unordered_map<Node*, Node*> lcaMapping1 = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(tree1, speciesTree, geneSpeciesMapping1);

            string t2str = Util::GetFileContent(parts[1]);
            Node* pre_tree2 = NewickLex::ParseNewickString(t2str);
            Node* tree2 = pre_tree2;
            if (pre_tree2->GetNbChildren() == 1)
            {
                tree2 = pre_tree2->GetChild(0);
            }

            unordered_map<Node*, Node*> geneSpeciesMapping2 = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(tree2, speciesTree, "___", 1);
            unordered_map<Node*, Node*> lcaMapping2 = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(tree2, speciesTree, geneSpeciesMapping2);

            vector<Node*> nads1 = GeneSpeciesTreeUtil::Instance()->GetNADNodes(tree1, speciesTree, lcaMapping1);

            vector<Node*> nads2 = GeneSpeciesTreeUtil::Instance()->GetNADNodes(tree2, speciesTree, lcaMapping2);



            //same # of leaves, but less nads = interesting
            if (geneSpeciesMapping1.size() == geneSpeciesMapping2.size() && nads1.size() < nads2.size())
            {
                cout<<parts[0]<<":"<<parts[1]<<" --> "<<nads2.size()<<" - "<<nads1.size()<<" = "<<(nads2.size() - nads1.size())<<endl;
                cout<<geneSpeciesMapping1.size()<<" "<<geneSpeciesMapping2.size()<<endl;
                treePairsOfInterest[parts[0]] = parts[1];

                string v70corrected_file = Util::ReplaceAll(parts[1], "fishtrees", "fishtrees_corrected_v2");
                /*string v70corrected_newick = Util::GetFileContent(v70corrected_file);

                Node* pre_v70corrected = NewickLex::ParseNewickString(v70corrected_newick);
                Node* v70_corrected = pre_v70corrected;
                if (v70_corrected->GetNbChildren() == 1)
                    v70_corrected = v70_corrected->GetChild(0);
                unordered_map<Node*, Node*> geneSpeciesMapping70_corrected = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(v70_corrected, speciesTree, "___", 1);
                unordered_map<Node*, Node*> lcaMapping70_corrected = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(v70_corrected, speciesTree, geneSpeciesMapping70_corrected);
                vector<Node*> nads70Corrected = GeneSpeciesTreeUtil::Instance()->GetNADNodes(v70_corrected, speciesTree, lcaMapping70_corrected);

                v70corrected_newick = NewickLex::ToNewickString(v70_corrected);
                delete pre_v70corrected;*/

                string tree1newick = NewickLex::ToNewickString(tree1);
                string tree2newick = NewickLex::ToNewickString(tree2);

                string f74 = parts[0];
                string f70 = parts[1];

                string pairFile = "/u/lafonman/Projects/nad_data/v70_vs_v74/" + Util::GetPathFilename(f70) + "__vs__" + Util::GetPathFilename(f74);

                string pairContent = "V70FILE=" + f70 + "\n" + "V70NADS=" + Util::ToString((int)nads2.size()) + "\n" + "V70TREE=" + tree2newick + "\n" +
                                      "V70CORRFILE=" + v70corrected_file + "\n" +  //"V70CORRNADS=" + Util::ToString((int)nads70Corrected.size()) + "\n" + "V70CORRTREE=" + v70corrected_newick + "\n" +
                                     "V74FILE=" + f74 + "\n" + "V74NADS=" + Util::ToString((int)nads1.size()) + "\n" + "V74TREE=" + tree1newick + "\n";

                cout<<"WRITING "<<pairFile<<endl;
                Util::WriteFileContent(pairFile, pairContent);

            }

            delete pre_tree1;
            delete pre_tree2;
        }
    }

    delete speciesTree;

    return treePairsOfInterest;
}







map<string, map<string, double> > loadDistances(string filename)
{
    map<string, map<string, double> > distances;


    std::ifstream ifs( filename );
    std::string content( (std::istreambuf_iterator<char>(ifs) ),
                         (std::istreambuf_iterator<char>()    ));


    map<int, vector<string> > splitLines;
    vector<string> lines = Util::Split(content, "\n");

    //first line is irrelevant - it has no genes
    for (int i = 1; i < lines.size(); i++)
    {
        vector<string> split = Util::Split(lines[i], " ", false);
        splitLines[i] = split;
    }


    for (int i = 1; i < lines.size(); i++)
    {
        if (lines[i] != "")
        {
            vector<string> split = splitLines[i];
            string lineGene = split[0];

            map<string, double> geneDists;
            if (distances.find(lineGene) != distances.end())
            {
                geneDists = distances[lineGene];
            }

            for (int j = 1; j < split.size(); j++)
            {
                vector<string> otherSplit = splitLines[j];
                string otherGene = otherSplit[0];

                geneDists[otherGene] = Util::ToDouble(split[j]);
            }

            distances[lineGene] = geneDists;
        }

    }




    return distances;
}



void TestPolySolverDistances_Easy()
{
    //string s = "((a,b), (c,d))";
    string s = "((((a, (c,d)), b), e), (f,g))";
    Node* speciesTree = NewickLex::ParseNewickString(s, true);

    //string g = "(a,b,c,(a,d)ad)";
    string g = "(a, c, b, e, (d,(a,b))dab);";
    Node* geneTree = NewickLex::ParseNewickString(g, true);

    unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree);
    unordered_map<Node*, Node*> lcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, geneSpeciesMapping);

    PolySolverDistance solver;
    vector<Node*> v;
    for (int i = 0; i < geneTree->GetNbChildren(); i++)
    {
        v.push_back(geneTree->GetChild(i));
    }

    map<Node*, map<Node*, double> > distances;
    vector<Node*> leaves = geneTree->GetLeafVector();
    for (int i = 0; i < leaves.size(); i++)
    {
        for (int j = i + 1; j < leaves.size(); j++)
        {
            distances[leaves[i]][leaves[j]] = 10;
            distances[leaves[j]][leaves[i]] = 10;
        }
    }



    Node* solved = solver.SolvePolytomy(v, speciesTree, lcaMapping, distances);

    cout<<NewickLex::ToNewickString(solved);

    delete solved;
    delete geneTree;
    delete speciesTree;

}


void CorrectPolytomyByDistances(Node* speciesTree, Node* geneTree, string distfilename, string outfilename)
{

    if (Util::FileExists(distfilename))
    {
        map<string, map<string, double> > rawDists = loadDistances(distfilename);


        cout<<"DISTANCES LOADED"<<endl;


        unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree, ";;", 1);
        unordered_map<Node*, Node*> lcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, geneSpeciesMapping);


        cout<<"MAPPING DONE"<<endl;

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

        cout<<"DISTANCE PARSING DONE"<<endl;

        PolySolverDistance solver;
        Node* solved = solver.SolvePolytomy(leaves, speciesTree, lcaMapping, distances);

        if (solved)
        {
            //GeneSpeciesTreeUtil::Instance()->RelabelGenes(solved, ";;", "__");
            //GeneSpeciesTreeUtil::Instance()->RelabelGenesByIndex(solved, ";;", 0);

            if (outfilename != "")
            {
                ofstream outfile;
                outfile.open (outfilename);

                //outfile<<"dlscore="<<((PSDNodeInfo*)solved->nodeInfo)->dlScore<<endl;

                outfile<<NewickLex::ToNewickString(solved);
                outfile.close();
            }
            else
            {
                cout<<endl<<NewickLex::ToNewickString(solved)<<endl;
            }

            delete solved;
        }



        /*Node* solved2 = PolySolver::Instance()->SolvePolytomies(geneTree, speciesTree,lcaMapping);

        if (solved2)
        {
            cout<<NewickLex::ToNewickString(solved2);
            delete solved2;
        }*/




        rawDists.clear();
        geneSpeciesMapping.clear();
        lcaMapping.clear();
        distances.clear();
        leaves.clear();
    }
    else
    {
        cout<<"DISTANCES TREE DOESN'T EXIST"<<endl;
    }



}




void TestPolySolverDistances()
{

    std::ifstream sifs( "/u/lafonman/Projects/lyon_data/speciestree.newick" );
    std::string spcontent( (std::istreambuf_iterator<char>(sifs) ),
                         (std::istreambuf_iterator<char>()    ));
    Node* speciesTree = NewickLex::ParseNewickString(spcontent, true);
    sifs.close();

    bool doOut = false;

    std::ifstream gifs( "/u/lafonman/Projects/lyon_data/star_trees.trees" );
    std::string gcontent( (std::istreambuf_iterator<char>(gifs) ),
                         (std::istreambuf_iterator<char>()    ));
    vector<string> glines = Util::Split(gcontent, "\n", false);
    gifs.close();

    //for (int noFamily = 1; noFamily <= glines.size(); noFamily++)
    for (int noFamily = 1; noFamily <= 1; noFamily++)
    {
        string filename = "/u/lafonman/Projects/lyon_data/resolved/resolved_famille_" + Util::ToString(noFamily) + ".trees";

        if (!doOut || !Util::FileExists(filename))
        {

            Node* geneTree = NewickLex::ParseNewickString(glines[noFamily - 1]);

            cout<<"FAMILY "<<noFamily<<"/"<<glines.size()<<endl;

            string distfilename = "/u/lafonman/Projects/lyon_data/distances/dists/famille_" + Util::ToString(noFamily) + ".dist";
            if (Util::FileExists(distfilename))
            {
                map<string, map<string, double> > rawDists = loadDistances(distfilename);


                cout<<"DISTANCES LOADED"<<endl;


                unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree, ";;", 1);
                unordered_map<Node*, Node*> lcaMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(geneTree, speciesTree, geneSpeciesMapping);


                cout<<"MAPPING DONE"<<endl;

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
                            cout<<"DISTANCE PROBLEM !"<<endl;
                        }

                    }
                }

                cout<<"DISTANCES CONVERTED"<<endl;

                PolySolverDistance solver;
                Node* solved = solver.SolvePolytomy(leaves, speciesTree, lcaMapping, distances);

                if (solved)
                {
                    GeneSpeciesTreeUtil::Instance()->RelabelGenes(solved, ";;", "__");
                    //GeneSpeciesTreeUtil::Instance()->RelabelGenesByIndex(solved, ";;", 0);


                    if (doOut)
                    {
                        ofstream outfile;
                        outfile.open (filename);
                        outfile<<"dlscore="<<((PSDNodeInfo*)solved->nodeInfo)->dlScore<<endl;

                        outfile<<NewickLex::ToNewickString(solved);
                        outfile.close();
                    }
                    else
                    {
                        cout<<endl<<NewickLex::ToNewickString(solved)<<endl;
                    }

                    delete solved;
                }



                /*Node* solved2 = PolySolver::Instance()->SolvePolytomies(geneTree, speciesTree,lcaMapping);

                if (solved2)
                {
                    cout<<NewickLex::ToNewickString(solved2);
                    delete solved2;
                }*/

                delete geneTree;


                rawDists.clear();
                geneSpeciesMapping.clear();
                lcaMapping.clear();
                distances.clear();
                leaves.clear();
            }
            else
            {
                cout<<"DISTANCES "<<noFamily<<" DOESN'T EXIST"<<endl;
            }
        }
        else
        {
            cout<<"FAMILY "<<noFamily<<" ALREADY DONE"<<endl;
        }

    }

    delete speciesTree;


    /*string stemp = "(((((((((((((((((((((((G0000238:0.254104,G0000128:0.138004):0.0729985,((((((((G0000233:0.0312839,G0000375:0.0339745):0.0183137,G0000252:0.0483493):0.0123972,G0000076:0.0575857):0.0660037,G0000301:0.114321):0.00206085,((G0000266:0.123225,(((G0000140:0.187436,G0000251:0.0896724):0.0155173,G0000329:0.116305):0.0052462,G0000376:0.0958872):0.0104403):0.0101655,(G0000248:0.0412211,G0000208:0.034207):0.0521993):0.0158816):0.00614127,((((G0000306:0.00179904,G0000015:0.00284494):0.00393032,G0000330:0.00541568):0.016964,G0000243:0):0.0764587,G0000122:0.166105):0.00769098):0.0151767,((G0000130:0.000778958,G0000358:0.00155567):0.00105457,G0000011:0.0036293):0.140156):0.0170132,G0000319:0.185534):0.00108844):0.0175918,(G0000230:0.206489,((G0000099:0.0931847,G0000062:0.112208):0.0913747,((G0000009:0.144446,(G0000250:0.150835,G0000388:0.0904893):0.0208382):0.0168665,G0000035:0.10944):0.0171593):0.0306789):0.0138467):0.0121691,(G0000036:0.0987043,G0000169:0.129969):0.0924758):0.0212983,(G0000217:0.164318,G0000203:0):0.116857):0.0156302,G0000054:0.298452):0.0136007,(((G0000172:0.169204,G0000316:0.162494):0.0263685,G0000195:0.149796):0.00937932,G0000359:0.274194):0.0673321):0.0459707,((((G0000234:0.161002,(G0000187:0.147452,(((G0000100:0.00166174,G0000027:0.000689359):0.0212046,G0000020:0.0273562):0.0175188,G0000270:0.0296382):0.119326):0.0111113):0.0731233,G0000150:0.333152):0.0157082,G0000188:0.328664):0.0109479,G0000098:0.33978):0.0378279):0.00640623,((((((((((((G0000239:0.0323665,G0000087:0.0253543):0.0335093,((G0000218:0,G0000241:5.66764e-05):0.0666888,(G0000390:0.0497128,(G0000115:0.016331,G0000013:0):0.0150839):0.00789842):0.00480799):0.00711252,((G0000379:0.0622064,((G0000074:0.0218011,G0000211:0.00935357):0.0466934,G0000364:0.00423561):0.00583292):0.00683702,G0000227:0.039692):0.00485544):0.00716675,((((G0000309:0.0376255,(((G0000092:0.0151195,G0000089:0.0798491):0.00444591,G0000277:0.00373353):0.000739372,G0000374:0.0039851):0.0418162):0.011639,(G0000348:0.0263819,(G0000258:0,G0000151:0.0124086):0.0300047):0.0132706):0.0032652,G0000284:0.0265949):0.00638542,(G0000032:0.0195105,G0000025:0.0165541):0.0444836):0.00410855):0.00658047,G0000165:0.0593693):0.00620732,((G0000189:0.0718947,G0000339:0.0529691):0.0321596,G0000213:0.104534):0.0125843):0.0156414,(G0000026:0.124463,G0000119:0.144746):0.0412624):0.0116114,(G0000181:0.110879,G0000322:0.108412):0.0021633):0.0198142,(G0000090:0.0367866,G0000385:0.03976):0.110796):0.021221,((G0000236:0.115972,((G0000345:0.0759825,G0000201:0.0713972):0.00937515,(G0000215:0.0632544,G0000377:0.0771187):0.00832306):0.0324308):0.0482615,(G0000343:0.147102,(((G0000344:0.0206868,G0000176:0.0252953):0.0237803,((G0000137:0.0116006,(G0000069:5e-07,G0000334:5e-07):0.0368811):0.0230429,G0000166:0.0318714):0.0276466):0.0409955,G0000071:0.109764):0.0546842):0.0161282):0.0152849):0.0669916,(((G0000001:0.145887,(((G0000075:0.0161956,G0000355:0.0276125):0.0538864,G0000059:0.106642):0.0141668,(G0000353:0.0108203,G0000051:0.0662981):0.073258):0.039811):0.0265179,G0000053:0.239039):0.0234587,((((((((((G0000179:0.0528768,(G0000149:0.059944,G0000012:0.033062):0.00871846):0.00513302,G0000175:0.0486453):0.0152828,G0000327:0.0656791):0.0134622,G0000275:0.0587703):0.0116866,G0000216:0.0429223):0.00948067,G0000080:0.101843):0.00698863,(G0000378:0.0840052,((((((G0000143:0.0119849,G0000242:0):0.00149461,G0000116:0.00325137):0.000823821,G0000351:0.00393437):0.00459557,G0000298:0.00315928):0.00668744,G0000207:0.0205994):0.00891754,G0000336:0.026472):0.0500824):0.0193103):0.00864064,G0000219:0.123648):0.0893406,G0000105:0.251094):0.129511,(((((G0000147:0.0957905,(G0000382:0.133427,G0000373:0.08135):0.0372777):0.0103639,(G0000022:0.0782341,G0000273:0.0560289):0.0926416):0.0182005,G0000321:0.197084):0.0209834,G0000293:0.117406):0.0227216,G0000292:0.249332):0.189295):0.0239755):0.0297606):0.0157476,(((((((((((((G0000308:0,G0000157:0.0236339):0.0635372,((G0000008:5e-07,G0000338:5e-07):0.0524399,(G0000064:0.0563956,G0000372:0.0619171):0.00594599):0.0104854):0.00568274,(((((G0000124:0.00466094,G0000107:0.00716895):0.0439466,G0000312:0.0763911):0.0163762,(G0000314:0.00202576,G0000135:0.000325341):0.0570286):0.0475607,G0000131:0.0788245):0.00546017,G0000253:0.145196):0.0129583):0.00355531,((((((G0000305:0.0333976,G0000393:0.022635):0.00982589,G0000094:0.055891):0.00700321,(G0000347:0.00497359,G0000315:0.0068563):0.0407704):0.0227938,(G0000263:0.0106923,G0000383:0.00832615):0.0765791):0.00525226,((G0000342:0.00131528,G0000366:0.00339431):0.0384992,((G0000256:0.0602861,G0000056:0.0357752):0.145839,(G0000155:0,G0000118:0.00264524):0.0109236):0.0497321):0.00879628):0.00747669,((G0000341:0.0798716,G0000186:0.0607256):0.00298241,((((G0000267:0.00275783,G0000363:0):0.0151728,(G0000224:5e-07,G0000246:5e-07):0.0109419):0.00495283,(G0000061:0.00656785,G0000280:0.0100468):0.00297472):0.00496461,G0000065:0.0364552):0.0348886):0.0220159):0.011582):0.00483815,(G0000170:0.0754326,G0000245:0.0981123):0.0391722):0.00557885,G0000024:0.0877647):0.0136428,(G0000042:0.012402,G0000057:0.0168048):0.097808):0.00461992,(((((((((((((G0000307:0.0296867,((((G0000126:0.00414267,G0000297:0.00530627):0.0014775,G0000068:0.00588924):0.0110415,G0000332:0):0.0189658,G0000389:0.0300846):0.0135213):0.0423784,G0000088:0.0988066):0.00776214,G0000354:0.0768725):0.00545261,(((G0000038:0.0856525,G0000299:0.0588753):0.0214461,G0000369:0.0938086):0.00774934,G0000399:0.0630708):0.0149017):0.00464544,G0000226:0.0536006):0.00868628,G0000395:0.107277):0.0114101,(G0000030:0.132089,((G0000102:0.0580691,G0000084:0.0658453):0.0556828,G0000259:0.018384):0.123502):0.0301456):0.00982659,(G0000380:0.0643543,G0000365:0):0.105416):0.0014685,G0000014:0.0441129):0.00612246,G0000214:0.0680296):0.00323342,(G0000396:0.0634655,G0000333:0.0944818):0.0560379):0.0135846,G0000173:0.109621):0.0151327,G0000247:0.113611):0.0295244):0.0198373,(G0000177:0.075743,G0000199:0.0773728):0.15383):0.0108924,((G0000002:0.138166,G0000291:0.100234):0.0207133,G0000136:0.108256):0.0448408):0.0116072,(G0000006:0.0473773,G0000204:0.208545):0.118956):0.00636387,(G0000368:0.154378,G0000210:0.0931481):0.0484865):0.0351512,(((G0000220:0.0221878,((G0000324:0.0187206,G0000325:0.00755593):0.0507103,(G0000350:0.0312571,G0000160:0.0324026):0.058946):0.0212395):0.0594598,G0000106:0.382974):0.0354578,(G0000193:0.0729688,G0000200:0):0.197208):0.0157199):0.0357093):0.017328):0.182013,((((((((((((G0000302:0.0688022,G0000295:0.0382618):0.0531652,G0000328:0.096841):0.0134979,(G0000047:0.0852456,G0000095:0.0704402):0.0359007):0.0221091,(((G0000043:0.0136215,(G0000180:0.00827061,G0000163:0.00777278):0.00255339):0.00272335,G0000156:0.0099659):0.0350702,(G0000361:0.0130459,G0000112:0.000684437):0.04076):0.0876201):0.0239656,((G0000268:0.00420036,G0000037:0.00953):0.00511186,G0000168:0.00523362):0.130108):0.0342027,(((G0000129:0.159461,G0000222:0.0916949):0.0710347,((G0000349:0.0519375,G0000337:0.0525155):0.0118504,G0000290:0.0899785):0.0213969):0.0421443,((G0000078:0.0263061,(G0000079:0,G0000153:0.00804381):0.0239895):0.0315189,G0000229:0):0.0414779):0.0245535):0.0104014,(G0000120:0.14324,G0000276:0.220112):0.0287547):0.0258542,G0000183:0.139017):0.0435562,((G0000340:0.141181,(G0000132:0.0875588,((G0000067:0.0266466,G0000274:0.025098):0.064247,G0000285:0.109645):0.00958785):0.056239):0.059836,(G0000261:0.20822,(G0000326:0.1561,G0000021:0.155797):0.171323):0.00878397):0.0212091):0.00871352,G0000161:0.195792):0.0619015,(G0000264:0.376526,(((G0000387:0.102143,G0000282:0.114563):0.0317743,(G0000142:0.107075,(G0000196:0.0348208,G0000117:0.0241254):0.111119):0.0380187):0.0693245,G0000190:0.30136):0.054272):0.0206824):0.0177608,((G0000007:0.0767187,G0000110:0.0462049):0.00160214,G0000318:0.0831313):0.13104):0.0163301):0.0209031,(((G0000231:0.0497182,(((G0000237:0.0516788,(G0000004:0.0190328,G0000255:0.0277478):0.0177557):0.00614502,(G0000300:0.035764,(G0000141:0.0132259,G0000279:0.0121416):0.0154246):0.0218044):0.00339947,(G0000097:0.0282539,G0000055:0.041564):0.0203941):0.00611402):0.00715904,((G0000304:0.035537,G0000357:0.00318123):0.014323,((G0000133:0.015907,G0000311:0.0118098):0.0341107,G0000202:0.0526125):0.00588487):0.0374798):0.0175516,((G0000049:0.118143,G0000391:0.0659817):0.0118227,G0000052:0.147553):0.00960379):0.0584512):0.0368121,((G0000323:0.0849497,G0000104:0.107922):0.00479717,G0000397:0.0935021):0.07795):0.00844732,(((G0000303:0.0529079,G0000096:0.0411912):0.0196889,G0000034:0.0710547):0.0454727,G0000278:0.092802):0.0186217):0.0141345,((((((((((G0000232:0.0407478,G0000185:0.0255938):0.00151421,G0000085:0.0504971):0.0120567,(((G0000254:0.047336,G0000272:0.0397181):0.00304586,G0000018:0.0376271):0.00769446,G0000063:0.0496725):0.00725901):0.000377601,(G0000269:0.0382153,G0000171:0.015128):0.0260251):0.00191631,(G0000182:0.166485,(G0000031:0.015493,G0000154:0.168883):0.0409872):0.023102):0.00592851,((((G0000265:0.0596114,G0000394:0.0822491):0.00308956,(G0000225:0.0228228,G0000194:0.0312334):0.0875524):0.00494924,(G0000070:0.0540335,((G0000144:0.0530003,G0000367:0.0514528):0.0270125,G0000103:0.0558046):0.0115033):0.00707241):0.00538808,((((G0000223:0.0074437,((G0000244:5e-07,G0000017:5e-07):0.000319523,G0000289:0.00422889):0.00721474):0.000398085,(G0000240:0.00793168,G0000271:0.00349275):0.00147732):0.00953415,G0000123:0.0196335):0.0248222,G0000050:0.0646316):0.0172963):0.0048855):0.00304317,G0000039:0.0772309):0.00892112,(((((G0000045:0.0208771,(G0000005:0.00753226,G0000072:0.00975958):0.0440358):0.0327746,(((((G0000073:0.0048342,G0000197:0):0.00691329,G0000083:0.00671658):0.0011335,G0000205:0.00226745):0.00286697,G0000384:0.0173784):0.00183075,G0000257:0.00793465):0.0133735):0.0185879,(G0000060:0.129211,G0000392:0.0451232):0.0461435):0.0245091,((G0000346:0.0801503,(G0000262:0.0558961,(G0000191:0.0315442,G0000081:0):0.0749134):0.00467256):0.00665298,(G0000082:0.0738212,G0000086:0.0954328):0.060213):0.00319987):0.0222745,((G0000174:0.0335326,G0000101:0.0229647):0.0177369,G0000287:0.0533442):0.0987275):0.0136543):0.00104954,G0000138:0.0672671):0.00618238,((((G0000121:0.0351704,G0000029:0.0326856):0.0228512,G0000228:0.0322139):0.00906625,((G0000041:0.0240071,G0000164:0.0108019):0.0698096,G0000145:0.0601663):0.00329216):0.0483143,G0000286:0.0765282):0.0193162):0.0253941):0.0120678,((G0000184:0.0922566,G0000016:0.0749396):0.0549199,G0000209:0.218096):0.0399826):0.0121898,((G0000091:0,G0000108:0.0294835):0.131897,((G0000125:0.00225812,G0000249:0.00241844):0.0819815,(G0000371:0.101304,(G0000212:0,G0000111:0.0418368):0.014519):0.0293001):0.0245049):0.0284952):0.00577495,((G0000040:0.0581406,G0000198:0.135194):0.0135651,G0000356:0.10786):0.0328128):0.00195213,((G0000046:0.0840683,(G0000023:0.0333032,G0000288:0.0260632):0.0377257):0.0156898,G0000019:0.15576):0.0106312):0.00428702,(((((((G0000127:0.00732703,G0000148:0.00812144):0.0512521,(G0000320:0.00107893,G0000296:0.00119207):0.0179016):0.0133706,G0000033:0.0179105):0.0200632,(G0000381:0.00180253,G0000206:0.000468479):0.047508):0.0509057,((G0000400:0.00510334,G0000162:0.00121225):0.0360144,G0000114:0.0441275):0.0397061):0.00666291,G0000192:0.278976):0.0137388,G0000317:0.108078):0.00936001):0.00454211,(((((G0000048:0.00373776,G0000109:0.0099926):0.000536628,(G0000386:0.000313056,G0000352:0.00195795):0.0109936):0.0250722,((G0000260:0.0324704,(G0000310:0.0411993,G0000028:0.0234001):0.016536):0.0476199,(((G0000066:0.00332612,G0000159:0.00350764):0.0010512,G0000113:0.004128):0.00709636,G0000167:0.0087801):0.0200207):0.0112931):0.0251157,G0000139:0.0571912):0.00825669,(G0000003:0.0555191,G0000158:0.046332):0.00393582):0.00376923):0.00655434,((((G0000178:0.0416317,G0000134:0.0135932):0.0511039,G0000360:0.0677763):0.00914353,(((G0000221:0,G0000370:0.00757352):0.0132516,G0000331:0.0559057):0.0477504,(G0000058:0.0506587,(G0000010:0.00756026,G0000281:0.0154653):0.0425074):0.00537587):0.0211926):0.00371898,(G0000077:0.0712864,G0000362:0.0383978):0.0429632):0.00454819):0.00346914,((G0000235:0.0215668,G0000398:0.0324895):0.0108403,(((G0000313:0.0291546,G0000294:0.0137142):0.0267902,G0000093:0.0770604):0.00142076,G0000146:0.0743064):0.00497842):0.013824):0.001636,(G0000044:0.00769062,G0000335:0.00603974):0.082168,(G0000152:0.00702507,G0000283:0.00901832):0.0352884);";
    Node* geneTemp = NewickLex::ParseNewickString(stemp);
    GeneSpeciesTreeUtil::Instance()->RelabelGenesByIndex(geneTemp, ":", 0);
    cout<<endl<<NewickLex::ToNewickString(geneTemp)<<endl;
    delete geneTemp;




    stemp = "(((((((((((((((((((((((G0000351__pan_troglodytes, G0000297__homo_sapiens), G0000289__gorilla_gorilla), G0000271__pongo_abelii), (G0000060__nomascus_leucogenys, G0000109__nomascus_leucogenys)), G0000045__macaca_mulatta), G0000005__callithrix_jacchus), ((((((G0000017__pan_troglodytes, (((G0000066__homo_sapiens, G0000073__homo_sapiens), G0000116__homo_sapiens), G0000197__homo_sapiens)), (G0000028__gorilla_gorilla, G0000065__gorilla_gorilla)), G0000048__pongo_abelii), G0000143__nomascus_leucogenys), G0000207__macaca_mulatta), G0000072__callithrix_jacchus)), G0000191__tarsius_syrichta), (((((((G0000089__pan_troglodytes, G0000224__homo_sapiens), G0000113__gorilla_gorilla), G0000298__pongo_abelii), G0000167__nomascus_leucogenys), G0000240__macaca_mulatta), G0000123__callithrix_jacchus), G0000284__tarsius_syrichta)), (G0000226__microcebus_murinus, (G0000070__otolemur_garnettii, G0000158__otolemur_garnettii))), (((((G0000101__mus_musculus, G0000312__rattus_norvegicus), G0000380__dipodomys_ordii), G0000346__ictidomys_tridecemlineatus), G0000103__cavia_porcellus), (G0000088__oryctolagus_cuniculus, G0000144__ochotona_princeps))), (((((((G0000235__ailuropoda_melanoleuca, G0000354__mustela_putorius_furo), G0000294__canis_familiaris), G0000093__felis_catus), G0000077__equus_caballus), (G0000342__myotis_lucifugus, G0000327__pteropus_vampyrus)), ((((((G0000010__tursiops_truncatus, ((((((G0000008__bos_taurus, G0000038__bos_taurus), G0000041__bos_taurus), G0000058__bos_taurus), G0000074__bos_taurus), G0000149__bos_taurus), G0000154__bos_taurus)), (G0000012__tursiops_truncatus, G0000164__bos_taurus)), (G0000157__tursiops_truncatus, G0000211__bos_taurus)), G0000085__sus_scrofa), ((G0000185__tursiops_truncatus, G0000221__bos_taurus), G0000134__sus_scrofa)), G0000179__vicugna_pacos)), ((((G0000026__erinaceus_europaeus, G0000111__erinaceus_europaeus), G0000119__erinaceus_europaeus), G0000210__erinaceus_europaeus), G0000080__sorex_araneus))), (((((G0000025__loxodonta_africana, G0000032__loxodonta_africana), G0000171__procavia_capensis), (G0000040__loxodonta_africana, G0000182__procavia_capensis)), G0000042__echinops_telfairi), (G0000031__dasypus_novemcinctus, G0000016__choloepus_hoffmanni))), ((((G0000007__macropus_eugenii, G0000105__sarcophilus_harrisii), (G0000096__macropus_eugenii, G0000110__sarcophilus_harrisii)), G0000002__monodelphis_domestica), ((G0000136__macropus_eugenii, G0000199__sarcophilus_harrisii), G0000034__monodelphis_domestica))), ((((((((((((G0000092__pan_troglodytes, G0000244__homo_sapiens), G0000126__gorilla_gorilla), G0000068__pongo_abelii), G0000223__nomascus_leucogenys), G0000267__macaca_mulatta), G0000260__callithrix_jacchus), G0000081__tarsius_syrichta), (G0000263__microcebus_murinus, G0000173__otolemur_garnettii)), (((((G0000135__mus_musculus, G0000320__rattus_norvegicus), G0000392__dipodomys_ordii), G0000395__ictidomys_tridecemlineatus), G0000253__cavia_porcellus), (G0000189__oryctolagus_cuniculus, G0000209__ochotona_princeps))), (((((((G0000286__ailuropoda_melanoleuca, G0000390__mustela_putorius_furo), G0000313__canis_familiaris), G0000094__felis_catus), G0000181__equus_caballus), (G0000366__myotis_lucifugus, G0000369__pteropus_vampyrus)), (((G0000281__tursiops_truncatus, G0000232__bos_taurus), G0000175__sus_scrofa), G0000227__vicugna_pacos)), (G0000212__erinaceus_europaeus, G0000125__sorex_araneus))), (((G0000138__loxodonta_africana, G0000198__procavia_capensis), G0000057__echinops_telfairi), (G0000046__dasypus_novemcinctus, G0000023__choloepus_hoffmanni))), ((G0000177__macropus_eugenii, G0000291__sarcophilus_harrisii), G0000104__monodelphis_domestica))), G0000343__ornithorhynchus_anatinus), (((((((G0000141__gallus_gallus, G0000255__meleagris_gallopavo), G0000106__anas_platyrhynchos), G0000237__ficedula_albicollis), G0000001__pelodiscus_sinensis), ((((G0000304__gallus_gallus, G0000279__meleagris_gallopavo), G0000166__anas_platyrhynchos), G0000344__ficedula_albicollis), G0000049__pelodiscus_sinensis)), ((((G0000311__gallus_gallus, G0000325__meleagris_gallopavo), G0000202__anas_platyrhynchos), (G0000353__ficedula_albicollis, G0000051__taeniopygia_guttata)), G0000071__pelodiscus_sinensis)), G0000052__anolis_carolinensis)), (((((((((G0000020__xenopus_tropicalis, G0000027__xenopus_tropicalis), G0000100__xenopus_tropicalis), G0000117__xenopus_tropicalis), G0000142__xenopus_tropicalis), G0000150__xenopus_tropicalis), G0000187__xenopus_tropicalis), G0000196__xenopus_tropicalis), G0000234__xenopus_tropicalis), G0000264__xenopus_tropicalis)), (((((((((((((((G0000159__pan_troglodytes, G0000246__homo_sapiens), G0000203__gorilla_gorilla), G0000257__pongo_abelii), G0000280__nomascus_leucogenys), G0000307__macaca_mulatta), G0000336__callithrix_jacchus), G0000139__tarsius_syrichta), (G0000348__microcebus_murinus, G0000262__otolemur_garnettii)), ((((G0000192__mus_musculus, G0000385__rattus_norvegicus), G0000394__dipodomys_ordii), G0000400__ictidomys_tridecemlineatus), (G0000367__oryctolagus_cuniculus, G0000339__ochotona_princeps))), (((((((G0000305__ailuropoda_melanoleuca, G0000393__mustela_putorius_furo), G0000315__canis_familiaris), G0000214__felis_catus), G0000362__equus_caballus), ((((G0000014__myotis_lucifugus, (((G0000056__pteropus_vampyrus, G0000087__pteropus_vampyrus), G0000118__pteropus_vampyrus), G0000152__pteropus_vampyrus)), (G0000044__myotis_lucifugus, G0000155__pteropus_vampyrus)), (G0000145__myotis_lucifugus, G0000239__pteropus_vampyrus)), (G0000165__myotis_lucifugus, G0000256__pteropus_vampyrus))), (((G0000299__tursiops_truncatus, G0000331__bos_taurus), G0000178__sus_scrofa), G0000360__vicugna_pacos)), (G0000219__erinaceus_europaeus, G0000151__sorex_araneus))), (((G0000170__loxodonta_africana, G0000245__procavia_capensis), G0000091__echinops_telfairi), (G0000184__dasypus_novemcinctus, G0000039__choloepus_hoffmanni))), ((G0000215__macropus_eugenii, G0000303__sarcophilus_harrisii), G0000201__monodelphis_domestica)), G0000006__ornithorhynchus_anatinus), (((((G0000324__gallus_gallus, G0000355__meleagris_gallopavo), G0000220__anas_platyrhynchos), (G0000097__ficedula_albicollis, G0000055__taeniopygia_guttata)), G0000193__pelodiscus_sinensis), G0000053__anolis_carolinensis)), G0000270__xenopus_tropicalis)), G0000188__latimeria_chalumnae), (((((((((((((G0000285__xiphophorus_maculatus, G0000266__oryzias_latipes), G0000156__gasterosteus_aculeatus), G0000295__oreochromis_niloticus), (G0000036__takifugu_rubripes, G0000122__tetraodon_nigroviridis)), G0000217__gadus_morhua), (((((G0000011__danio_rerio, G0000037__danio_rerio), G0000054__danio_rerio), G0000120__danio_rerio), G0000128__danio_rerio), G0000130__danio_rerio)), ((((((G0000293__xiphophorus_maculatus, G0000326__oryzias_latipes), G0000163__gasterosteus_aculeatus), G0000302__oreochromis_niloticus), (G0000062__takifugu_rubripes, G0000169__tetraodon_nigroviridis)), G0000229__gadus_morhua), G0000168__danio_rerio)), ((((((G0000140__xiphophorus_maculatus, G0000095__oryzias_latipes), G0000043__gasterosteus_aculeatus), G0000252__oreochromis_niloticus), (G0000129__takifugu_rubripes, G0000222__tetraodon_nigroviridis)), G0000238__gadus_morhua), G0000172__danio_rerio)), ((((((G0000250__xiphophorus_maculatus, G0000208__oryzias_latipes), G0000112__gasterosteus_aculeatus), G0000274__oreochromis_niloticus), G0000273__tetraodon_nigroviridis), G0000319__gadus_morhua), G0000230__danio_rerio)), (((((G0000316__xiphophorus_maculatus, G0000349__oryzias_latipes), G0000180__gasterosteus_aculeatus), G0000375__oreochromis_niloticus), G0000321__gadus_morhua), G0000261__danio_rerio)), ((((G0000376__xiphophorus_maculatus, G0000361__gasterosteus_aculeatus), G0000195__oreochromis_niloticus), G0000078__gadus_morhua), G0000268__danio_rerio)), ((((((G0000021__xiphophorus_maculatus, G0000132__xiphophorus_maculatus), G0000047__oryzias_latipes), (G0000015__gasterosteus_aculeatus, G0000035__gasterosteus_aculeatus)), G0000233__oreochromis_niloticus), G0000079__gadus_morhua), G0000276__danio_rerio)), (((((G0000328__xiphophorus_maculatus, G0000373__oryzias_latipes), G0000243__gasterosteus_aculeatus), G0000382__oreochromis_niloticus), G0000153__gadus_morhua), G0000292__danio_rerio))), (((((((((((((((((G0000205__pan_troglodytes, G0000277__homo_sapiens), G0000242__gorilla_gorilla), G0000083__pongo_abelii), G0000332__nomascus_leucogenys), G0000310__macaca_mulatta), G0000341__callithrix_jacchus), G0000186__tarsius_syrichta), (G0000383__microcebus_murinus, G0000309__otolemur_garnettii)), ((((((((((G0000381__mus_musculus, (((G0000033__rattus_norvegicus, G0000102__rattus_norvegicus), G0000107__rattus_norvegicus), G0000124__rattus_norvegicus)), (G0000206__mus_musculus, G0000127__rattus_norvegicus)), (G0000225__mus_musculus, G0000148__rattus_norvegicus)), (G0000259__mus_musculus, G0000174__rattus_norvegicus)), (G0000314__mus_musculus, G0000194__rattus_norvegicus)), (G0000084__mus_musculus, G0000287__rattus_norvegicus)), G0000317__dipodomys_ordii), (((G0000114__ictidomys_tridecemlineatus, G0000131__ictidomys_tridecemlineatus), G0000162__ictidomys_tridecemlineatus), G0000213__ictidomys_tridecemlineatus)), G0000019__cavia_porcellus), (G0000003__oryctolagus_cuniculus, G0000024__ochotona_princeps))), (((((((G0000322__ailuropoda_melanoleuca, G0000398__mustela_putorius_furo), G0000347__canis_familiaris), G0000275__felis_catus), G0000399__equus_caballus), (G0000254__myotis_lucifugus, G0000272__pteropus_vampyrus)), (((G0000308__tursiops_truncatus, G0000338__bos_taurus), G0000379__sus_scrofa), G0000372__vicugna_pacos)), (G0000368__erinaceus_europaeus, G0000249__sorex_araneus))), (((G0000269__loxodonta_africana, G0000333__procavia_capensis), G0000108__echinops_telfairi), G0000288__choloepus_hoffmanni)), ((G0000377__macropus_eugenii, G0000323__sarcophilus_harrisii), G0000236__monodelphis_domestica)), G0000204__ornithorhynchus_anatinus), ((((G0000334__gallus_gallus, G0000357__meleagris_gallopavo), G0000300__anas_platyrhynchos), (G0000160__ficedula_albicollis, G0000176__taeniopygia_guttata)), G0000200__pelodiscus_sinensis)), G0000282__xenopus_tropicalis), G0000098__latimeria_chalumnae), ((((((G0000251__xiphophorus_maculatus, G0000248__oryzias_latipes), G0000147__gasterosteus_aculeatus), G0000290__oreochromis_niloticus), (G0000022__takifugu_rubripes, G0000099__tetraodon_nigroviridis)), G0000183__gadus_morhua), G0000340__danio_rerio))), (((((((((((((((((G0000352__pan_troglodytes, G0000386__homo_sapiens), G0000374__gorilla_gorilla), G0000061__pongo_abelii), G0000384__nomascus_leucogenys), G0000363__macaca_mulatta), G0000389__callithrix_jacchus), G0000050__tarsius_syrichta), G0000378__otolemur_garnettii), (((((G0000090__mus_musculus, G0000296__rattus_norvegicus), G0000365__dipodomys_ordii), G0000265__ictidomys_tridecemlineatus), G0000030__cavia_porcellus), (G0000082__oryctolagus_cuniculus, G0000086__ochotona_princeps))), ((((((((G0000013__ailuropoda_melanoleuca, G0000146__mustela_putorius_furo), ((G0000121__canis_familiaris, G0000216__canis_familiaris), G0000218__canis_familiaris)), G0000029__felis_catus), (((G0000115__ailuropoda_melanoleuca, G0000228__mustela_putorius_furo), G0000241__canis_familiaris), G0000063__felis_catus)), G0000018__equus_caballus), (G0000335__myotis_lucifugus, G0000283__pteropus_vampyrus)), ((G0000364__tursiops_truncatus, G0000370__bos_taurus), G0000064__sus_scrofa)), (G0000371__erinaceus_europaeus, G0000258__sorex_araneus))), ((G0000396__loxodonta_africana, G0000356__procavia_capensis), G0000247__echinops_telfairi)), ((G0000397__macropus_eugenii, G0000345__sarcophilus_harrisii), G0000318__monodelphis_domestica)), G0000278__ornithorhynchus_anatinus), ((((((G0000004__gallus_gallus, G0000069__gallus_gallus), G0000133__meleagris_gallopavo), (G0000075__gallus_gallus, G0000137__meleagris_gallopavo)), G0000059__anas_platyrhynchos), (G0000231__ficedula_albicollis, G0000350__taeniopygia_guttata)), G0000391__pelodiscus_sinensis)), G0000387__xenopus_tropicalis), G0000190__latimeria_chalumnae), (((((((G0000329__xiphophorus_maculatus, G0000388__oryzias_latipes), G0000301__gasterosteus_aculeatus), G0000009__oreochromis_niloticus), ((G0000337__xiphophorus_maculatus, G0000306__gasterosteus_aculeatus), G0000067__oreochromis_niloticus)), ((G0000359__xiphophorus_maculatus, G0000330__gasterosteus_aculeatus), G0000076__oreochromis_niloticus)), G0000161__gadus_morhua), G0000358__danio_rerio)));";
    geneTemp = NewickLex::ParseNewickString(stemp);
    GeneSpeciesTreeUtil::Instance()->RelabelGenesByIndex(geneTemp, "__", 0);
    cout<<endl<<NewickLex::ToNewickString(geneTemp)<<endl;
    delete geneTemp;*/
}


void TestPolySolverNAD()
{
    /*string s = "((a,b), (c,d))";
    Node* speciesTree = NewickLex::ParseNewickString(s, true);

    string g = "(a,b,c,(a,d)ad)";
    Node* geneTree = NewickLex::ParseNewickString(g, true);

    unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree);

    PolySolverNAD solver;
    Node* solved = solver.SolvePolytomies(geneTree, speciesTree, geneSpeciesMapping);

    cout<<NewickLex::ToNewickString(solved);

    delete solved;
    delete geneTree;
    delete speciesTree;*/

    string s = "((((a, (c,d)), b), e), (f,g))";
    Node* speciesTree = NewickLex::ParseNewickString(s, true);

    /*string g = "(a, c, b, e, (d,(a,b))dab)";
    Node* geneTree = NewickLex::ParseNewickString(g, true);

    unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree);

    PolySolverNAD solver;
    Node* solved = solver.SolvePolytomies(geneTree, speciesTree, geneSpeciesMapping);*/

    string g = "((((((d, (a,b)), g),a),b),c), (f, (f, e)));";
    Node* geneTree = NewickLex::ParseNewickString(g, true);

    unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTree, speciesTree);


    PolySolverNAD solver;
    PolySolverCorrectionInfo info = solver.CorrectHighestNAD(geneTree, speciesTree, geneSpeciesMapping);

    Node* solved = info.correction;

    if (solved)
        cout<<NewickLex::ToNewickString(solved);

    if (solved)
        delete solved;
    delete geneTree;
    delete speciesTree;

}



void CorrectFolderNADs(string folderIn, string folderOut)
{
    string speciesNewick = "((((((Xiphophorus maculatus, Oryzias latipes), Gasterosteus aculeatus),Oreochromis niloticus), (Takifugu rubripes,Tetraodon nigroviridis)), Gadus morhua), (Astyanax mexicanus,Danio rerio))";
    Node* speciesTree = NewickLex::ParseNewickString(speciesNewick, true);

    tinydir_dir dir;
    tinydir_open(&dir, folderIn.c_str());

    int cpt = 0, cptCorrected = 0;
    int cptDidNothing = 0, cptKilledOne = 0, cptKilledMore = 0, cptCreatedMore = 0;

    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);


        if (file.is_dir)
        {
            ;
        }
        else
        {
            string filename = string(file.name);    //for some reason, calling file.name later causes crashes
            string fullname = folderIn + filename;


            if (fullname.find(".newick") != string::npos)
            {

                std::ifstream ifs( fullname );
                std::string content( (std::istreambuf_iterator<char>(ifs) ),
                                     (std::istreambuf_iterator<char>()    ));
                //cout<<fullname<<endl<<content<<endl;


                Node* geneTree = NewickLex::ParseNewickString(content);
                Node* realGeneTree = geneTree;

                //HACK ALERT !  The following is due to a bug in which the tree might have a single child root.
                //PolySolver expects its input gene trees to be binary
                if (realGeneTree->GetNbChildren() == 1)
                {
                    realGeneTree = realGeneTree->GetChild(0);
                }

                unordered_map<Node*, Node*> geneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(realGeneTree, speciesTree, "___", 1);

                PolySolverNAD solver;
                PolySolverCorrectionInfo info = solver.CorrectHighestNAD(realGeneTree, speciesTree, geneSpeciesMapping);

                Node* solved = info.correction;
                if (solved)
                {
                    Node* realSolved = solved;
                    if (solved->GetNbChildren() == 1)
                        realSolved = solved->GetChild(0);

                    //TODO : lca mappings were already computed by solver
                    //cout<<NewickLex::ToNewickString(realSolved);
                    //return;

                    unordered_map<Node*, Node*> beforeLCAMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(realGeneTree, speciesTree, geneSpeciesMapping);
                    vector<Node*> nadsBefore = GeneSpeciesTreeUtil::Instance()->GetNADNodes(realGeneTree, speciesTree, beforeLCAMapping);


                    unordered_map<Node*, Node*> afterGeneSpeciesMapping = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(realSolved, speciesTree, "___", 1);
                    unordered_map<Node*, Node*> afterLCAMapping = GeneSpeciesTreeUtil::Instance()->GetLCAMapping(realSolved, speciesTree, afterGeneSpeciesMapping);
                    vector<Node*> nadsAfter = GeneSpeciesTreeUtil::Instance()->GetNADNodes(realSolved, speciesTree, afterLCAMapping);

                    cout<<"BEFORE : "<<nadsBefore.size()<<", AFTER : "<<nadsAfter.size()<<endl;

                    if (nadsBefore.size() == nadsAfter.size())
                        cptDidNothing++;
                    if (nadsBefore.size() == nadsAfter.size() + 1)
                        cptKilledOne++;
                    if (nadsBefore.size() > nadsAfter.size() + 1)
                        cptKilledMore++;
                    if (nadsBefore.size() < nadsAfter.size())
                        cptCreatedMore++;

                    cptCorrected++;


                    string baseDir = "/u/lafonman/Projects/PolytomySolver/data/";
                    string treePairsDir = baseDir + "treepairs/";
                    string fastaDir = baseDir + "fasta/";
                    string alignementsDir = baseDir + "alignments/";
                    string conseloutDir = baseDir + "conselout/";
                    string treeID = Util::GetSubstringAfter(Util::ReplaceAll(fullname, ".newick", ""), "/");

                    //AUTester tester;
                    //tester.RunAUTest(realGeneTree, solved, treeID, "___", 0, 1, treePairsDir, fastaDir, alignementsDir, conseloutDir);


                    cout<<"WRITING "<<folderOut<<filename<<endl;
                    string outtreefile = folderOut + filename;
                    string newick = NewickLex::ToNewickString(solved);
                    string beforeNewick = NewickLex::ToNewickString(realGeneTree);


                    string txt = "ORIGINALFILE=" + fullname + "\n" +
                            "ORIGINALTREE=" + beforeNewick + "\n" +
                            "NADSBEFORE=" + Util::ToString((int)nadsBefore.size()) + "\n" +
                            "NADSAFTER=" + Util::ToString((int)nadsAfter.size()) + "\n" +
                            "POLYSIZE1=" + Util::ToString(info.firstPolySize) + "\n" +
                            "POLYSIZE2=" + Util::ToString(info.secondPolySize) + "\n" +
                            "CORRECTEDCLADE=";
                    for (int i = 0; i < info.nadCladeGenes.size();i++)
                    {
                        if (i != 0)
                            txt += ";";
                        txt += info.nadCladeGenes[i];
                    }
                    txt = txt + "\n" + "CORRECTEDTREE=" + newick;


                    Util::WriteFileContent(outtreefile, txt);

                    delete solved;
                }

                delete geneTree;

                cpt++;
            }
        }


        tinydir_next(&dir);
    }

    tinydir_close(&dir);

    cout<<"FILES="<<cpt<<endl<<"CORRECTED="<<cptCorrected<<endl<<
          "KILLED ONE="<<cptKilledOne<<endl<<"KILLED MORE="<<cptKilledMore<<endl<<
          "DID NOTHING="<<cptDidNothing<<endl<<"HARMFUL="<<cptCreatedMore<<endl;

    delete speciesTree;
}




void EvaluateRFDistFile()
{
    /*tinydir_dir dir;
    tinydir_open(&dir, "/u/lafonman/Projects/all_protein_trees_70/fishtrees_corrected_v2/");
    int cnt = 0;
    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);


        if (!file.is_dir)
        {
            string fullname = "/u/lafonman/Projects/all_protein_trees_70/fishtrees_corrected_v2/" + string(file.name);


            string poly1 = "";
            string poly2 = "";
            string corrclade = "";
            string nadsBefore = "";
            string nadsAfter = "";
            string v70corrcontent = Util::GetFileContent(fullname);
            vector<string> v70lines = Util::Split(v70corrcontent, "\n", false);

            for (int v70l = 0; v70l < v70lines.size(); v70l++)
            {
                vector<string> kv = Util::Split(v70lines[v70l], "=", false);

                string key = kv[0];
                string val = kv[1];

                if (key == "POLYSIZE1")
                    poly1 = val;
                if (key == "POLYSIZE2")
                    poly2 = val;
                if (key == "CORRECTEDCLADE")
                    corrclade = val;
                if (key == "NADSBEFORE")
                    nadsBefore = val;
                if (key == "NADSAFTER")
                    nadsAfter = val;


            }

            if (Util::ToInt(poly1) > 3 || Util::ToInt(poly2) > 2)
            {
                cout<<"P1="<<poly1<<" P2="<<poly2<<" NB="<<nadsBefore<<" NA="<<nadsAfter<<" F="<<fullname<<endl;
                cnt++;
            }

        }
        tinydir_next(&dir);
    }
    tinydir_close(&dir);

    cout<<cnt<<endl;*/



    string content = Util::GetFileContent("/u/lafonman/Projects/nad_data/rfdists_v3.txt");

    vector<string> lines = Util::Split(content, "\n", false);

    for (int i = 0; i < lines.size(); i++)
    {
        string rfdist = "";
        string v70corrfile = "";
        string nadsBefore = "";
        string nadsAfter = "";
        vector<string> params = Util::Split(lines[i], ";;", false);

        for (int p = 0; p < params.size(); p++)
        {
            vector<string> kv = Util::Split(params[p], "=");
            string key = kv[0];
            string val = kv[1];

            if (key == "RF")
                rfdist = val;
            if (key == "V70CORRFILE")
                v70corrfile = val;


        }


        string poly1 = "";
        string poly2 = "";

        if (v70corrfile != "")
        {
            string v70corrcontent = Util::GetFileContent(v70corrfile);
            vector<string> v70lines = Util::Split(v70corrcontent, "\n", false);

            for (int v70l = 0; v70l < v70lines.size(); v70l++)
            {
                vector<string> kv = Util::Split(v70lines[v70l], "=", false);

                string key = kv[0];
                string val = kv[1];

                if (key == "POLYSIZE1")
                    poly1 = val;
                if (key == "POLYSIZE2")
                    poly2 = val;
                if (key == "NADSBEFORE")
                    nadsBefore = val;
                if (key == "NADSAFTER")
                    nadsAfter = val;

            }
        }

        if (Util::ToInt(poly1) > 2 || Util::ToInt(poly2) > 2)
        {
            //if (rfdist.find("0/") == 0)
                cout<<i<<" : RF="<<rfdist<<" P1="<<poly1<<" P2="<<poly2<<" NB="<<nadsBefore<<" NA="<<nadsAfter<<endl;
        }
    }
}


int main_test(int argc, char *argv[])
{
    //TestPolySolverNAD();
    //CorrectFolderNADs("/u/lafonman/Projects/all_protein_trees/fishtrees/");
    //TestPolySolverDistances();

    //FindCorrespondingTrees("/u/lafonman/Projects/all_protein_trees/fishtrees_nomexicanus/", "/u/lafonman/Projects/all_protein_trees_70/fishtrees/");

    //FindTreesInWhichEnsemblCorrectedNADs();
    //CorrectFolderNADs("/u/lafonman/Projects/all_protein_trees_70/fishtrees/", "/u/lafonman/Projects/all_protein_trees_70/fishtrees_corrected_v2/");
    EvaluateRFDistFile();


    return 0;







    //THIS PART BELOW TESTS POLYTOMYSOLVER

    string s = "((a,b), (c,d))";
    Node* speciesTree = NewickLex::ParseNewickString(s, true);




    string g;
    Node* geneTree;
    unordered_map<Node*, Node*> geneSpeciesMapping;


    //TEST 1 : solve a Newick gene tree
    //we map all extant genes to the species with the same label

    g = "(((a, a, b, c)), b, (a, d, d))";
    geneTree = NewickLex::ParseNewickString(g);
    geneSpeciesMapping = PolySolver::Instance()->GetGeneSpeciesMappingByPrefix(geneTree, speciesTree);


    Node* res = PolySolver::Instance()->SolvePolytomies(geneTree, speciesTree, geneSpeciesMapping);
    string rstr = NewickLex::ToNewickString(res);
    cout<<rstr<<endl;
    delete res;


    delete geneTree;



    //TEST 2 : solve a manually created gene tree

    //creates the (((a, a, b, c)), b, (a, d, d)) tree
    geneTree = new Node(false);

    Node* g1 = geneTree->AddChild();
    Node* g1child = g1->AddChild();
    g1child->SetLabel("a");
    g1child = g1->AddChild();
    g1child->SetLabel("a");
    g1child = g1->AddChild();
    g1child->SetLabel("b");
    g1child = g1->AddChild();
    g1child->SetLabel("c");

    Node* g2 = geneTree->AddChild();
    g2->SetLabel("b");

    Node* g3 = geneTree->AddChild();
    Node* g3child = g3->AddChild();
    g3child->SetLabel("a");
    g3child = g3->AddChild();
    g3child->SetLabel("d");
    g3child = g3->AddChild();
    g3child->SetLabel("d");

    geneSpeciesMapping.clear();
    geneSpeciesMapping = PolySolver::Instance()->GetGeneSpeciesMappingByPrefix(geneTree, speciesTree);

    Node* res2 = PolySolver::Instance()->SolvePolytomies(geneTree, speciesTree, geneSpeciesMapping);
    string rstr2 = NewickLex::ToNewickString(res2);
    cout<<rstr2<<endl;
    delete res2;

    delete geneTree;



    //TEST 3 : another test
    g = "(a, a, a, a, b, (b, b, c, c, d, d, (c, d, (a, a, b, b, c, c, d, d) )), c, (a, b, (c, c, c, d, a, b, a), (c, d, (a, b)) ) )";
    geneTree = NewickLex::ParseNewickString(g);

    geneSpeciesMapping.clear();
    geneSpeciesMapping = PolySolver::Instance()->GetGeneSpeciesMappingByPrefix(geneTree, speciesTree);

    //as of Sept 2013, result should give
    //((((((((a, b), (c, d)), ((a, b), (c, d))), (c, d)), (b, (c, d))), (b, (c, d))), (((((a, a), b), (((c, c), c), d)), ((a, b), (c, d))), (a, b))), (((((a, a), a), a), b), c));

    //or dup/loss equivalent
    Node* res3 = PolySolver::Instance()->SolvePolytomies(geneTree, speciesTree, geneSpeciesMapping);
    string rstr3 = NewickLex::ToNewickString(res3);
    cout<<rstr3<<endl;
    delete res3;

    delete geneTree;





    //TEST 4 : an example from Eric Tannier
    g = "(((((((YEPES2_4_PE2050:1.0E-8,YERPA_1_PE1481:1.0E-8)r90,(YEPES3_1_PE2171:1.0E-8,YERPN_3_PE1591:1.0E-8)r91)r89,YERPG_2_PE2164:1.0E-8)r88,YERPP_3_PE981:1.0E-8)r86,YERPY_1_PE3056:0.11933584)1.0,((YERPG_1_PE79:1.0E-8,YEPES1_4_PE10:1.0E-8)r81,YERPP_1_PE9:1.0E-8)r80)0.905,((YEPES2_3_PE10:1.0E-8,YERPA_2_PE9:1.0E-8)r71,YERPN_1_PE9:1.0E-8)r70)n24;";
    geneTree = NewickLex::ParseNewickString(g);


    string stemp;
    PolySolver::Instance()->RestrictGeneTreeByBranchSupport(geneTree, 2.8);

    //output the restricted version
    stemp = NewickLex::ToNewickString(geneTree, false);
    cout<<stemp<<endl;

    delete speciesTree;
    s = "((((((YERPA:2,YEPES2:2):1,(YERPN:1,YEPES3:1):2):1,YERPG:4):1,YEPES1:5):1,YERPP:6):4,(YERPY:9,(YERP3:8,(YERPB:7,YEPSE1:7):1):1):1);";
    speciesTree = NewickLex::ParseNewickString(s, true);

    geneSpeciesMapping.clear();
    geneSpeciesMapping = PolySolver::Instance()->GetGeneSpeciesMappingByPrefix(geneTree, speciesTree);

    //as of Sept 2013, result should give

    //or dup/loss equivalent
    Node* res4 = PolySolver::Instance()->SolvePolytomies(geneTree, speciesTree, geneSpeciesMapping);
    string rstr4 = NewickLex::ToNewickString(res4);
    cout<<rstr4<<endl;
    delete res4;

    delete geneTree;





    delete speciesTree;



    return -1;
}

#endif
