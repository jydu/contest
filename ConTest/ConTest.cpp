//
// File: ConTest.cpp
// Created by: Julien Dutheil
// Created on: Fri Sept 8 17:42 2006
//
//
//Copyright or © or Copr. Julien Dutheil
//
//Julien.Dutheil@univ-montp2.fr
//
//This software is a computer program whose purpose is to detect positions
//within a set of aligned sequence that are evolutionarily constrained.
//
//This software is governed by the CeCILL license under French law and
//abiding by the rules of distribution of free software.  You can  use, 
//modify and/ or redistribute the software under the terms of the CeCILL
//license as circulated by CEA, CNRS and INRIA at the following URL
//"http://www.cecill.info". 
//
//As a counterpart to the access to the source code and  rights to copy,
//modify and redistribute granted by the license, users are provided only
//with a limited warranty  and the software's author,  the holder of the
//economic rights,  and the successive licensors  have only  limited
//liability. 
//
//In this respect, the user's attention is drawn to the risks associated
//with loading,  using,  modifying and/or developing or reproducing the
//software by the user in light of its specific status of free software,
//that may mean  that it is complicated to manipulate,  and  that  also
//therefore means  that it is reserved for developers  and  experienced
//professionals having in-depth computer knowledge. Users are therefore
//encouraged to load and test the software's suitability as regards their
//requirements in conditions enabling the security of their systems and/or 
//data to be ensured and,  more generally, to use and operate it in the 
//same conditions as regards security. 
//
//The fact that you are presently reading this means that you have had
//knowledge of the CeCILL license and that you accept its terms.


// From the STL:

#include <iostream>
using namespace std;
#include <cstdlib>
#include <stdexcept>

// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/AlphabetIndex/SimpleIndexDistance.h>

// From bpp-phyl:
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/HomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Simulation/MutationProcess.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Mapping/SubstitutionRegister.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "contest parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options,").endLine();
  (*ApplicationTools::message << "  and the example option file in this package.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

/******************************************************************************/

vector<double> computeNorms(const ProbabilisticSubstitutionMapping& mapping)
{
  size_t nbVectors = mapping.getNumberOfSites();
  vector<double> vect(nbVectors);
  for (size_t i = 0; i < nbVectors; i++)
    vect[i] = SubstitutionMappingTools::computeNormForSite(mapping, i);
  return vect;
}

/******************************************************************************/

void doCalculations(ofstream& out,
    DRHomogeneousTreeLikelihood* tl,
    DRTreeParsimonyScore* ts,
    const SiteContainer* sites,
    vector<string>& properties,
    vector<AlphabetIndex2*>& weights,
    bool outputSites)
{
  TotalSubstitutionRegister* reg = new TotalSubstitutionRegister(tl->getSubstitutionModel()); 
  WeightedSubstitutionCount* subCount = new UniformizationSubstitutionCount(tl->getSubstitutionModel(), reg); 
  vector<unsigned int> scores;
  if (ts)
    scores = ts->getScoreForEachSite();

  //Here it starts:

  //Compute posterior rates:
  vector<double> pr = tl->getPosteriorRateOfEachSite();
  //Compute vectors:
  auto_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(*tl, *subCount, false));
  vector<double> normsUnweighted = computeNorms(*mapping);

  vector< vector<double> > norms;
  for (size_t i = 0; i < properties.size(); i++)
  {
    subCount->setWeights(weights[i], false);
    mapping.reset(SubstitutionMappingTools::computeSubstitutionVectors(*tl, *subCount, false));
    norms.push_back(computeNorms(*mapping));
  }
   
  size_t nbSites = sites->getNumberOfSites();
  size_t nbSeqs = sites->getNumberOfSequences();
  for (size_t j = 0; j < nbSites; j++)
  {
    //For each site:
    const Site* site = &sites->getSite(j);
    if (outputSites) out << "[" << site->getPosition() << "]\t";
    if (ts) out << scores[j] << "\t";
    out << pr[j] << "\t" << normsUnweighted[j];
    for (size_t i = 0; i < properties.size(); i++)
    {
      //For each property...

      //1) get the norm of the weighted vector:
      out << "\t" << norms[i][j];

      //This only works for index1 ppt...
      SimpleIndexDistance* sid = dynamic_cast<SimpleIndexDistance*>(weights[i]);
      if (sid) {
        //2) compute observed statistics:
        double min = -log(0.), max = log(0.), sum = 0;
        double n = 0;
        for (size_t k = 0; k < nbSeqs; k++)
        {
          try
          {
            double value = sid->getAlphabetIndex1().getIndex(site->getValue(k));
            if (value < min) min = value;
            if (value > max) max = value;
            sum += value;
            n++;
          }
          catch (BadIntException& ex) {} //Generic character ignored.
        }
        //3) write results:
        out << "\t" << min << "\t" << max << "\t" << (sum / n);
      } else {
        //3) write results:
        out << "\tNA\tNA\tNA";
      }
    }
    out << endl;
  }
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  cout << endl;
  cout << endl;
  cout << "***********************************************************" << endl;
  cout << "* This is ConTest      version 1.2.0     date: 03/10/2017 *" << endl;
  cout << "*           A program to detect constrained sites.        *" << endl;
  cout << "***********************************************************" << endl;
  cout << endl;

  try
  {
  
  // **************************
  // * Retrieving parameters: *
  // **************************
  
  if (argc == 1)
  { // No argument, show display some help and leave.
    help();
    exit(-1);
  }

  BppApplication contest(argc, argv, "ConTest");
  contest.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(contest.getParams(), "", false);

  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, contest.getParams());
  
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, contest.getParams());
  delete allSites;

  ApplicationTools::displayResult("# of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("# of sites", TextTools::toString(sites->getNumberOfSites()));
  
  Tree* tmpTree = PhylogeneticsApplicationTools::getTree(contest.getParams());
  TreeTemplate<Node>* tree = new TreeTemplate<Node>(*tmpTree);
  delete tmpTree;
  ApplicationTools::displayResult("# number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  ApplicationTools::displayResult("# number of sons at root", TextTools::toString(tree->getRootNode()->getNumberOfSons()));
  
  SubstitutionModel* model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, 0, sites, contest.getParams());
  
  DiscreteDistribution* rDist = PhylogeneticsApplicationTools::getRateDistribution(contest.getParams());

  HomogeneousSequenceSimulator * simulator = new HomogeneousSequenceSimulator(model, rDist, tree);
  bool continuousSim = ApplicationTools::getBooleanParameter("simulations.continuous", contest.getParams(), false, "", true, false);
  ApplicationTools::displayResult("Rate distribution for simulations", (continuousSim ? "continuous" : "discrete"));
  simulator->enableContinuousRates(continuousSim);
  
  //Reading properties:
  unsigned int nbProperties = ApplicationTools::getParameter<unsigned int>("properties.number", contest.getParams(), 1, "", true, false);
  vector<string> propertiesName;
  vector<string> propertiesDesc;
  for (unsigned int i = 0; i < nbProperties; i++)
  {
    string name = ApplicationTools::getStringParameter("properties.name", contest.getParams(), "None", TextTools::toString(i+1), true, true);
    string desc = ApplicationTools::getStringParameter("properties.type", contest.getParams(), "None", TextTools::toString(i+1), true, true);
    if (name == "None")
      throw Exception("Error! A valid name should be provided for property " + TextTools::toString(i));
    if (name == "None")
      throw Exception("Error! A valid description should be provided for property " + TextTools::toString(i));
    propertiesName.push_back(name);
    propertiesDesc.push_back(desc);
  }
  vector<AlphabetIndex2*> weights;
  for (size_t i = 0; i < propertiesDesc.size(); i++)
  {
    AlphabetIndex2* weight = SequenceApplicationTools::getAlphabetIndex2(alphabet, propertiesDesc[i], "", false);
    ApplicationTools::displayResult("Loading property:", propertiesName[i]);
    weights.push_back(weight);
  }
   
  // Analysing real data:
  string path = ApplicationTools::getAFilePath("output.file", contest.getParams(), true, false);
  ofstream out(path.c_str(), ios::out);
  out << "Group\tParc\tpr\tnorm";
  //For each property, compute the corresponding norm,
  //and range and mean for observed indices.
  for (size_t i = 0; i < propertiesName.size(); i++)
  {
    string p = propertiesName[i];
    out << "\tnorm." << p << "\tmin." << p << "\tmax." << p << "\tmean." << p;
  }
  out << endl;
    
  DRHomogeneousTreeLikelihood* tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, false, false);
  tl->initialize();
  DRTreeParsimonyScore* ts = new DRTreeParsimonyScore(*tree, *sites, false);
  ApplicationTools::displayTask("Computing observed rates, norms and stats");
  doCalculations(out, tl, ts, sites, propertiesName, weights, true);
  out.close();
  delete tl;
  delete ts;
  ApplicationTools::displayTaskDone();
  
   
  // Performs simulations:
  
  ApplicationTools::displayTask("Computing null distribution", true);
  string simpath = ApplicationTools::getAFilePath("null.output.file", contest.getParams(), true, false);
  ofstream simout(simpath.c_str(), ios::out);
  simout << "pr\tnorm";
  //For each property, compute the corresponding norm,
  //and range and mean for observed indices.
  for (size_t i = 0; i < propertiesName.size(); i++)
  {
    string p = propertiesName[i]; 
    simout << "\tnorm." << p << "\tmin." << p << "\tmax." << p << "\tmean." << p;
  }
  simout << endl;
  
  size_t nbRepCPU = ApplicationTools::getParameter<size_t>("null.nb_rep_CPU", contest.getParams(), 10);
  size_t nbRepRAM = ApplicationTools::getParameter<size_t>("null.nb_rep_RAM", contest.getParams(), 100);

  OutputStream* tmp = ApplicationTools::warning;
  ApplicationTools::warning = 0;
  for (size_t i = 0; i < nbRepCPU; i++)
  {
    //Generate data set:
    auto_ptr<SiteContainer> simSites(simulator->simulate(nbRepRAM));
    auto_ptr<DRHomogeneousTreeLikelihood> simTl(new DRHomogeneousTreeLikelihood(*tree, *simSites, model, rDist, false, false));
    simTl->initialize();
    doCalculations(simout, simTl.get(), 0, simSites.get(), propertiesName, weights, false);
    ApplicationTools::displayGauge(i, nbRepCPU-1, '=');
  }
  ApplicationTools::warning = tmp;

  simout.close();
  ApplicationTools::displayTaskDone();
  
  delete alphabet;
  delete tree;
  delete model;
  delete rDist;
  delete simulator;

  contest.done();
  
  } 
  catch(Exception e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

