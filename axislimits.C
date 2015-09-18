#include <iostream>
#include <assert.h>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "SetBranchAddress.C"

Double_t Z1Masscut = -1;
Double_t Z2Masscut = -1;
Double_t ZZMassmin = 100;
Double_t ZZMassmax = 140;

Double_t pi = TMath::Pi();
Int_t maxevents = -1;
Double_t discriminantconstant = 1;
TF1 *cfunction = 0;
Double_t *cfunctionparameters = 0;

enum Statistic {Minimum, Maximum, Average, RMS};

using namespace std;

Int_t binsScatterPlotx = 1000;
Int_t binsScatterPloty = 1000;
Int_t binsHistogram = 100;
Int_t runNumberBins = 30;
Int_t binsProfileResolution = 30;    //for everything but runNumber and nHits
                                     //(nHits gets a bin for each integer between the minimum and the maximum)

Double_t findStatistic(Statistic what,Int_t nFiles,TString *files,TString var)
{
    cout << "Don't go here" << endl;
    assert(0);
    Float_t x = 0, x1 = 0, x2 = 0, x3 = 0, wt = 1;
    x1 = x1; x2 = x2; x3 = x3;
    Float_t ZZMass = 0;
    Float_t Z1Mass = 999;
    Float_t Z2Mass = 999;

    Double_t totalweight = 0, totallength = 0;
    Double_t result = 0;
    if (what == Minimum) result = 1e100;
    if (what == Maximum) result = -1e100;

    Double_t average = 0;
    if (what == RMS)
        average = findStatistic(Average,nFiles,files,var);

    for (Int_t j = 0; j < nFiles; j++)
    {
        branches.clear();
        variables.clear();
        if (files[j] == "") continue;
        TFile *f = TFile::Open(files[j]);
        TTree *tree = (TTree*)f->Get("SelectedTree");
        if (tree == 0)
            tree = (TTree*)f->Get("GenTree");
        Int_t length = TMath::Min(tree->GetEntries(),maxevents > 0 ? (Long64_t)maxevents : tree->GetEntries());

        if (var == "discriminant")
        {
            SetBranchAddress(tree,"p0plus_VAJHU_NEW","p0plus_VAJHU",&x1);
            SetBranchAddress(tree,"bkg_VAMCFM_NEW","bkg_VAMCFM",&x2);
        }
        else
            SetBranchAddress(tree,var,TString(var).Prepend("Gen"),&x);
        SetBranchAddress(tree,"mH","ZZMass","GenZZMass",&ZZMass);
        if (Z1Masscut > 0)
            SetBranchAddress(tree,"mZ1","Z1Mass","GenZ1Mass",&Z1Mass);
        if (Z2Masscut > 0)
            SetBranchAddress(tree,"mZ2","Z2Mass","GenZ2Mass",&Z2Mass);

        if (var != "wt" && (what == Average || what == RMS))
            SetBranchAddress(tree,"wt","MC_weight",&wt);

        for (Int_t i = 0; i<length; i++)
        {
            GetEntry(tree,i);
            if (var == "discriminant")
                x = x1 / (x1 + discriminantconstant * x2);
            if (ZZMass < ZZMassmin || ZZMass > ZZMassmax || Z1Mass < Z1Masscut || Z2Mass < Z2Masscut) continue;

            totalweight += wt;
            if (wt > 0)
                totallength++;
            if (wt < 0)
                totallength--;

            if (what == Minimum && x < result)
                result = x;
            if (what == Maximum && x > result)
                result = x;
            if (what == Average)
                result += x * wt;
            if (what == RMS)
                result += (x - average) * (x - average) * wt * wt;
        }
        delete f;         //automatically closes the file
    }
    if (what == Average) result /= totalweight;
    if (what == RMS)  result = sqrt(result / (((totallength - 1) / totallength) * totalweight));  
                   //http://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
    return result;
}

Double_t findAverage(Int_t nFiles,TString *files,TString var)
{
    return findStatistic(Average,nFiles,files,var);
}

Double_t findMin(Int_t nFiles,TString *files,TString var)
{
    return findStatistic(Minimum,nFiles,files,var);
}

Double_t findMax(Int_t nFiles,TString *files,TString var)
{
    return findStatistic(Maximum,nFiles,files,var);
}

Double_t findRMS(Int_t nFiles,TString *files,TString var)
{
    return findStatistic(RMS,nFiles,files,var);
}


//These functions are for 1 file

Double_t findStatistic(Statistic what,TString file,TString var)
{
    return findStatistic(what,1,&file,var);
}

Double_t findAverage(TString file,TString var)
{
    return findStatistic(Average,file,var);
}

Double_t findMin(TString file,TString var)
{
    return findStatistic(Minimum,file,var);
}

Double_t findMax(TString file,TString var)
{
    return findStatistic(Maximum,file,var);
}

Double_t findRMS(TString file,TString var)
{
    return findStatistic(RMS,file,var);
}




//This puts the axis limits that should be used for jhuGenPlot in min and max.
//Default axis limits are defined for pt, qoverpt, dxy, dz, theta, eta, and phi.
//For run number and nHits, the minimum and maximum are used.
//For any other variable, average +/- 5*rms are used.
//To use this instead of the default values, just comment out the part that says [else] if (var == "?") {min = ?; max = ?;}

void axislimits(Int_t nFiles,TString *files,TString var,Double_t &min,Double_t &max)
{
    if (var.Contains("phi") || var.Contains("Phi"))
    {
        max = pi;
        min = -pi;
    }
    else if (var.Contains("cos"))
    {
        max = 1;
        min = -1;
    }
    else if (var == "dEta_VBF")
    {
        max = 12;
        min = -12;
    }
    else if (var == "GenDRjet")
    {
        max = 12;
        min = 0;
    }
    else if (var == "mJJ_VBF" || var == "GenDijetMass")
    {
        max = 4000;
        min = 0;
    }
    else if (var == "pTH" || var == "GenHPt")
    {
        max = 1200;
        min = 0;
    }
    else if (var == "mH" || var == "GenHMass")
    {
        max = ZZMassmax;
        min = ZZMassmin;
    }
    else if (var == "mZ1" || var == "GenZ1Mass")
    {
        max = /*92*/200;
        min = TMath::Max(4.,Z1Masscut);
    }
    else if (var == "mZ2" || var == "GenZ2Mass")
    {
        max = /*40*/200;
        min = TMath::Max(4.,Z2Masscut);
    }
    else if (var == "discriminant")
    {
        min = 0;
        max = 1.1;
    }
    else
    {
        cout << "No axis limits for " << var << ".  Using average +/- 5 * rms." << endl;
        Double_t average = findAverage(nFiles,files,var);
        Double_t rms = findRMS (nFiles,files,var);
        max = TMath::Min(average + 5 * rms,findMax(nFiles,files,var));
        min = TMath::Max(average - 5 * rms,findMin(nFiles,files,var));
    }
}
