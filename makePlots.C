#include "misalignmentDependence.C"
#include "TString.h"
#include "TSystem.h"

const Int_t xsize = 1;
const Int_t ysize = 18;

TString xvariables[xsize] = {""};
TString yvariables[ysize] = {"mH", "mZ1", "mZ2", "costhetastar_VBF", "costheta1_VBF", "costheta2_VBF", "Phi_VBF", "Phi1_VBF",
                                                 "costhetastar_ZZ4l", "costheta1_ZZ4l", "costheta2_ZZ4l", "Phi_ZZ4l", "Phi1_ZZ4l",
                             "mJJ_VBF", "dPhi_VBF", "dEta_VBF", "sqrtminusq2", "pTH"
                            };

//********************************************************************
//general functions, to allow any choice of variables
//it makes plots for each pair of variables if matrix[x][y] is true
//the order of the variables in the matrix is the same as shown above,
// in xvariables and yvariables.
//********************************************************************

void makePlots(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString directory,
               Bool_t matrix[xsize][ysize])
{
    stufftodelete->SetOwner(true);
/*
    for (Int_t i = 0, totaltime = 0; i < nFiles; i++)
    {
        TFile *f = 0;
        bool exists = false;
        if (files[i] == "") exists = true;

        for (int j = 1; j <= 60*24 && !exists; j++, totaltime++)  //wait up to 1 day for the validation to be finished
        {
            f = TFile::Open(files[i]);
            if (f != 0)
                exists = f->IsOpen();
            delete f;
            if (exists) continue;
            gSystem->Sleep(60000);
            cout << "It's been ";
            if (j >= 60)
                cout << j/60 << " hour";
            if (j >= 120)
                cout << "s";
            if (j % 60 != 0 && j >= 60)
                cout << " and ";
            if (j % 60 != 0)
                cout << j%60 << " minute";
            if (j % 60 >= 2)
                cout << "s";
            cout << endl;
        }
        if (!exists) return;
        if (i == nFiles - 1 && totaltime > nFiles)
            gSystem->Sleep(60000);
    }
*/
    TString directorytomake = directory;
    gSystem->mkdir(directorytomake,true);
    if (misalignment != "")
    {
        directorytomake.Append("/fits");
        gSystem->mkdir(directorytomake);
    }

    for (Int_t x = 0; x < xsize; x++)
    {
        for (Int_t y = 0; y < ysize; y++)
        {
            if (false) continue;        //this line is to make it easier to do e.g. all plots involving Delta eta
                                        //(replace false with yvariables[y] != "eta")

            if (!matrix[x][y]) continue;

            if (xvariables[x] == "" && yvariables[y] == "") continue;

            Int_t nPlots = nFiles+4;                     //scatterplot for each (if you uncomment it), profile, resolution, and fits for each.
            vector<TString> s;

            TString slashstring = "";
            if (directory.Last('/') != directory.Length() - 1) slashstring = "/";

            vector<TString> plotnames;
            for (Int_t i = 0; i < nFiles; i++)
            {
                plotnames.push_back(names[i]);   //this is plotnames[i]
                plotnames[i].ReplaceAll(" ","");
            }

            plotnames.push_back("");             //this is plotnames[nFiles], but gets changed
            if (xvariables[x] == "")
                plotnames[nFiles] = "hist";
            else
                plotnames[nFiles] = "profile";

            plotnames.push_back("resolution");   //this is plotnames[nFiles+1]

            plotnames.push_back("");             //this is plotnames[nFiles+2]
            plotnames.push_back("");             //this is plotnames[nFiles+3]

            if (plotnames[nFiles] == "profile")
            {
                plotnames[nFiles+2] = ".profile";
                plotnames[nFiles+2].Prepend(misalignment);
                plotnames[nFiles+3] = ".resolution";
                plotnames[nFiles+3].Prepend(misalignment);
                plotnames[nFiles+2].Prepend("fits/");
                plotnames[nFiles+3].Prepend("fits/");
            }
            else if (fitHist(misalignment))
            {
               plotnames[nFiles+2] = ".hist";
               plotnames[nFiles+2].Prepend(misalignment);
               plotnames[nFiles+2].Prepend("fits/");
            }
            else
            {
                plotnames[nFiles+2] = "profile.";
                plotnames[nFiles+2].Append(misalignment);
                plotnames[nFiles+3] = "resolution.";
                plotnames[nFiles+3].Append(misalignment);
            }

            TString xvarstring = xvariables[x];
            if (xvariables[x] != "" && yvariables[y] != "") xvarstring.Append(".");

            TString yvarstring = yvariables[y];

            for (Int_t i = 0; i < nPlots; i++)
            {
                stringstream ss;
                ss << directory << slashstring << plotnames[i] << "."
                   << xvarstring << yvarstring << ".pngepsroot";
                s.push_back(ss.str());
            }

            Int_t i;
            for (i = 0; i < nFiles; i++)
            {
                if (xvariables[x] == "" || yvariables[y] == "") continue;
                //uncomment this section to make scatterplots
                /*
                jhuGenPlot(files[i],xvariables[x],yvariables[y],false,false,s[i]);
                stufftodelete->Clear();
                for ( ; gROOT->GetListOfCanvases()->GetEntries() > 0; )
                    deleteCanvas( gROOT->GetListOfCanvases()->Last());
                for ( ; gROOT->GetListOfFiles()->GetEntries() > 0; )
                    delete (TFile*)gROOT->GetListOfFiles()->Last();
                */
            }

            if (xvariables[x] != "" && yvariables[y] != "")
            {
                //make profile
                TCanvas *c1 = jhuGenPlot(nFiles,files,names,xvariables[x],yvariables[y],false,s[i]);
                if (misalignmentDependence(c1,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                           true,false,s[i+2]))
                {
                    s[i+2].ReplaceAll(".png",".parameter.png");
                    misalignmentDependence(c1,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                               false,false,s[i+2]);
                }
                stufftodelete->Clear();
                for ( ; gROOT->GetListOfCanvases()->GetEntries() > 0; )
                    deleteCanvas( gROOT->GetListOfCanvases()->Last());
                for ( ; gROOT->GetListOfFiles()->GetEntries() > 0; )
                    delete (TFile*)gROOT->GetListOfFiles()->Last();

                //make resolution plot
                TCanvas *c2 = jhuGenPlot(nFiles,files,names,xvariables[x],yvariables[y],true ,s[i+1]);
                if (misalignmentDependence(c2,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                           true,true,s[i+3]))
                {
                    s[i+3].ReplaceAll(".png",".parameter.png");
                    misalignmentDependence(c2,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                               false,true,s[i+3]);
                }
                stufftodelete->Clear();
                for ( ; gROOT->GetListOfCanvases()->GetEntries() > 0; )
                    deleteCanvas( gROOT->GetListOfCanvases()->Last());
                for ( ; gROOT->GetListOfFiles()->GetEntries() > 0; )
                    delete (TFile*)gROOT->GetListOfFiles()->Last();
            }
            else
            {
                //make histogram
                TCanvas *c1 = jhuGenPlot(nFiles,files,names,xvariables[x],yvariables[y],false,s[i]);
                if (misalignmentDependence(c1,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                           true,false,s[i+2]))
                {
                    if (fitHist(misalignment))
                    {
                        s[i+2].ReplaceAll(".png",".parameter.png");
                        misalignmentDependence(c1,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                               false,false,s[i+2]);
                    }
                    else
                        misalignmentDependence(c1,nFiles,names,misalignment,values,phases,xvariables[x],yvariables[y],
                                               true,true,s[i+3]);
                }
                stufftodelete->Clear();
                for ( ; gROOT->GetListOfCanvases()->GetEntries() > 0; )
                    deleteCanvas( gROOT->GetListOfCanvases()->Last());
                for ( ; gROOT->GetListOfFiles()->GetEntries() > 0; )
                    delete (TFile*)gROOT->GetListOfFiles()->Last();
            }
            cout << y + ysize * x + 1 << "/" << xsize*ysize << endl;
        }
    }
}

void makePlots(Int_t nFiles,TString *files,TString *names,TString directory, Bool_t matrix[xsize][ysize])
{
    makePlots(nFiles,files,names,"",(Double_t*)0,(Double_t*)0,directory,
              matrix);
}

void makePlots(TString file,TString misalignment,Double_t *values,Double_t *phases,TString directory,Bool_t matrix[xsize][ysize])
{
    int n = file.CountChar(',') + 1;
    TString *files = new TString[n];
    TString *names = new TString[n];
    setTDRStyle();
    vector<Color_t> tempcolors = colors;
    vector<Style_t> tempstyles = styles;
    for (int i = 0; i < n; i++)
    {
        TString thisfile = nPart(i+1,file,",");
        int numberofpipes = thisfile.CountChar('|');
        if (numberofpipes >= 0 && nPart(numberofpipes+1,thisfile,"|").IsDigit())
        {
            if (numberofpipes >= 1 && nPart(numberofpipes,thisfile,"|").IsDigit())
            {
                colors[i] = nPart(numberofpipes,thisfile,"|").Atoi();
                styles[i] = nPart(numberofpipes+1,thisfile,"|").Atoi();
                thisfile.Remove(thisfile.Length() - nPart(numberofpipes,thisfile,"|").Length() - nPart(numberofpipes+1,thisfile,"|").Length() - 2);
            }
            else
            {
                colors[i] = nPart(numberofpipes + 1,thisfile,"|").Atoi();
                thisfile.Remove(thisfile.Length() - nPart(numberofpipes+1,thisfile,"|").Length() - 2);
            }
        }
        files[i] = nPart(1,thisfile,"=",true);
        names[i] = nPart(2,thisfile,"=",false);
    }
    if (n == 1 && names[0] == "")
        names[0] = "scatterplot";     //With 1 file there's no legend, so this is only used in the filename of the scatterplots, if made
    makePlots(n,files,names,misalignment,values,phases,directory,matrix);
    delete[] files;
    delete[] names;
    colors = tempcolors;
    styles = tempstyles;
}

void makePlots(TString file,TString directory,Bool_t matrix[xsize][ysize])
{
    makePlots(file,"",(Double_t*)0,(Double_t*)0,directory,matrix);
}

//***************************************************************************
//functions to make plots for 1 row, column, or cell of the matrix
//examples:
//   xvar = "nHits", yvar = "ptrel" - makes plots of nHits vs Delta_pt/pt_org
//   xvar = "all",   yvar = "pt"    - makes all plots involving Delta_pt
//                                    (not Delta_pt/pt_org)
//   xvar = "",      yvar = "all"   - makes all histograms of Delta_???
//                                    (including Delta_pt/pt_org)
//***************************************************************************

void makePlots(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString directory,
               TString xvar,TString yvar)
{
    if (yvar == "")
    {
        yvar = xvar;
        xvar = "";
    }
    Bool_t matrix[xsize][ysize];
    for (int x = 0; x < xsize; x++)
        for (int y = 0; y < ysize; y++)
        {
            bool xmatch = (xvar == "all" || xvar == xvariables[x]);
            bool ymatch = (yvar == "all" || yvar == yvariables[y]);
            matrix[x][y] = (xmatch && ymatch);
        }
    makePlots(nFiles,files,names,misalignment,values,phases,directory,matrix);
}

void makePlots(Int_t nFiles,TString *files,TString *names,TString directory,
               TString xvar,TString yvar)
{
    makePlots(nFiles,files,names,"",(Double_t*)0,(Double_t*)0,directory,
              xvar,yvar);
}

void makePlots(TString file,TString misalignment,Double_t *values,Double_t *phases,TString directory,
               TString xvar,TString yvar)
{
    int n = file.CountChar(',') + 1;
    TString *files = new TString[n];
    TString *names = new TString[n];
    setTDRStyle();
    vector<Color_t> tempcolors = colors;
    vector<Style_t> tempstyles = styles;
    for (int i = 0; i < n; i++)
    {
        TString thisfile = nPart(i+1,file,",");
        int numberofpipes = thisfile.CountChar('|');
        if (numberofpipes >= 0 && nPart(numberofpipes+1,thisfile,"|").IsDigit())
        {
            if (numberofpipes >= 1 && nPart(numberofpipes,thisfile,"|").IsDigit())
            {
                colors[i] = nPart(numberofpipes,thisfile,"|").Atoi();
                styles[i] = nPart(numberofpipes+1,thisfile,"|").Atoi();
                thisfile.Remove(thisfile.Length() - nPart(numberofpipes,thisfile,"|").Length() - nPart(numberofpipes+1,thisfile,"|").Length() - 2);
            }
            else
            {
                colors[i] = nPart(numberofpipes + 1,thisfile,"|").Atoi();
                thisfile.Remove(thisfile.Length() - nPart(numberofpipes+1,thisfile,"|").Length() - 2);
            }
        }
        files[i] = nPart(1,thisfile,"=",true);
        names[i] = nPart(2,thisfile,"=",false);
    }
    if (n == 1 && names[0] == "")
        names[0] = "scatterplot";     //With 1 file there's no legend, so this is only used in the filename of the scatterplots, if made
    makePlots(n,files,names,misalignment,values,phases,directory,xvar,yvar);
    delete[] files;
    delete[] names;
    colors = tempcolors;
    styles = tempstyles;
}

void makePlots(TString file,TString directory,TString xvar,TString yvar)
{
    makePlots(file,"",(Double_t*)0,(Double_t*)0,directory,xvar,yvar);
}

//***************************
//functions to make all plots
//***************************

void makePlots(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString directory)
{
    makePlots(nFiles,files,names,misalignment,values,phases,directory,"all","all");
}

void makePlots(Int_t nFiles,TString *files,TString *names,TString directory)
{
    makePlots(nFiles,files,names,"",(Double_t*)0,(Double_t*)0,directory);
}

void makePlots(TString file,TString misalignment,Double_t *values,Double_t *phases,TString directory)
{
    int n = file.CountChar(',') + 1;
    TString *files = new TString[n];
    TString *names = new TString[n];
    setTDRStyle();
    vector<Color_t> tempcolors = colors;
    vector<Style_t> tempstyles = styles;
    for (int i = 0; i < n; i++)
    {
        TString thisfile = nPart(i+1,file,",");
        int numberofpipes = thisfile.CountChar('|');
        if (numberofpipes >= 0 && nPart(numberofpipes+1,thisfile,"|").IsDigit())
        {
            if (numberofpipes >= 1 && nPart(numberofpipes,thisfile,"|").IsDigit())
            {
                colors[i] = nPart(numberofpipes,thisfile,"|").Atoi();
                styles[i] = nPart(numberofpipes+1,thisfile,"|").Atoi();
                thisfile.Remove(thisfile.Length() - nPart(numberofpipes,thisfile,"|").Length() - nPart(numberofpipes+1,thisfile,"|").Length() - 2);
            }
            else
            {
                colors[i] = nPart(numberofpipes + 1,thisfile,"|").Atoi();
                thisfile.Remove(thisfile.Length() - nPart(numberofpipes+1,thisfile,"|").Length() - 2);
            }
        }
        files[i] = nPart(1,thisfile,"=",true);
        names[i] = nPart(2,thisfile,"=",false);
    }
    if (n == 1 && names[0] == "")
        names[0] = "scatterplot";     //With 1 file there's no legend, so this is only used in the filename of the scatterplots, if made
    makePlots(n,files,names,misalignment,values,phases,directory);
    delete[] files;
    delete[] names;
    colors = tempcolors;
    styles = tempstyles;
}

void makePlots(TString file,TString directory)
{
    makePlots(file,"",(Double_t*)0,(Double_t*)0,directory);
}
