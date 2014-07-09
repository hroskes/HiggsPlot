#include "tdrstyle.C"
#include "axislimits.C"
#include "placeLegend.C"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TText.h"
#include <iostream>

using namespace std;

TList *stufftodelete = new TList();

void deleteCanvas(TObject *canvas);
void placeholder(TString saveas = "",Bool_t wide = false);
void saveplot(TCanvas *c1,TString saveas);

bool logscale = false;

TCanvas *jhuGenPlot(Int_t nFiles,TString *files,TString *names,TString xvar,TString yvar,
                    Bool_t resolution,
                    TString saveas = "")
{
    stufftodelete->SetOwner(true);
    cout << xvar << " " << yvar << endl;
    if (xvar == "" && yvar == "")
        return 0;

    PlotType type;
    if (xvar == "")       type = Histogram;
    else if (yvar == "") {type = Histogram;    yvar = xvar; xvar = "";}
    else if (resolution)  type = Resolution;
    else if (nFiles < 1)  type = ScatterPlot;
    else                  type = Profile;
    if (nFiles < 1) nFiles = 1;

    const Int_t n = nFiles;
    
    setTDRStyle();
    gStyle->SetOptStat(0);        //for histograms, the mean and rms are included in the legend if nFiles >= 2
                                  //if nFiles == 1, there is no legend, so they're in the statbox
    //if (type == Histogram && nFiles == 1)
        //gStyle->SetOptStat(1110);
    //for a scatterplot, this is needed to show the z axis scale
    //for histograms, this is needed so that 10^-? on the right is not cut off
    //if (type == ScatterPlot || (type == Histogram))
    {
        gStyle->SetCanvasDefW(/*678*/1000);
        gStyle->SetPadRightMargin(0.115);
    }
    //else
    //{
    //    gStyle->SetCanvasDefW(600);
    //    gStyle->SetPadRightMargin(0.04);
    //}

    Int_t nBinsScatterPlotx = binsScatterPlotx;
    Int_t nBinsScatterPloty = binsScatterPloty;
    Int_t nBinsHistogram = binsHistogram;
    Int_t nBinsProfileResolution = binsProfileResolution;

    vector<TH1*> p;
    Int_t lengths[n];
    Float_t totalweights[n];

    Double_t xmin = -1,xmax = 1,ymin = -1,ymax = 1;
    if (type == Profile || type == ScatterPlot || type == Resolution)
        axislimits(nFiles,files,xvar,xmin,xmax);
    if (type == Profile || type == ScatterPlot || type == Histogram || type == Resolution)
        axislimits(nFiles,files,yvar,ymin,ymax);

    TString meansrmss[n];
    Bool_t used[n];

    for (Int_t i = 0; i < n; i++)
    {
        totalweights[i] = 0;
        stringstream sid;
        sid << "p" << i;
        TString id = sid.str();

        //for a profile or resolution, it fills a histogram, q[j], for each bin, then gets the mean and width from there.
        vector<TH1D*> q;

        if (type == ScatterPlot)
            p.push_back(new TH2F(id,"",nBinsScatterPlotx,xmin,xmax,nBinsScatterPloty,ymin,ymax));
        if (type == Histogram)
            p.push_back(new TH1D(id,"",nBinsHistogram,ymin,ymax));
        if (type == Resolution || type == Profile)
        {
            p.push_back(new TH1D(id,"",nBinsProfileResolution,xmin,xmax));
            for (Int_t j = 0; j < nBinsProfileResolution; j++)
            {
                stringstream sid2;
                sid2 << "q" << i << j;
                TString id2 = sid2.str();
                q.push_back(new TH1D(id2,"",nBinsHistogram,ymin,ymax));
            }
        }
        stufftodelete->Add(p[i]);
        p[i]->SetBit(kCanDelete,true);

        used[i] = true;
        if (files[i] == "")
        {
            used[i] = false;
            p[i]->SetLineColor(kWhite);
            p[i]->SetMarkerColor(kWhite);
            for (unsigned int j = 0; j < q.size(); j++)
                delete q[j];
            continue;
        }

        TFile *f = TFile::Open(files[i]);
        TTree *tree = (TTree*)f->Get("SelectedTree");
        if (tree == 0)
            tree = (TTree*)f->Get("GenTree");

        branches.clear();
        variables.clear();

        lengths[i] = TMath::Min(tree->GetEntries(),maxevents > 0 ? (Long64_t)maxevents : tree->GetEntries());

        Float_t x  = 0, y  = 0,
                x1 = 0, y1 = 0,
                x2 = 0, y2 = 0,
                x3 = 0, y3 = 0;
        x1 = x1; x2 = x2; x3 = x3; y1 = y1; y2 = y2; y3 = y3;
        Float_t wt = 1;
        Float_t ZZMass = 0;
        Float_t Z1Mass = 999;
        Float_t Z2Mass = 999;

        if ((type == Profile || type == ScatterPlot || type == Resolution) && xvar != yvar)
        {
            if (xvar == "discriminant")
            {
                SetBranchAddress(tree,"p0plus_VAJHU_NEW","p0plus_VAJHU",&x1);
                SetBranchAddress(tree,"bkg_VAMCFM_NEW","bkg_VAMCFM",&x2);
            }
            else
                SetBranchAddress(tree,xvar,TString(xvar).Prepend("Gen"),&x);
        }
        if (yvar == "discriminant")
        {
            SetBranchAddress(tree,"p0plus_VAJHU_NEW","p0plus_VAJHU",&y1);
            SetBranchAddress(tree,"bkg_VAMCFM_NEW","bkg_VAMCFM",&y2);
        }
        else
            SetBranchAddress(tree,yvar,TString(yvar).Prepend("Gen"),&y);
        SetBranchAddress(tree,"ZZMass","GenZZMass",&ZZMass);
        if (Z1Masscut > 0)
            SetBranchAddress(tree,"Z1Mass","GenZ1Mass",&Z1Mass);
        if (Z2Masscut > 0)
            SetBranchAddress(tree,"Z2Mass","GenZ2Mass",&Z2Mass);

        if (xvar != "wt" && yvar != "wt")
            SetBranchAddress(tree,"wt","MC_weight",&wt);

        Int_t notincluded = 0;
        cout << lengths[i] << endl;

        for (Int_t j = 0; j<lengths[i]; j++)
        {
            GetEntry(tree,j);

            if (cfunction != 0)
            {
                double mzz = ZZMass;
                discriminantconstant = cfunction->EvalPar(&mzz,cfunctionparameters);
            }

            if (xvar == "discriminant")
            {
                x = x1 / (x1 + discriminantconstant * x2);
            }
            if (yvar == "discriminant")
            {
                y = y1 / (y1 + discriminantconstant * y2);
            }
            if (ZZMass < ZZMassmin || ZZMass > ZZMassmax || Z1Mass < Z1Masscut || Z2Mass < Z2Masscut)
            {
                notincluded++;
                continue;
            }

            totalweights[i] += wt;

            if (ymin <= y && y < ymax && xmin <= x && x < xmax)
            {
                if (type == Histogram)
                    p[i]->Fill(y,wt);
                if (type == ScatterPlot)
                    ((TH2*)p[i])->Fill(x,y,wt);
                if (type == Resolution || type == Profile)
                {
                    int which = (p[i]->Fill(x,0)) - 1;
                    //get which q[j] by filling p[i] with nothing.  (TH1::Fill returns the bin number)
                    //p[i]'s actual contents are set later.
                    if (which >= 0 && (unsigned)which < q.size()) q[which]->Fill(y,wt);
                }
            }

            if (lengths[i] < 10 ? true : 
                (((j+1)/(int)(pow(10,(int)(log10(lengths[i]))-1)))*(int)(pow(10,(int)(log10(lengths[i]))-1)) == j + 1 || j + 1 == lengths[i]))
            //print when j+1 is a multiple of 10^x, where 10^x has 1 less digit than lengths[i]
            // and when it's finished
            //For example, if lengths[i] = 123456, it will print this when j+1 = 10000, 20000, ..., 120000, 123456
            //So it will print between 10 and 100 times: 10 when lengths[i] = 10^x and 100 when lengths[i] = 10^x - 1
            {
                cout << j + 1 << "/" << lengths[i] << ": "; 
                if (type == Profile || type == ScatterPlot || type == Resolution)
                    cout << x << ", " << y << endl;
                if (type == Histogram)
                    cout << y << ", " << ZZMass << endl;
            }
        }
        lengths[i] -= notincluded;

        meansrmss[i] = "";
        if (type == Histogram)
        {
            cout << "Average = " << p[i]->GetMean() << endl;
            cout << "RMS     = " << p[i]->GetRMS()  << endl;
            stringstream meanrms;
            meanrms.precision(3);
            meanrms << "#mu=" << p[i]->GetMean() << ", #sigma=" << p[i]->GetRMS();
            meansrmss[i] = meanrms.str();
        }

        if (type == Resolution)
        {
            for (Int_t j = 0; j < nBinsProfileResolution; j++)
            {
                p[i]->SetBinContent(j+1,q[j]->GetRMS());
                p[i]->SetBinError  (j+1,q[j]->GetRMSError());
                delete q[j];
            }
        }

        if (type == Profile)
        {
            for (Int_t j = 0; j < nBinsProfileResolution; j++)
            {
                p[i]->SetBinContent(j+1,q[j]->GetMean());
                p[i]->SetBinError  (j+1,q[j]->GetMeanError());
                delete q[j];
            }
        }

        setAxisLabels(p[i],type,xvar,yvar);

        p[i]->SetLineColor(colors[i]);
        if (type == Resolution || type == Profile)
        {
            p[i]->SetMarkerColor(colors[i]);
            p[i]->SetMarkerStyle(20+i);
        }
        else
        {
            p[i]->SetMarkerColor(kWhite);
            p[i]->SetMarkerStyle(1);
        }
    }

    TH1 *firstp = 0;
    for (int i = 0; i < n; i++)
    {
        if (used[i])
        {
            firstp = p[i];
            break;
        }
    }

    TCanvas *c1 = TCanvas::MakeDefCanvas();

    if (logscale)
        c1->SetLogy();

    TH1 *maxp = firstp;
    if (type == ScatterPlot)
    {
        p[0]->Scale(1.0/totalweights[0]);
        firstp->Draw("COLZ");
    }
    else if (type == Resolution || type == Profile)
    {
        vector<TGraphErrors*> g;
        TMultiGraph *list = new TMultiGraph();
        for (Int_t i = 0, ii = 0; i < n; i++, ii++)
        {
            if (!used[i])
            {
                ii--;
                continue;
            }
            g.push_back(new TGraphErrors(p[i]));
            for (Int_t j = 0; j < g[ii]->GetN(); j++)
            {
                if (g[ii]->GetY()[j] == 0 && g[ii]->GetEY()[j] == 0)
                {
                    g[ii]->RemovePoint(j);
                    j--;
                }
            }
            list->Add(g[ii]);
        }
        list->Draw("AP");
        Double_t yaxismax = list->GetYaxis()->GetXmax();
        Double_t yaxismin = list->GetYaxis()->GetXmin();
        delete list;       //automatically deletes g[i]
        if (yaxismin > 0)
        {
            yaxismax += yaxismin;
            yaxismin = 0;
        }
        firstp->GetYaxis()->SetRangeUser(yaxismin,yaxismax);
    }
    else if (type == Histogram)
    {
        Bool_t allthesame = true;
        for (Int_t i = 1; i < n && allthesame; i++)
        {
            if (totalweights[i] != totalweights[0])
                allthesame = false;
        }
        if (/*!allthesame*/true)
        {
            for (Int_t i = 0; i < n; i++)
            {
                cout << p[i]->Integral() << endl;
                p[i]->Scale(1.0/totalweights[i]);     //This does NOT include events that are out of the run number range (minrun and maxrun).
                                                      //It DOES include events that are out of the histogram range.
                cout << p[i]->Integral() << endl;
            }
        }
        maxp = (TH1D*)firstp->Clone("maxp");
        stufftodelete->Add(maxp);
        maxp->SetBit(kCanDelete,true);
        maxp->SetLineColor(kWhite);
        for (Int_t i = 1; i <= maxp->GetNbinsX(); i++)
        {
            for (Int_t j = 0; j < n; j++)
            {
                if (!used[j])
                    continue;
                maxp->SetBinContent(i,TMath::Max(maxp->GetBinContent(i),p[j]->GetBinContent(i)));
            }
        }
        maxp->Draw();
        if (!logscale)
            maxp->SetMinimum(0);
    }

    TLegend *legend = new TLegend(.6,.7,.9,.9,"","br");
    stufftodelete->Add(legend);
    legend->SetBit(kCanDelete,true);
    if (n == 1 && !used[0])
    {
        deleteCanvas(c1);
        stufftodelete->Clear();
        return 0;
    }
    for (Int_t i = 0; i < n; i++)
    {
        if (!used[i])
            continue;
        if (type == Resolution || type == Profile)
        {
            if (p[i] == firstp)
                p[i]->Draw("P");
            else
                p[i]->Draw("same P");
            legend->AddEntry(p[i],names[i],"pl");
        }
        else if (type == Histogram)
        {
            p[i]->Draw("same");
            legend->AddEntry(p[i],names[i],"l");
            legend->AddEntry((TObject*)0,meansrmss[i],"");
        }
    }
    if (type != ScatterPlot)
    {
        if (legend->GetListOfPrimitives()->At(0) == 0)
        {
            stufftodelete->Clear();
            deleteCanvas(c1);
            return 0;
        }

        
        c1->Update();
        Double_t x1min  = .98*gPad->GetUxmin() + .02*gPad->GetUxmax();
        Double_t x2max  = .02*gPad->GetUxmin() + .98*gPad->GetUxmax();
        Double_t y1min  = .98*gPad->GetUymin() + .02*gPad->GetUymax();
        Double_t y2max  = .02*gPad->GetUymin() + .98*gPad->GetUymax();
        Double_t width  = .4*(x2max-x1min);
        Double_t height = (1./20)*legend->GetListOfPrimitives()->GetEntries()*(y2max-y1min);
        if (type == Histogram)
        {
            width *= 2;
            height /= 2;
            legend->SetNColumns(2);
        }
        Double_t newy2max = placeLegend(legend,width,height,x1min,y1min,x2max,y2max);
        if (!logscale)
            maxp->GetYaxis()->SetRangeUser(gPad->GetUymin(),(newy2max-.02*gPad->GetUymin())/.98);
                
        legend->SetFillStyle(0);
        legend->Draw();
    }

    if (saveas != "")
        saveplot(c1,saveas);

    return c1;
}


//make a 1D histogram of Delta_yvar

TCanvas *jhuGenPlot(Int_t nFiles,TString *files,TString *names,TString var,
                    TString saveas = "")
{
    return jhuGenPlot(nFiles,files,names,"",var,false,saveas);
}



//For 1 file

TCanvas *jhuGenPlot(TString file,TString xvar,TString yvar,Bool_t profile,
                    Bool_t resolution,
                    TString saveas = "")
{
    Int_t nFiles = 0;
    if (profile)                       //it interprets nFiles < 1 as 1 file, make a scatterplot
        nFiles = 1;
    TString *files = &file;
    TString name = "";
    TString *names = &name;
    return jhuGenPlot(nFiles,files,names,xvar,yvar,resolution,saveas);
}

//make a 1D histogram of Delta_yvar

TCanvas *jhuGenPlot(TString file,TString var,
                    TString saveas = "")
{
    Int_t nFiles = 1;
    TString *files = &file;
    TString name = "";
    TString *names = &name;
    return jhuGenPlot(nFiles,files,names,var,saveas);
}

void placeholder(TString saveas,Bool_t wide)
{
    setTDRStyle();
    if (wide)
        gStyle->SetCanvasDefW(678);
    else
        gStyle->SetCanvasDefW(600);
    TText line1(.5,.6,"This is a placeholder so that when there are");
    TText line2(.5,.4,"4 plots per line it lines up nicely");
    line1.SetTextAlign(22);
    line2.SetTextAlign(22);
    TCanvas *c1 = TCanvas::MakeDefCanvas();
    line1.Draw();
    line2.Draw();
    if (saveas != "")
        saveplot(c1,saveas);
    deleteCanvas(c1);
}

void saveplot(TCanvas *c1,TString saveas)
{
    if (saveas == "")
        return;
    saveas.ReplaceAll("#","");
    TString saveas2 = saveas,
            saveas3 = saveas,
            saveas4;
    saveas2.ReplaceAll(".pngepsroot","");
    saveas3.Remove(saveas3.Length()-11);
    if (saveas2 == saveas3)
    {
        stringstream s1,s2,s3;
        s1 << saveas2 << ".png";
        s2 << saveas2 << ".eps";
        s3 << saveas2 << ".root";
        saveas2 = s1.str();
        saveas3 = s2.str();
        saveas4 = s3.str();
        c1->SaveAs(saveas2);
        c1->SaveAs(saveas3);
        c1->SaveAs(saveas4);
        return;
    }
    else
    {
        c1->SaveAs(saveas);
    }
}

void deleteCanvas(TObject *canvas)
{
    if (canvas == 0) return;
    if (!canvas->InheritsFrom("TCanvas"))
    {
        delete canvas;
        return;
    }
    TCanvas *c1 = (TCanvas*)canvas;
    TList *list = c1->GetListOfPrimitives();
    list->SetOwner(true);
    list->Clear();
    delete c1;
}

