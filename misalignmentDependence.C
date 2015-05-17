#include "TGraph2DErrors.h"
#include "jhuGenPlot.C"

using namespace std;

Bool_t fitHist(TString misalignment);

//Note: this misalignmentDependence is adapted from track splitting.  I have not changed the name
//      but it's for general fitting, not only misalignments.


//This can do three different things:
// (1) if xvar == "", it will plot the mean (if !resolution) or width (if resolution) of Delta_yvar as a function
//     of the misalignment values, as given in values.  misalignment (e.g. sagitta, elliptical) will be used as the
//     x axis label.
// (2) if xvar != "", it will fit the profile/resolution to a function.  If parameter > 0, it will plot the parameter given by parameter as
//     a function of the misalignment.  parametername is used as the y axis label.  You can put a semicolon in parametername
//     to separate the name from the units.  Functionname describes the funciton, and is put in brackets in the y axis label.
//     For example, to fit to Delta_pt = [0]*(eta_org-[1]), you could use functionname = "#Deltap_{T} = A(#eta_{org}-B)",
//     parameter = 0, and parametername = "A;GeV".
// (3) if parameter < 0, it will draw the profile/resolution along with the fitted functions.
//     The parameter of interest is still indicated by parameter, which is transformed to -parameter - 1.
//     For example, -1 --> 0, -2 --> 1, -3 --> 2, ...
//     This parameter's value and error will be in the legend.  You still need to enter parametername and functionname,
//     because they will be used for labels.

//The best way to run misalignmentDependence is through makePlots.  If you want to run misalignmentDependence directly,
//the LAST function, all the way at the bottom of this file, is probably the most practical to use (for all three of these).


// The first function takes a canvas as its argument.  This canvas needs to have been produced with jhuGenPlot using
// the same values of xvar, yvar, and resolution or something strange could happen.

void misalignmentDependence(TCanvas *c1old,
                            Int_t nFiles,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TF1 *function,Int_t parameter,TString parametername = "",TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    cout << saveas << endl;
    if (c1old == 0) return;
    c1old = (TCanvas*)c1old->Clone("c1old");
    if (misalignment == "" || yvar == "") return;
    Bool_t drawfits = (parameter < 0);
    if (parameter < 0)
        parameter = -parameter - 1;   //-1 --> 0, -2 --> 1, -3 --> 2, ...
    TString yaxislabel = nPart(1,parametername);
    TString parameterunits = nPart(2,parametername);
    if (parameterunits != "")
        yaxislabel.Append(" (").Append(parameterunits).Append(")");
    TList *list = c1old->GetListOfPrimitives();
    int n = list->GetEntries() - 2 - (xvar == "");

    setTDRStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetFitFormat("5.4g");
    gStyle->SetFuncColor(2);
    gStyle->SetFuncStyle(1);
    gStyle->SetFuncWidth(1);
    if (!drawfits)
    {
        gStyle->SetCanvasDefW(678);
        gStyle->SetPadRightMargin(0.115);
    }

    TH1 **p = new TH1*[n];
    TF1 **f = new TF1*[n];
    Bool_t used[n];
    for (Int_t i = 0; i < n; i++)
    {
        stringstream s0;
        s0 << "p" << i;
        TString pname = s0.str();
        p[i] = (TH1*)list->/*At(i+1+(xvar == ""))*/FindObject(pname);
        used[i] = (p[i] != 0);
        if (used[i])
            p[i]->SetDirectory(0);
        if (xvar == "" && function == 0)
            continue;
        stringstream s;
        s << function->GetName() << i;
        TString newname = s.str();
        f[i] = (TF1*)function->Clone(newname);
        stufftodelete->Add(f[i]);
    }

    Double_t *result = new Double_t[nFiles];
    Double_t *error  = new Double_t[nFiles];
    if (xvar == "" && function == 0)
    {
        yaxislabel = axislabel(yvar,'y',resolution);
        for (Int_t i = 0; i < nFiles; i++)
        {
            if (!used[i]) continue;
            if (!resolution)
            {
                result[i] = p[i]->GetMean();
                error[i]  = p[i]->GetMeanError();
            }
            else
            {
                result[i] = p[i]->GetRMS();
                error[i]  = p[i]->GetRMSError();
            }
            cout << result[i] << " +/- " << error[i] << endl;
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            if (!used[i]) continue;
            f[i]->SetLineColor(colors[i]);
            f[i]->SetLineStyle(styles[i]);
            f[i]->SetLineWidth(1);
            p[i]->SetMarkerColor(colors[i]);
            p[i]->SetMarkerStyle(20+i);
            p[i]->SetLineColor(colors[i]);
            p[i]->SetLineStyle(styles[i]);
            p[i]->Fit(f[i],"IM");
            error[i]  = f[i]->GetParError (parameter);
            if (function->GetName() == TString("sine"))
            {
                if (f[i]->GetParameter(0) < 0)
                {
                    f[i]->SetParameter(0,-f[i]->GetParameter(0));
                    f[i]->SetParameter(2,f[i]->GetParameter(2)+pi);
                }
                while(f[i]->GetParameter(2) >= 2*pi)
                    f[i]->SetParameter(2,f[i]->GetParameter(2)-2*pi);
                while(f[i]->GetParameter(2) < 0)
                    f[i]->SetParameter(2,f[i]->GetParameter(2)+2*pi);
            }
            result[i] = f[i]->GetParameter(parameter);
        }
    }


    TCanvas *c1 = TCanvas::MakeDefCanvas();

    if (drawfits && !(xvar == "" && function == 0) && yvar != "")
    {
        TString legendtitle = "[";
        legendtitle.Append(functionname);
        legendtitle.Append("]");
        TLegend *legend = new TLegend(.7,.7,.9,.9,legendtitle,"br");
        stufftodelete->Add(legend);
        TString drawoption = "";
        TH1 *maxp = (TH1*)list->FindObject("maxp");
        if (maxp != 0)
        {
            maxp->Draw();
            drawoption = "same";
        }
        for (int i = 0; i < n; i++)
        {
            if (!used[i]) continue;
            p[i]->Draw(drawoption);
            f[i]->Draw("same");
            drawoption = "same";

            stringstream s;
            s.precision(3);
            s << nPart(1,parametername) << " = " <<  result[i] << " #pm " << error[i];
            if (parameterunits != "") s << " " << parameterunits;
            TString str = s.str();
            legend->AddEntry(p[i],names[i],"pl");
            legend->AddEntry(f[i],str,"l");
        }
        c1->Update();
        Double_t x1min  = .98*gPad->GetUxmin() + .02*gPad->GetUxmax();
        Double_t x2max  = .02*gPad->GetUxmin() + .98*gPad->GetUxmax();
        Double_t y1min  = .98*gPad->GetUymin() + .02*gPad->GetUymax();
        Double_t y2max  = .02*gPad->GetUymin() + .98*gPad->GetUymax();
        Double_t width  = .4*(x2max-x1min);
        Double_t height = (1./20)*legend->GetListOfPrimitives()->GetEntries()*(y2max-y1min);
        width *= 2;
        height /= 2;
        legend->SetNColumns(2);

        Double_t newy2max = placeLegend(legend,width,height,x1min,y1min,x2max,y2max);
        p[0]->GetYaxis()->SetRangeUser(gPad->GetUymin(),(newy2max-.02*gPad->GetUymin())/.98);
        if (maxp != 0)
            maxp->GetYaxis()->SetRangeUser(gPad->GetUymin(),(newy2max-.02*gPad->GetUymin())/.98);

        legend->SetFillStyle(0);
        legend->Draw();
    }
    else
    {
        if (values == 0) return;

        Bool_t phasesmatter = false;
        if (misalignment == "elliptical" || misalignment == "sagitta" || misalignment == "skew")
        {
            if (phases == 0)
            {
                cout << "This misalignment has a phase, but you didn't supply the phases!" << endl
                     << "Can't produce plots depending on the misalignment value." << endl;
                return;
            }
            int firstnonzero = -1;
            for (Int_t i = 0; i < nFiles; i++)
            {
                if (values[i] == 0) continue;                    //if the amplitude is 0 the phase is arbitrary
                if (firstnonzero == -1) firstnonzero = i;
                if (phases[i] != phases[firstnonzero])
                    phasesmatter = true;
            }
        }

        if (!phasesmatter)
        {
            TGraphErrors *g = new TGraphErrors(nFiles,values,result,(Double_t*)0,error);
            g->SetName("");
            stufftodelete->Add(g);

            TString xaxislabel = "#epsilon_{";
            xaxislabel.Append(misalignment);
            xaxislabel.Append("}");
            g->GetXaxis()->SetTitle(xaxislabel);
            if (xvar != "")
            {
                yaxislabel.Append("   [");
                yaxislabel.Append(functionname);
                yaxislabel.Append("]");
            }
            g->GetYaxis()->SetTitle(yaxislabel);

            g->SetMarkerColor(colors[0]);
            g->SetMarkerStyle(20);

            g->Draw("AP");
            Double_t yaxismax = g->GetYaxis()->GetXmax();
            Double_t yaxismin = g->GetYaxis()->GetXmin();
            if (yaxismin > 0)
            {
                yaxismax += yaxismin;
                yaxismin = 0;
            }
            g->GetYaxis()->SetRangeUser(yaxismin,yaxismax);
            g->Draw("AP");
        }
        else
        {
            double *xvalues = new double[nFiles];
            double *yvalues = new double[nFiles];      //these are not physically x and y (except in the case of skew)
            for (int i = 0; i < nFiles; i++)
            {
                xvalues[i] = values[i] * cos(phases[i]);
                yvalues[i] = values[i] * sin(phases[i]);
            }
            TGraph2DErrors *g = new TGraph2DErrors(nFiles,xvalues,yvalues,result,(Double_t*)0,(Double_t*)0,error);
            g->SetName("");
            stufftodelete->Add(g);
            delete[] xvalues;        //A TGraph2DErrors has its own copy of xvalues and yvalues, so it's ok to delete these copies.
            delete[] yvalues;
            
            TString xaxislabel = "#epsilon_{";
            xaxislabel.Append(misalignment);
            xaxislabel.Append("}cos(#delta)");
            TString realyaxislabel = xaxislabel;
            realyaxislabel.ReplaceAll("cos(#delta)","sin(#delta)");
            g->GetXaxis()->SetTitle(xaxislabel);
            g->GetYaxis()->SetTitle(realyaxislabel);
            TString zaxislabel = /*"fake"*/yaxislabel;         //yaxislabel is defined earlier
            if (xvar != "")
            {
                zaxislabel.Append("   [");
                zaxislabel.Append(functionname);
                zaxislabel.Append("]");
            }
            g->GetZaxis()->SetTitle(zaxislabel);
            g->SetMarkerStyle(20);
            g->Draw("pcolerr");
        }
    }

    if (saveas != "")
    {
        saveplot(c1,saveas);
        delete[] p;
        delete[] f;
        delete[] result;
        delete[] error;
        delete c1old;
    }
}


//This version allows you to show multiple parameters.  It runs the previous version multiple times, once for each parameter.
//saveas will be modified to indicate which parameter is being used each time.

void misalignmentDependence(TCanvas *c1old,
                            Int_t nFiles,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TF1 *function,Int_t nParameters,Int_t *parameters,TString *parameternames,TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    for (int i = 0; i < nParameters; i++)
    {
        TString saveasi = saveas;
        TString insert = nPart(1,parameternames[i]);
        insert.Prepend(".");
        saveasi.Insert(saveasi.Last('.'),insert);    //insert the parameter name before the file extension
        misalignmentDependence(c1old,
                               nFiles,names,misalignment,values,phases,xvar,yvar,
                               function,parameters[i],parameternames[i],functionname,
                               resolution,
                               saveasi);
    }
}


//This version does not take a canvas as its argument.  It runs jhuGenPlot to produce the canvas.

void misalignmentDependence(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TF1 *function,Int_t parameter,TString parametername = "",TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    misalignmentDependence(jhuGenPlot(nFiles,files,names,xvar,yvar,resolution,""),
                           nFiles,names,misalignment,values,phases,xvar,yvar,
                           function,parameter,parametername,functionname,
                           resolution,saveas);
}

void misalignmentDependence(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TF1 *function,Int_t nParameters,Int_t *parameters,TString *parameternames,TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    for (int i = 0; i < nParameters; i++)
    {
        TString saveasi = saveas;
        TString insert = nPart(1,parameternames[i]);
        insert.Prepend(".");
        saveasi.Insert(saveasi.Last('.'),insert);    //insert the parameter name before the file extension
        misalignmentDependence(nFiles,files,names,misalignment,values,phases,xvar,yvar,
                               function,parameters[i],parameternames[i],functionname,
                               resolution,
                               saveasi);
    }
}


// This version allows you to use a string for the function.  It creates a TF1 using this string and uses this TF1

void misalignmentDependence(TCanvas *c1old,
                            Int_t nFiles,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TString function,Int_t parameter,TString parametername = "",TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    TF1 *f = new TF1("func",function.Data());
    misalignmentDependence(c1old,nFiles,names,misalignment,values,phases,xvar,yvar,f,parameter,parametername,functionname,resolution,saveas);
    delete f;
}

void misalignmentDependence(TCanvas *c1old,
                            Int_t nFiles,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TString function,Int_t nParameters,Int_t *parameters,TString *parameternames,TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    for (int i = 0; i < nParameters; i++)
    {
        TString saveasi = saveas;
        TString insert = nPart(1,parameternames[i]);
        insert.Prepend(".");
        saveasi.Insert(saveasi.Last('.'),insert);    //insert the parameter name before the file extension
        misalignmentDependence(c1old,
                               nFiles,names,misalignment,values,phases,xvar,yvar,
                               function,parameters[i],parameternames[i],functionname,
                               resolution,
                               saveasi);
    }
}


void misalignmentDependence(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TString function,Int_t parameter,TString parametername = "",TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    TF1 *f = new TF1(TString("func"),function.Data());
    misalignmentDependence(nFiles,files,names,misalignment,values,phases,xvar,yvar,f,parameter,parametername,functionname,resolution,saveas);
    delete f;
}

void misalignmentDependence(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                            TString function,Int_t nParameters,Int_t *parameters,TString *parameternames,TString functionname = "",
                            Bool_t resolution = false,
                            TString saveas = "")
{
    for (int i = 0; i < nParameters; i++)
    {
        TString saveasi = saveas;
        TString insert = nPart(1,parameternames[i]);
        insert.Prepend(".");
        saveasi.Insert(saveasi.Last('.'),insert);    //insert the parameter name before the file extension
        misalignmentDependence(nFiles,files,names,misalignment,values,phases,xvar,yvar,
                               function,parameters[i],parameternames[i],functionname,
                               resolution,
                               saveasi);
    }
}




//This version does not take a function as its argument.  It automatically determines what function, parameter,
//functionname, and parametername to use based on misalignment, xvar, yvar, and resolution.
//However, you have to manually put into the function which plots to fit to what shapes.
//I have put in several fits using elliptical and sagitta misalignments.
//The 2012A data, using the prompt geometry, is a nice example if you want to see an elliptical misalignment.
//If drawfits is true, it draws the fits; otherwise it plots the parameter as a function of misalignment as given by values.

//If the combination of misalignment, xvar, yvar, resolution, has a default function to use, it returns true,
// otherwise it returns false.

//This is the version called by makeThesePlots.C

Bool_t misalignmentDependence(TCanvas *c1old,
                              Int_t nFiles,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                              Bool_t drawfits = true,
                              Bool_t resolution = false,
                              TString saveas = "")
{
    TF1 *f = 0;
    TString functionname = "";

    //if only one parameter is of interest
    TString parametername = "";
    Int_t parameter = 9999;

    //if multiple parameters are of interest
    Int_t nParameters = -1;
    TString *parameternames = 0;
    Int_t *parameters = 0;

    if (misalignment == "breitwigner")
    {
        if (xvar == "" && yvar == "mH")
        {
            f = new TF1("breitwigner","[2]*TMath::BreitWigner(x,[0],[1])");
            f->SetParameters(91,2);
            f->SetParLimits(0,70,105);
            f->SetParLimits(1,0,10);
            nParameters = 2;
            Int_t tempParameters[2] = {0,1};
            TString tempParameterNames[2] = {"M;GeV","#gamma;GeV"};
            parameters = tempParameters;
            parameternames = tempParameterNames;
            functionname = "Breit Wigner";
        }
    }

    if (xvar == "" && functionname == "")
    {
        if (c1old == 0 || misalignment == "" || values == 0) return false;
        misalignmentDependence(c1old,nFiles,names,misalignment,values,phases,xvar,yvar,(TF1*)0,0,"","",resolution,saveas);
        return true;
    }

    if (functionname == "") return false;
    if (drawfits)
    {
        parameter = -parameter-1;
        for (int i = 0; i < nParameters; i++)
            parameters[i] = -parameters[i]-1;
    }
    if (nParameters > 0)
        misalignmentDependence(c1old,nFiles,names,misalignment,values,phases,xvar,yvar,
                               f,nParameters,parameters,parameternames,functionname,resolution,saveas);
    else
        misalignmentDependence(c1old,nFiles,names,misalignment,values,phases,xvar,yvar,
                               f,parameter,parametername,functionname,resolution,saveas);
    delete f;
    return true;
    
}

//This is the most practically useful version.  It does not take a canvas, but produces it automatically and then determines what
//function to fit it to.

Bool_t misalignmentDependence(Int_t nFiles,TString *files,TString *names,TString misalignment,Double_t *values,Double_t *phases,TString xvar,TString yvar,
                              Bool_t drawfits = true,
                              Bool_t resolution = false,
                              TString saveas = "")
{
    return misalignmentDependence(jhuGenPlot(nFiles,files,names,xvar,yvar,resolution,""),
                                  nFiles,names,misalignment,values,phases,xvar,yvar,
                                  drawfits,resolution,saveas);
}

//return true if any histograms are fit for this misalignment, false otherwise
Bool_t fitHist(TString misalignment)
{
    if (misalignment == "breitwigner")
        return true;
    return false;
}
