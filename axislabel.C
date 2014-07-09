#include "TString.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include <sstream>
using namespace std;

enum PlotType {ScatterPlot,Profile,Histogram,Resolution};
//ScatterPlot:  make a scatterplot of yvar vs. xvar
//Profile:      make a profile of yvar vs. xvar
//Histogram:    make a histogram of yvar
//Resolution:   make a plot of (width of yvar) vs xvar

//This just puts the variable name into a fancy format, with greek letters and subscripts
TString fancyname(TString variable)
{
    if (variable.Contains("Mass"))
    {
        variable.ReplaceAll("Mass","");
        variable.Prepend("m");
    }
    variable.ReplaceAll("star","*");
    variable.ReplaceAll("theta","#theta");
    variable.ReplaceAll("phi","#Phi");
    variable.ReplaceAll("hel","");
    variable.ReplaceAll("Z1","_{Z1}");
    variable.ReplaceAll("Z2","_{Z2}");
    variable.ReplaceAll("ZZ","_{4l}");
    variable.ReplaceAll("Phi*_{Z1}","Phi_{1}");
    return variable;
}

//this gives the units, to be put in the axis label
TString units(TString variable)
{
    if (variable.Contains("Mass"))
        return "GeV";
    return "";
}


//this gives the full axis label, including units.
TString axislabel(TString variable, Char_t axis, Bool_t resolution = false)
{
    stringstream s;
    if (resolution && axis == 'y')
        s << "#sigma(";
    s << fancyname(variable);
    if (resolution && axis == 'y')
        s << ")";
    if (units(variable) != "")
        s << " (" << units(variable) << ")";
    return s.str();
}

void setAxisLabels(TH1 *p, PlotType type,TString xvar,TString yvar)
{
    if (type == Histogram)
        p->SetXTitle(axislabel(yvar,'y',false));
    if (type == ScatterPlot || type == Profile || type == Resolution)
        p->SetXTitle(axislabel(xvar,'x'));

    if (type == ScatterPlot || type == Profile)
        p->SetYTitle(axislabel(yvar,'y',false));
    if (type == Resolution)
        p->SetYTitle(axislabel(yvar,'y',true));
}

void setAxisLabels(TMultiGraph *p, PlotType type,TString xvar,TString yvar)
{
    if (type == Histogram)
        p->GetXaxis()->SetTitle(axislabel(yvar,'y',false));
    if (type == ScatterPlot || type == Profile || type == Resolution)
        p->GetXaxis()->SetTitle(axislabel(xvar,'x'));

    if (type == ScatterPlot || type == Profile)
        p->GetYaxis()->SetTitle(axislabel(yvar,'y',false));
    if (type == Resolution)
        p->GetYaxis()->SetTitle(axislabel(yvar,'y',true));
}


//This divides a string at semicolons
//e.g. nPart(1,"a;b;c;d;e") = a
//if part <= 0 or part > (number of semicolons + 1), it returns ""
//It's used in misalignmentDependence.C
//It's similar to TString::Tokenize but doesn't require deleting anything
TString nPart(Int_t part,TString string,TString delimit = ";")
{
    if (part <= 0) return "";
    for (int i = 1; i < part; i++)    //part-1 times
    {
        if (string.Index(delimit) < 0) return "";
        string.Replace(0,string.Index(delimit)+1,"",0);
    }
    if (string.Index(delimit) >= 0)
        string.Remove(string.Index(delimit));
    return string;
}
