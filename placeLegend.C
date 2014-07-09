#include "TLegend.h"
#include "TLegendEntry.h"
#include "TObject.h"
#include "TClass.h"
#include "TH1.h"
#include "TList.h"
#include "TProfile.h"
#include "TF1.h"
#include <iostream>

using namespace std;

//This function places a legend in a canvas in the place where it fits.
//It searches through the canvas (using a grid of legendGrid*legendGrid) looking for a place where it fits
//The legend has a margin of margin*width on the right and left, and margin*height on the top and bottom.
//It decides whether the legend fits (in fitsHere) by looping through the entries in the legend.
//It can handle objects of several classes
//For a TH1: if the markercolor is 0 (white), it assumes that the histogram will be drawn as a histogram
//           otherwise it assumes that it will be drawn with error bars (Draw("P"))
//If it has no way of handling a class, it will place the legend around everything else and give an error message.
//More classes can be added into fitsHere, at the bottom.
//It is safe to use a legend with (TObject*)0x0 in it - that entry will just be ignored


//If the legend doesn't fit anywhere, more space will be added on top, multiplying the height of the canvas by 1+increaseby
// and then the process will repeat
//The height of the legend will also be multiplied by 1+increaseby, so that the visible height of the legend remains the same.
//The actual canvas' height is not automatically multiplied - you have to do this manually.
//The new y2max is returned by placelegend.  y1min is not changed.


//l is the legend
//width is the legend's width
//height is the legend's height
//x1min, y1min, x2max, and y2max are the corners of the canvas.

Int_t legendGrid = 100;
Double_t margin = .1;
Double_t increaseby = .1;

Bool_t fitsHere(TLegend *l,Double_t x1, Double_t y1, Double_t x2, Double_t y2);

Double_t placeLegend(TLegend *l, Double_t width, Double_t height, Double_t x1min, Double_t y1min, Double_t x2max, Double_t y2max)
{
    for (int i = legendGrid; i >= 0; i--)
    {
        for (int j = legendGrid; j >= 0; j--)
        {
            Double_t x1 = x1min * (1-(double)i/legendGrid) + (x2max - width)  * (double)i/legendGrid - margin*width;
            Double_t y1 = y1min * (1-(double)j/legendGrid) + (y2max - height) * (double)j/legendGrid - margin*height;
            Double_t x2 = x1 + (1+2*margin) * width;
            Double_t y2 = y1 + (1+2*margin) * height;
            if (fitsHere(l,x1,y1,x2,y2))
            {
                x1 += margin*width;
                y1 += margin*height;
                x2 -= margin*width;
                y2 -= margin*height;
                l->SetX1(x1);
                l->SetY1(y1);
                l->SetX2(x2);
                l->SetY2(y2);
                return y2max;
            }
        }
    }
    Double_t newy2max = y2max + increaseby * (y2max-y1min);
    Double_t newheight = height * (newy2max - y1min) / (y2max - y1min);
    return placeLegend(l,width,newheight,x1min,y1min,x2max,newy2max);
}

Bool_t fitsHere(TLegend *l,Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
    Bool_t fits = true;
    TList *list = l->GetListOfPrimitives();
    for (Int_t k = 0; list->At(k) != 0 && fits; k++)
    {
        TObject *obj = ((TLegendEntry*)(list->At(k)))->GetObject();
        if (obj == 0) continue;
        TClass *cl = obj->IsA();

        //Histogram, drawn as a histogram
        if (cl->InheritsFrom("TH1") && !cl->InheritsFrom("TH2") && !cl->InheritsFrom("TH3")
         && cl != TProfile::Class() && ((TH1*)obj)->GetMarkerColor() == kWhite)
        {
            Int_t where = 0;
            TH1 *h = (TH1*)obj;
            for (Int_t i = 1; i <= h->GetNbinsX() && fits; i++)
            {
                if (h->GetBinLowEdge(i) + h->GetBinWidth(i) < x1) continue;   //to the left of the legend
                if (h->GetBinLowEdge(i)                     > x2) continue;   //to the right of the legend
                if (h->GetBinContent(i) > y1 && h->GetBinContent(i) < y2) fits = false;   //inside the legend
                if (h->GetBinContent(i) < y1)
                {
                    if (where == 0) where = -1;             //below the legend
                    if (where == 1) fits = false;           //a previous bin was above it so there's a vertical line through it
                }
                if (h->GetBinContent(i) > y2)
                {
                    if (where == 0) where = 1;              //above the legend
                    if (where == -1) fits = false;          //a previous bin was below it so there's a vertical line through it
                }
            }
            continue;
        }
        //Histogram, drawn with Draw("P")
        else if (cl->InheritsFrom("TH1") && !cl->InheritsFrom("TH2") && !cl->InheritsFrom("TH3")
              && cl != TProfile::Class())
        //Probably TProfile would be the same but I haven't tested it
        {
            TH1 *h = (TH1*)obj;
            for (Int_t i = 1; i <= h->GetNbinsX() && fits; i++)
            {
                if (h->GetBinLowEdge(i) + h->GetBinWidth(i)/2 < x1) continue;
                if (h->GetBinLowEdge(i)                       > x2) continue;
                if (h->GetBinContent(i) > y1 && h->GetBinContent(i) < y2) fits = false;
                if (h->GetBinContent(i) + h->GetBinError(i) > y2 && h->GetBinContent(i) - h->GetBinError(i) < y2) fits = false;
                if (h->GetBinContent(i) + h->GetBinError(i) > y1 && h->GetBinContent(i) - h->GetBinError(i) < y1) fits = false;
            }
        }
        else if (cl->InheritsFrom("TF1") && !cl->InheritsFrom("TF2"))
        {
            //TF1 *f = (TF1*)obj;
            //Double_t max = f->GetMaximum(x1,x2);
            //Double_t min = f->GetMinimum(x1,x2);
            //if (min < y2 && max > y1) fits = false;
        }
        // else if (cl->InheritsFrom(...... add more objects here
        else
        {
            cout << "Don't know how to place the legend around objects of type " << obj->ClassName() << "." << endl
                 << "Add this class into placeLegend.C if you want it to work properly." << endl
                 << "The legend will still be placed around any other objects." << endl;
        }
    }
    return fits;
}
