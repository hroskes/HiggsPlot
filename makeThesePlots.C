#include "-----directory of plotting code-----/makePlots.C"

const int nFiles = -----number of files-----;
TString files[nFiles] = {
    "-----file1-----",
    "-----file2-----",
    ...
};
TString names[nFiles] = {
    "-----name1-----",
    "-----name2-----",         //these are used in the legend
    ...
};

TString plotdirectory = "-----directory to save plots-----";

/*******************************
To make ALL plots (313 in total)
*******************************/

void makeThesePlots()
{
    makePlots(nFiles,files,names,plotdirectory);
}

/**************************************************************************
to make all plots involving a single x or y variable, or both

   xvar = "nHits", yvar = "ptrel" - makes plots of nHits vs Delta_pt/pt_org
                                    (4 total - profile and resolution,
                                     of Delta_pt/pt_org and its pull
                                     distribution)
   xvar = "all",   yvar = "pt"    - makes all plots involving Delta_pt
                                    (not Delta_pt/pt_org)
                                    (38 plots total:
                                     histogram and pull distribution, and
                                     their mean and width as a function
                                     of the 9 x variables)
   xvar = "",      yvar = "all"   - makes all histograms of all y variables
                                    (including Delta_pt/pt_org)
                                    (16 plots total - 8 y variables,
                                     regular and pull histograms)
**************************************************************************/

void makeThesePlots(TString xvar,TString yvar)
{
    makePlots(nFiles,files,names,plotdirectory,xvar,yvar);
}

/*************************************************************************************
To make a custom selection of plots
The matrix should look like this - use true to make that plot, and false not to make it.
        y var   pt rel  pt      eta     phi     dz      dxy     theta   q/pt    (hist)
x var     +---------------------------------------------------------------------------
pt        |
eta       |
phi       |
dz        |
dxy       |
theta     |
q/pt      |
run #     |
nHits     |
(hist)    |

*************************************************************************************/


void makeThesePlots(Bool_t matrix[xsize][ysize])
{
    makePlots(nFiles,files,names,plotdirectory,matrix);
}

