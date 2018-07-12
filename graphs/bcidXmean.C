/******************************************************************************
 *                         TILECAL SIMULATOR
 *
 * This code tests the signal pileup.
 *
 *
 * Bernardo S. Peralva    <bernardo@iprj.uerj.br>
 * Guilherme I. Gonçalves <ggoncalves@iprj.uerj.br>
 *
 * Copyright (C) 2018 Bernardo & Guilherme

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 *   http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

#include "RConfig.h"
#include "TMath.h"
#include "TGraph.h"
#include "../utils/matrix.C"

void bcidXmean()
{
    // read output data
    TMatrixD data;
    readmatrix("../out/wiener_vs_of2.dat", data);

    // get X axis length
    Int_t nrows = data.GetNrows();

    // define data series
    TVectorD bcid(nrows),
             woMean(nrows),
             wgMean(nrows),
             of2Mean(nrows);

    // get data series from output data
    for (Int_t i = 0; i < nrows ; i++)
    {
        bcid[i]    = data[i][0];
        woMean[i]  = data[i][1];
        wgMean[i]  = data[i][3];
        of2Mean[i] = data[i][5];
    }

    // set style
    TStyle *defStyle = new TStyle("Modern", "Modern Style");
    defStyle->SetTitleOffset(1, "xyz");
    defStyle->SetTitleSize(0.045, "xyz");
    defStyle->SetLabelSize(0.04, "xyz");
    defStyle->SetLegendTextSize(0.04);
    gROOT->SetStyle("Modern");

    // create canvas and multigraph
    TCanvas* c1 = new TCanvas("bcidXmean", "BCID x Mean");
    TMultiGraph *mg = new TMultiGraph();

    // create the 1st TGraph
    TGraph *gr1 = new TGraph(nrows, bcid.GetMatrixArray(), woMean.GetMatrixArray());
    gr1->SetLineWidth(2);
    gr1->SetLineColor(kRed);
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerColor(kRed);
    gr1->SetTitle("Wiener-Hopf Otimizado");

    // create the 2nd TGraph
    TGraph *gr2 = new TGraph(nrows, bcid.GetMatrixArray(), wgMean.GetMatrixArray());
    gr2->SetLineWidth(2);
    gr2->SetLineColor(kBlue);
    gr2->SetMarkerStyle(8);
    gr2->SetMarkerColor(kBlue);
    gr2->SetTitle("Wiener-Hopf Generalizado");

    // create the 3rd TGraph
    TGraph *gr3 = new TGraph(nrows, bcid.GetMatrixArray(), of2Mean.GetMatrixArray());
    gr3->SetLineWidth(2);
    gr3->SetMarkerStyle(8);
    gr3->SetTitle("OF");

    // put the graphs in the multigraph
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);

    // draw the multigraph
    mg->Draw("ACP");

    // set axis title
    mg->GetXaxis()->SetTitle("#acute{I}ndice da colis#tilde{a}o no trem de colis#tilde{o}es");
    mg->GetYaxis()->SetTitle("M#acute{e}dia do erro de estimac#tilde{a}o (GeV)");
    mg->GetXaxis()->CenterTitle();
    mg->GetYaxis()->CenterTitle();

    // build legend
    c1->BuildLegend(0.57,0.99,0.99,0.85);

    c1->Print("../out/bcidXmean.eps");
}
