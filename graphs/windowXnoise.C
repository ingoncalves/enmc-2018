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

void windowXnoise()
{
    const Int_t WINDOW_SIZE(7);
    const Int_t NBCID(40);

    // read noise data
    TMatrixD NOISES_TEST;
    readmatrix("../data/tile_e4mu200_test.dat", NOISES_TEST);

    TMatrixD series(NBCID, WINDOW_SIZE);

    // initialize histograms
    for (Int_t bcid = 1; bcid <= NBCID; bcid++)
    {
        // filter noises by bcid
        TMatrixD NOISES_BCID;
        filtermatrix(NOISES_TEST, NOISES_BCID, [&](const TMatrixDRow_const & _row)
        {
            return _row[0] == bcid;
        });

        // initialize histograms by each window position
        TH1D histograms[WINDOW_SIZE];
        for (Int_t i = 0; i < WINDOW_SIZE; i++)
        {
            histograms[i] = TH1D("", "", 100, -2000.0, 2000.0);
        }

        // fill histograms
        for (Int_t i = 0; i < NOISES_BCID.GetNrows(); i++)
        {
            for (Int_t j = 0; j < WINDOW_SIZE; j++)
            {
                (histograms[j]).Fill(NOISES_BCID[i][j + 3]);
            }
        }

        // get mean by each window position
        for (Int_t i = 0; i < WINDOW_SIZE; i++)
        {
            series[bcid - 1][i] = (histograms[i]).GetMean();
        }

    }

    Double_t xSerie[WINDOW_SIZE] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    Int_t color;

    // set style
    TStyle *defStyle = new TStyle("Modern", "Modern Style");
    defStyle->SetTitleOffset(1, "xyz");
    defStyle->SetTitleSize(0.045, "xyz");
    defStyle->SetLabelSize(0.04, "xyz");
    defStyle->SetLegendTextSize(0.04);
    gROOT->SetStyle("Modern");

    ///////////////////////////////////////////////////////////////////////////
    // BCID 1 to 7
    TCanvas* c1      = new TCanvas("windowXnoise_1to7", "Window x Noise for BCID 1 to 7");
    TMultiGraph *mg1 = new TMultiGraph();

    color = kBlue - 10;
    for (Int_t i = 1; i <= 7; i++)
    {
        TString title;
        title.Form("BCID %d", i);
        TGraph* g = new TGraph(WINDOW_SIZE, xSerie, series[i - 1].GetPtr());
        g->SetLineWidth(2);
        g->SetLineColor(color);
        g->SetMarkerStyle(8);
        g->SetMarkerColor(color++);
        g->SetTitle(title.Data());
        mg1->Add(g);
    }

    mg1->Draw("ACP");
    mg1->GetXaxis()->SetTitle("#acute{I}ndice da amostra temporal na janela de leitura");
    mg1->GetYaxis()->SetTitle("Amplitude (contagens de ADC)");
    mg1->GetXaxis()->CenterTitle();
    mg1->GetYaxis()->CenterTitle();
    c1->BuildLegend(0.75, 0.59, 0.89, 0.25);
    c1->Print("../out/windowXnoise_1to7.eps");
    ///////////////////////////////////////////////////////////////////////////



    ///////////////////////////////////////////////////////////////////////////
    // BCID 11 to 18
    TCanvas* c2      = new TCanvas("windowXnoise_11to18", "Window x Noise for BCID 12 to 18");
    TMultiGraph *mg2 = new TMultiGraph();
    color = kRed - 1;

    for (Int_t i = 12; i <= 18; i++)
    {
        TString title;
        title.Form("BCID %d", i);
        TGraph* g = new TGraph(WINDOW_SIZE, xSerie, series[i - 1].GetPtr());
        g->SetLineWidth(2);
        g->SetLineColor(color);
        g->SetMarkerStyle(8);
        g->SetMarkerColor(color--);
        g->SetTitle(title.Data());
        mg2->Add(g);
    }

    mg2->Draw("ACP");
    mg2->GetXaxis()->SetTitle("#acute{I}ndice da amostra temporal na janela de leitura");
    mg2->GetYaxis()->SetTitle("Amplitude (contagens de ADC)");
    mg2->GetXaxis()->CenterTitle();
    mg2->GetYaxis()->CenterTitle();
    c2->BuildLegend(0.15, 0.59, 0.29, 0.25);
    c2->Print("../out/windowXnoise_11to18.eps");
    ///////////////////////////////////////////////////////////////////////////
}
