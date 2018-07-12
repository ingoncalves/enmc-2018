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

#include <iostream>
#include "RConfig.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TRandom.h"
#include "../lib/signal.C"
#include "../lib/shaper.C"

void pileup()
{
    const Int_t    WINDOW_SIZE(7);
    const Double_t PEDESTAL(50.0);
    const Double_t SAMPLING_RATE(25);

    TRandom generator(0);

    // read eletronic pulse shaper file
    TVectorD shaper;
    Double_t shaperResolution;
    UInt_t shaperZeroIndex;
    readShaperFromFile("../data/pulsehi_physics.dat",
                       shaperResolution,
                       shaperZeroIndex,
                       shaper);

    // X axis
    Double_t x[WINDOW_SIZE];
    for (Int_t i = 0; i < WINDOW_SIZE; i++)
    {
        x[i] = (i - Int_t(WINDOW_SIZE) / 2) * SAMPLING_RATE;
    }

    TVectorD signal(WINDOW_SIZE), pileup(WINDOW_SIZE), resulting(WINDOW_SIZE);
    Double_t amplitude, phase, pileupAmplitude, pileupPhase, pileupLag;

    // pileup lag
    pileupLag = (5 - Int_t(WINDOW_SIZE) / 2) * SAMPLING_RATE;

    // generate the signal and normalize it
    generateSignal(WINDOW_SIZE, 0.0, 1.0, 0.0, 0.0, PEDESTAL, 0.0, shaper,
                   shaperResolution, shaperZeroIndex, amplitude, phase, signal);
    signal *= 1.0 / signal.Max();

    // generate the pileup and normalize it
    generateSignal(WINDOW_SIZE, 0.0, 1.0, 0.0, 0.0, PEDESTAL, pileupLag, shaper,
                   shaperResolution, shaperZeroIndex, pileupAmplitude, pileupPhase, pileup);
    pileup *= 1.0 / pileup.Max();

    // generate the resulting sinal
    resulting = signal + pileup;

    // set style
    TStyle *defStyle = new TStyle("Modern", "Modern Style");
    defStyle->SetTitleOffset(1, "xyz");
    defStyle->SetTitleSize(0.045, "xyz");
    defStyle->SetLabelSize(0.04, "xyz");
    defStyle->SetLegendTextSize(0.04);
    gROOT->SetStyle("Modern");

    // create canvas and multigraph
    TCanvas* c1 = new TCanvas("pileup", "Pileup");
    TMultiGraph *mg = new TMultiGraph();

    // create the 1st TGraph
    TGraph *gr1 = new TGraph(WINDOW_SIZE, x, signal.GetMatrixArray());
    gr1->SetLineWidth(2);
    gr1->SetLineColor(kBlack);
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerColor(kBlack);
    gr1->SetTitle("sinal de interesse");

    // create the 2nd TGraph
    TGraph *gr2 = new TGraph(WINDOW_SIZE, x, pileup.GetMatrixArray());
    gr2->SetLineWidth(2);
    gr2->SetLineColor(kRed);
    gr2->SetMarkerStyle(8);
    gr2->SetMarkerColor(kRed);
    gr2->SetTitle("sinal empilhado");

    // create the 3rd TGraph
    TGraph *gr3 = new TGraph(WINDOW_SIZE, x, resulting.GetMatrixArray());
    gr3->SetLineWidth(2);
    gr3->SetLineColor(kMagenta);
    gr3->SetMarkerStyle(8);
    gr3->SetMarkerColor(kMagenta);
    gr3->SetTitle("sinal resultante");

    // put the graphs in the multigraph
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);

    // draw the multigraph
    mg->Draw("ACP");

    // set axis title
    mg->GetXaxis()->SetTitle("Tempo (s)");
    mg->GetYaxis()->SetTitle("Amplitude (u. a.)");
    mg->GetXaxis()->CenterTitle();
    mg->GetYaxis()->CenterTitle();


    // build legend
    c1->BuildLegend(0.13, 0.69, 0.43, 0.83);

    c1->Print("../out/pileup.eps");

}
