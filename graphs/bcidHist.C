/******************************************************************************
 *                         TILECAL SIMULATOR
 *
 * This code implements a simulator for TileCal (The ATLAS Tile Calorimeter)
 * with a pileup scenario.
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
#include <fstream>
#include <string>
#include "RConfig.h"
#include "TMath.h"
#include "TGraph.h"
#include "TRandom.h"

#include "../lib/signal.C"
#include "../lib/shaper.C"
#include "../utils/matrix.C"
#include "../utils/statistics.C"
#include "../utils/units.C"

void bcidHist()
{
    const Int_t WINDOW_SIZE(7);
    const Int_t BCID(2);

    // get OF2 weights
    TVectorD weightsOF2;
    readvector("../data/of2_weights.dat", weightsOF2);

    // get General Wiener weights
    TVectorD weightsWG;
    readvector("../data/wg_weights.dat", weightsWG);

    // get Optimal Wiener weights
    TMatrixD weightsWO;
    readmatrix("../data/wo_weights.dat", weightsWO);

    // get pileup noises
    TMatrixD NOISES_TEST;
    readmatrix("../data/tile_e4mu200_test.dat", NOISES_TEST);

    // read eletronic pulse shaper file
    TVectorD shaper;
    Double_t shaperResolution;
    UInt_t shaperZeroIndex;
    readShaperFromFile("../data/pulsehi_physics.dat",
                       shaperResolution,
                       shaperZeroIndex,
                       shaper);

    // define histograms
    TH1* ampWO = new TH1D( "WO", "Wiener-Hopf Otimizado", 100, -10.0, 10.0);
    TH1* ampWG = new TH1D( "WG", "Wiener-Hopf Generalizado", 100, -10.0, 10.0);
    TH1* ampOF2 = new TH1D( "OF", "OF", 100, -10.0, 10.0);

    // filter noises by bcid
    TMatrixD NOISES_BCID;
    filtermatrix(NOISES_TEST, NOISES_BCID, [&](const TMatrixDRow_const & _row)
    {
        return _row[0] == BCID;
    });

    Int_t nsamples   = NOISES_BCID.GetNrows();

    std::cout << "bcid     : " << BCID << std::endl;
    std::cout << "samples  : " << nsamples << std::endl;

    // generate pulses
    for (int i = 0; i < nsamples; i++)
    {
        // signal and desired amplitude
        TVectorD pulse;
        Double_t d, phase;
        generateSignal(WINDOW_SIZE, 0, 0, 0, 0, 0, 0, shaper, shaperResolution, shaperZeroIndex, d, phase, pulse);

        // build the signal
        TVectorD signal(WINDOW_SIZE);
        for (int j = 0; j < WINDOW_SIZE; j++)
        {
            signal[j] = NOISES_BCID[i][j + 3] + pulse[j];
        }

        // estimations
        Double_t aproxWO  = 0.0;
        Double_t aproxWG  = 0.0;
        Double_t aproxOF2 = 0.0;

        // inner product signal * weights
        for (int k = 0; k < WINDOW_SIZE; k++)
        {
            aproxWO  += weightsWO[BCID - 1][k + 1] * signal[k];
            aproxWG  += weightsWG[k] * signal[k];
            aproxOF2 += weightsOF2[k] * signal[k];
        }

        // sum to bias
        aproxWO += weightsWO[BCID - 1][WINDOW_SIZE + 1];
        aproxWG += weightsWG[WINDOW_SIZE];

        // calculate estimation error
        Double_t errorWO = adc2gev(aproxWO - d);
        Double_t errorWG = adc2gev(aproxWG - d);
        Double_t errorOF2 = adc2gev(aproxOF2 - d);

        // histogram the amplitude value
        ampWO->Fill(errorWO);
        ampWG->Fill(errorWG);
        ampOF2->Fill(errorOF2);

        // clear variables
        pulse.Clear();
        signal.Clear();
    }

    // set style
    TStyle *defStyle = new TStyle("Modern", "Modern Style");
    defStyle->SetTitleOffset(1.1, "xyz");
    defStyle->SetTitleSize(0.045, "xyz");
    defStyle->SetLabelSize(0.04, "xyz");
    defStyle->SetLegendTextSize(0.04);
    gROOT->SetStyle("Modern");

    // draw histograms
    TCanvas* c1 = new TCanvas("bcidHist", "BCID Histogram");
    THStack *hs = new THStack("hs", "");

    ampWO->SetLineWidth(2);
    ampWO->SetLineColor(kRed);
    hs->Add(ampWO);

    ampWG->SetLineWidth(2);
    ampWG->SetLineColor(kBlue);
    hs->Add(ampWG);

    ampOF2->SetLineWidth(2);
    ampOF2->SetLineColor(kBlack);
    hs->Add(ampOF2);

    hs->Draw("nostack");
    hs->GetXaxis()->SetTitle("Erro (GeV)");
    hs->GetYaxis()->SetTitle("Eventos");
    hs->GetXaxis()->CenterTitle();
    hs->GetYaxis()->CenterTitle();

    // build legend
    c1->BuildLegend(0.57, 0.89, 0.99, 0.75);

    c1->Print("../out/bcidHist.eps");
}
