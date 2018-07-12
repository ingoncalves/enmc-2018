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

#include "./lib/wiener.C"
#include "./lib/signal.C"
#include "./lib/shaper.C"
#include "./utils/matrix.C"
#include "./utils/statistics.C"
#include "./utils/units.C"

/*
 * The main function. This procedure compares three methods efficiency:
 *  - OF (Optimal Filter)
 *  - Optimal Wiener Hopf
 *  - General Wiener Hopf
 *  The result is written in the file "out/wiener_vs_of2.dat".
 */
int main()
{
    const unsigned WINDOW_SIZE(7);

    // get OF2 weights
    TVectorD weightsOF2;
    readvector("./data/of2_weights.dat", weightsOF2);

    // get General Wiener weights
    TVectorD weightsWG;
    readvector("./data/wg_weights.dat", weightsWG);

    // get Optimal Wiener weights
    TMatrixD weightsWO;
    readmatrix("./data/wo_weights.dat", weightsWO);

    // get noises dataset
    TMatrixD NOISES_TEST;
    readmatrix("./data/tile_e4mu200_test.dat", NOISES_TEST);

    // read eletronic pulse shaper file
    TVectorD shaper;
    Double_t shaperResolution;
    UInt_t shaperZeroIndex;
    readShaperFromFile("./data/pulsehi_physics.dat",
                       shaperResolution,
                       shaperZeroIndex,
                       shaper);

    Int_t totalsamples = 0;

    // output file
    const std::string FILENAME("./out/wiener_vs_of2.dat");
    std::ofstream output(FILENAME);

    for (Int_t bcid = 1; bcid <= 40; bcid++)
    {
        // filter noises by bcid
        TMatrixD NOISES_BCID;
        filtermatrix(NOISES_TEST, NOISES_BCID, [&](const TMatrixDRow_const & _row)
        {
            return _row[0] == bcid;
        });

        Int_t    nsamples = NOISES_BCID.GetNrows();
        TVectorD ampWO(nsamples);
        TVectorD ampWG(nsamples);
        TVectorD ampOF2(nsamples);
        totalsamples += nsamples;

        ampWO.Zero();
        ampWG.Zero();
        ampOF2.Zero();

        std::cout << "bcid " << bcid << ": " << nsamples << " samples" << std::endl;

        // generate pulses
        for (int i = 0; i < nsamples; i++)
        {
            // signal and desired amplitude
            TVectorD pulse;
            Double_t d, phase;
            generateSignal(WINDOW_SIZE, 0, 0, 0, 0, 0, 0, shaper, shaperResolution, shaperZeroIndex, d, phase, pulse);

            // sum the noise with the known pulse
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
                aproxWO  += weightsWO[bcid - 1][k + 1] * signal[k];
                aproxWG  += weightsWG[k] * signal[k];
                aproxOF2 += weightsOF2[k] * signal[k];
            }

            // sum to bias
            aproxWO += weightsWO[bcid - 1][WINDOW_SIZE + 1];
            aproxWG += weightsWG[WINDOW_SIZE];

            // store the amplitude value in GeV
            ampWO[i]  = adc2gev(aproxWO - d);
            ampWG[i]  = adc2gev(aproxWG - d);
            ampOF2[i] = adc2gev(aproxOF2 - d);

            pulse.Clear();
            signal.Clear();
        }

        // write results in file
        output << bcid << " ";
        output << mean(ampWO) << " ";
        output << rms(ampWO) << " ";
        output << mean(ampWG) << " ";
        output << rms(ampWG) << " ";
        output << mean(ampOF2) << " ";
        output << rms(ampOF2) << std::endl;

        // clear variables
        NOISES_BCID.Clear();
        ampWO.Clear();
        ampWG.Clear();
        ampOF2.Clear();

    }

    std::cout << "total samples: " << totalsamples << std::endl;

    // fecha o arquivo
    output.close();

    return 0;
}
