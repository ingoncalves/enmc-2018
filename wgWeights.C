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

/*
 * This procedure calculates the General Wiener-Hopf weights, and
 * writes them to the file "data/wg_weights.dat".
 */
void wgWeights()
{
    const unsigned WINDOW_SIZE(7);

    TMatrixD NOISES_TRAIN;
    readmatrix("./data/tile_e4mu200_train.dat", NOISES_TRAIN);

    // read eletronic pulse shaper file
    TVectorD shaper;
    Double_t shaperResolution;
    UInt_t shaperZeroIndex;
    readShaperFromFile("./data/pulsehi_physics.dat",
                       shaperResolution,
                       shaperZeroIndex,
                       shaper);

    // get noise samples
    Int_t nsamples = NOISES_TRAIN.GetNrows();

    // variables
    TVectorD d(nsamples);
    TVectorD weights(WINDOW_SIZE + 1);
    TMatrixD X(nsamples, WINDOW_SIZE + 1);

    d.Zero();
    weights.Zero();
    X.Zero();

    // define o arquivo de saida
    const std::string FILENAME("./data/wg_weights.dat");
    std::ofstream output(FILENAME);

    // generate Wiener params
    for (Int_t i = 0; i < nsamples; i++)
    {
        // signal and desired amplitude
        TVectorD pulse;
        Double_t amplitude, phase;
        generateSignal(WINDOW_SIZE, 0, 0, 0, 0, 0, 0, shaper, shaperResolution, shaperZeroIndex, amplitude, phase, pulse);

        // sum the noise with the known pulse
        for (Int_t j = 0; j < WINDOW_SIZE; j++)
        {
            X[i][j] = NOISES_TRAIN[i][j + 3] + pulse[j];
        }

        // additional element
        X[i][7] = 1;

        // desired amplitude
        d[i] = amplitude;

        pulse.Clear();
    }

    // get Wiener weights
    wiener(X, d, weights);

    // write weights in file
    for (Int_t i = 0; i < WINDOW_SIZE + 1; i++)
        output << weights[i] << std::endl;

    // close the file
    output.close();
}
