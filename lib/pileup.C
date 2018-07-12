/******************************************************************************
 *                         TILECAL SIMULATOR
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
#include "TRandom.h"

/*
 * This function generates a pileup noise vector.
 * For each sample position, there is a probability of _prob x 100%
 * to generate a pileup signal, which will be added to a resultant
 * pileup vector.
 *
 * @param _size Signal samples length
 * @param _phaseMean Phase mean
 * @param _phaseStddev Phase standard deviation
 * @param _defMean Deformation mean
 * @param _defStddev Deformation standard deviation
 * @param _ped Pedestal
 * @param _lag Signal lag in nanoseconds
 * @param _shaper Vector with eletronic signal shaper
 * @param _shaperResolution Shaper resolution in nanoseconds
 * @param _shaperZeroIndex Shaper index of time series zero
 * @param _prob Probability of generate a pileup in a sample position [0.0;1.0]
 * @param _signal Random signal vector generated
 */
void generatePileup(
    const UInt_t&   _size,
    const Double_t& _phaseMean,
    const Double_t& _phaseStddev,
    const Double_t& _defMean,
    const Double_t& _defStddev,
    const Double_t& _ped,
    const TVectorD& _shaper,
    const Double_t& _shaperResolution,
    const UInt_t&   _shaperZeroIndex,
    const UInt_t&   _samplingRate,
    const Double_t& _prob,
    TVectorD &      _pileup)
{
    _pileup.ResizeTo(_size);

    TRandom generator(0);
    TVectorD noise(_size);
    Double_t pileupLag, pileupAmplitude, pileupPhase;

    for (int j = 0; j < _size; j++)
    {
        if (generator.Uniform(0.0, 1.0) < _prob)
        {
            noise.Clear();
            pileupLag = (j - Int_t(_size) / 2) * _samplingRate;
            generateSignal(_size, _phaseMean, _phaseStddev, _defMean, _defStddev, _ped,
                           pileupLag, // lag in nanoseconds
                           _shaper, _shaperResolution, _shaperZeroIndex, pileupAmplitude, pileupPhase, noise);
            _pileup += noise;
        }
    }

}
