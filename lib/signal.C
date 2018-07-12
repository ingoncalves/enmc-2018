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
 * This function generates a random signal vector.
 * The signal is composed by the multiplication of a random amplitude
 * value with a signal shape. After that, a pedestal is added to the signal.
 * Both amplitude and phase are random values. The amplitude follows the uniform
 * distribution, whose range is between 0 to 1023. The phase follows the normal
 * distribution, whose mean and standard deviation are given by parameter.
 * There is a deformation of each signal sample, which is a random value also,
 * following the normal distribution, whose mean and standard deviation
 * are given by parameter.
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
 * @param _amplitude Random amplitude value generated
 * @param _phase Random phase value generated
 * @param _signal Random signal vector generated
 */
void generateSignal(
    const UInt_t&   _size,
    const Double_t& _phaseMean,
    const Double_t& _phaseStddev,
    const Double_t& _defMean,
    const Double_t& _defStddev,
    const Double_t& _ped,
    const Double_t& _lag,
    const TVectorD& _shaper,
    const Double_t& _shaperResolution,
    const UInt_t&   _shaperZeroIndex,
    Double_t&       _amplitude,
    Double_t&       _phase,
    TVectorD &      _signal)
{
    _signal.ResizeTo(_size);
    TRandom generator(0);

    // random amplitude between [0,1023] - uniform
    _amplitude = generator.Integer(1024);
    
    // random amplitude with exp distribution
    //_amplitude = generator.Exp(300);

    // random phase - normal
    _phase = generator.Gaus(_phaseMean, _phaseStddev);

    const Double_t SAMPLING_RATE(25);

    for (int i = 0; i < _size; i++)
    {
        // random deformation - normal
        Double_t deformation = generator.Gaus(_defMean, _defStddev);
        Int_t shaperIndex = Int_t(_shaperZeroIndex)
                            - Int_t(_lag / _shaperResolution)
                            + ( i - Int_t(_size) / 2) * ( SAMPLING_RATE / _shaperResolution )
                            + round(_phase / _shaperResolution);

        if (shaperIndex < 0 || shaperIndex > _shaper.GetNoElements() - 1) shaperIndex = 0;

        _signal[i] = _amplitude * _shaper[shaperIndex] + _ped + deformation;
    }

}
