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
#include "TRandom.h"

/*
 * This function generates a random noise vector, following
 * the normal distribution, whose mean and standard deviation
 * are given by parameter.
 * The value of the parameter _ped will be added to
 * each noise sample.
 *
 * @param _size Noise samples length
 * @param _mean Noise mean
 * @param _stddev Noise standard deviation
 * @param _ped Pedestal
 * @param _noise Random noise vector generated
 */
void generateNoise(
    const UInt_t&   _size,
    const Double_t& _mean,
    const Double_t& _stddev,
    const Double_t& _ped,
    TVectorD &      _noise)
{
    _noise.ResizeTo(_size);

    // generate noise
    TRandom generator(0);
    for (int i = 0; i < _size; i++)
        _noise[i] = _ped + generator.Gaus(_mean, _stddev);
}
