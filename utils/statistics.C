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

#include "TMath.h"
#include "RConfig.h"

Double_t mean(const TVectorD& _data)
{
    if (_data.GetNrows() == 0)
        return 0.0;

    Double_t sum = 0.0;
    for (Int_t i = 0; i < _data.GetNoElements(); ++i)
    {
        sum += _data[i];
    }
    return sum / _data.GetNoElements();
}

Double_t rms(const TVectorD& _data)
{
    if (_data.GetNrows() == 0)
        return 0.0;

    Double_t sum = 0.0;
    Double_t m   = mean(_data);
    for (Int_t i = 0; i < _data.GetNoElements(); ++i)
    {
        sum += (_data[i] - m) * (_data[i] - m);
    }
    return TMath::Sqrt(sum / (_data.GetNoElements() - 1.0));
}
