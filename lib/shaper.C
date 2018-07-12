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

#include <vector>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "RConfig.h"

/*
 * This function generates the eletronic signal shaper vector,
 * reading the shaper values from a time series in a external file.
 *
 * The shaper file has two columns. The first is the time in ns and
 * the second is the shape weight.
 *
 * External shaper file example:
 *
 *      -75.5 0.00000000
 *      -75.0 0.00002304
 *      -74.5 0.00005178
 *      -74.0 0.00008592
 *      ...
 *      -1.5 0.99758100
 *      -1.0 0.99892900
 *      -0.5 0.99973300
 *       0.0 1.00000000
 *       0.5 0.99973500
 *       1.0 0.99894400
 *       1.5 0.99763200
 *       ...
 *       123.5 0.00196603
 *       124.0 0.00191204
 *       124.5 0.00185470
 *
 * @param _path External file path
 * @param _resolution Generated shaper resolution
 * @param _zeroIndex Time series zero index
 * @param _shaper Output shaper vector
 */
void readShaperFromFile(
    const char* _path,
    Double_t&   _resolution,
    UInt_t&     _zeroIndex,
    TVectorD &  _shaper)
{
    std::vector<Double_t> times;
    std::vector<Double_t> weights;
    std::ifstream file;

    file.open(_path);
    if (file)
    {
        unsigned i(0);
        Double_t a, b;
        while (file >> a >> b) //loop on the input operation, not eof
        {
            times.push_back(a);
            weights.push_back(b);
            if (a == 0.0) _zeroIndex = i;
            i++;
        }
    }
    else
    {
        throw std::invalid_argument( "invalid shaper path" );
    }

    // read eletronic pulse shaper file
    int size = weights.size();
    _resolution = size > 2 ? times[1] - times[0] : times[0];

    _shaper.ResizeTo(size);
    for (int j = 0; j < size; j++)
    {
        _shaper[j] = weights[j];
    }

    weights.clear();
    times.clear();
}
