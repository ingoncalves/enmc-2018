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

void wiener(
    const TMatrixD& _X,
    const TVectorD& _d,
    TVectorD& _weights)
{
    const unsigned M = _X.GetNrows();
    const unsigned N = _X.GetNcols();

    Double_t sum;

    // matrix R
    TMatrixD R(N, N);
    R.Zero();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum = 0.0;
            for (int k = 0; k < M; k++)
            {
                sum += _X[k][i] * _X[k][j];
            }
            R[i][j] = sum / M;
        }
    }

    // vector p
    TVectorD p(N);
    p.Zero();

    // calculando p
    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (int k = 0; k < M; k++)
        {
            sum += _X[k][i] * _d[k];
        }
        p[i] = sum / M;
    }

    // solve the linear system
    R.Invert();
    _weights = R * p;

    R.Clear();
    p.Clear();
}
