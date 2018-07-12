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

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include "TMath.h"
#include "RConfig.h"

void readmatrix(
    const char* _path,
    TMatrixD& _out)
{
    std::vector<std::vector<Double_t>> m;

    std::ifstream file;
    std::string line;

    file.open(_path);
    if (!file.is_open())
    {
        throw std::invalid_argument( "invalid shaper path" );
    }

    while (!std::getline(file, line, '\n').eof())
    {
        istringstream reader(line);

        std::vector<Double_t> lineData;
        std::string::const_iterator i = line.begin();

        while (!reader.eof())
        {
            Double_t val;
            reader >> val;

            if (reader.fail())
                break;

            lineData.push_back(val);
        }

        m.push_back(lineData);
    }

    file.close();

    Int_t rows = m.size();
    Int_t cols = m[0].size();

    _out.ResizeTo(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            _out[i][j] = m[i][j];
}


void readvector(
    const char* _path,
    TVectorD& _out)
{
    std::vector<Double_t> v;

    std::ifstream file;
    std::string line;

    file.open(_path);
    if (!file.is_open())
    {
        throw std::invalid_argument( "invalid shaper path" );
    }

    while (!std::getline(file, line, '\n').eof())
    {
        istringstream reader(line);

        std::string::const_iterator i = line.begin();

        while (!reader.eof())
        {
            Double_t val;
            reader >> val;

            if (reader.fail())
                break;

            v.push_back(val);
        }
    }

    file.close();

    _out.ResizeTo(v.size());
    for (int i = 0; i < v.size(); i++)
        _out[i] = v[i];
}


template<typename FilterFunction>
void filtermatrix(
    const TMatrixD& _input,
    TMatrixD& _output,
    const FilterFunction _filterFn)
{
    std::vector<Double_t> aux;
    Int_t nrows = 0;

    for (Int_t i = 0; i < _input.GetNrows(); i++)
    {
        if (_filterFn(_input[i]))
        {
            for (Int_t j = 0; j < _input.GetNcols(); j++)
            {
                aux.push_back(_input[i][j]);
            }
            nrows++;
        }
    }

    if (nrows > 0)
    {
        _output.Use(nrows, _input.GetNcols(), aux.data());
    }
    else
    {
        _output.ResizeTo(nrows, _input.GetNcols());
    }
}
