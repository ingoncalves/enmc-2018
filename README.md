# ENERGY ESTIMATION BASED ON WIENER-HOPF FILTERING FOR THE ATLAS TILE CALORIMETER

In high-energy calorimeter systems operating in high-luminosity conditions, the signal pile-up effect is observed. Therefore, the deterministic approximation typically employed for energy estimation, becomes unrealistic in severe signal pile-up conditions. This work presents an approach for energy estimation in the ATLAS experiment at LHC, based on the Wiener-Hopf Filter theory where the data modeling is not required. Simulation data sets were used for design and performance evaluation considering severe signal pile-up conditions. The results showed that the proposed approach of the Wiener-Hopf Filter estimator outperforms the current method used for energy estimation in the ATLAS tile calorimeter, with severe pile-up condictions.

#### Keywords
Signal reconstruction, Wiener-Hopf filters, Optimal filters, High-energy calorimetry.

#### Authors

- Guilherme I. Gonçalves <ggoncalves@iprj.uerj.br>
- Bernardo S. Peralva <bernardo@iprj.uerj.br>
- Luciano M. A. Filho <luciano.andrade@ufjf.edu.br>
- Augusto S. Cerqueira <augusto.santiago@ufjf.edu.br>
- José M. Seixas <seixas@lps.ufrj.br>


XXI ENMC - Encontro Nacional de Modelagem Computacional</br>
8 a 11 de Outubro de 2018 - Instituto Federal Fluminense - Búzios, RJ - Brasil

Read the full article [here](./enmc2018.pdf).

## Replicate results

To replicate our results, follow these steps:

### Requirements

Download and install ROOT https://root.cern.ch/

### Generate Wiener-Hopf Weights

You can calculate **Optimal Wiener-Hopf** weights by:

    root woWeights.C

and also, the **General Wiener-Hopf** wights by:

    root wgWeights.C

### Comparing OF with Wiener-Hopf

You can replicate results of the comparison between OF and Wiener-Hopf methods by:

    root main.C

### Generating figures

You can generate figures of merit by:

    cd ./graphs
    root pileup.C
    root windowXnoise.C
    root bcidHist.C
    root bcidXmean.C
    root bcidXrms.C
