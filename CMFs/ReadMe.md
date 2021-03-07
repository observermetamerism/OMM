*** File description
1. UN_World_Population_2019.xlsx: Population Distribution reported by the United Nation, 2019. For more information: https://population.un.org/wpp/DataQuery/

2. MonteCarlndObsLMS_N=10000.mat: 10,000 individual observers' LMS conefundametals created based on the 2019 UN World Population report. The cone fundamentals range from 390nm to 780nm with 5 nm steps. 

3. MonteCarIndObsAge_N=10000.mat: The 10,000 individual observers' age

4. MonteCarIndObsParams_N=10000.mat: The 10,000 individual observers' physiological parameters. For more information: https://www.rit.edu/cos/colorscience/re_AsanoObserverFunctions.php

5. CatLMS_N=100.mat: 100 Categorical observers (LMS) created from MonteCarlndObsLMS_N=10000.mat using the method proposed by Yuta Asano & Mark Fairchild. For more information: Asano, Y, Fairchild, MD. Categorical observers for metamerism. Color Res Appl. 2020; 45: 576– 585. https://doi.org/10.1002/col.22493

6. Cat_2deg_CMFs_N=100.mat: 100 2-deg Categorical observers (xyz-like CMFs) converted from CatLMS_N=100.mat. For more information: https://www.rit.edu/cos/colorscience/re_AsanoObserverFunctions.php

*** Note that OMM_N_10 & OMM_N_50 are first 10 and 50 observers out of Cat_2deg_CMFs_N=100.mat.