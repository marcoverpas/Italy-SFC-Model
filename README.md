# Empirical SFC macroeconomic model for Italy using [Bimets](https://github.com/cran/bimets) package

**Authors**: Rosa Canelli and Marco Veronese Passarella

**Last change**: 9 January 2024

**Description**: This repository contains two types of code. The file named `Italy_SFC_model.r` replicates the baseline and the experiments discussed in: Canelli, R., Fontana, G., Realfonzo, R. and Veronese Passarella, M. (2022) [Is the Italian government debt sustainable? Scenarios after the Covid-19 shock](https://academic.oup.com/cje/article-abstract/46/3/581/6584486), *Cambridge Journal of Economics*. Notice that an early version of it was used to produce the simulations discussed in: Canelli, R., Fontana, G., Realfonzo, R. and Veronese Passarella, M. (2021) [Are EU Policies Effective to Tackle the Covid-19 Crisis? The Case of Italy](https://www.tandfonline.com/doi/full/10.1080/09538259.2021.1876477), *Review of Political Economy*. The other files allow replicating the baseline and the experiments discussed in: Canelli, R., Fontana, G., Realfonzo, R. and Veronese Passarella, M. (2024) [Energy crisis, economic growth and public finance in Italy](https://www.sciencedirect.com/journal/energy-economics), *Energy Economics*. A detailed description of the latter is provided below.


### 1_Model_and_data

The first step to create the model is to define the structure of the economy, that is, the level of aggregation implied by the balance sheet and the transactions-flow matrix. [...]

$$
\begin{array}{l|r|r|r|r|r|r|r}
\hline
  & Households & Firms & Government & Banks & ECB & Foreign & Total\\
\hline
Cash and reserves & 200683 & 0 & 0 & 10817 & -211500 & 0 & 0\\
\hline
Deposits & 1428434 & 0 & 0 & -1428434 & 0 & 0 & 0\\
\hline
Securities & 233263 & 0 & -2678397 & 1366294 & 868289 & 210551 & 0\\
\hline
Loans & -763488 & -871902 & 0 & 1635390 & 0 & 0 & 0\\
\hline
Shares & 1372850 & -1372850 & 0 & 0 & 0 & 0 & 0\\
\hline
Other net FA & 1583746 & 284629 & 323282 & -1563895 & -783662 & 155900 & 0\\
\hline
Net financial wealth & 4055488 & -1960123 & -2355115 & 20172 & -126873 & 366451 & 0\\
\hline
Column total & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\hline
\end{array}
$$


The second step is to upload data. This code takes the observed series from a Dropbox folder, containing Eurostat data for Italy over the period 1995-2021.

```R
#Upload data: time series for transactions-flow matrix and balance sheet
DataEE <- read.csv("https://www.dropbox.com/scl/fi/58un32aub0asx1tu1t53c/Data_Aalborg.csv?rlkey=642w2hwa8qo069nmn1rc848i2&dl=1") 
```

Running this chunk of code allows creating the dataset used by the model to estimate the coefficients and attribute initial values to the variables.

| Year | CONS    | INV    | GOV    | X      | IM     | Y       | TAX         | INVnet      | TR     | WB     | INTh  | INTlh       | INTgh       | INTm         | INTf  | INTg        | INTgb        | INTb        | INTcb        | INTrow | Ff          | DIV         | FUf    | Fb          | Fcb          | YD      | OPh          | OPf         | OPg         | OPb         | OPcb   | OProw        | Yrow        | DEF1        | Hh     | Hbd   | Hs      | Mh      | Ms       | Bh     | Deb         | Bb      | Bcb    | Bs           | Brow       | Lf       | Lh      | Ls      | Eh      | Es       | OAh     | OAf    | OAg        | OAb      | OAcb    | OArow      | Vh      | NVh     | Vf       | Vg       | Vb      | Vcb     | Vrow    | K           | Rlh         | Rl          | Rb          | Rstar       | Prod        | ProdR       | Delta | Kappa       | KappaN      | Py      | Pc    | Exr   | Exr42 | Exr42R    | IMen      | Pen         | Pen2        | Pim         | Pim_en      | Py_row      | RHO         | w           | gid          | gkn         | gy           | Rbb          | Rbrow       | Rbh         | Re          | gw           | CPI          | PI          | PIen         | Nd      | Ns      |
|------|---------|--------|--------|--------|--------|---------|-------------|-------------|--------|--------|-------|-------------|-------------|--------------|-------|-------------|--------------|-------------|--------------|--------|-------------|-------------|--------|-------------|--------------|---------|--------------|-------------|-------------|-------------|--------|--------------|-------------|-------------|--------|-------|---------|---------|----------|--------|-------------|---------|--------|--------------|------------|----------|---------|---------|---------|----------|---------|--------|------------|----------|---------|------------|---------|---------|----------|----------|---------|---------|---------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------|-------------|-------------|---------|-------|-------|-------|-----------|-----------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|--------------|-------------|--------------|--------------|-------------|-------------|-------------|--------------|--------------|-------------|--------------|---------|---------|
| 1995 | 581704  | 198068 | 172298 | 243824 | 207651 | 988243  | 304502.9179 | 61690.466   | 32837  | 353786 | 94539 | 22058.98524 | 202553.2366 | -85955.25137 | 32051 | 101685.9999 | -29707.98536 | 110357.2512 | -85955.25137 | 14796  | 594208      | 435525      | 158683 | 110357.2512 | -85955.25137 | 709948  | -12593.3333  | -8198       | -26684.6667 | 108193.2514 | 9950   | -70667.25137 | 11869868.24 | 88273.33333 | 42951  | 3915  | -46866  | 642843  | -642843  | 407396 | 1179962.142 | 390724  | 106668 | -1179962.142 | 275174.142 | -480014  | -166659 | 646673  | 337925  | -337925  | 276065  | 104061 | 239854.142 | -483688  | -80593  | -55699.142 | 1707180 | 1540521 | -713878  | -940108  | -85219  | -20791  | 219475  | 2964729     | 0.13236     | 0.066770969 | 0.086177341 | 0.086875    | 49.59938769 | 75.25206367 | 0.046 | 0.333333333 | 0.34940601  | 65.911  | 65.5  | 1.308 | 77.5  | 92.40989  | 15433.334 | 49.50860467 | 56.1921773  | 64.63870703 | 51.80703164 | 55.27815729 | 0.006090134 | 17.75633015 | 0.037188395  | 0.023283628 | 0.051381348  | -0.076033173 | 0.053769587 | 0.497190048 | 1.288821484 | 0.074269715  | 0.045801527  | 0.045075936 | 0            | 19924.5 | 22560.6 |
| 1996 | 611812  | 202429 | 184246 | 247921 | 200535 | 1045873 | 339106.6832 | 66051.466   | 32914  | 378499 | 96012 | 21362.90615 | 196733.19   | -79358.28384 | 28312 | 108261.0004 | -21732.90575 | 107300.2842 | -79358.28384 | 12619  | 677120      | 514076      | 163044 | 107300.2842 | -79358.28384 | 763826  | -25868.601   | 38058       | -49285.399  | -15041      | 976    | 51161        | 12432671.32 | 65672.601   | 43739  | 4488  | -48227  | 666811  | -666811  | 490407 | 1245634.743 | 439727  | 97545  | -1245634.743 | 217955.743 | -482938  | -176991 | 659929  | 370549  | -370549  | 298020  | 100224 | 190568.743 | -537593  | -69133  | 17913.257  | 1869526 | 1692535 | -753263  | -1055066 | -100260 | -19815  | 235869  | 3030780.466 | 0.128183333 | 0.058624503 | 0.091749554 | 0.084375    | 52.20568342 | 75.79002268 | 0.046 | 0.345083721 | 0.369781837 | 68.882  | 68.5  | 1.27  | 86.1  | 101.71875 | 19682.28  | 51.77617435 | 58.72509228 | 68.18104965 | 54.64617037 | 58.3075213  | 0.006730543 | 18.8931151  | 0.022017691  | 0.022279091 | 0.058315617  | -0.055622142 | 0.045858233 | 0.482904079 | 1.521272472 | 0.064021391  | 0.045801527  | 0.045075936 | 0.045801527  | 20033.7 | 22737.8 |
| 1997 | 645631  | 213028 | 193160 | 263802 | 223262 | 1092359 | 348240.7178 | 73612.09856 | 35143  | 397066 | 79460 | 18600.83965 | 153836.1258 | -55775.28619 | 28260 | 94403.00059 | -10613.83906 | 92022.28678 | -55775.28619 | 6956   | 701730      | 567744      | 133986 | 92022.28678 | -55775.28619 | 914897  | 91702.431    | 34697       | -29706.431  | -41340      | 407    | -55760       | 13182538.9  | 30240.569   | 46923  | 4901  | -51824  | 638201  | -638201  | 559168 | 1275875.312 | 441022  | 92006  | -1275875.312 | 183679.312 | -493212  | -188351 | 681563  | 448317  | -448317  | 457543  | 109224 | 160862.312 | -630885  | -59590  | -37154.312 | 2150152 | 1961801 | -832305  | -1115013 | -141600 | -19408  | 146525  | 3104392.565 | 0.105094833 | 0.057297876 | 0.075787064 | 0.063958333 | 54.29705441 | 76.8524924  | 0.046 | 0.351875279 | 0.377800508 | 70.651  | 69.8  | 1.134 | 87.5  | 100.6284  | 20382.45  | 52.75878788 | 60.23324663 | 71.82719165 | 57.56850286 | 61.42565315 | 0.007679399 | 19.73665636 | 0.052359099  | 0.024288166 | 0.044447079  | -0.024137338 | 0.031914736 | 0.313690722 | 1.532169834 | 0.044648077  | 0.018978102  | 0.0256816   | 0.018978102  | 20118.2 | 22865.7 |
| 1998 | 679941  | 224747 | 199624 | 274042 | 239497 | 1138857 | 363999.2758 | 81944.94203 | 58835  | 395312 | 62467 | 16279.38098 | 122672.1812 | -43925.80025 | 15700 | 85175.00053 | -1030.380451 | 74874.80078 | -43925.80025 | 7459   | 835871      | 745987      | 89884  | 74874.80078 | -43925.80025 | 958284  | -15192.525   | 108026      | -43185.475  | -124795     | -1730  | 76877        | 13879204.32 | 23560.525   | 50305  | 5521  | -55826  | 609760  | -609760  | 513956 | 1299435.837 | 577525  | 78127  | -1299435.837 | 129827.837 | -504062  | -205507 | 709569  | 562925  | -562925  | 708705  | 99819  | 117676.837 | -949250  | -43439  | 66488.163  | 2445651 | 2240144 | -967168  | -1181759 | -266395 | -21138  | 196316  | 3186337.507 | 0.086431083 | 0.031146962 | 0.066758091 | 0.047916667 | 56.0861343  | 77.52271562 | 0.046 | 0.357418823 | 0.384542359 | 72.348  | 71.2  | 1.121 | 90.4  | 102.12797 | 12666.42  | 81.00496309 | 61.68001765 | 74.8634512  | 60.00202299 | 64.02222168 | 0.009054382 | 19.4682229  | 0.055011548  | 0.02639645  | 0.042566592  | -0.002336347 | 0.040608819 | 0.219383408 | 1.663972145 | -0.013600757 | 0.020057307  | 0.024019476 | 0.535383324  | 20305.5 | 23111.2 |
| 1999 | 710248  | 237770 | 206206 | 272643 | 251716 | 1175151 | 356021.3367 | 91198.4747  | 62587  | 410340 | 54993 | 13040.49806 | 115548.0802 | -47514.58213 | 13073 | 71724.00058 | -5035.497487 | 68592.58271 | -47514.58213 | 8726   | 560877      | 558925      | 1952   | 68592.58271 | -47514.58213 | 926019  | 126602.754   | -190861     | 58330.246   | 12396       | -1685  | -4783        | 14675566.04 | 32010.246   | 56333  | 5795  | -62128  | 610642  | -610642  | 438251 | 1331446.083 | 610011  | 75770  | -1331446.083 | 207414.083 | -572790  | -239433 | 812223  | 750109  | -750109  | 840013  | 119913 | 176007.083 | -1071386 | -36465  | -28082.083 | 2695348 | 2455915 | -1202986 | -1155439 | -253999 | -22823  | 179332  | 3277535.981 | 0.06345525  | 0.022823373 | 0.055196262 | 0.027083333 | 57.17411294 | 77.83026537 | 0.046 | 0.358547094 | 0.386592589 | 73.46   | 72.3  | 1.066 | 91.8  | 101.55784 | 16561.2   | 92.65119048 | 62.62804911 | 77.66434745 | 62.24690268 | 66.41751068 | 0.009490012 | 19.96409441 | 0.057945156  | 0.028621725 | 0.031868795  | -0.008719099 | 0.067212088 | 0.224820958 | 0.992894258 | 0.025470815  | 0.015449438  | 0.015370155 | 0.143771776  | 20553.9 | 23279.3 |
| 2000 | 751730  | 259037 | 220195 | 318172 | 307621 | 1241513 | 375224.8494 | 108270.3449 | 66311  | 427487 | 53814 | 16814.76105 | 108571.6979 | -37942.93688 | 17088 | 72578.99948 | -8496.761572 | 63348.93636 | -37942.93688 | 10447  | 780481      | 560823      | 219658 | 63348.93636 | -37942.93688 | 925612  | 129052.913   | -16457      | 14867.087   | -241870     | -274   | 114681       | 17137543.71 | 21803.087   | 59394  | 6644  | -66038  | 631780  | -631780  | 503331 | 1353249.17  | 557289  | 85497  | -1353249.17  | 207132.17  | -641875  | -273340 | 915215  | 747248  | -747248  | 961384  | 146758 | 190874.17  | -1343237 | -42556  | 86776.83   | 2903137 | 2629797 | -1242365 | -1162375 | -495869 | -23097  | 293909  | 3385806.326 | 0.070227417 | 0.026622006 | 0.054511407 | 0.040416667 | 59.46968826 | 79.53044862 | 0.046 | 0.366681635 | 0.397059351 | 74.776  | 74.2  | 0.924 | 88    | 96.48973  | 34898.31  | 94.79664835 | 63.75       | 79.69396408 | 63.8736124  | 68.15321166 | 0.010516319 | 20.47704585 | 0.08944358   | 0.033034068 | 0.056471041  | -0.013928866 | 0.050367843 | 0.24773862  | 0.747655341 | 0.0256937    | 0.026279391  | 0.017914511 | 0.02315629   | 20876.4 | 23418.8 |
| 2001 | 775878  | 271825 | 239959 | 334498 | 318022 | 1304138 | 343517.9719 | 116077.909  | 68186  | 451041 | 52025 | 19920.33585 | 99117.42021 | -27172.08436 | 16530 | 75157.99956 | -7151.336287 | 56471.08392 | -27172.08436 | 10364  | 821681      | 614287      | 207394 | 56471.08392 | -27172.08436 | 720897  | -177595.112  | -14886      | 854.112     | 192428      | -1733  | 932          | 18063016.94 | 66957.112   | 49167  | 8156  | -57323  | 690239  | -690239  | 566966 | 1420206.282 | 534466  | 86879  | -1420206.282 | 231895.282 | -688808  | -303134 | 991942  | 736025  | -736025  | 835553  | 118037 | 191728.282 | -1147766 | -54386  | 56833.718  | 2877950 | 2574816 | -1306796 | -1228478 | -303441 | -24830  | 288729  | 3501884.235 | 0.0728775   | 0.023997979 | 0.055538922 | 0.0425      | 61.21133041 | 79.44880318 | 0.046 | 0.372410369 | 0.40375049  | 77.045  | 75.9  | 0.896 | 89.9  | 97.08403  | 34288.28  | 98.23146825 | 60.81       | 83.42047245 | 66.86035743 | 71.34007175 | 0.011816197 | 21.17016733 | 0.049367465  | 0.034283682 | 0.050442484  | -0.012832366 | 0.050035685 | 0.19692294  | 0.822065767 | 0.033848705  | 0.022911051  | 0.030343961 | 0.036233558  | 21305.5 | 23573.9 |
| 2002 | 798351  | 291267 | 250419 | 329564 | 319342 | 1350259 | 391298.2531 | 130180.3252 | 71440  | 470108 | 46028 | 19814.37917 | 81591.92619 | -15749.54702 | 17230 | 70159.00011 | -9526.379064 | 43267.54713 | -15749.54702 | 13843  | 811498      | 594283      | 217215 | 43267.54713 | -15749.54702 | 858604  | 24775.706    | -51423      | -21476.706  | -59291      | -933   | 108348       | 17984691    | 16469.294   | 53408  | 8970  | -62378  | 716998  | -716998  | 637330 | 1436675.576 | 511446  | 65401  | -1436675.576 | 222498.576 | -716643  | -336487 | 1053130 | 777177  | -777177  | 786643  | 112972 | 170251.576 | -1219280 | -28786  | 178199.424 | 2971556 | 2635069 | -1380848 | -1266424 | -362732 | -25763  | 400698  | 3632064.56  | 0.065365083 | 0.024042654 | 0.04940057  | 0.032083333 | 62.23624958 | 78.21867053 | 0.046 | 0.347       | 0.404172649 | 79.567  | 77.9  | 0.946 | 92.2  | 99.27912  | 32653.125 | 105.7765564 | 63.15       | 86.39226109 | 69.24220501 | 73.88150563 | 0.012510495 | 21.66825684 | 0.071523958  | 0.037174366 | 0.035365122  | -0.017824107 | 0.059695048 | 0.143909734 | 0.807422302 | 0.023527897  | 0.026350461  | 0.032734116 | 0.076809278  | 21695.7 | 23902.1 |
| 2003 | 828304  | 295388 | 264649 | 324994 | 318642 | 1394693 | 393258.1745 | 128313.0302 | 74195  | 490889 | 47510 | 19618.3698  | 90145.08292 | -23016.71312 | 14382 | 66123.00043 | -14327.36938 | 42689.71355 | -23016.71312 | 13322  | 990303      | 724781      | 265522 | 42689.71355 | -23016.71312 | 890185  | -96621.539   | 100881      | 29264.539   | -78778      | -279   | 45533        | 16579854.78 | 34725.539   | 64267  | 8975  | -73242  | 754034  | -754034  | 642341 | 1471401.115 | 511490  | 68819  | -1471401.115 | 248751.115 | -773097  | -371352 | 1144449 | 747473  | -747473  | 860187  | 109856 | 199516.115 | -1352390 | -21619  | 204449.885 | 3068302 | 2696950 | -1410714 | -1271885 | -441510 | -26042  | 453201  | 3760377.591 | 0.0583035   | 0.018603099 | 0.046025005 | 0.0225      | 63.4081062  | 77.25912151 | 0.046 | 0.370891743 | 0.402510012 | 82.072  | 80.1  | 1.131 | 98.4  | 105.94827 | 34191.43  | 105.3707108 | 64.38       | 89.12955539 | 71.43610862 | 76.22240308 | 0.011902646 | 22.31770135 | 0.01414853   | 0.03532785  | 0.032907761  | -0.028013455 | 0.059874541 | 0.141441769 | 0.932581638 | 0.029972162  | 0.028241335  | 0.031482901 | -0.00383682  | 21995.5 | 24141.3 |
| 2004 | 856707  | 308799 | 278255 | 348699 | 340141 | 1452319 | 380777.1062 | 135821.6308 | 74916  | 509540 | 46295 | 20447.26004 | 84747.52034 | -18005.2603  | 14173 | 64586.99987 | -14956.26017 | 37669.26017 | -18005.2603  | 12801  | 1049131     | 882027      | 167104 | 37669.26017 | -18005.2603  | 1041694 | -127976.154  | 120525      | -14715.846  | -22910      | 694    | 44383        | 16589586.29 | 54986.154   | 76070  | 8639  | -84709  | 784680  | -784680  | 718263 | 1526387.269 | 481986  | 74525  | -1526387.269 | 251613.269 | -806836  | -421412 | 1228248 | 857892  | -857892  | 866444  | 112319 | 184800.269 | -1398613 | -15164  | 250213.731 | 3303349 | 2881937 | -1552409 | -1341587 | -464420 | -25348  | 501827  | 3896199.221 | 0.055061667 | 0.017566147 | 0.043894897 | 0.02        | 64.93945681 | 77.06758223 | 0.046 | 0.372752757 | 0.404838856 | 84.263  | 81.9  | 1.244 | 100.4 | 107.67591 | 38846.224 | 96.41309005 | 68.91       | 91.86997485 | 73.632517   | 78.56597314 | 0.011009584 | 22.78373472 | 0.045401303  | 0.036119147 | 0.041318054  | -0.029240572 | 0.051461076 | 0.131935406 | 1.180011853 | 0.020881782  | 0.02247191   | 0.026696072 | -0.085010537 | 22364.2 | 24277.5 |
| 2005 | 885949  | 316194 | 292969 | 367451 | 368928 | 1493635 | 388702.5627 | 136968.8358 | 76130  | 531929 | 47503 | 22393.55976 | 90820.76384 | -20924.20408 | 15074 | 64506.99964 | -16365.56012 | 42026.20372 | -20924.20408 | 10976  | 970860      | 915219      | 55641  | 42026.20372 | -20924.20408 | 1219425 | -4679.641    | 24228       | 9482.641    | 14648       | -8930  | -34749       | 17493109.5  | 65827.641   | 89146  | 8737  | -97883  | 822325  | -822325  | 735715 | 1592214.91  | 499382  | 82215  | -1592214.91  | 274902.91  | -876585  | -479078 | 1355663 | 1057758 | -1057758 | 989547  | 121381 | 194282.91  | -1491229 | -18610  | 204628.09  | 3694491 | 3215413 | -1812962 | -1397932 | -449772 | -34278  | 479531  | 4033168.057 | 0.05313935  | 0.017196279 | 0.042261228 | 0.020208333 | 66.79493773 | 77.70648192 | 0.046 | 0.370337903 | 0.401841653 | 85.958  | 83.7  | 1.244 | 99.3  | 105.95688 | 47940.26  | 75.45       | 76.69       | 94.2        | 75.5        | 80.55857946 | 0.010624753 | 23.78771549 | 0.023947616  | 0.035154474 | 0.028448295  | -0.03395443  | 0.043622501 | 0.126444998 | 1.066823097 | 0.04406568   | 0.021978022  | 0.02011559  | -0.217429916 | 22361.5 | 24238.8 |
| 2006 | 921983  | 340933 | 302522 | 406392 | 419143 | 1552687 | 412935.4616 | 155407.2694 | 84673  | 555830 | 50027 | 26932.7072  | 101583.2882 | -24623.58098 | 19164 | 65618.99958 | -24421.70762 | 46298.58056 | -24623.58098 | 13081  | 681954      | 645212      | 36742  | 46298.58056 | -24623.58098 | 1196782 | 227676.881   | -295739     | 49413.119   | -108750     | -3773  | 131172       | 18391491.59 | 64502.119   | 103968 | 10322 | -114290 | 879188  | -879188  | 735374 | 1656717.029 | 487660  | 93789  | -1656717.029 | 339894.029 | -954392  | -541251 | 1495643 | 1289471 | -1289471 | 1023462 | 126710 | 243696.029 | -1672959 | -17550  | 296640.971 | 4031463 | 3490212 | -2117153 | -1413021 | -558522 | -38051  | 636535  | 4188575.327 | 0.056217792 | 0.020079799 | 0.041212401 | 0.027916667 | 68.3704904  | 77.88491115 | 0.046 | 0.370695733 | 0.403542447 | 87.784  | 85.6  | 1.256 | 99.3  | 105.4757  | 61681.2   | 95.70833333 | 78.98       | 98.4        | 95.7        | 82.64410275 | 0.011740379 | 24.47522886 | 0.078239941  | 0.038532307 | 0.039535763  | -0.04890386  | 0.047584073 | 0.138074238 | 0.609980733 | 0.028902034  | 0.022700119  | 0.021242933 | 0.26850011   | 22709.9 | 24364.2 |
| 2007 | 953816  | 359238 | 307266 | 441837 | 447317 | 1614840 | 465321.4558 | 166563.535  | 78465  | 576340 | 56865 | 34286.6271  | 118836.814  | -27685.1869  | 30044 | 73006.99985 | -45126.62725 | 46889.18675 | -27685.1869  | 26982  | 998664      | 657999      | 340665 | 46889.18675 | -27685.1869  | 679025  | -272211.731  | -9792       | 26685.731   | 144725      | -6743  | 117336       | 18332622.93 | 21101.731   | 112983 | 11364 | -124347 | 915020  | -915020  | 762538 | 1677818.76  | 462908  | 93824  | -1677818.76  | 358548.76  | -1071324 | -602022 | 1673346 | 1199296 | -1199296 | 827606  | 134894 | 270381.76  | -1646395 | -14271  | 427784.24  | 3817443 | 3215421 | -2135726 | -1407437 | -413797 | -44794  | 786333  | 4355138.862 | 0.063347    | 0.028043804 | 0.044067272 | 0.038541667 | 70.68339892 | 78.57203081 | 0.046 | 0.370789555 | 0.40412414  | 89.96   | 87.3  | 1.371 | 101   | 106.39722 | 59734.24  | 94.025      | 86.23       | 98.9        | 94          | 85.39684094 | 0.012419401 | 25.22706282 | 0.053690901  | 0.039766155 | 0.040029317  | -0.092537069 | 0.079383566 | 0.161600511 | 0.510286001 | 0.030718159  | 0.019859813  | 0.024788116 | -0.017588158 | 22846.1 | 24327   |
| 2008 | 973266  | 356666 | 320333 | 439932 | 452498 | 1637699 | 457702.286  | 156329.6124 | 81770  | 595572 | 65938 | 41162.11041 | 147428.9744 | -40328.86402 | 39097 | 76688.00002 | -57908.11039 | 62679.86404 | -40328.86402 | 27496  | 915292      | 562142      | 353150 | 62679.86404 | -40328.86402 | 760338  | -150061.578  | -87738      | -10330.422  | 290011      | -4201  | -37680       | 17674323.9  | 61417.578   | 123257 | 13153 | -136410 | 974991  | -974991  | 809949 | 1739236.338 | 472953  | 105039 | -1739236.338 | 351295.338 | -1134965 | -624988 | 1759953 | 1141864 | -1141864 | 577420  | 137587 | 260051.338 | -1394854 | -17624  | 437419.662 | 3627481 | 3002493 | -2139242 | -1479185 | -123786 | -48995  | 788715  | 4511468.474 | 0.0683731   | 0.034447758 | 0.045706963 | 0.038541667 | 71.08346246 | 77.16398444 | 0.046 | 0.363007967 | 0.394170122 | 92.12   | 90.4  | 1.471 | 103.3 | 107.93314 | 76410     | 122.7666667 | 83.27       | 103.6       | 122.8       | 88.16565515 | 0.013490381 | 25.85048895 | -0.007159599 | 0.035895437 | 0.014155582  | -0.12509637  | 0.07668692  | 0.193339839 | 0.468726653 | 0.024712593  | 0.035509737  | 0.024010671 | 0.305681114  | 23039.1 | 24703.4 |
| 2009 | 953184  | 307703 | 326155 | 353292 | 363078 | 1577256 | 410194.439  | 100175.4502 | 99915  | 590431 | 45050 | 29728.59587 | 92382.19295 | -17603.59708 | 17698 | 66364.99993 | -29036.59594 | 35993.59701 | -17603.59708 | 20623  | 1017050     | 564863      | 452187 | 35993.59701 | -17603.59708 | 905176  | -20882.158   | 47923       | -17766.842  | 65986       | -11416 | -63844       | 17578903.82 | 99844.158   | 129192 | 11664 | -140856 | 994220  | -994220  | 779292 | 1839080.496 | 609381  | 119299 | -1839080.496 | 331108.496 | -1146281 | -654677 | 1800958 | 997363  | -997363  | 709095  | 148886 | 242284.496 | -1485583 | -38854  | 424171.504 | 3609162 | 2954485 | -1994758 | -1596796 | -57800  | -60411  | 755280  | 4611643.924 | 0.047566667 | 0.015439495 | 0.038157551 | 0.012291667 | 69.6378712  | 74.34701087 | 0.046 | 0.342015998 | 0.366467855 | 93.666  | 91.1  | 1.395 | 105.1 | 109.31352 | 52081.4   | 88.21666667 | 84.81       | 94.6        | 88.2        | 90.05561801 | 0.01173181  | 26.06828437 | -0.137279696 | 0.022204622 | -0.036907271 | -0.061394253 | 0.05870559  | 0.114059272 | 0.494685006 | 0.008425195  | 0.007743363  | 0.016782458 | -0.281428184 | 22649.4 | 24556   |
| 2010 | 978454  | 331599 | 331166 | 404013 | 433952 | 1611280 | 421807.5657 | 119463.3795 | 107096 | 599034 | 34172 | 26394.39438 | 59409.22375 | 1157.170637  | 10465 | 66268.00032 | -14186.39406 | 21515.82969 | 1157.170637  | 19888  | 1013316     | 613378      | 399938 | 21515.82969 | 1157.170637  | 880005  | -73383.264   | 11535       | 71182.264   | 4682        | -22786 | 8770         | 18852954.48 | 81565.264   | 133437 | 11172 | -144609 | 988877  | -988877  | 732051 | 1920645.76  | 638344  | 151353 | -1920645.76  | 398897.76  | -1154988 | -686301 | 1841289 | 915991  | -915991  | 771981  | 144560 | 313466.76  | -1555046 | -89941  | 414979.24  | 3542337 | 2856036 | -1926419 | -1607179 | -53118  | -83197  | 813877  | 4731107.304 | 0.040316667 | 0.0090607   | 0.036033224 | 0.01        | 71.70039693 | 76.21620721 | 0.046 | 0.340571434 | 0.366240927 | 94.075  | 92.6  | 1.326 | 101   | 104.33731 | 66130.02  | 105.6083333 | 94.19       | 99.5        | 105.6       | 91.06593797 | 0.011297664 | 26.65643189 | 0.077659301  | 0.025904728 | 0.021571641  | -0.023280007 | 0.060064904 | 0.076234869 | 0.614999754 | 0.022561804  | 0.016465423  | 0.004366579 | 0.197147175  | 22472.4 | 24528.1 |
| 2011 | 1007661 | 337468 | 326718 | 443061 | 466151 | 1648757 | 454145.6216 | 119837.064  | 110671 | 608088 | 38915 | 31564.12683 | 66948.11742 | 3531.009403  | 12701 | 73204       | -20791.12682 | 19942.9906  | 3531.009403  | 23516  | 1174074     | 792739      | 381335 | 19942.9906  | 3531.009403  | 917706  | -198504.369  | 146106      | 105275.369  | 111366      | -12726 | -151517      | 18926146.73 | 52916.369   | 141838 | 11535 | -153373 | 984230  | -984230  | 743323 | 1973562.129 | 523412  | 176406 | -1973562.129 | 530421.129 | -1175094 | -704239 | 1879333 | 844950  | -844950  | 755979  | 137492 | 418742.129 | -1371802 | -118956 | 178544.871 | 3470320 | 2766081 | -1882552 | -1554820 | 58248   | -95923  | 708966  | 4850944.368 | 0.045991667 | 0.010808497 | 0.038114264 | 0.0125      | 73.13020337 | 76.50643222 | 0.046 | 0.339883716 | 0.365296473 | 95.587  | 95.3  | 1.392 | 101.6 | 104.48254 | 76271.13  | 133.8916667 | 103.32      | 107.7       | 133.9       | 92.94986771 | 0.011719822 | 26.97159078 | 0.017699088  | 0.025329602 | 0.023259148  | -0.032570412 | 0.058952449 | 0.091452805 | 0.865444093 | 0.011822996  | 0.029157667  | 0.016072283 | 0.267813462  | 22545.5 | 24606.8 |
| 2012 | 995710  | 288966 | 321754 | 460981 | 443052 | 1624359 | 436848.3013 | 65822.55909 | 120511 | 600264 | 37724 | 36778.88178 | 69509.18801 | 4993.693761  | 12686 | 80829.00101 | -16705.88076 | 27765.30725 | 4993.693761  | 23032  | 1032254     | 808402      | 223852 | 27765.30725 | 4993.693761  | 1280778 | 122959.994   | 20845       | -161315.994 | -83650      | -3494  | 104655       | 20616647.33 | 81252.006   | 137756 | 13380 | -151136 | 1046373 | -1046373 | 736578 | 2054814.135 | 746424  | 192181 | -2054814.135 | 379631.135 | -1161701 | -697100 | 1858801 | 912080  | -912080  | 915462  | 126115 | 257426.135 | -1597634 | -140462 | 439092.865 | 3748249 | 3051149 | -1947666 | -1797388 | -25402  | -99417  | 818724  | 4916766.927 | 0.052225    | 0.010920194 | 0.040955894 | 0.00875     | 72.15909697 | 74.3402396  | 0.046 | 0.330371365 | 0.351000189 | 97.066  | 98.4  | 1.285 | 99    | 102.51107 | 83664.24  | 150.1916667 | 108.71      | 111         | 150.2       | 93.74142715 | 0.012787027 | 26.66560051 | -0.143723257 | 0.01356902  | -0.014797814 | -0.031917267 | 0.043422101 | 0.093511418 | 0.95674537  | -0.01134491  | 0.032528856  | 0.015472815 | 0.121740213  | 22510.8 | 25201.8 |
| 2013 | 982190  | 272433 | 319441 | 461783 | 423095 | 1612752 | 455823.3594 | 46261.72137 | 131490 | 594339 | 37854 | 35859.98583 | 85903.61056 | -12189.62473 | 11661 | 74784.99964 | -23152.9862  | 36557.62437 | -12189.62473 | 24224  | 910531      | 727101      | 183430 | 36557.62437 | -12189.62473 | 1111389 | 39870.735    | -96221      | -5962.735   | -118788     | 30741  | 150360       | 20727022.91 | 82082.265   | 139992 | 13043 | -153035 | 1076752 | -1076752 | 643794 | 2136896.4   | 965725  | 187388 | -2136896.4   | 339989.4   | -1106952 | -684676 | 1791628 | 1036255 | -1036255 | 968231  | 106538 | 251463.4   | -1837834 | -103029 | 614630.6   | 3865024 | 3180348 | -2036669 | -1885433 | -144190 | -68676  | 954620  | 4963028.648 | 0.051441667 | 0.010534332 | 0.036395019 | 0.005416667 | 72.85619418 | 74.20675716 | 0.046 | 0.324953192 | 0.343826695 | 98.18   | 99.7  | 1.328 | 101.7 | 104.41029 | 72200.4   | 141.8083333 | 108.7       | 108.3       | 141.8       | 95.63143633 | 0.012113281 | 26.8493095  | -0.057214344 | 0.009408972 | -0.007145588 | -0.031018545 | 0.063809308 | 0.116625273 | 0.797189939 | 0.006889362  | 0.013211382  | 0.011476727 | -0.055817567 | 22136.1 | 25204.8 |
| 2014 | 986310  | 275995 | 317979 | 473719 | 426597 | 1627406 | 487469.4046 | 47695.68219 | 142072 | 595991 | 35086 | 33320.89867 | 91045.6278  | -22638.72913 | 8050  | 71390.99944 | -21486.89923 | 42522.72857 | -22638.72913 | 24471  | 1090953     | 814731      | 276222 | 42522.72857 | -22638.72913 | 1135780 | -7153.324    | 67588       | -144554.676 | -39772      | -9189  | 133081       | 21619705.53 | 66611.324   | 141353 | 12882 | -154235 | 1110770 | -1110770 | 536907 | 2203507.724 | 1273284 | 202850 | -2203507.724 | 190466.724 | -1082212 | -679562 | 1761774 | 1057339 | -1057339 | 1163011 | 103109 | 106908.724 | -2121132 | -126480 | 874583.276 | 4009380 | 3329818 | -2036442 | -2096599 | -183962 | -77865  | 1065050 | 5010724.33  | 0.048666667 | 0.007438469 | 0.033408732 | 0.001583333 | 73.27027149 | 73.95285635 | 0.046 | 0.324784581 | 0.343716797 | 99.077  | 99.9  | 1.329 | 103.4 | 104.6746  | 51756.01  | 130.5583333 | 101.85      | 104.9       | 130.6       | 97.87056077 | 0.01159736  | 26.83314574 | 0.013074774  | 0.009610197 | 0.009086332  | -0.022249501 | 0.071975773 | 0.141420435 | 0.786226363 | -0.000602018 | 0.002006018  | 0.00913628  | -0.079332432 | 22211   | 25447   |
| 2015 | 1005936 | 283186 | 316344 | 491905 | 442016 | 1655355 | 506679.773  | 52692.68081 | 141366 | 609506 | 28901 | 28060.81389 | 76995.17867 | -20033.36478 | 5335  | 65123.99918 | -13475.8147  | 39953.36396 | -20033.36478 | 21638  | 928480      | 694646      | 233834 | 39953.36396 | -20033.36478 | 1134133 | 126440.409   | -112034     | -22738.409  | 9458        | 951    | -2077        | 25354512.33 | 36187.591   | 144595 | 11845 | -156440 | 1137500 | -1137500 | 420995 | 2239695.315 | 1361780 | 289731 | -2239695.315 | 167189.315 | -1059755 | -679918 | 1739673 | 1154348 | -1154348 | 1280495 | 128309 | 84170.315  | -2150302 | -210205 | 867532.685 | 4137933 | 3458015 | -2085794 | -2155525 | -174504 | -76914  | 1034722 | 5063417.011 | 0.0412925   | 0.005034182 | 0.029554695 | 0.0005      | 73.91660601 | 73.91660601 | 0.046 | 0.326924485 | 0.346291842 | 100     | 100   | 1.11  | 100   | 100       | 44458.08  | 100.0083333 | 95.9        | 100         | 100         | 100         | 0.010413187 | 27.21628585 | 0.02605482   | 0.010515981 | 0.017173957  | -0.010583511 | 0.113605146 | 0.143405056 | 0.656975672 | 0.014278613  | 0.001001001  | 0.009315987 | -0.233995021 | 22394.9 | 25428.2 |
| 2016 | 1019579 | 297798 | 322650 | 497339 | 441578 | 1695788 | 492490.8411 | 64880.81749 | 145038 | 625000 | 22308 | 23771.06648 | 53680.81605 | -7601.749576 | 5205  | 63428.00057 | 1081.934091  | 37659.75014 | -7601.749576 | 16267  | 1119880     | 752840      | 367040 | 37659.75014 | -7601.749576 | 1023768 | -66586.909   | 54297       | 41327.909   | 150612      | -9643  | -170007      | 25801943.97 | 46226.909   | 150624 | 11441 | -162065 | 1182519 | -1182519 | 363801 | 2285922.224 | 1315743 | 417139 | -2285922.224 | 189239.224 | -1039509 | -686157 | 1725666 | 1126524 | -1126524 | 1324893 | 149481 | 125498.224 | -1894223 | -341631 | 635981.776 | 4148361 | 3462204 | -2016552 | -2160424 | -23892  | -86557  | 825221  | 5128297.829 | 0.034961667 | 0.005007172 | 0.028319924 | 0.0000833   | 74.77711781 | 73.93865348 | 0.046 | 0.330672683 | 0.351058495 | 101.134 | 99.9  | 1.107 | 102.1 | 100.84031 | 36762.5   | 80.675      | 98.35       | 96          | 80.7        | 102.3973588 | 0.009675109 | 27.55987106 | 0.051598596  | 0.012813643 | 0.024425576  | 0.0007945    | 0.097296888 | 0.127509391 | 0.65217768  | 0.012624251  | -0.001       | 0.01134     | -0.193317223 | 22677.9 | 25689.9 |
| 2017 | 1046342 | 313526 | 327002 | 533720 | 483997 | 1736593 | 493832.7575 | 77624.29988 | 142440 | 642230 | 21320 | 20590.42798 | 46722.17343 | -4811.745452 | 3572  | 62428.001   | 8671.573029  | 37645.74646 | -4811.745452 | 11846  | 1006081     | 689357      | 316724 | 37645.74646 | -4811.745452 | 1218330 | 179170.011   | -84710      | 39424.989   | -57156      | 1275   | -78004       | 26323897.53 | 42848.989   | 155752 | 12001 | -167753 | 1206658 | -1206658 | 330574 | 2328771.213 | 1248226 | 523682 | -2328771.213 | 226289.213 | -1009435 | -697152 | 1706587 | 1157381 | -1157381 | 1480979 | 153462 | 164923.213 | -1841204 | -441211 | 483050.787 | 4331344 | 3634192 | -2013354 | -2163848 | -81048  | -85282  | 709340  | 5205922.128 | 0.030008333 | 0.003538613 | 0.027309766 | 0           | 75.6682106  | 74.27991891 | 0.046 | 0.333580287 | 0.354957562 | 101.869 | 101.3 | 1.13  | 103.4 | 101.39634 | 52193.31  | 95.65       | 101.25      | 99.4        | 95.7        | 103.8605307 | 0.009945652 | 27.98375606 | 0.052814324  | 0.015136465 | 0.02406256   | 0.006590628  | 0.062598016 | 0.128427831 | 0.611932813 | 0.015380515  | 0.014014014  | 0.007267586 | 0.18562132   | 22950.1 | 25857   |
| 2018 | 1066167 | 328194 | 334454 | 555394 | 512818 | 1771391 | 480055.6522 | 88721.58209 | 141572 | 665490 | 16570 | 18683.6736  | 29951.57742 | 5302.096182  | 1934  | 61354.9994  | 16425.3258   | 31740.90322 | 5302.096182  | 9676   | 1326322     | 907171      | 419151 | 31740.90322 | 5302.096182  | 863136  | -419352.251  | 222355      | 99226.251   | 220437      | -3082  | -119584      | 26726264.67 | 52023.251   | 162421 | 12232 | -174653 | 1224330 | -1224330 | 313141 | 2380794.464 | 1202244 | 541877 | -2380794.464 | 323532.464 | -1018031 | -711907 | 1729938 | 1075545 | -1075545 | 1367631 | 171179 | 264149.464 | -1580695 | -455588 | 233323.536 | 4143068 | 3431161 | -1922397 | -2116645 | 139389  | -88364  | 556856  | 5294643.711 | 0.0268      | 0.001899746 | 0.026346512 | 0           | 76.53979104 | 74.34151252 | 0.046 | 0.334562833 | 0.356671486 | 102.957 | 102.5 | 1.181 | 106.7 | 103.37155 | 59646.3   | 115.1333333 | 106.79      | 102.4       | 115.1       | 106.2717832 | 0.00999077  | 28.75506624 | 0.046783999  | 0.017042434 | 0.020038086  | 0.013158936  | 0.04275944  | 0.090604758 | 0.783813627 | 0.027562782  | 0.011846002  | 0.010680384 | 0.203694023  | 23143.4 | 25898.9 |
| 2019 | 1074680 | 327703 | 334512 | 567784 | 508031 | 1796648 | 506154.8658 | 84149.38931 | 146882 | 681104 | 14113 | 18527.37968 | 29505.85369 | 3134.525981  | 3429  | 56871.99979 | 12795.62011  | 31617.4738  | 3134.525981  | 11436  | 1082293     | 869739      | 212554 | 31617.4738  | 3134.525981  | 1426904 | 189603.392   | -29822      | -124910.392 | -98786      | -18379 | 82294        | 28774461.9  | 28976.608   | 165890 | 12654 | -178544 | 1286064 | -1286064 | 289695 | 2409771.072 | 1370663 | 573547 | -2409771.072 | 175866.072 | -995791  | -728244 | 1724035 | 1228692 | -1228692 | 1541288 | 186937 | 139239.072 | -1780685 | -501746 | 414966.928 | 4511629 | 3783385 | -2037546 | -2270532 | 40603   | -106743 | 590833  | 5378793.1   | 0.026025    | 0.003443494 | 0.023887824 | 0           | 77.17560137 | 74.26229167 | 0.046 | 0.334024374 | 0.3556951   | 103.923 | 103.2 | 1.12  | 105.8 | 101.08213 | 55150.68  | 108.25      | 107.02      | 101.5       | 108.3       | 108.7763154 | 0.009839324 | 29.25704467 | -0.001496066 | 0.015893305 | 0.014258286  | 0.010643114  | 0.035347303 | 0.094225457 | 0.808649568 | 0.017457043  | 0.006829268  | 0.009382558 | -0.059785755 | 23280   | 25861.5 |
| 2020 | 963726  | 293530 | 343580 | 488941 | 429156 | 1660621 | 428608.9298 | 46105.51741 | 202196 | 642972 | 10558 | 16955.9478  | 19250.6187  | 8263.329099  | 3779  | 54290.99891 | 16112.05111  | 28583.66981 | 8263.329099  | 10665  | 1126688     | 801632      | 325056 | 28583.66981 | 8263.329099  | 1149580 | -107752.74   | 112818      | -91533.26   | 130876      | -14959 | -29449       | 26993213.3  | 163194.74   | 185433 | 11474 | -196907 | 1371058 | -1371058 | 263335 | 2572965.812 | 1443952 | 764972 | -2572965.812 | 100706.812 | -1041503 | -736292 | 1777795 | 1234134 | -1234134 | 1651571 | 269617 | 47705.812  | -1690684 | -689767 | 411557.188 | 4705531 | 3969239 | -2006020 | -2525260 | 171479  | -121702 | 512264  | 5424898.617 | 0.023283333 | 0.00362841  | 0.022529526 | 0           | 72.73822717 | 68.88024466 | 0.046 | 0.306110974 | 0.323621459 | 105.601 | 103   | 1.142 | 108.1 | 101.91719 | 31741.38  | 78.45       | 99.6        | 96.3        | 78.5        | 111.9985301 | 0.00836872  | 28.16334576 | -0.1042804   | 0.008571722 | -0.075711547 | 0.011754933  | 0.060642737 | 0.066451332 | 0.652427134 | -0.037382413 | -0.001937984 | 0.01614657  | -0.275288684 | 22830.1 | 25140.6 |
| 2021 | 1030124 | 357215 | 352718 | 582192 | 540198 | 1782051 | 483365.7672 | 107669.6636 | 188601 | 692915 | 10905 | 14916.66234 | 12622.01724 | 13199.6451   | 2326  | 60678.00031 | 25091.33797  | 29134.35521 | 13199.6451   | 9765   | 1141970.282 | 738858.2821 | 403112 | 29134.35521 | 13199.6451   | 1116373 | -60674.87011 | 55160.28211 | 275576.588  | -151307     | -5171  | -113584      | 29441497.75 | 105431.588  | 200683 | 10817 | -211500 | 1428434 | -1428434 | 233263 | 2678397.4   | 1366294 | 868289 | -2678397.4   | 210551.4   | -871902  | -763488 | 1635390 | 1372850 | -1372850 | 1583746 | 284629 | 323282.4   | -1563895 | -783662 | 155899.6   | 4818976 | 4055488 | -1960123 | -2355115 | 20172   | -126873 | 366451  | 5532568.281 | 0.020259167 | 0.002667731 | 0.023582902 | 0           | 79.33378445 | 74.72335354 | 0.046 | 0.322101944 | 0.344334175 | 106.17  | 105   | 1.183 | 109.3 | 101.53931 | 64858.995 | 112.4833333 | 125.39      | 105         | 85.59190031 | 114.2006256 | 0.007572628 | 30.84736029 | 0.216962491  | 0.019847314 | 0.073123247  | 0.01737685   | 0.096964642 | 0.047931408 | 0.598685623 | 0.095301693  | 0.019417476  | 0.005388207 | 0.433821967  | 22462.7 | 24829.5 |

[work in progress] 🛠️



