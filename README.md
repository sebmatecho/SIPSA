# SIPSA-index: Project overview

This repository contains code for the paper 
> [**Proposal for a SIPSA index and its relation to food inflation for the Colombian case: empirical evidence**](https://revistas.usantotomas.edu.co/index.php/estadistica/article/view/5635/5549).

SIPSA is actually a branch of DANE (Departamento Administrativo Nacional de Estadística - National Administrative Department of Statistics, public agency in Colombia dealing with national statistics) in charge of collecting food prices for every product sold in marketplaces located in every city ranging from medium to big sized, across the country. The main aim was to utilize such information as a proxy for understanding and thus, predicting, food inflation (wich, for colombian market, represents ~25% of overall inflation). Several conclusions: 

- Construct an index based on information available directly from DANE's webpage. 
- Using such index, an improvement of 40% is observed when attempting to forecast food inflation (based on train/test samples split) using a Box-Jenkins approach (SARIMAX model). 
- The proposed index displays, in the short run, a similar behaviour as food inflation, leading to an improvement of ~20% in joint forecast (based on train/test samples split) using VAR models
- The proposed index displays, in the long run, a similar behaviour as food inflation, leading to a equilibrium in the sense of Johansen, leading to the conclusion that index is worthy tool to be considered when approaching food inflation's forecasting. 

# Data Wrangling and index construction 

SIPSA weekly reports would display the price of each product at each marketplace clustered by food type: i) vegetables and greens, ii) fresh fruit, iii) tubers, roots and bananas, iv) grains and cereals, v) fish, vi) eggs and dairy products, vii) meats, and viii) processed products. Each dataset cotains ~4500 rows. Weekly data raging from late 2012 to mid 2018 was considered (310 datasets). Then, weekly, median price for each product and median price for clustered products were considered. Finally, a unique weekly value was obtained considering the participation of such clusters in food inflation (i), 9.4% ii) 5%, iii) 4.1%, iv) 18.6%, v) 21.4%, vi) 23.7, vii) 3% and viii) 15.87$, respectively). Data files were obtained directly from DANE's official webpage.

Once the proposed index was computed, it was compared with correspondant food inflation. A correlation of 0.6 was found. In the following figure, the red line displays food variation and blue line displays the proposed SIPSA index. 

![Tux, the Linux mascot](/Images/ambas.png)

# Model bulding 

Three approaches were used in order to study the relation between food inflation and the proposed SIPSA index: i) SARIMAX model (for individual forecasting purposes), ii) VAR model (for joint forecasting purposes) and iii) VEC model (for equilibrium equation purposes). 


## SARIMAX model 

After an exhaustive search process, a <img src="https://render.githubusercontent.com/render/math?math=\text{SARIMAX}(10,1,13)\times (1,1,0)_{12}"> was established as the best model. Such model is indexed by the following parameters. 

| Parameter      | <img src="https://render.githubusercontent.com/render/math?math=\phi_1">  | <img src="https://render.githubusercontent.com/render/math?math=\phi_{10}">  | <img src="https://render.githubusercontent.com/render/math?math=\theta_1"> | <img src="https://render.githubusercontent.com/render/math?math=\phi_{13}"> | <img src="https://render.githubusercontent.com/render/math?math=\Phi_13">  | <img src="https://render.githubusercontent.com/render/math?math=\text{SIPSA}_t">  | <img src="https://render.githubusercontent.com/render/math?math=\text{SIPSA}_{t-1}"> |
| :------------: | :-------------: | :--------------: | :-------------: | :----------------: | :-----------: | :-------------------: | :----------------------: |
| Estimation     | -0\.466         | -0\.261          | -0\.492         | 0\.323             | -0\.457       | 0\.300                | 0\.157                   |
| Standard error | 0\.1147         | 0\.1115          | 0\.1357         | 0\.1246            | 0\.1130       | 0\.0425               | 0\.0434                  |
| t-value        | -4\.06          | -2\.34           | -3\.62          | 2\.60              | -4\.04        | 7\.06                 | 3\.63                    |

This model displays no evidencie of autocorrelation (p-value for Ljung-Box test of 0.1686 was obtained) as well of no normality on its residuals (p-value for Jarque-Bera test <0.01 was found). Additionaly, visual inspection of residuals shows good adequacy. BIC was used as general goodness of fit measure, for this model, we report a BIC of -120. 

![Tux, the Linux mascot](/Images/valarimax.png)


## VAR model 

After an exhaustive search process, a SVAR(1) was established as the best model. Such model is indexed by the following parameters. 
<img src="https://render.githubusercontent.com/render/math?math=\text{Inf}_t = 0.26475 -0.00195t+ 0.47304 \text{Inf}_{t-1} +0.06235 \text{SISPA}_{t-1} + \sum_{i=1}^{11}\text{SDinf}_i \bm{1}_i,">
<img src="https://render.githubusercontent.com/render/math?math=\text{SIPSA}_t = 0.34206-0.00464t+ 0.44231\text{Inf}_{t-1} -0.0003\text{SISPA}_{t-1} + \sum_{i=1}^{11}\text{SDSIPSA}_i \bm{1}_i ">

Where

|  i                | 1        | 2       | 3       | 4       | 5       | 6       | 7       | 8       | 9       | 10      | 11      |
| :---------------------: | :------: | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: |
|<img src="https://render.githubusercontent.com/render/math?math=\text{SDinf}_i">   | -1\.1881 | -1\.055 | -0\.941 | -1\.384 | -1\.417 | -0\.991 | -1\.567 | -1\.131 | -1\.428 | -1\.147 | -0\.877 |
|<img src="https://render.githubusercontent.com/render/math?math=\text{SDSIPSA}_i">  | -0\.930  | -0\.293 | -0\.590 | -1\.755 | -1\.210 | -0\.060 | -1\.597 | -0\.895 | -0\.655 | -0\.535 | -0\.200 |

Once the model was obtained, stability was verified by means of the impulse-response approach and residuals were checked by means of Portmanteau, Arch and Mutivariate normality test at usual 5% of significance.

![Tux, the Linux mascot](/Images/fluctuacion.png)

 At the end, the proposed model aims to provide a reasonable pattern understanding in the short run relation between the proposed SIPSA index and food inflation. This is clearly determined by the variance descomposition. SIPSA proposed index explains ~40% of the variance in the forecast of food inflation. 
 
 ![Tux, the Linux mascot](/Images/fevdsipsa.eps)
 
 ## VEC model
 
 

