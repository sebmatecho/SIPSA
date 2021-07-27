# SIPSA-index: Project overview

This repository contains code for the paper 
> [**Proposal for a SIPSA index and its relation to food inflation for the Colombian case: empirical evidence**](https://revistas.usantotomas.edu.co/index.php/estadistica/article/view/5635/5549).

SIPSA is actually a branch of DANE (Departamento Administrativo Nacional de Estadística - National Administrative Department of Statistics, public agency in Colombia dealing with national statistics) in charge of collecting food prices for every product sold in marketplaces located in every city ranging from medium to big sized, across the country. The main aim was to utilize such information as a proxy for understanding and thus, predicting, food inflation (wich, for colombian market, represents ~25% of overall inflation). Several conclusions: 

- Construct an index based on information available directly from DANE's webpage. 
- Using such index, an improvement of 40% is observed when attempting to forecast food inflation (based on train/test samples split) using a Box-Jenkins approach (SARIMAX model). 
- The proposed index displays, in the short run, a similar behaviour as food inflation, leading to an improvement of ~20% in joint forecast (based on train/test samples split) using VAR models
- The proposed index displays, in the long run, a similar behaviour as food inflation, leading to a equilibrium in the sense of Johansen, leading to the conclusion that index is worthy tool to be considered when approaching food inflation's forecasting. 

# Data Wrangling and index construction 

SIPSA weekly reports would display the price of each product at each marketplace clustered by food type: i) vegetables and greens, ii) fresh fruit, iii) tubers, roots and bananas, iv) grains and cereals, v) fish, vi) eggs and dairy products, vii) meats, and viii) processed products. Each dataset cotains ~4500 rows. Weekly data raging from late 2012 to mid 2018 was considered (280 datasets). Then, weekly, median price for each product and median price for clustered products were considered. Finally, a unique weekly value was obtained considering the participation of such clusters in food inflation (i), 9.4% ii) 5%, iii) 4.1%, iv) 18.6%, v) 21.4%, vi) 23.7, vii) 3% and viii) 15.87$, respectively). Data files were obtained directly from DANE's official webpage.

Once the proposed index was computed, it was compared with correspondant food inflation. A correlation of 0.6 was found. 


![Tux, the Linux mascot](/Images/ambas.png)

# Model bulding 

# Model Performance




