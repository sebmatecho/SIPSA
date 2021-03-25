# SIPSA-index: Project overview

This repository contains code for the paper 
> [**Proposal for a SIPSA index and its relation to food inflation for the Colombian case: empirical evidence**](https://revistas.usantotomas.edu.co/index.php/estadistica/article/view/5635/5549).

SIPSA is actually a branch of DANE (Departamento Administrativo Nacional de Estadística - National Administrative Department of Statistics) in charge of collecting food prices   for every product sold in marketplaces located in every city, ranging from medium to big, in the country. Our aim was to utilize such information as a proxy for understanding and thus, predicting, food inflation (wich, for colombian market, represents ~25% of overall inflation). Several conclusions: 

- Construct an index based on information available directly from DANE's webpage. 
- Using such index, an improvement of 40% is observed when attempting to forecast food inflation (based on train/test samples split) using a Box-Jenkins approach (SARIMAX model). 
- The proposed index displays, in the short run, a similar behaviour as food inflation, leading to an improvement of ~20% in joint forecast (based on train/test samples split) using VAR models
- The proposed index displays, in the long run, a similar behaviour as food inflation, leading to a equilibrium in the sense of Johansen, leading to the conclusion that index is worthy tool to be considered when approaching food inflation's forecasting. 

# Data Wrangling and index construction 

Data files were obtained directly from DANE's official webpage

    ![Tux, the Linux mascot](/Images/ambas.eps)

# Model bulding 

# Model Performance




