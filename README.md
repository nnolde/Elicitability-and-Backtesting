# Elicitability-and-Backtesting
R code to accompany paper "Elicitability and Backtesting: Perspectives for Banking Regulation" by N. Nolde and J.F. Ziegel (2016)

Files to reproduce the simulation study:

- Rcode_sim_RMforecasting.R:  Computation of one-step ahead risk measure forecasts
- Rcode_sim_TradBacktesting.R: Traditional backtesting analysis
- Rcode_sim_ComBacktesting.R: Comparative backtesting analysis
- Rcode_sim_RMdataprep.R: Preparation of forecasting outputs for backtesting analyses
- Sim3norm3.RDATA,Sim3std3.RDATA,Sim3sstd3.RDATA: Stored risk measure forecasts 
- simdat.RDATA: simulated time series
- Rfns.R: contains all relevant functions

Files to reproduce the data analysis (NASDAQ Composite index)
- RcodeDataAnalysis.R: relevant R code
- nasdaq2016.csv: data file

Files to reproduce small out-of-sample size simulation study (Section D of the Online Supplement)
- Rcode_sim_SmallSizeEffect.R
