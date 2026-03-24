Modelling heatwaves in central France: a case-study in extremal dependence
H. C. Winter and J. A. Tawn
Appl. Statist., 65 (2016), 345 -- 365


Data set:
Data are from EUROPEAN CLIMATE ASSESSMENT & DATASET (ECA&D)

Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment. Int. J. of Climatol., 22, 1441-1453.
 Data and metadata available at http://www.ecad.eu 

Data contained in a text file with 4 columns: (1) Source identifier; (2) Date; (3) Maximum temperature; (4) Quality code. For more information see the header.


Programs and use (all code is written in R):
load_data.R - Loads the data file, produces some diagnostic plots and obtains time-lagged data
fit_marg_dep_mods.R - Fit the marginal model and different dependence models from the manuscript.
fit_marg_dep_mods_funcs.R - Functions required to fit marginal and dependence models.
dur_dists_from_dep.R - Apply two different types of tail chain simulation technique outlined in manuscript Appendix with the dependence parameters estimated from fit_marg_dep_mods.R.
dur_dists_from_dep_funcs.R - Functions required to obtain tail chains.


Hugo Winter
Department of Mathematics and Statistics
Lancaster University
Lancaster
LA1 4YF
UK

E-mail: hugocwinter@gmail.com
