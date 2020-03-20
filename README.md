# Florida-Niche
Predict the probability of occurrence of 23 coral species, including the critically endangered Acropora cervicornis, using observations at 985 sites from 2011–2015. 

Description of Documents and Code Files
	Used for replication and extension of manuscript results
	
	
	
FOLDER: trainingData
	
	CSV files containing the raw coral data and species information used in the main scripts. 
	This data is merged with the environmental data to create the points necessary to train
	and test the models.
	
	
	
FOLDER: predictorData

	TIF/TXT files containing the raster environmental data used in the main scripts. This data
	is merged with the raw coral data to create the points necessary to train and test the models.
	
	
	
FOLDER: polygons

	Shapefiles containing the extents within which the data is prepared and run as well as the
	polygons for splitting the area of interest into regions for data collation.
	
	
FOLDER: kmz

	KMZ files for each valid species created from the results found in the results folder.
	
	
	
NicheSpaceModel_bulkScript.R
	
	This is the primary script. This performs all of the necessary actions for loading the coral and
	environmental data and preparing it for model development. A model for each species is then developed
	in a loop. The model development attempts to find the optimal combination of predictors by training
	models and then checking performance until an optimal combination is found. The resulting models
	are saved in the results folder. Additional all data and all data predictor CSV files are produced
	once all of the models are developed.



NicheSpaceModel_plots.R

	Produces the top4PD jpeg in the results folder which plots the top predictors using 4 of the most
	relevant species in an attempt to understand the trending.
	
	
	
NicheSpaceModel_kmz.R

	Produces the kmz files found in the kmz folder from the results in the results folder for each
	valid species.



Disturbance R code.R

	Uses the ACER results and different cutoff values to analyze the landscape metrics and spatial
	information of the resulting polygons. The cutoffs are used to simulate different levels of
	disturbances.



brt.functions.R

	Additional functions imported to run within the other scripts.
