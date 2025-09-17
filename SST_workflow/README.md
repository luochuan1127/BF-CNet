**Overview**

Due to file size limitations, the three *.mat data files required for the SST_workflow folder have been uploaded to Google Drive.
Please download them from the link below and place them into the SST_workflow directory before running the project step by step as described in README.md.
Google Drive:
https://drive.google.com/drive/folders/17IY4e7kf1v0a5ZPjYIJrl3AaxAsaTJ3L?usp=sharing

This project demonstrates the end-to-end workflow for SST (Sea Surface Temperature) data, including dataset creation, model training, and evaluation. This version is organized as the simplest directly executable release, making it easy for reviewers and users to run without additional setup, including:



Dataset creation from CMEMS data.

Model training using UFEM, which is built upon the BayesNF framework (DOI: 10.1038/s41467-024-51477-5).

Evaluation of the trained model, with deterministic metrics and uncertainty quantification.



The code provides a complete pipeline for preparing SST data, training the model, and assessing its performance under cloud-masked conditions.



**Workflow**

1\. Dataset Creation

 	Run code01\_dataset\_creation.m.

 	This script calls SST\_dataset.mat and creates the cloud-masked datasets:

 		train\_CMEMS\_SST.csv

 		validation\_CMEMS\_SST.csv

 	Outputs are stored in the input/ folder.

2\. Model Training

 	Run Model\_SST.ipynb.

 	This notebook trains the model, saves the trained model to the model/ folder, and generates reconstruction results stored in the output/ folder.

3\. Best Quantile Selection (SSIM-based)

 	Run code02\_best\_quantile.m.

 	This script selects the optimal quantile based on SSIM evaluation.

4\. Result Display and Deterministic Evaluation

 	Run code03\_result\_display.m.

 	This script computes deterministic evaluation metrics and produces visualizations:

 		Comparison plots showing cloud-gap filling performance

 		Uncertainty distribution maps

5\. Uncertainty Evaluation

 	Run code04\_uncertainty\_evaluation.m.

 	This script performs uncertainty quantification and evaluation.



**Evaluation Metrics**



The evaluation includes both deterministic and probabilistic metrics:

1\. Deterministic Metrics

 	MAE (Mean Absolute Error)

 	RMSE (Root Mean Squared Error)

 	SSIM (Structural Similarity Index)

 	PSNR (Peak Signal-to-Noise Ratio)

2\. Uncertainty Metrics

 	PICP (Prediction Interval Coverage Probability)

 	ECE (Expected Calibration Error)



**File Structure**



1\. code01\_dataset\_creation.m → Build training/validation datasets.

2\. Model\_SST.ipynb → Train the model and generate reconstruction results.

3\. code02\_best\_quantile.m → Select the best quantile based on SSIM.

4\. code03\_result\_display.m → Compute metrics \& visualize results.

5\. code04\_uncertainty\_evaluation.m → Perform uncertainty evaluation.



input/ → Contains input datasets.

model/ → Stores trained models.

output/ → Stores reconstruction results and visualizations.


