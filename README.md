# ProtBENCH: A Comprehensive Benchmark for Machine Learning-based Prediction of Drug/Compound-Target Interactions
In this study/repository, we performed a rigorous benchmark analysis to investigate the representation capability of various protein featurization techniques together with different datasets (and their structured splits), modeling approaches and evaluation strategies for DTI prediction modelling. Below, we summarized the major contributions of our study to the literature:
* We proposed protein family-specific challenging benchmark datasets with high coverage on both compound and protein spaces, which can be used as reliable, reference/gold-standard datasets for fair evaluation of models at multiple difficulty levels on DTI modelling tasks.
* We employed a network-based strategy for splitting train/test folds of these protein family-specific datasets, which ensures that train and test folds are totally dissimilar from each other with a minimum loss of data points. Hence, it can aid researchers in designing more powerful and robust DTI prediction models that have a real translational value.
* We extended the scope of our study by involving the state-of-the-art protein representation learning methods and discussed their potential in DTI prediction, in which the studies regarding their usage in DTI prediction modelling are limited.
As the first near comprehensive benchmark study of protein representation methods in computational drug discovery and repurposing, we hope this study will assist ongoing work in the context of representing biomolecules for high performance drug-target interaction prediction.

The study is summarized in the schematic workflow below.

![ProtBENCH_workflow](https://user-images.githubusercontent.com/13165170/167027598-d4b13b56-7266-4c22-897b-d847021fe462.png)

Please refer to our pre-print for more information:

Atas, H. & Doğan, T. (2022). How to Approach Machine Learning-based Prediction of Drug/ Compound-Target Interactions. *bioRxiv*, 2022.05.01.490207.

DOI/URL: [10.1101/2022.05.01.490207](https://biorxiv.org/cgi/content/short/2022.05.01.490207)

## Programming Environment and Files
**Descriptions of folders and files:**

*	**datasets** folder includes benchmark datasets constructed by applying extensive filtering operations in different scales.
    * **small_scale** folder contains bioactivity data of compound-centric datasets in binary classification format (i.e., 1:active, 0:inactive). Each folder named with center compound ChEMBL id and name (e.g., ChEMBL44_Genistein folder)  involves 5-fold train/test splits of protein samples having experimentally measured bioactivity data for clusters of corresponding center compound. Protein representations to be used for the construction of target feature-based classification models are present in the **feature_vectors** folder.
    * **medium_scale** folder contains train and test datapoints of modified Davis kinase benchmark dataset. Both protein representations and ECFP4 compound fingerprints are involved in the **feature_vectors** folder.
    * **large_scale** folder contains protein family-specific bioactivity datasets. For each protein family, bioactivity samples are distributed into train and test sets by applying three different splitting strategy, and involved in separate folders named **fully_dissimilar_split**, **dissimilar_compound_split**, and **random_split**. Protein representations and compound fingerprints used for all split datasets are stored in the **feature_vectors** folder.
*	**scripts** folder includes script files required for the construction of models on different dataset scales. It also involves “score_metrics.py” script file utilized for performance calculation of PCM models on medium- and large-scale datasets.
*	**results** folder contains test performance results of all models on each scale. 

**Dependencies:**

* Python 3.6
* Scikit-learn 0.22
* Pandas 1.1.5
* Numpy 1.19.5

**Step-by-step operation:**
1. Install dependencies.
2. Clone this repository.
3. Download datasets from [here](https://drive.google.com/drive/folders/1lID4vX9c4hm7hWTqzsVFj4ZmMvjmL2EN) and uncompress the “datasets.zip” file. Place the uncompressed folder in the cloned repository at the same level as the **results** and **scripts** folders. 
4. Set the location of **scripts** folder as the current working directory, and run the corresponding script to build models at the data scale of your interest. Instead of creating all models of the selected scale at once, you can easily edit the script file according to your purpose before running it. You need to uncomment "save model" and "save predictions" code blocks in the script file to save the constructed model and its predicted outputs, respectively.

**Example commands to run the scripts for building models of interest:**
   
For generating models on small-scale, run the "small_scale_target_feature_based_svm_models.py" script for SVM algorithm and "small_scale_target_feature_based_rf_models.py" for RF algorithm by defining the folder name of the dataset (i.e., "ChEMBL44_Genistein", "ChEMBL50_Quercetin", "ChEMBL83_Tamoxifen", "ChEMBL91_Miconazole", "ChEMBL104_Clotrimazole", "ChEMBL633_Amiodarone", "ChEMBL808_Econazole", "ChEMBL116438_Curcumin", and "ChEMBL295698_Levoketoconazole") as parameter, e.g.:
   
```
python small_scale_target_feature_based_svm_models.py ChEMBL44_Genistein
```
   
For generating models on medium-scale, run the "medium_scale_pcm_rf_models.py" script without defining any parameter, e.g.:
   
```
python medium_scale_pcm_rf_models.py
```
   
For generating models on large-scale, run the "large_scale_pcm_rf_models.py" script by defining the split name (i.e., "fully_dissimilar_split", "dissimilar_compound_split", and "random_split") and protein family of interest (i.e., "epigenetic-regulators", "ion-channels", "membrane-receptors", "transcription-factors", "transporters", "hydrolases", "oxidoreductases", "proteases", "transferases", and "other-enzymes") as parameters, respectively, e.g.:
   
```
python large_scale_pcm_rf_models.py fully_dissimilar_split proteases
```
   
## License

Copyright (C) 2022 HUBioDataLab

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.


