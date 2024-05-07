# ProtBENCH: A Comprehensive Benchmark for Machine Learning-based Prediction of Drug/Compound-Target Interactions
In this study/repository, we performed a rigorous benchmark analysis to investigate the representation capability of various protein featurization techniques together with different datasets (and their structured splits), modeling approaches and evaluation strategies for DTI prediction modelling. Below, we summarized the major contributions of our study to the literature:
* We proposed protein family-specific challenging benchmark datasets with high coverage on both compound and protein spaces, which can be used as reliable, reference/gold-standard datasets for fair evaluation of models at multiple difficulty levels on DTI modelling tasks.
* We employed a network-based strategy for splitting train/test folds of these protein family-specific datasets, which ensures that train and test folds are totally dissimilar from each other with a minimum loss of data points. Hence, it can aid researchers in designing more powerful and robust DTI prediction models that have a real translational value.
* We extended the scope of our study by involving the state-of-the-art protein representation learning methods and discussed their potential in DTI prediction, in which the studies regarding their usage in DTI prediction modelling are limited.
As the first near comprehensive benchmark study of protein representation methods in computational drug discovery and repurposing, we hope this study will assist ongoing work in the context of representing biomolecules for high performance drug-target interaction prediction.

The study is summarized in the schematic workflow below.

![ProtBENCH_workflow](https://user-images.githubusercontent.com/13165170/167027598-d4b13b56-7266-4c22-897b-d847021fe462.png)

Please refer to our article for more information:

Atas Guvenilir, H., Doğan, T. (2023). How to approach machine learning-based prediction of drug/compound–target interactions. *Journal of Cheminformatics*, 15:16.

DOI/URL: [10.1186/s13321-023-00689-w](https://doi.org/10.1186/s13321-023-00689-w)

## Programming Environment and Files
**Descriptions of folders and files:**

*	**datasets** folder includes benchmark datasets constructed by applying extensive filtering operations in different scales (this folder will be created together with its sub-folders and files when the “datasets.zip” file is downloaded from [here](https://drive.google.com/file/d/1zVOyFIEOo33yeF3vFE8paz5pS5H5Z99N/view?usp=sharing) and uncompressed, as explained below).
    * **small_scale** folder contains bioactivity data of compound-centric datasets in binary classification format (i.e., 1:active, 0:inactive). Each folder named with center compound ChEMBL id and name (e.g., ChEMBL44_Genistein folder)  involves 5-fold train/test splits of protein samples having experimentally measured bioactivity data for clusters of corresponding center compound. Protein representations to be used for the construction of target feature-based classification models are present in the **feature_vectors** folder.
    * **medium_scale** folder contains train and test datapoints of modified Davis kinase benchmark dataset. Both protein representations and ECFP4 compound fingerprints are involved in the **feature_vectors** folder.
    * **large_scale** folder contains protein family-specific bioactivity datasets. For each protein family, bioactivity samples are distributed into train and test sets by applying three different splitting strategy, and involved in separate folders named **fully_dissimilar_split**, **dissimilar_compound_split**, and **random_split**. Protein representations and compound fingerprints used for all split datasets are stored in the **feature_vectors** folder.
*	**scripts** folder includes script files required for the construction of models on different dataset scales. It also involves “score_metrics.py” script file utilized for performance calculation of PCM models on medium- and large-scale datasets.
*	**results** folder contains test performance results of all models on each scale. 
*	**initial_files** folder involves protein-protein and compound-compound pairwise similarity files as well as initial bioactivity datasets of each protein family required for "network_based_train_test_splitting.py" script file as input.
*	**train_test_splitting** folder contains "network_based_train_test_splitting.py" script file to reproduce fully_dissimilar_split, dissimilar_compound_split, and random_split based train/test bioactivity files of large-scale protein family-specific datasets. It can be also used to split your own bioactivity datasets by supplying input files with required formats and adjusting partition parameters. 

**Dependencies:**

* Python 3.6.13
* Scikit-learn 0.22.1
* Pandas 1.1.5
* Numpy 1.19.5
* Python-louvain 0.13 (only required if the user runs "network_based_train_test_splitting.py" script)

**Step-by-step operation:**
1. We highly recommend you to use conda platform for installing dependencies properly. After installation of appropriate [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) version for your operating system, create and activate conda environment, and then install dependencies in the activated environment, as below:
```
conda create -n protbench python=3.6.13
conda activate protbench
conda install -c anaconda scikit-learn=0.22.1
conda install -c anaconda pandas=1.1.5
pip install numpy==1.19.5
pip install python-louvain==0.13
```
2. Clone this repository.

To run models directly;

3. Download datasets from [here](https://drive.google.com/file/d/1zVOyFIEOo33yeF3vFE8paz5pS5H5Z99N/view?usp=sharing) and uncompress the “datasets.zip” file. Place the uncompressed folder in the cloned repository at the same level as the **results** and **scripts** folders. 
4. Set the location of **scripts** folder as the current working directory, and run the corresponding script to build models at the data scale of your interest. Instead of creating all models of the selected scale at once, you can easily edit the script file according to your purpose before running it. You need to uncomment "save model" and "save predictions" code blocks in the script file to save the constructed model and its predicted outputs, respectively.

To reproduce train/test splits;

3. Download input files from [here](https://drive.google.com/file/d/1d9ggqtmz9x0cBIdGAmR0URn2FpqX_zcH/view?usp=share_link) and uncompress the “initial_files.zip” file. Place the uncompressed folder into the **train_test_splitting** folder in the cloned repository.
4. Create an empty folder with the name "final_files" and sub-folders "fully_dissimilar_split", "dissimilar_compound_split", and "random_split" in it to save the output.
5. Run the "network_based_train_test_splitting.py" script by uncommenting the line of the splitting method and the family of interest. Please note that the fully-dissimilar-split datasets you create may slightly differ from the ones used in this study due to differences in partition parameters and further filtering operations used for the removal of data points with missing protein representations.

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


