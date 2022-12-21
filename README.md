# Reproduction for the methodology of "Are Factor Regression and Sparse Regression Adequate?"

# 1. **Real Data Analysis** <br />

## Data  <br />

### 1.1.  Data Description <br />
FRED-MD is a large macroeconomic database designed for the empirical analysis of “big data.” The datasets of monthly and quarterly observations mimic the coverage of datasets already used in the literature, but they add three appealing features. They are updated in real-time through the FRED database. They are publicly accessible, facilitating the replication of empirical work.  <br />
### 1.2. Data Availability <br />
 The original dataset can be found in https://research.stlouisfed.org/econ/mccracken/fred-databases/. <br />
 In our .zip file, the original Dataset is provided in folder ``Code_Reproduction/Real_Data/Pre_process/current.csv." <br />
 The description of the data is provied in "Code_Reproduction/Real_Data/Pre_process/Description_of_Variables.pdf."  <br /> This description describes the meaning of all columns of this datset.

## Codes Description and Implementation Details:
 There are in total four folders under the folder: ``Code_Reproduction/Real_Data."  <br /> 1. Pre_process, <br />2. Prediction_Table_3 <br />  3. PCR_Adequate_Table4, <br /> 4. Sparse_Adequate_Table4. <br />
### Folder `Code_Reproduction/Real_Data/Pre_process':
 It contains the Codes for real data pre-processing. The document `current.csv' in that folder is our downloaded orignal dataset from the aforementioned website. <br />
 
 The code **realdata_read.R** contains the codes for pre-processing the data. The line 4 in "realdata_read.R", we have a variable called: "choose_time". If one lets choose_time=1 and run the codes, the output is the processed dataset **realdata_115_126.csv**. If one lets choose_time=2, and run the codes, the output is the processed dataset **realdata_189_126.csv**.<br />
 
 Finally, for document "Description_of_Variables.pdf" inside that document, it contains the detailed description of all variables in the orignal dataset "current.csv". <br />

### Folder `Code_Reproduction/Real_Data/Prediction_Table3':
This folder contains the codes for reproducing the prediction results of Table 3 in our paper. <br />

*For code file **realdata_prediction_115.R**, it first reads the data "realdata_115_126.csv" in that folder. On the line 54 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE" line during "08.2010-02.2020" in Table 3. If we set j=81 and run the codes, we reproduce the line "GS5" line during "08.2010-02.2020" in Table 3. <br />

*For code file **realdata_prediction_189.R**, it first reads the data "realdata_189_126.csv" in that folder. On the line 54 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE" line during "02.1992-10.2007" in Table 3. If we set j=81 and run the codes, we reproduce the line "GS5" line during "02.1992-10.2007" in Table 3. <br />

For both .R files, the output are printed in line 149-161 and 189 for different methods :Lasso (sparse linear regression), PCR (latent factor regression), Ridge (Ridge regression), El-Net (Elastic Net), RF (Random Forest), and FarmSelect (Factor adjusted Lasso), respectively.

### Folder `Code_Reproduction/Real_Data/PCR_Adequate_Table4':
This folder contains the codes for reproducing the hypothesis testing results of column "LA_factor" of Table 4 in our paper. <br />

*For code file **Real_data_PCR_115.R**, it first reads the data "realdata_115_126.csv" in that folder. On the line 17 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "LA_factor," during "08.2010-02.2020" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "LA_factor", during "08.2010-02.2020" in Table 4. <br />


*For code file **Real_data_PCR_189.R**, it first reads the data "realdata_189_126.csv" in that folder. On the line 17 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "LA_factor," during "02.1992-10.2007" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "LA_factor", during "02.1992-10.2007" in Table 4. <br />

For both .R files, the output are printed p-values of the test on whether PCR is adequate.

### Folder `Code_Reproduction/Real_Data/Sparse_Adequate_Table4':
This folder contains the codes for reproducing the hypothesis testing results of column "SP_Linear" of Table 4 in our paper.

*For code file **Real_data_Sparse_115.R**, it first reads the data "realdata_115_126.csv" in that folder. On the line 94 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "SP_Linear," during "08.2010-02.2020" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "SP_Linear", during "08.2010-02.2020" in Table 4. <br />


*For code file **Real_data_Sparse_189.R**, it first reads the data "realdata_189_126.csv" in that folder. On the line 94 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "SP_Linear," during "02.1992-10.2007" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "SP_Linear", during "02.1992-10.2007" in Table 4. <br />
 
For both .R files, the output are printed p-values of the test on whether Sparse Regression is adequate.
