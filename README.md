# Reproduction for the methodology of "Are Factor Regression and Sparse Regression Adequate?"

# 1. **Real Data Analysis** <br />

## Data  <br />

### 1.1.  Data Description <br />
FRED-MD is a large macroeconomic database designed for the empirical analysis of “big data.” The datasets of monthly and quarterly observations mimic the coverage of datasets already used in the literature, but they add three appealing features. They are updated in real-time through the FRED database. They are publicly accessible, facilitating the replication of empirical work.  <br />
### 1.2. Data Availability <br />
 The original dataset can be found in https://research.stlouisfed.org/econ/mccracken/fred-databases/. <br />
 In our .zip file, the original Dataset is provided in folder ``Code_Reproduction/Real_Data/Pre_process/current.csv." <br />
 The description of the data is provied in "Code_Reproduction/Real_Data/Pre_process/Description_of_Variables.pdf."  <br /> This description describes the meaning of all columns of this datset.

## Codes Description and Implementation <br />
 There are in total four folders under the folder: ``Code_Reproduction/Real_Data."  <br /> 1. Pre_process, 2. Prediction_Table 3, 3. PCR_Adequate_Table4, 4. Sparse_Adequate_Table4. <br />
### For Folder `Code_Reproduction/Real_Data/Pre_process':
 It contains the Codes for real data pre-processing. The document `current.csv' in that folder is our downloaded orignal dataset from the aforementioned website. <br />
 
 The code "realdata_read.R" contains the codes for pre-processing the data. The line 4 in "realdata_read.R", we have a variable called: "choose_time". If one lets choose_time=1 and run the codes, the output is the processed dataset "realdata_115_126.csv". If one lets choose_time=2, and run the codes, the output is the processed dataset "realdata_189_126.csv."<br />
 
 Finally, for document "Description_of_Variables.pdf" inside that document, it contains the detailed description of all variables in the orignal dataset "current.csv". <br />


 
