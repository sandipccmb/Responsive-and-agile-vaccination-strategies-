Responsive and agile vaccination strategies against COVID-19 in India.

Authors: Sandip Mandal; Nimalan Arinaminpathy; Balram Bhargava; Samiran Panda
 
This repository contains codes and data used to simulate and analyze COVID-19 vaccination strategies in India.
1.	Impact of rapidly responsive vaccination strategies 
2.	Variation of impact of responsive vaccination and NPIs with threshold TPR 
Note: The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat), with plots also being constructed in MATLAB. 

 OS System requirements
The codes developed here are tested on Windows operating system (Windows 10 Home: 64-bit). However as Matlab and Octave are available for most operating systems, codes should run on Mac OSX and Linux as well.

Installation guide
MATLAB
Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. The current codes were developed and tested on MATLAB R2018b.

Codes & their functionality

Setup_model.m
This model has six state variables which are further divided into three age groups, two risk group (co-morbid or not) and two groups according to their vaccination status. All the variables and parameter values are assigned in this script. 

Make_model2.m
This is a function which specifies the transmission model in matrix form, capturing the linear and non-linear parts separately.

goveqs_basis3.m
Using the matrix formulation constructed in Make_model2, this computes the time derivatives for each state variable given values for those states.

 linspecer.m
This function is used to plot multiple lines with distinguishable and nice looking colors.

jbfill.m
This routine will shade the area of a 2-D plot between two user defined vectors. 
(Ref. John Bockstege (2021). Shade area between two curves (https://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves), MATLAB Central File Exchange. Retrieved April 2, 2021.)

get_init.m
This function assigns the initial conditions for simulation, depending on vaccination coverage, and existing seroprevalence. 

get_aggregators.m
Function input ind_sets is a cell array, each element a vector showing the set of indices that need to be aggregated over.

get_address.m
We make use of Matlab ‘structures’ (the equivalent of ‘lists’ in R) to assist in book-keeping on the indices corresponding to different state variables. This function assists in constructing those indices. For example, the index named i.A.v0.r1.a3 is the index corresponding to the asymptomatic, non-vaccinated, with co-morbidity and among elderly population. 

Alloc_parameters.m
For a given parameter vector x, this function allocates the corresponding values to each of the parameters being sampled.

Instructions for use
In the above titled article, there is one figure (Figure 1) in the main text and one figure (Figure S2) in the supplementary document. To generate these figures and output table (Table S3) select the appropriate script and run it in MATLAB. Save the data file as instructed at the bottom of each script and then run the file named as ‘Figure X.m’ to get the figure. Parameters should be adjustment for infection preventing vaccine and symptomatic disease preventing vaccine. 

Setup_model.m
Model variables and parameter values are assigned in this script. To setup the model for infection preventing vaccine choose the parameter value p.c1 = 0.6 (when vaccine efficacy is 60%) and p.c3 = 0. For disease preventing vaccine p.c3 = 0.6 (when vaccine efficacy is 60%) and p.c1 = 0;

Find_Pn.m
To find the prevalence of non-COVID symptoms run this code for a particular value of R0. This is estimated from the fact that the TPR during the peak epidemic of the 1st wave in India was about 10%.  

Simulate1.m
After setting up the model, run this code to get data for the plot “Figure 1B” of this article. Assign R0 = 2.0, Pn = 5298. Vaccination is done uniformly among 75% of the adult population (>18y) in one month. Save the output as ‘simulate_sev.mat’.

Simulate2.m
To get the data for Figure S2 of the above mentioned article, (i.e. to see the change in impact by varying threshold TPR) run this code. This figure requires four data set - for two different R0 values (R0 = 2 and 2.5) with two types of vaccination (disease preventing and infection preventing vaccine). First setup the model as per the nature of vaccination (Setup_model.m) and then for each run of this code save the data as below:
simulate2_sev.mat: (for R0 = 2, disease preventing vaccine)     
simulate2p5_sev.mat: (for R0 = 2.5, disease preventing vaccine)     
simulate2_sus.mat: (for R0 = 2, infection preventing vaccine)     
simulate2p5_sus.mat: (for R0 = 2.5, infection preventing vaccine)     

Figure1.m
Two panel figure plotted from the data stored as ‘districts.xlsx’ and ‘simulate_sev.mat’.

Figure1 A: illustrates the TPR and cases in Delhi from 17 December 2020 to 27 March 2021.  
Figure1 B: illustrates daily deaths per 100,000 population for three scenarios: (i) no vaccination, no restrictions, (ii) vaccination, no restrictions, (iii) vaccination, with 25% reduction in transmission at same time as vaccination is being deployed (and then restored to 100% straight after). This is shown for a disease-preventing vaccine.


FigureS2.m
Two panel figure plotted from the data stored as ‘simulate2_sev.mat’, ‘simulate2p5_sev.mat’, ‘simulate2_sus.mat’, ‘simulate2p5_sus.mat’.

Figure S2 A: Illustrate deaths averted vs threshold TPR, for a disease-preventing vaccine for R0 = 2 and R0 = 2.5.
Figure S2 B: Similar as Figure S2 A, but for an infection-preventing vaccine. 


TableS3.m
For each of the scenario choose appropriate data set and find the value of percent reduction of cumulative cases (here named as z or z1) or deaths (defined as x or x1). First column of each of these vector is corresponding to threshold TPR = 0.5%.  



MAT files
simulate2_sev.mat 
simulate2p5_sev.mat 
simulate2_sus.mat 
simulate2p5_sus.mat 

Excel file
districts.xlsx

