# ADAPT-DSR
This repository houses the code for the AMRx-Liverpool seed project 'Personalised antimicrobial susceptibility reporting'  

The objectives of the project are:  
1. To develop a prototype user interface that will choose how to report antimicrobial susceptibility test results based on patients' past medical history
2. To write up and disseminate the results in a peer-reviewed academic journal

The project components are:
1. A model that predicts the binary outcome of whether a patient with UTI has a **complicated UTI** or **uncomplicated UTI**
2. A model that predicts the binary outcome of whether a patient with bacteraemia has has **UTI** or **does not have UTI**
3. A user interface that uses the above predictions to advise how susceptibility results on urine and blood cultures should be reported

### Processing Steps
filter_cultures.py
label_and_join_data.py
process_cultures.py
add_demographics_from_admissions_data.py
add_demographics_from_patients_data.py
