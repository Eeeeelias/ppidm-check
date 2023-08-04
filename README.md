# PPIDM-check

## What's it about?

PPIDM-check aims to replicate & validate the findings of the PPIDM method proposed by Alborzi et al. (2021). PPIDM-check does this by conducting another overlap analysis with experimentally resolved DDIs. 
In addition to that, PPIDM-check rewires the PPIDM-computed DDI network to evaluate it's robustness to false information. As a third and last evaluation method, we constructed an extended version of the PPI/DDI 
graph used by NEASE (Louadi et al., 2021) and compared the updated findings with existing literature. 

To then make use of the predictions, we built another extended DDI graph that can be integrated into DIGGER (Louadi et al., 2021).

---

## Replication

### 1. ppidm_run
Go into the ppidm_run folder and uncomment the functions that you want to run. This will take a long time. Comments 
about what the functions do are included in the `main.py` file. Rudimentary p-value analysis can be done using the 
accompanying jupyter notebook.

### 2. ppidm_validation
These are three separate files that run independently of each other. However, you will need to make pickle files of 
certain files (Namely: `train_set`, `random_train` and `negative_train`). The code to make these files is usually in the 
respective files but commented out. Uncomment those lines to write them to a file. Furthermore, to get each of the 
interaction files needed for the did comparison, look at `PPIDM_pvalue_analysis.ipynb` and execute the code specified 
there.

### 3. digger_extension
This part contains the code to extend the current DIGGER interactions but also some validation using NEASE. To get the
extended Graph run `ddi_network.py` (this will take a while). After that, you can execute `digger_data_explorer` to 
create the actual .pkl for the Graph. Note that you'll have to uncomment the actual dumping of the data.  
For the NEASE valdiation check out the other files in this directory. If you want to run this, be aware that we've 
modified the NEASE source package to make it easier to work with. The base functionality hasn't been changed however!

---

## Other information

This project is part of my bachelors thesis.
