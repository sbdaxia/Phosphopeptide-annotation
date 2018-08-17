# Metabolic Enzyme Phosphosite Functional Annotation

The overarching question: **how to perform high-throughput functional annotation of phopho-sites of proteins in the phospho-proteomics dataset**

## Experimental concepts
- Using metabolites as functional outputs, using metabolome as the high-throughput readouts
- Using 10 tissue sites as different perturbations
- Using changes in metabolites to predict the function of up-steam or down-stream enzymatic functions
- Identify enzymatic changes that correlate with phospho-proteomics changes
- Use phospho-proteomics change to then predict functions
- Look for motif homology and identify potential kinase that regulate function of certain phosphorylation
- Experimentally validate this in cell or animal models, confirm functional output of certain phosphorylation

## Development
### Overlap metabolic enzyme list with phospho-peptide list
- Metabolic enzyme list used
  - [**nature10350-s2.xls**](https://www.nature.com/articles/nature10350)
  - **proteins.csv** from [SMPDB dataset](https://smpdb.ca) 
- Phospho-peptide list used 
  - phosho-mouse_single-phosphopeptides_filtered.xls
  - Generated from HMS Gygi lab
