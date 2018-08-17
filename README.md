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
- Procedure
	- Enterez ID for all enzymes in the metabolic enzyme list were fed through Uniprot to obtain their corresponding Uniprot ID and their Uniprot annotation, generated file
      - **Human Uniprot.xlsx**
        ```
        #Transform Gene ID from masterlist to Human Uniprot ID
        nature10350_s2$`UniprotID_human` <- 0
        for(i in 1:2752) { 
          for (j in 1:2727) {
            if (nature10350_s2$`Entrez Gene ID`[i] == Human_Uniprot$`Entrez Gene ID`[j]) {
              nature10350_s2$`UniprotID_human`[i] <- Human_Uniprot$Entry[j]
            }
          }
        }
        ```
  - then Uniprot annotation (human) were changed from human to mouse and fed back into Uniprot to obtain the corresponding mouse enzyme Uniprot ID, generated file
    - **Mouse Uniprot.xlsx**
      ```
      #Get the name of the proteins from Human Annotation, change them to mouse annotation, retreive corresponding Mouse Uniprot ID
      Protein_name <- data.frame(do.call('rbind', strsplit(as.character(Human_Uniprot$`Entry name`), '_', fixed = TRUE)))
      Human_Uniprot$Gene <- as.character(Protein_name$X1)
      Mouse_Uniprot$Gene <- as.character(data.frame(do.call('rbind', strsplit(as.character(Mouse_Uniprot$`Entry name`), '_', fixed = TRUE)))[,1])
      Human_Uniprot$Mouse_UniproID <- 1
      for (x in 1:2474){
        for (z in 1:2727) {
          if (Mouse_Uniprot$Gene[x] == Human_Uniprot$Gene[z]) {
            Human_Uniprot$Mouse_UniproID[z] <- Mouse_Uniprot$Entry[x]
          } 
        }
      }

      # Label the masterlist with mouse Uniprot ID
      for (a in 1:2752) {
        for (b in 1:2727) {
          if (masterlist$UniprotID_human[a] == Human_Uniprot$Entry[b]) {
            masterlist$UniprotID_Mouse[a] <- Human_Uniprot$Mouse_UniproID[b]
          }
        }
      }
      ```
  - enzymes with mouse Uniprot ID was used to fetch all overlapped detected phosphopeptides from the phosphopeptide list to generate file
    - **overlapped_phosphopeptides.csv**
  - a list of metabolic enzymes that were detected in phosphoproteomics and described in the master metabolic enzyme list were compiled, generated file, generated 508 proteins
    - **overlapped_metabolic_enzymes.xlsx**
      ```
      #Print out all overlapped proteins between Phospho-peptides and metabolic masterlist
      overlap <- c(NA)
      for (c in 1:21304) {
        for (d in 1:2752) {
          if (phosho_mouse_single_phosphopeptides_filtered$UniprotID[c] == masterlist$UniprotID_Mouse[d]) {
            overlap <- c(overlap, c)
          }
        }
      }
      overlapped_phosphopeptides <- phosho_mouse_single_phosphopeptides_filtered[overlap,]
      ```
  - extracted all proteins categorized as metabolic from the master protein list of SMPDB, **proteins.csv** and their corresponding Uniprot ID
  - list all proteins that were detected in the phosphor-peptide dataset and their Uniprot ID
  - intersect the Unitprot ID from SMPDB dataset and the phosphor-peptide dataset, 141 proteins were identified in both datasets,
  - intersect the 141 proteins with the 508 proteins identified before, 30 new proteins were detected in the new list of metabolic proteins, their Uniprot IDs were exported, generated file
    - **added_enzyme_ID.csv**
      ```
      # since SMPDB has their ownlist of metabolic enzymes, this is to extract the metabolic enzymes from the list
      metab_enzyme_index <- c()
      for (i in 1: dim(proteins)[1]) {
        if (proteins$`Pathway Subject`[i] == "Metabolic") {
          metab_enzyme_index <- c(metab_enzyme_index, i)
        }
      }
      metabolic_enzymes <- proteins[metab_enzyme_index,]
      metabolic_enzymes_ID <- as.character(metabolic_enzymes$`Uniprot ID`)

      # extract all proteins that were detected in the phospho-proteomics dataset
      detected_protein <- as.data.frame(table(phosho_mouse_single_phosphopeptides_filtered$UniprotID))
      detected_protein_ID <- as.character(detected_protein$Var1)

      # intersect the detected protein with all metabolism related protein protein from SMPDB
      library(dplyr)
      new_overlap_enzymes <- intersect(metabolic_enzymes_ID, detected_protein_ID)

      # find the difference between the new_overlap_enzymes and the overlapped_metabolic_enzymes.csv to see 
      # what are the enzymes that are extracted from SMPDB that have not been covered in our list
      previous_enzyme_ID <- as.character(overlapped_metabolic_enzymes$Entry)
      previous_enzymes <- intersect(new_overlap_enzymes, previous_enzyme_ID)
      added_enzyme_SMPDB <- setdiff(new_overlap_enzymes, previous_enzyme_ID)
      ```
  - these Uniprot IDs were fed into Uniprot.org to extract protein information, those information was appended to **overlapped_metabolic_enzymes.xlsx**
  
### Extract metabolites that are up/downstream of all overlapped enzymes
- Database used
  - Small molecule pathway database ([SMPDB](http://smpdb.ca))
- Procedure
  - Uniprot ID for all overlapped enzymes were used to search through SMPDB to extract pathways that involves target enzyme
    - Since the Uniprot ID was for mouse proteins, only mouse metabolic pathways were extracted
  - substrates/substrate HMDB ID/Product/Product HMDB ID were recorded manually, substrates from different pathways were grouped together and products from different pathways were grouped together
    - for multi-step reactions (the same enzyme catalyze multiple steps of a pathway), all intermediate metabolites were recorded twice, once as substrate and once as metabolite
    - for reversible reactions, reversible was marked in the comment section
    - for transporters, metabolites that are being transported are only recorded once in the product section
- all information was recorded in file
  - **overlapped_metabolic_enzymes.xlsx**
 
### Overlap metabolomics data with metabolites for detected enzymes
- this is to remove enzymes that have no detectable metabolite
- Database used
  - Human Metabolite Database ([HMDB](http://www.hmdb.ca)): 
- Metabolomics dataset 
  - **Tissue Metabolomics.xlsx**
- Procedure
  - names of metabolites that are detected in the metabolomics dataset are entered into the HMDB dataset to extract their corresponding HMDB ID, IDs are logged and generated file 
    - **Tissue Metabolomics-HMDB.xlsx**
⁃	for all metabolites connected to enzymes in the overlapped_metabolic_enzymes.xlsx file, each was matched with the metabolomics dataset to remove any undetected metabolites
⁃	Then enzymes with no detected substrates and products were removed from the list
⁃	enzymes with only ATP/ADP/NAD/NADH/NADP/NADPH/NA as substrate and product were removed from the list
⁃	generated file with 87 proteins
⁃	enzyme_metabolite_match.csv

  
