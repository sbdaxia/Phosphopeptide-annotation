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
	```
	#remove metabolites with no HMDB ID
	emptyrow <- which(is.na(HMDB_ID$`HMDB ID`))
	HMDB_ID <- HMDB_ID[-emptyrow,]

	#remove proteins in the enzyme list that had no detected product from SMPDB
	emptyrow_enzyme <- which(is.na(enzyme_metabolite_match$`Product HMDB`))
	enzyme_metabolite_match <- enzyme_metabolite_match[-emptyrow_enzyme,]

	#fill all empty substrate HMDB cells with 9999999
	enzyme_metabolite_match[is.na(enzyme_metabolite_match$`Substrate HMDB`), "Substrate HMDB"] <- "9999999"
	enzyme_metabolite_match[is.na(enzyme_metabolite_match$`Substrate`), "Substrate"] <- "9999999"

	#remove empty columns
	enzyme_metabolite_match$X13 = NULL
	enzyme_metabolite_match$X14 = NULL
	```
- for all metabolites connected to enzymes in the **overlapped_metabolic_enzymes.xlsx** file, each was matched with the metabolomics dataset to remove any undetected metabolites
	```
	#overlapp substrates in the enzyme list with metabolite list to remove undetected substrates
	for (n in 1:dim(enzyme_metabolite_match)[1]) {
	  substrate <- strsplit(enzyme_metabolite_match$Substrate[n], ", ")[[1]]
	  substrateid <- strsplit(enzyme_metabolite_match$`Substrate HMDB`[n], ", ")[[1]]
	  newsubstrate <- c()
	  newsubstrateid <- c()
	  for (i in 1:length(substrateid)) {
	    if (sum(grepl(substrateid[i], HMDB_ID$`HMDB ID`)) > 0) {
	      newsubstrateid <- c(newsubstrateid, substrateid[i])
	      newsubstrate <- c(newsubstrate, substrate[i])
	    }
	  }
	  enzyme_metabolite_match$Substrate[n] <- paste(newsubstrate, collapse =", ")
	  enzyme_metabolite_match$`Substrate HMDB`[n] <- paste(newsubstrateid, collapse = ", ")
	}

	#overlapp products in the enzyme list with metabolite list to remove undetected products
	for (n in 1:dim(enzyme_metabolite_match)[1]) {
	  product <- strsplit(enzyme_metabolite_match$Product[n], ", ")[[1]]
	  productid <- strsplit(enzyme_metabolite_match$`Product HMDB`[n], ", ")[[1]]
	  newproduct <- c()
	  newproductid <- c()
	  for (i in 1:length(productid)) {
	    if (sum(grepl(productid[i], HMDB_ID$`HMDB ID`)) > 0) {
	      newproductid <- c(newproductid, productid[i])
	      newproduct <- c(newproduct, product[i])
	    }
	  }
	  enzyme_metabolite_match$Product[n] <- paste(newproduct, collapse =", ")
	  enzyme_metabolite_match$`Product HMDB`[n] <- paste(newproductid, collapse = ", ")
	}
	```
- Then enzymes with no detected substrates and products were removed from the list
	```
	#remove enzymes that do not have detected substrates and detected products
	emptylist <- c()
	for (n in 1:dim(enzyme_metabolite_match)[1]) {
	  if ((enzyme_metabolite_match$`Substrate HMDB`[n] == '') && (enzyme_metabolite_match$`Product HMDB`[n] == '')) {
	    emptylist <- c(emptylist, n)
	  }
	}

	enzyme_metabolite_match <- enzyme_metabolite_match[-emptylist, ]
	```
- enzymes with only ATP/ADP/NAD/NADH/NADP/NADPH/NA as substrate and product were removed from the list
	```
	# remove emzymes that only have ATP, ADP, NAD, NADH, NADPH, NADP as substrate and metabolite
	# ATP = 0000538
	# ADP = 0001341
	# NAD = 0000902
	# NADH = not detected in metabolomics 
	# NADPH = not detected in metabolomics 
	# NADP = not detected in metabolomcis
	useless_substrate <- c()
	for (i in 1: dim(enzyme_metabolite_match)[1]) {
	  if (enzyme_metabolite_match$`Substrate HMDB`[i] == "0000538" | enzyme_metabolite_match$`Substrate HMDB`[i] == "0001341" | enzyme_metabolite_match$`Substrate HMDB`[i] == "0000902" | is.na(enzyme_metabolite_match$`Substrate HMDB`[i])) {
	    useless_substrate <- c(useless_substrate, i)
	  }
	}

	useless_product <- c()
	for (n in 1: dim(enzyme_metabolite_match)[1]) {
	  if (enzyme_metabolite_match$`Product HMDB`[n] == "0000538" | enzyme_metabolite_match$`Product HMDB`[n] == "0001341" | enzyme_metabolite_match$`Product HMDB`[n] == "0000902" | is.na(enzyme_metabolite_match$`Product HMDB`[n])) {
	    useless_product <- c(useless_product, n)
	  }
	}

	useless_proteins <- intersect(useless_substrate, useless_product)

	enzyme_metabolite_match <- enzyme_metabolite_match[-useless_proteins,]

	write.csv(enzyme_metabolite_match, "enzyme_metabolite_match.csv")
	```
- generated file with 87 proteins
	- **enzyme_metabolite_match.csv**
	
### Overlap metabolomics-phosphoproteomics-metabolic enzymes
- Procedure
	- enzymes in the enzyme_metabolite_match.csv were used to extract phosphopeptides that were detected of these enzymes, generated file
		- **phosphopeptide_enzyme_metabolite_overlap.csv**
	```
	#extract peptides that belong to enzymes that have detectable metabolites
	overlap_index <- c()
	for (i in 1: dim(phosho_mouse_single_phosphopeptides_filtered)[1]) {
	  if (sum(grepl(phosho_mouse_single_phosphopeptides_filtered$UniprotID[i], enzyme_metabolite_match$Entry)) > 0) {
	    overlap_index <- c(overlap_index, i)
	  }
	}

	phosphopeptide_enzyme_metabolite_overlap <- phosho_mouse_single_phosphopeptides_filtered[overlap_index,]

	write.csv(phosphopeptide_enzyme_metabolite_overlap, "phosphopeptide_enzyme_metabolite_overlap.csv")
	```

### Z-scoring overlapped phosphorus-peptide quantification and metabolomics quantification
- Starting datasets
	- **phosphopeptide_enzyme_metabolite_overlap.csv**
	- **Tissue Metabolomics.xlsx**
#### Z-scoring the phospho-peptide dataset
- quantification from the three replicate runs were combined by adding up for each corresponding tissue site
	- standard z-scoring was performed across tissue sites and z-scores were recorded in the new dataset that contains only the overall quantification and z-scores, generated
		- **zscore_phosphopeptide.csv**
```
library(mosaic)

zscore_phosphopeptide <- phosphopeptide_enzyme_metabolite_overlap[,2:7]
zscore_phosphopeptide <- cbind(zscore_phosphopeptide, phosphopeptide_enzyme_metabolite_overlap[,15:18])
zscore_phosphopeptide$kidney <- NA
zscore_phosphopeptide$liver<- NA
zscore_phosphopeptide$spleen<- NA
zscore_phosphopeptide$lung<- NA
zscore_phosphopeptide$brain<- NA
zscore_phosphopeptide$heart<- NA
zscore_phosphopeptide$brown_adipose_tissue<- NA
zscore_phosphopeptide$pancreas<- NA
zscore_phosphopeptide$white_adipose_tissue<- NA
zscore_phosphopeptide$testes<- NA

for (i in 1: dim(phosphopeptide_enzyme_metabolite_overlap)[1]) {
  zscore_phosphopeptide[i,11:20] <- phosphopeptide_enzyme_metabolite_overlap[i, 23:32] + phosphopeptide_enzyme_metabolite_overlap[i, 34:43] + phosphopeptide_enzyme_metabolite_overlap[i,45:54]
}

tablecolnames <- c("kidney_zscore", "liver_zscore", "spleen_zscore", "lung_zscore", "brain_zscore", "heart_zscore", "brown_adipose_tissue_zscore", "pancreas_zscore", "white_adipose_tissue_zscore", "testes_zscore")

zscore1 <-zscore(as.numeric(zscore_phosphopeptide[1, 11:20]))
zscore_table <- as.data.frame(rbind(zscore1))
for (n in 2: dim(zscore_phosphopeptide)[1]) {
  organzscore <- zscore(as.numeric(zscore_phosphopeptide[n, 11:20]))
 zscore_table <- as.data.frame(rbind(zscore_table, organzscore))
}

names(zscore_table) <- tablecolnames
zscore_phosphopeptide <- cbind(zscore_phosphopeptide,zscore_table)
write.csv(zscore_phosphopeptide, "zscore_phosphopeptide.csv")
```
#### Z-scoring the metabolomics dataset
- metabolomics quantification normalized to protein amount from the Tissue Metabolomics.xlsx was extracted and transposed into a csv file so the metabolites are the row names and tissue samples were the column names, this is done in excel since transposing dataframe in R is quite tricky
- metabolites that are doubled were split up, with each carry the same quantification, HMDB annotation and ID were put in the same dataset, generated file
	- **metabolomitc_HMDB_annotated.csv**
- from metabolomitc_HMDB_annotated.csv, all quantifications from one of the six replicates of each organ were added up together, zscoring was performed across all tissue sites and z-scores were recorded, generated file 
	- zscore_metabolomics.csv
```
library(mosaic)

# extracts the names and HMDB annotations from the metabolomics dataset
zscore_metabolomics <- metabolomitc_HMDB_annotated[, 1:3]

#take and six replicate measurement for each organ site and add them all together for the final zscore calculation
Brain <- metabolomitc_HMDB_annotated[,4] + metabolomitc_HMDB_annotated[,5] + metabolomitc_HMDB_annotated[,6] + metabolomitc_HMDB_annotated[,7] + metabolomitc_HMDB_annotated[,8] + metabolomitc_HMDB_annotated[,9]
Pancreas <- metabolomitc_HMDB_annotated[,10] + metabolomitc_HMDB_annotated[,11] + metabolomitc_HMDB_annotated[,12] + metabolomitc_HMDB_annotated[,13] + metabolomitc_HMDB_annotated[,14] + metabolomitc_HMDB_annotated[,15]
Lung <- metabolomitc_HMDB_annotated[,16] + metabolomitc_HMDB_annotated[,17] + metabolomitc_HMDB_annotated[,18] + metabolomitc_HMDB_annotated[,19] + metabolomitc_HMDB_annotated[,20] + metabolomitc_HMDB_annotated[,21]
Liver <- metabolomitc_HMDB_annotated[,22] + metabolomitc_HMDB_annotated[,23] + metabolomitc_HMDB_annotated[,24] + metabolomitc_HMDB_annotated[,25] + metabolomitc_HMDB_annotated[,26] + metabolomitc_HMDB_annotated[,27]
SkM <- metabolomitc_HMDB_annotated[,28] + metabolomitc_HMDB_annotated[,29] + metabolomitc_HMDB_annotated[,30] + metabolomitc_HMDB_annotated[,31] + metabolomitc_HMDB_annotated[,32] + metabolomitc_HMDB_annotated[,33]
Kidney <- metabolomitc_HMDB_annotated[,34] + metabolomitc_HMDB_annotated[,35] + metabolomitc_HMDB_annotated[,36] + metabolomitc_HMDB_annotated[,37] + metabolomitc_HMDB_annotated[,38] + metabolomitc_HMDB_annotated[,39]
organlist <- c(Brain, Pancreas, Lung, Liver, SkM, Kidney)
zscore_metabolomics <- cbind(zscore_metabolomics, organlist)
columnnames <- c("Brain", "Pancreas", "Lung", "Liver", "SkM", "Kidney")
names(zscore_metabolomics)[4:9] <- columnnames

#calculate the zscores for each metabolite
zscorelist <- data.frame()
for (i in 1:dim(zscore_metabolomics)[1]) {
  met_zscore <- zscore(as.numeric(zscore_metabolomics[i, 4:9]))
  zscorelist <- rbind(zscorelist, met_zscore)
}
zscore_col_names <- c("Brain_zscore", "Pancreas_zscore", "Lung_zscore", "Liver_zscore", "SkM_zscore", "Kidney_zscore")
names(zscorelist) <- zscore_col_names

#combine calculated zscores with the quantified metabolomics dataset
zscore_metabolomics <- cbind(zscore_metabolomics, zscorelist)

write.csv(zscore_metabolomics, "zscore_metabolomics.csv")
```



# Authors
Bing Shui
