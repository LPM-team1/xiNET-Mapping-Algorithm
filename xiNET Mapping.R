

###############################################################
#  Crosslink Analyzer for Xinet _ Proteome Disocover / Plink  #
###############################################################


#Activation of libraries
library(tidyverse)
library(dplyr)
library(Biostrings)
library(data.table)
library(stringr)
library(seqinr)
library(readr)
library(purrr)


#Directory of Archives:
setwd('M:/ACADEMIC/YourPath/NPQ_Curves_XiNET')

#Directory of database:
Database<-'YourDataBase.fasta'

#Directory of files:
File<-'Input_Crosslink'

#Output for Files:
#Final archive localization:
Output_file<-'XL-Analysis'

#Name of final files:
Name_Output<-'Xinet_crosslink'

#------------------------------------------#
#Is a PLink file?: 'Yes'/'No'(PD)
PLink <- 'No'

#Filter at beginning by Score? 'Yes'/'No'
Filtering<-'Yes'

Min_Score<-40

#Filter by cross links repeated in replicas by treatments or in total? 'Treatment'/'Total'
Replicas_in<-'Treatment'

#How many replicas ?
Minimun_replicas<-2

##All in one:
#Name of files per condition:
Subjected_File<-'NPQv1_SoHL_08-08-25_CSMsCSVcomma_consensus1_matched_psms'
Format<-'.csv'

# Name of files per condition:
treatment_groups <- list(
  "NPQ_00-05" = c("Xpl1_010533", "Xpl1_010586", "Xpl1_010587"),
  "NPQ_05-10" = c("Xpl1_010535", "Xpl1_010577", "Xpl1_010578"),
  "NPQ_10-15" = c("Xpl1_010536", "Xpl1_010580", "Xpl1_010581"),
  "NPQ_15-20" = c("Xpl1_010538", "Xpl1_010583", "Xpl1_010584"))

#End of the files names:
File_type<-".raw"

###########################################################################################
#Start of the code:

# Define the main directory name with the specified path
main_dir <- file.path(getwd(), Output_file)
subdirs <- c("Crosslink_Peptides",
             "Crosslink_Peptides_Reasigned",
             "Xinet_by_Exlusive&Common")
# Create directories if it does not exist
if (!dir.exists(main_dir)) {
  dir.create(main_dir, recursive = TRUE)
}
for (subdir in subdirs) {
  subdir_path <- file.path(main_dir, subdir)
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path)
  }
}
# Confirmation message
cat("Directories created successfully.\n")
#Processing names and tretaments:
Replica_Names <- sub("\\.raw$", "", unlist(treatment_groups))
treatment_groups <- treatment_groups[sapply(treatment_groups, length) >= Minimun_replicas]
Experimental_Conditions <- rep(names(treatment_groups), times = sapply(treatment_groups, length))
Unique_Conditions <- names(treatment_groups)

# Convert Plink files
if (PLink == 'Yes') {
  Min_Score <- 1 / as.numeric(Min_Score)
  ruta <- paste0(getwd(),
                 .Platform$file.sep,
                 File,
                 .Platform$file.sep,
                 Subjected_File,
                 Format)
  PLink_OrigTable <- read.csv(ruta)
  TitlesFist <- names(PLink_OrigTable)
  TitlesSecond <- PLink_OrigTable [1, ]
  titlesFile <- paste(TitlesFist, TitlesSecond, sep = "|")
  PLink_OrigTable <- PLink_OrigTable[-1, ]
  
  colnames(PLink_OrigTable) <- c(
    "Peptide_Order",
    "Peptide|Spectrum_Order",
    "Peptide_Mass|Title",
    "Modifications|Charge",
    "Proteins|Precursor_Mass",
    "Protein_Type|Evalue",
    "Score",
    "Precursor_Mass_Error(Da)",
    "Precursor_Mass_Error(ppm)",
    "Alpha_Evalue",
    "Beta_Evalue"
  )
  Conv_Table <- PLink_OrigTable %>%
    mutate(group_id = cumsum(!is.na(`Peptide_Order`)))
  
  # Step 2: Extract the rows with group-level metadata
  group_info <- Conv_Table %>%
    filter(!is.na(`Peptide_Order`)) %>%
    select(
      group_id,
      Group_number = `Peptide_Order`,
      Group_Peptide = `Peptide|Spectrum_Order`,
      Group_Peptide_Mass = `Peptide_Mass|Title`,
      Group_Modifications = `Modifications|Charge`,
      Group_Proteins = `Proteins|Precursor_Mass`,
      Group_Protein_Type = `Protein_Type|Evalue`
    ) %>%
    mutate(
      `Modifications A` = "",
      `Modifications B` = "",
      `Sequence A` = str_extract(`Group_Peptide`, "^[^()]+"),
      `Crosslinker Position A` = as.numeric(str_extract(
        `Group_Peptide`, "(?<=\\()[0-9]+(?=\\))"
      )),
      `Sequence B` = str_extract(`Group_Peptide`, "(?<=-)[^-()]+"),
      `Crosslinker Position B` = as.numeric(
        str_extract_all(`Group_Peptide`, "(?<=\\()[0-9]+(?=\\))") %>%
          map_chr( ~ if (length(.) >= 2)
            .[2]
            else
              NA_character_) %>%
          as.numeric()
      ),
      `Protein Accession A` = str_extract(Group_Proteins, "^[^(]+(?= \\()"),
      `Leading Protein Position A` = as.numeric(str_extract(
        Group_Proteins, "(?<= \\()[0-9]+(?=\\))"
      )),
      `Protein Accession B` = str_extract(Group_Proteins, "(?<=\\)-)[^(]+(?= \\()"),
      `Leading Protein Position B` = as.numeric(
        str_extract_all(Group_Proteins, "(?<= \\()[0-9]+(?=\\))") %>%
          map_chr( ~ if (length(.) >= 2)
            .[2]
            else
              NA_character_) %>%
          as.numeric()
      )
    )
  
  # Separate and assign information from peptide modifications
  for (i in 1:nrow(group_info)) {
    if (!(group_info$Group_Modifications[i] %in% c("Null", "null", "NA", "", NA))) {
      mods <- unlist(strsplit(group_info$Group_Modifications[i], ";"))
      Modif_A <- ""
      Modif_B <- ""
      for (mod in mods) {
        Modif_number <- as.numeric(str_extract(mod, "(?<=\\]\\()[0-9]+(?=\\))"))
        if (Modif_number > nchar(group_info$`Sequence A`[i])) {
          Modif_number <- Modif_number - nchar(group_info$`Sequence A`[i]) - 3
          Modif_string <- paste0(as.character(str_extract(mod , "^[^\\]]*\\]\\(")), Modif_number, ")")
          Modif_B <- ifelse(
            Modif_B == "" | is.na(Modif_B),
            Modif_string,
            paste(Modif_B, " ; ", Modif_string)
          )
        } else {
          Modif_A <- ifelse(Modif_A == "" | is.na(Modif_A),
                            mod,
                            paste0(Modif_A, " ; ", mod))
        }
      }
      group_info$`Modifications A`[i] <- Modif_A
      group_info$`Modifications B`[i] <- Modif_B
    } else {
      group_info$`Modifications A`[i] <- ""
      group_info$`Modifications B`[i] <- ""
    }
  }
  
  
  # Step 3: Join group metadata onto the full dataset
  Conv_Table <- Conv_Table %>%
    left_join(group_info, by = "group_id") %>%
    filter(is.na(`Peptide_Order`)) %>%
    mutate(
      combined_id = paste0(Group_number, ".", `Peptide|Spectrum_Order`, " | ", group_id),
      File_ID = paste0(sub("\\..*$", "", `Peptide_Mass|Title`), ".raw"),
      `Sequence A` = str_extract(`Group_Peptide`, "^[^()]+"),
      `Crosslinker Position A` = as.numeric(str_extract(
        `Group_Peptide`, "(?<=\\()[0-9]+(?=\\))"
      )),
      `Sequence B` = str_extract(`Group_Peptide`, "(?<=-)[^-()]+"),
      `Crosslinker Position B` = as.numeric(
        str_extract_all(`Group_Peptide`, "(?<=\\()[0-9]+(?=\\))") %>%
          map_chr( ~ if (length(.) >= 2)
            .[2]
            else
              NA_character_) %>%
          as.numeric()
      ),
      `Protein Accession A` = str_extract(Group_Proteins, "^[^(]+(?= \\()"),
      `Leading Protein Position A` = as.numeric(str_extract(
        Group_Proteins, "(?<= \\()[0-9]+(?=\\))"
      )),
      `Protein Accession B` = str_extract(Group_Proteins, "(?<=\\)-)[^(]+(?= \\()"),
      `Leading Protein Position B` = as.numeric(
        str_extract_all(Group_Proteins, "(?<= \\()[0-9]+(?=\\))") %>%
          map_chr( ~ if (length(.) >= 2)
            .[2]
            else
              NA_character_) %>%
          as.numeric()
      )
    )
  
  
  #Create final table
  OrigTable <- Conv_Table %>%
    transmute(
      Checked = "Unknown",
      Sequence = `Group_Peptide`,
      Crosslinker = "Unknown",
      `Crosslink Type` = `Group_Protein_Type`,
      `Crosslink Strategy` = "Unknown",
      `Identified By` = "Unknown",
      `# Proteins` = "Unknown",
      `# Identified MS2 Scans` = 1,
      `XlinkX Score` = 1 / as.numeric(`Score`),
      `Delta XlinkX Score` = "Unknown",
      `First Scan` = `Peptide_Mass|Title`,
      `RT [min]` = "NA",
      `m/z [Da]` = as.numeric(`Proteins|Precursor_Mass`) / as.numeric(`Modifications|Charge`),
      `Charge` = `Modifications|Charge`,
      `MH+ [Da]` = `Proteins|Precursor_Mass`,
      `DeltaM [ppm]` = `Precursor_Mass_Error(ppm)`,
      
      # Plink Info:
      `DeltaM (Da)` = `Precursor_Mass_Error(Da)`,
      `Alpha_Evalue` = `Alpha_Evalue`,
      `Beta_Evalue` = `Beta_Evalue`,
      `Group Peptide Mass` = `Group_Peptide_Mass`,
      `Group_XL.Number_XL | Group` = `combined_id`,
      `Reporter Ion Score` = "Unknown",
      `Comp. Voltage` = "Unknown",
      `Accessions` = paste0(`Protein Accession A`, " - ", `Protein Accession B`),
      `Sequence A` = paste0(`Sequence A`),
      `Modifications A` = `Modifications A`,
      `Crosslinker Position A` = as.numeric(`Crosslinker Position A`),
      `Protein Accession A` = `Protein Accession A`,
      `Leading Protein Position A` = `Leading Protein Position A`,
      `Sequence B` = paste0(`Sequence B`),
      `Modifications B` = `Modifications B`,
      `Crosslinker Position B` = as.numeric(`Crosslinker Position B`),
      `Protein Accession B` = `Protein Accession B`,
      `Leading Protein Position B` = `Leading Protein Position B`,
      `All Scans` = "Unknown",
      `Spectrum File` = `File_ID`,
      `File ID` = "Unknown",
      `Validator Q-value` = "Unknown",
      `Quan Info` = "Unknown",
      `Precursor Abundance` = "Unknown",
      `Apex RT [min]` = "Unknown"
    ) %>%
    # Remove all columns where all values == "Unknown"
    select(where( ~ !all(.x == "Unknown")))
} else {
  ruta <- paste0(getwd(),
                 .Platform$file.sep,
                 File,
                 .Platform$file.sep,
                 Subjected_File,
                 Format)
  OrigTable <- readr::read_csv(ruta)
}

#Processing:
#order peptides:
OrigTable <- OrigTable %>%
  # Separar la columna Sequence en dos partes
  separate(
    Sequence,
    into = c("swap.Seq1", "swap.Seq2"),
    sep = "-",
    remove = FALSE
  ) %>%
  rowwise() %>%
  mutate(
    swap.ProtA = `Protein Accession A`,
    swap.SeqA  = `Sequence A`,
    swap.ModA  = `Modifications A`,
    swap.PosA  = `Crosslinker Position A`,
    swap.LeadA = `Leading Protein Position A`,
    swap   = replace_na(swap.Seq1 > swap.Seq2, FALSE)
  ) %>%
  #Do a swap on columns if true:
  mutate(
    # Accesions
    `Protein Accession A` = if_else(swap, `Protein Accession B`, `Protein Accession A`),
    `Protein Accession B` = if_else(swap, swap.ProtA, `Protein Accession B`),
    # Swap of Sequences
    `Sequence A` = if_else(swap, `Sequence B`, `Sequence A`),
    `Sequence B` = if_else(swap, swap.SeqA, `Sequence B`),
    # Midififations
    `Modifications A` = if_else(swap, `Modifications B`, `Modifications A`),
    `Modifications B` = if_else(swap, swap.ModA, `Modifications B`),
    # Xl position
    `Crosslinker Position A` = if_else(swap, `Crosslinker Position B`, `Crosslinker Position A`),
    `Crosslinker Position B` = if_else(swap, swap.PosA, `Crosslinker Position B`),
    # Protein Positions
    `Leading Protein Position A` = if_else(swap, `Leading Protein Position B`, `Leading Protein Position A`),
    `Leading Protein Position B` = if_else(swap, swap.LeadA, `Leading Protein Position B`),
    # rebuild Sequence column
    Sequence = paste(
      gsub("\\[|\\]", "", `Sequence A`),
      gsub("\\[|\\]", "", `Sequence B`),
      sep = "-"
    ),
    Accessions = paste(`Protein Accession A`, `Protein Accession B`, sep = "-")
  ) %>%
  ungroup() %>%
  select(-starts_with("swap"))




#Separate conditions, generate global files per condition, and clean of non important data:
ModifiedTable <- OrigTable %>%
  filter(if (Filtering == 'Yes')
    `XlinkX Score` >= Min_Score
    else
      TRUE) %>%
  group_by(`Spectrum File`) %>%
  filter(any(startsWith(
    as.character(`Spectrum File`), Replica_Names
  ))) %>%
  ungroup() %>%
  mutate(
    `Total Precursor Abundances` = ifelse(
      exists("Precursor Abundance", where = .),
      rowSums(select(., starts_with(
        "Precursor Abundance"
      )), na.rm = TRUE),
      NA
    ),
    `Total CSMs` = ifelse(
      exists("# Identified MS2 Scans", where = .),
      rowSums(select(., starts_with(
        "# Identified MS2 Scans"
      )), na.rm = TRUE),
      NA
    ),
    `Experimental group` = ifelse(
      sub(paste0('\\', File_type, '.*$'), "", `Spectrum File`) %in% Replica_Names,
      Experimental_Conditions
      [match(sub(paste0('\\', File_type, '.*$'), "", `Spectrum File`), Replica_Names)],
      NA
    ),
    #Change position in XL peptides from 0 to 1: Identify 1 number per amino acid, included when N-term is XL.
    `Crosslinker Position A` = ifelse(`Crosslinker Position A` == 0, 1, `Crosslinker Position A`),
    `Crosslinker Position B` = ifelse(`Crosslinker Position B` == 0, 1, `Crosslinker Position B`)
  ) %>%
  
  group_by(`Spectrum File`,
           `Sequence`,
           `Crosslinker Position A`,
           `Crosslinker Position B`) %>%
  mutate(
    `Total Precursor Abundances` = ifelse(
      exists("Precursor Abundance", where = .),
      sum(`Precursor Abundance`, na.rm =
            TRUE),
      NA
    ),
    `Total CSMs` = ifelse(
      exists("# Identified MS2 Scans", where = .),
      sum(`# Identified MS2 Scans`, na.rm = TRUE),
      NA
    )
  ) %>%
  slice_max(`XlinkX Score`, with_ties = FALSE) %>%
  ungroup()

#Change `XlinkX Score` for `Max.XlinkX Score`:
names(ModifiedTable)[names(ModifiedTable) == "XlinkX Score"] <- "Max. XlinkX Score"
#Extract data:
for (A in 1:length(ModifiedTable$Sequence)) {
  Peptide_string <- ModifiedTable$`Sequence A`[A]
  XL_position <- as.numeric(ModifiedTable$`Crosslinker Position A`[A])
  ModifiedTable$`Sequence A`[A] <- paste0(
    substring(Peptide_string, 1, XL_position - 1),
    "[",
    substring(Peptide_string, XL_position, XL_position),
    "]",
    substring(Peptide_string, XL_position + 1)
  )
  Peptide_string <- ModifiedTable$`Sequence B`[A]
  XL_position <- as.numeric(ModifiedTable$`Crosslinker Position B`[A])
  ModifiedTable$`Sequence B`[A] <- paste0(
    substring(Peptide_string, 1, XL_position - 1),
    "[",
    substring(Peptide_string, XL_position, XL_position),
    "]",
    substring(Peptide_string, XL_position + 1)
  )
}
name_out <- paste0('.Xlinks_All_Conditions_All_Rep_',
                   Sys.Date(),
                   "_",
                   Subjected_File,
                   ".csv")
rute_out_b <- paste0(
  getwd(),
  .Platform$file.sep,
  Output_file,
  .Platform$file.sep,
  "Crosslink_Peptides",
  .Platform$file.sep,
  name_out
)
write.csv(ModifiedTable, file = rute_out_b, row.names = FALSE)

# Spreading of XLinks data:
# 'Max. XlinkX Score':
spread_XlinkX_Score <- ModifiedTable %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Max. XlinkX Score`
  ) %>%
  spread(key = `Spectrum File`, value = `Max. XlinkX Score`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )), ~ paste0("XlinkX_Score_", .))
# 'Total CSMs':
spread_Total_CSMs <- ModifiedTable %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Total CSMs`
  ) %>%
  spread(key = `Spectrum File`, value = `Total CSMs`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )), ~ paste0("Total_CSM_", .))
# 'Total Precursor Abundances':
spread_Total_Precursor_Abundances <- ModifiedTable %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Total Precursor Abundances`
  ) %>%
  spread(key = `Spectrum File`, value = `Total Precursor Abundances`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )), ~ paste0("Total_Precursor_Abundances_", .))

# Join the spread data frames on 'Sequence':
spread_ModifiedTable <- ModifiedTable %>%
  group_by(`Sequence`, `Crosslinker Position A`, `Crosslinker Position B`) %>%
  slice_max(`Max. XlinkX Score`, with_ties = FALSE) %>%
  select(
    "Crosslink Type",
    "Sequence",
    "m/z [Da]",
    "Charge",
    "MH+ [Da]",
    "Sequence A",
    "Crosslinker Position A",
    "Modifications A",
    "Sequence B",
    "Crosslinker Position B",
    "Modifications B",
    "Accessions",
    "Protein Accession A",
    "Protein Accession B"
  ) %>%
  full_join(
    spread_XlinkX_Score,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  full_join(
    spread_Total_CSMs,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  full_join(
    spread_Total_Precursor_Abundances,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  #Delete '[' and ']' from sequence string:
  mutate(across(.cols = c("Sequence A", "Sequence B"),  ~ gsub("\\[|\\]", "", .)))

name_out <- paste0('Spread_by_All_Conditions_All_Rep_',
                   Sys.Date(),
                   "_",
                   Subjected_File,
                   ".csv")
rute_out_c <- paste0(
  getwd(),
  .Platform$file.sep,
  Output_file,
  .Platform$file.sep,
  "Crosslink_Peptides",
  .Platform$file.sep,
  name_out
)
write.csv(spread_ModifiedTable, file = rute_out_c, row.names = FALSE)

#Filter to get only all frequent cross links:
Nrep <- Minimun_replicas - 1
FrequentXL <- ModifiedTable %>%
  group_by(`Sequence`, `Crosslinker Position A`, `Crosslinker Position B`) %>%
  group_by(if (Replicas_in == 'Treatment')
    `Experimental group`
    else
      TRUE, .add = TRUE) %>%
  filter(n() > Nrep) %>%
  ungroup()

FrequentXL_Exp <- FrequentXL
name_out <- paste0(
  '.Frecuent_Crosslinks_(In at least ',
  Minimun_replicas,
  ' rep)_',
  Sys.Date(),
  "_",
  Subjected_File,
  ".csv"
)
rute_out_b <- paste0(
  getwd(),
  .Platform$file.sep,
  Output_file,
  .Platform$file.sep,
  "Crosslink_Peptides",
  .Platform$file.sep,
  name_out
)
write.csv(FrequentXL_Exp, file = rute_out_b, row.names = FALSE)

# Spreading of XLinks data:
# 'Max. XlinkX Score':
spread_XlinkX_Score <- FrequentXL %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Max. XlinkX Score`
  ) %>%
  spread(key = `Spectrum File`, value = `Max. XlinkX Score`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )),  ~ paste0("XlinkX_Score_", .))
# 'Total CSMs':
spread_Total_CSMs <- FrequentXL %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Total CSMs`
  ) %>%
  spread(key = `Spectrum File`, value = `Total CSMs`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )),  ~ paste0("Total_CSM_", .))
# 'Total Precursor Abundances':
spread_Total_Precursor_Abundances <- FrequentXL %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Total Precursor Abundances`
  ) %>%
  spread(key = `Spectrum File`, value = `Total Precursor Abundances`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )),  ~ paste0("Total_Precursor_Abundances_", .))

# Join the spread data frames on 'Sequence':
spread_ModifiedTable_min <- FrequentXL %>%
  group_by(`Sequence`, `Crosslinker Position A`, `Crosslinker Position B`) %>%
  slice_max(`Max. XlinkX Score`, with_ties = FALSE) %>%
  select(
    "Crosslink Type",
    "Sequence",
    "m/z [Da]",
    "Charge",
    "MH+ [Da]",
    "Sequence A",
    "Crosslinker Position A",
    "Modifications A",
    "Sequence B",
    "Crosslinker Position B",
    "Modifications B",
    "Accessions",
    "Protein Accession A",
    "Protein Accession B",
    "RT [min]"
  ) %>% #Added RT min
  full_join(
    spread_XlinkX_Score,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  full_join(
    spread_Total_CSMs,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  full_join(
    spread_Total_Precursor_Abundances,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  #Delete '[' and ']' from sequence string:
  mutate(across(.cols = c("Sequence A", "Sequence B"),  ~ gsub("\\[|\\]", "", .)))

name_out <- paste0(
  'Spread_by_Frecuent_Crosslinks_(In at least ',
  Minimun_replicas,
  ' rep)_',
  Sys.Date(),
  "_",
  Subjected_File,
  ".csv"
)
rute_out_c <- paste0(
  getwd(),
  .Platform$file.sep,
  Output_file,
  .Platform$file.sep,
  "Crosslink_Peptides",
  .Platform$file.sep,
  name_out
)
write.csv(spread_ModifiedTable_min,
          file = rute_out_c,
          row.names = FALSE)

#Filter to get only one frequent cross links in the experimental group:
FrequentXL <- ModifiedTable %>%
  group_by(`Sequence`, `Crosslinker Position A`, `Crosslinker Position B`) %>%
  group_by(if (Replicas_in == 'Treatment')
    `Experimental group`
    else
      TRUE, .add = TRUE) %>%
  filter(n() > Nrep) %>%
  slice_max(`Max. XlinkX Score`, with_ties = FALSE) %>%
  ungroup() %>%
  #Delete '[' and ']' from sequence string:
  mutate(across(.cols = c("Sequence A", "Sequence B"),  ~ gsub("\\[|\\]", "", .)))

#Get informative strands and set format of Xinet:
Xinet_Info <- data.frame(matrix(NA, nrow = nrow(FrequentXL), ncol = 14))

colnames(Xinet_Info) <- c(
  'Sequence_XL',
  'PepSeq1',
  'PepPos1',
  'PepSeq2',
  'PepPos2',
  'Protein1',
  'LinkPos1',
  'Protein2',
  'LinkPos2',
  'Score',
  'Protein1fasta',
  'Protein2fasta',
  'Spectrum File',
  'Experimental group'
)
Xinet_Info <- Xinet_Info %>%
  mutate(
    Sequence_XL = FrequentXL$`Sequence`,
    PepSeq1 = FrequentXL$`Sequence A`,
    LinkPos1 = as.numeric(FrequentXL$`Crosslinker Position A`),
    PepSeq2 = FrequentXL$`Sequence B`,
    LinkPos2 = as.numeric(FrequentXL$`Crosslinker Position B`),
    `Score` = as.numeric(FrequentXL$`Max. XlinkX Score`),
    `Spectrum File` = FrequentXL$`Spectrum File`,
    `Experimental group` = FrequentXL$`Experimental group`
  )

#- Database FASTA creation:
fasta_loc <- paste0(getwd(),
                    .Platform$file.sep,
                    'Input_Database',
                    .Platform$file.sep,
                    Database)
Seqnces <- read.fasta(
  file = fasta_loc,
  seqtype = "AA",
  as.string = TRUE,
  set.attributes = TRUE
)
SeqDatabse <- do.call(rbind, Seqnces)
SeqDatabse <- as.data.frame(SeqDatabse)
SeqDatabse <- SeqDatabse %>%
  mutate(Complete_Name = names(Seqnces)) %>%
  mutate(Accesion = gsub(".*\\|(.*?)\\|.*", "\\1", Complete_Name)) %>%
  select(Complete_Name, Accesion, Sequence = V1) %>%
  group_by(Complete_Name, Accesion, Sequence) %>%
  distinct() %>%
  ungroup() %>%
  arrange(Accesion)

#Blast search in database, completening of Xinet information
Xinet_All.XLs <- Xinet_Info %>%
  group_by(`Sequence_XL`, `LinkPos1`, `LinkPos2`) %>%
  slice_max(`Score`, with_ties = FALSE) %>%
  ungroup() %>%
  select(-c(`Spectrum File`, `Experimental group`))

SinglePeptStrands <- bind_rows(
  Xinet_All.XLs %>% select(
    Peptide = PepSeq1,
    LinkPos = LinkPos1,
    PepPos = PepPos1,
    Protein = Protein1,
    Protein_fasta = Protein1fasta
  ),
  Xinet_All.XLs %>% select(
    Peptide = PepSeq2,
    LinkPos = LinkPos2,
    PepPos = PepPos2,
    Protein = Protein2,
    Protein_fasta = Protein2fasta
  )
)

SinglePeptStrands <- SinglePeptStrands %>%
  mutate(
    Peptide = as.character(Peptide),
    LinkPos = as.numeric(LinkPos),
    PepPos = as.character(PepPos)
  ) %>%
  distinct(Peptide, LinkPos, .keep_all = TRUE)

#- Blasting:
print("Blasting of peptides is starting: Please do not stop the code.")
control <- NA
Cond_proc <- SinglePeptStrands
total_peptides <- nrow(Cond_proc)

for (H in 1:total_peptides) {
  pattern <- as.character(Cond_proc$Peptide[H])
  mtch_rslt <- vmatchPattern(pattern, as.character(SeqDatabse$Sequence), max.mismatch = 0)
  mtch_as.df <- as.data.frame(mtch_rslt)
  
  Cond_proc$Protein[H] <- SeqDatabse$Complete_Name[mtch_as.df[1, 1]]
  Cond_proc$Protein_fasta[H] <- SeqDatabse$Accesion[mtch_as.df[1, 1]]
  Cond_proc$PepPos[H] <- as.character(mtch_as.df$start[1])
  
  if (nrow(mtch_as.df) > 1) {
    for (G in 2:nrow(mtch_as.df)) {
      Variab <- paste(mtch_as.df$start[1:(G - 1)], SeqDatabse$Complete_Name[mtch_as.df$group[1:(G -
                                                                                                  1)]])
      if (!(paste(mtch_as.df$start[G], SeqDatabse$Complete_Name[mtch_as.df$group[G]]) %in% Variab)) {
        control <- 'Yes'
        vec_pep <- as.data.frame(Cond_proc[H, ])
        vec_pep$Protein[1] <- paste0(vec_pep$Protein[1], '; ', SeqDatabse$Complete_Name[mtch_as.df[G, 1]])
        vec_pep$Protein_fasta[1] <- paste0(vec_pep$Protein_fasta[1], '; ', SeqDatabse$Accesion[mtch_as.df[G, 1]])
        vec_pep$PepPos[1] <- paste0(vec_pep$PepPos[1], '; ', mtch_as.df$start[G])
        Cond_proc[H, ] <- vec_pep
      }
    }
  }
  
  # Progress every 50 peptides
  if (H %% 50 == 0) {
    percent <- round((H / total_peptides) * 100, 2)
    message(paste("Blasted", H, "peptides —", percent, "% completed"))
  }
}

print(paste0("End of blasting, ", total_peptides, " processed"))

#Warning:
if (is.na(control)) {
  print("Extra protein combinations were not created")
} else {
  print("EXTRA PROTEIN COMBINATIONS were introduced in some Cross-Links !!!")
  SinglePeptStrands <- NA
}

#Match in table with unique CrossLinks from the list of peptides (Cond_proc)
for (i in 1:nrow(Cond_proc)) {
  Peptide <- Cond_proc$Peptide[i]
  LinkPos <- Cond_proc$LinkPos[i]
  
  proteinID <- ifelse(
    is.na(Cond_proc$Protein[i]) ||
      Cond_proc$Protein[i] == "",
    "Error!",
    Cond_proc$Protein[i]
  )
  fastaID <- ifelse(
    is.na(Cond_proc$Protein_fasta[i]) ||
      Cond_proc$Protein_fasta[i] == "",
    "Error!",
    Cond_proc$Protein_fasta[i]
  )
  peptidePosition <- ifelse(is.na(Cond_proc$PepPos[i]) ||
                              Cond_proc$PepPos[i] == "",
                            "Error!",
                            Cond_proc$PepPos[i])
  
  # Matches to first peptide
  match_1 <- which(Xinet_All.XLs$PepSeq1 == Peptide &
                     Xinet_All.XLs$LinkPos1 == LinkPos)
  if (length(match_1) > 0) {
    Xinet_All.XLs$Protein1[match_1] <- proteinID
    Xinet_All.XLs$Protein1fasta[match_1] <- fastaID
    Xinet_All.XLs$PepPos1[match_1] <- peptidePosition
  }
  
  # Matches to second peptide
  match_2 <- which(Xinet_All.XLs$PepSeq2 == Peptide &
                     Xinet_All.XLs$LinkPos2 == LinkPos)
  if (length(match_2) > 0) {
    Xinet_All.XLs$Protein2[match_2] <- proteinID
    Xinet_All.XLs$Protein2fasta[match_2] <- fastaID
    Xinet_All.XLs$PepPos2[match_2] <- peptidePosition
  }
}





#Match in table with all CrossLinks from table with unique CrossLinks (Xinet_All.XLs)
idx <- match(Xinet_Info$Sequence_XL, Xinet_All.XLs$Sequence_XL)
# For all matched indices, assign values from Xinet_All.XLs
cols_to_update <- c("PepPos1",
                    "PepPos2",
                    "Protein1",
                    "Protein2",
                    "Protein1fasta",
                    "Protein2fasta")

for (col in cols_to_update) {
  Xinet_Info[[col]] <- Xinet_All.XLs[[col]][idx]
}


#Reasigning of protein names and IDs in crosslink tables:
#CSMs Table - frequents in at least in some replicas
Xinet_Info <- Xinet_Info %>%
  rowwise() %>%
  mutate(
    swap.ProtA = `Protein1`,
    swap.SeqA  = `PepSeq1`,
    swap.PosA  = `PepPos1`,
    swap.LinkA = `LinkPos1`,
    swap.FastaA = `Protein1fasta`,
    swap   = replace_na(`Protein1` > `Protein2`, FALSE)
  ) %>%
  #Do a swap on columns if true:
  mutate(
    # Accesions
    `Protein1` = if_else(swap, `Protein2`, `Protein1`),
    `Protein2` = if_else(swap, swap.ProtA, `Protein2`),
    # Swap of Sequences
    `PepSeq1` = if_else(swap, `PepSeq2`, `PepSeq1`),
    `PepSeq2` = if_else(swap, swap.SeqA, `PepSeq2`),
    # Xl position
    `PepPos1` = if_else(swap, `PepPos2`, `PepPos1`),
    `PepPos2` = if_else(swap, swap.PosA, `PepPos2`),
    # Protein Positions
    `LinkPos1` = if_else(swap, `LinkPos2`, `LinkPos1`),
    `LinkPos2` = if_else(swap, swap.LinkA, `LinkPos2`),
    # Protein FASTA
    `Protein1fasta` = if_else(swap, `Protein2fasta`, `Protein1fasta`),
    `Protein2fasta` = if_else(swap, swap.FastaA, `Protein2fasta`),
    # rebuild Sequence column
    Sequence_XL = paste(
      gsub("\\[|\\]", "", `PepSeq1`),
      gsub("\\[|\\]", "", `PepSeq2`),
      sep = "-"
    )
  ) %>%
  ungroup() %>%
  select(-starts_with("swap"))

# Generate dataframes for analysis:
Exp_Cond <- Unique_Conditions
Xinet_File <- List()
for (D in 1:length(Exp_Cond)) {
  Condition_procesed <- Xinet_Info %>%
    filter(`Experimental group` == Exp_Cond[D])
  Xinet_File[[paste0("XLink_Treatment_", Unique_Conditions[D])]] <- list(Condition_procesed)
}

#Reasigning of protein names and IDs in crosslink tables:
#CSMs Table - frequents in at least in some replicas
reasign_ModifiedTable_min <- FrequentXL_Exp %>%
  left_join(
    Xinet_All.XLs,
    by = c(
      "Sequence" = "Sequence_XL",
      "Crosslinker Position A" = "LinkPos1",
      "Crosslinker Position B" = "LinkPos2"
    )
  ) %>%
  mutate(
    `Protein Accession A` = Protein1,
    `Protein Accession B` = Protein2,
    `Leading Protein Position A` = PepPos1,
    `Leading Protein Position B` = PepPos2
  ) %>%
  rowwise() %>%
  mutate(
    swap.ProtA = `Protein Accession A`,
    swap.SeqA  = `Sequence A`,
    swap.ModA  = `Modifications A`,
    swap.PosA  = `Crosslinker Position A`,
    swap.LeadA = `Leading Protein Position A`,
    swap   = replace_na(`Protein Accession A` > `Protein Accession B`, FALSE)
  ) %>%
  #Do a swap on columns if true:
  mutate(
    # Accesions
    `Protein Accession A` = if_else(swap, `Protein Accession B`, `Protein Accession A`),
    `Protein Accession B` = if_else(swap, swap.ProtA, `Protein Accession B`),
    # Swap of Sequences
    `Sequence A` = if_else(swap, `Sequence B`, `Sequence A`),
    `Sequence B` = if_else(swap, swap.SeqA, `Sequence B`),
    # Midicifations
    `Modifications A` = if_else(swap, `Modifications B`, `Modifications A`),
    `Modifications B` = if_else(swap, swap.ModA, `Modifications B`),
    # Xl position
    `Crosslinker Position A` = if_else(swap, `Crosslinker Position B`, `Crosslinker Position A`),
    `Crosslinker Position B` = if_else(swap, swap.PosA, `Crosslinker Position B`),
    # Protein Positions
    `Leading Protein Position A` = if_else(swap, `Leading Protein Position B`, `Leading Protein Position A`),
    `Leading Protein Position B` = if_else(swap, swap.LeadA, `Leading Protein Position B`),
    # rebuild Sequence column
    Sequence = paste(
      gsub("\\[|\\]", "", `Sequence A`),
      gsub("\\[|\\]", "", `Sequence B`),
      sep = "-"
    )
  ) %>%
  ungroup() %>%
  select(-starts_with("swap")) %>%
  mutate(Accessions = paste(`Protein Accession A`, `Protein Accession B`, sep = " - ")) %>%
  relocate(
    `Protein Accession A`,
    `Leading Protein Position A`,
    `Protein Accession B`,
    `Leading Protein Position B`,
    Accessions,
    .after = `Protein Accession B`
  ) %>%
  select(-(PepSeq1:ncol(.)))

name_out <- paste0(
  '.Freq&Reasign_XLs_(In at least ',
  Minimun_replicas,
  ' rep)_',
  Sys.Date(),
  "_",
  Subjected_File,
  ".csv"
)
rute_out_b <- paste0(
  getwd(),
  .Platform$file.sep,
  Output_file,
  .Platform$file.sep,
  "Crosslink_Peptides_Reasigned",
  .Platform$file.sep,
  name_out
)
write.csv(reasign_ModifiedTable_min,
          file = rute_out_b,
          row.names = FALSE)

# Spreading of XLinks data from re-asigned tables:
# 'Max. XlinkX Score':
spread_XlinkX_Score <- reasign_ModifiedTable_min %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Max. XlinkX Score`
  ) %>%
  spread(key = `Spectrum File`, value = `Max. XlinkX Score`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )),  ~ paste0("XlinkX_Score_", .))
# 'Total CSMs':
spread_Total_CSMs <- reasign_ModifiedTable_min %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Total CSMs`
  ) %>%
  spread(key = `Spectrum File`, value = `Total CSMs`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )),  ~ paste0("Total_CSM_", .))
# 'Total Precursor Abundances':
spread_Total_Precursor_Abundances <- reasign_ModifiedTable_min %>%
  select(
    Sequence,
    `Crosslinker Position A`,
    `Crosslinker Position B`,
    `Spectrum File`,
    `Total Precursor Abundances`
  ) %>%
  spread(key = `Spectrum File`, value = `Total Precursor Abundances`) %>%
  rename_at(vars(-c(
    Sequence, `Crosslinker Position A`, `Crosslinker Position B`
  )),  ~ paste0("Total_Precursor_Abundances_", .))

# Join the spread data frames on 'Sequence':
spread_reasign_ModifiedTable_min <- reasign_ModifiedTable_min %>%
  group_by(`Sequence`, `Crosslinker Position A`, `Crosslinker Position B`) %>%
  slice_max(`Max. XlinkX Score`, with_ties = FALSE) %>%
  select(
    "Crosslink Type",
    "Sequence",
    "m/z [Da]",
    "Charge",
    "MH+ [Da]",
    "Sequence A",
    "Crosslinker Position A",
    "Modifications A",
    "Sequence B",
    "Crosslinker Position B",
    "Modifications B",
    "Accessions",
    "Protein Accession A",
    "Protein Accession B",
    "RT [min]"
  ) %>% #Added RT min
  full_join(
    spread_XlinkX_Score,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  full_join(
    spread_Total_CSMs,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  full_join(
    spread_Total_Precursor_Abundances,
    by = c("Sequence", "Crosslinker Position A", "Crosslinker Position B")
  ) %>%
  #Delete '[' and ']' from sequence string:
  mutate(across(.cols = c("Sequence A", "Sequence B"),  ~ gsub("\\[|\\]", "", .)))

name_out <- paste0(
  'Spread&Reasign_by_Freq_XLs_(In at least ',
  Minimun_replicas,
  ' rep)_',
  Sys.Date(),
  "_",
  Subjected_File,
  ".csv"
)
rute_out_c <- paste0(
  getwd(),
  .Platform$file.sep,
  Output_file,
  .Platform$file.sep,
  "Crosslink_Peptides_Reasigned",
  .Platform$file.sep,
  name_out
)
write.csv(spread_reasign_ModifiedTable_min,
          file = rute_out_c,
          row.names = FALSE)

#Look for crosslink present in all conditions:
All_Conditions <- data.frame()
for (Z in 1:length(Xinet_File)) {
  Xin <- as.data.frame(Xinet_File[Z])
  Xin <- Xin %>% mutate(across(everything(), as.character))
  All_Conditions <- bind_rows(All_Conditions, Xin)
  Xin <- Xin[, match('PepSeq1', names(Xin)):match('Protein2fasta', names(Xin))]
  Xin <- Xin %>%
    distinct()
  Name_Otpt <- paste0(Name_Output,
                      '-Treatment_',
                      Unique_Conditions[Z],
                      "_",
                      Sys.Date(),
                      ".csv")
  rute_out <- paste0(getwd(),
                     .Platform$file.sep,
                     Output_file,
                     .Platform$file.sep,
                     Name_Otpt)
  write.csv(Xin, file = rute_out, row.names = FALSE)
}

# Perform comparation by exclusive & common:
if (length(Exp_Cond) == 2) {
  Final <- All_Conditions %>%
    group_by(PepSeq1, LinkPos1, PepSeq2, LinkPos2) %>%
    mutate(group = ifelse(n() > 1, "Common", "Exclusive")) %>%
    ungroup() %>%
    select(
      PepSeq1,
      PepPos1,
      PepSeq2,
      PepPos2,
      Protein1,
      LinkPos1,
      Protein2,
      LinkPos2,
      Score,
      Protein1fasta,
      Protein2fasta,
      Experimental.group,
      group
    )
  
  #Save exclusive XLs in conditions:
  for (Z1 in 1:2) {
    Opt_Xl <- Final %>%
      filter(group == "Exclusive") %>%
      filter(Experimental.group == Exp_Cond[Z1]) %>%
      select(names(Final[, 1:11])) %>%
      distinct()
    Name_Otpt <- paste0(
      Name_Output,
      "_",
      "Exclusive_in.Treatment_",
      Unique_Conditions[Z1],
      "_",
      Sys.Date(),
      ".csv"
    )
    rute_out_c <- paste0(
      getwd(),
      .Platform$file.sep,
      Output_file,
      .Platform$file.sep,
      "Xinet_by_Exlusive&Common",
      .Platform$file.sep,
      Name_Otpt
    )
    write.csv(Opt_Xl, file = rute_out_c, row.names = FALSE)
  }
  
  #Save common XLs between conditions:
  Opt_Xl_com <- Final %>%
    filter(group == "Common") %>%
    select(names(Final[, 1:11])) %>%
    group_by(PepSeq1, LinkPos1, PepSeq2, LinkPos2) %>%
    slice_max(order_by = Score, with_ties = FALSE) %>%
    ungroup()
  Name_Otpt <- paste0(Name_Output, "_Common_in.Treatments_", Sys.Date(), ".csv")
  rute_out_c <- paste0(
    getwd(),
    .Platform$file.sep,
    Output_file,
    .Platform$file.sep,
    "Xinet_by_Exlusive&Common",
    .Platform$file.sep,
    Name_Otpt
  )
  write.csv(
    Opt_Xl_com,
    file = rute_out_c,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  #All XLs in the conditions:
  All <- Final %>%
    group_by(PepSeq1, LinkPos1, PepSeq2, LinkPos2) %>%
    slice_max(order_by = Score, with_ties = FALSE) %>%
    ungroup() %>%
    select(names(Final[, 1:11])) %>%
    distinct()
  Name_Otpt <- paste0(Name_Output, "_", '_all.Treatments', "_", Sys.Date(), ".csv")
  rute_out_c <- paste0(
    getwd(),
    .Platform$file.sep,
    Output_file,
    .Platform$file.sep,
    "Xinet_by_Exlusive&Common",
    .Platform$file.sep,
    Name_Otpt
  )
  write.csv(All, file = rute_out_c, row.names = FALSE)
}

#In case that exist multiple conditions or treatments:
if (length(Exp_Cond) > 2) {
  #Only common XLs between conditions:
  Opt_Xl_com <- All_Conditions %>%
    group_by(PepSeq1, LinkPos1, PepSeq2, LinkPos2) %>%
    filter(n() == length(Exp_Cond)) %>%
    slice_max(order_by = Score, with_ties = FALSE) %>%
    ungroup() %>%
    select(
      PepSeq1,
      PepPos1,
      PepSeq2,
      PepPos2,
      Protein1,
      LinkPos1,
      Protein2,
      LinkPos2,
      Score,
      Protein1fasta,
      Protein2fasta
    ) %>%
    distinct()
  Name_Otpt <- paste0(Name_Output, "_Common_in.Treatments_", Sys.Date(), ".csv")
  rute_out_c <- paste0(
    getwd(),
    .Platform$file.sep,
    Output_file,
    .Platform$file.sep,
    "Xinet_by_Exlusive&Common",
    .Platform$file.sep,
    Name_Otpt
  )
  write.csv(Opt_Xl_com, file = rute_out_c, row.names = FALSE)
  
  #All XLs in the conditions:
  All <- All_Conditions %>%
    group_by(PepSeq1, LinkPos1, PepSeq2, LinkPos2) %>%
    slice_max(order_by = Score, with_ties = FALSE) %>%
    ungroup() %>%
    select(
      PepSeq1,
      PepPos1,
      PepSeq2,
      PepPos2,
      Protein1,
      LinkPos1,
      Protein2,
      LinkPos2,
      Score,
      Protein1fasta,
      Protein2fasta
    ) %>%
    distinct()
  Name_Otpt <- paste0(Name_Output, "_", '_all.Treatments', "_", Sys.Date(), ".csv")
  rute_out_c <- paste0(
    getwd(),
    .Platform$file.sep,
    Output_file,
    .Platform$file.sep,
    "Xinet_by_Exlusive&Common",
    .Platform$file.sep,
    Name_Otpt
  )
  write.csv(All, file = rute_out_c, row.names = FALSE)
}
#End
print("The analysis has finished :D")
closeAllConnections()
rm(list = ls())
gc()
