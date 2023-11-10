# Read gene expression data from a tab-separated file
genesLUAD = read.table("genesLUAD.tsv", header = TRUE, sep = "\t")

# Filter data based on specific gene names: CAPN2, GAL3ST2, GPR27
data = genesLUAD[(genesLUAD[, 'gene_name'] == 'CAPN2') | (genesLUAD[, 'gene_name'] == 'GAL3ST2') | (genesLUAD[, 'gene_name'] == 'GPR27'), ]

# Display dimensions of the filtered data
dim(data)
# Output: 3 204

# Extract relevant columns for each gene
"patients" = colnames(data[4:201])
"CAPN2" = t(data[data['gene_name'] == 'CAPN2', 4:201])
"GAL3ST2" = t(data[data['gene_name'] == 'GAL3ST2', 4:201])
"GPR27" = t(data[data['gene_name'] == 'GPR27', 4:201])

# Create a data frame with extracted data
genes = data.frame(patients, CAPN2, GAL3ST2, GPR27)

# Read clinical and subtype data from tab-separated files
clinLUAD = read.table("clinLUAD.tsv", header = TRUE, sep = "\t")
subtypeLUAD = read.table("subtypeLUAD.tsv", header = TRUE, sep = "\t")

# Create a 'patient' column in clinLUAD to match with the 'patient' column in subtypeLUAD
clinLUAD$patient <- clinLUAD$bcr_patient_barcode

# Merge clinical and subtype data based on the 'patient' column
surv_data <- merge(x = subtypeLUAD, y = clinLUAD, by = "patient", all.x = TRUE)[, c("patient", "vital_status")]

# Extract the first 12 digits from the 'patients' column and replace '.' with '-'
genes$patient <- substr(genes$patients, 1, 12)
genes$patient <- gsub("\\.", "-", genes$patient)

# Remove duplicates from 'surv_data' and 'genes' based on the 'patient' column
length(unique(surv_data$patient))
surv_data <- surv_data[!duplicated(surv_data), ]
genes <- genes[!duplicated(genes), ]

# Load the dplyr library for data manipulation
library(dplyr)

# Perform a left join on multiple columns ('patient') between 'surv_data' and 'genes'
expr_and_surv_df <- surv_data %>% left_join(genes, by = c('patient'))

# Remove rows with missing values
expr_and_surv_df <- expr_and_surv_df[complete.cases(expr_and_surv_df), ]

# Remove duplicate rows based on the 'patient' column
expr_and_surv_df <- expr_and_surv_df[!duplicated(expr_and_surv_df$patient), ]

# Rename columns for better readability
colnames(expr_and_surv_df) = c("id_pat", "vital_status", "patients", "CAPN2", "GAL3ST2", "GPR27")

# Convert 'vital_status' to binary: "Alive" to 1, other values to 0
expr_and_surv_df$vital_status <- ifelse(expr_and_surv_df$vital_status == "Alive", 1, 0)

# Set row names to 'id_pat' column and remove unnecessary columns
rownames(expr_and_surv_df) <- expr_and_surv_df$id_pat
expr_and_surv_df <- subset(expr_and_surv_df, select = -c(id_pat, patients))

# Fit a logistic regression model with three genes: CAPN2, GAL3ST2, and GPR27
model = glm(vital_status ~ CAPN2 + GAL3ST2 + GPR27, data = expr_and_surv_df, family = binomial)

# Display the summary of the logistic regression model
summary(model)

# Fit a simplified model with only the CAPN2 gene
model2 = glm(vital_status ~ CAPN2, data = expr_and_surv_df, family = binomial)

# Display the summary of the simplified model
summary(model2)
