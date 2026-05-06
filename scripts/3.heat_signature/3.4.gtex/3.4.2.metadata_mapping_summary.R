library(readr)
library(dplyr)
library(stringr)

sub_expr <- read_csv("output/3.heat_signature/5.gtex/data/gene_tpm_2022-06-06_v10_adipose_subcutaneous_ssGSEA.csv")
vis_expr <- read_csv("output/3.heat_signature/5.gtex/data/gene_tpm_2022-06-06_v10_adipose_visceral_omentum_ssGSEA.csv")

meta_v8 <- read_tsv("local_data/gtex/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt")

# Adjust this depending on your column names.
# GTEx sample IDs usually look like GTEX-XXXX-...
# Donor IDs are often the first two fields, e.g. GTEX-XXXX.
get_donor <- function(x) str_extract(x, "^GTEX-[A-Z0-9]+")

sub_donors <- unique(sub_expr$subject)
vis_donors <- unique(vis_expr$subject)

meta_donor_col <- "SUBJID"  # change if needed
meta_donors <- unique(meta_v8[[meta_donor_col]])
meta_donors <- str_remove(meta_donors, "^GTEX-")

check_overlap <- function(expr_donors, label) {
  tibble(
    dataset = label,
    n_expression_donors = length(expr_donors),
    n_expression_donors_in_v8_meta = sum(expr_donors %in% meta_donors),
    n_expression_donors_missing_v8_meta = sum(!expr_donors %in% meta_donors)
  )
}

bind_rows(
  check_overlap(sub_donors, "V10 subcutaneous"),
  check_overlap(vis_donors, "V10 visceral_omentum")
)

setdiff(sub_donors, meta_donors)
setdiff(vis_donors, meta_donors)

write_csv(
  bind_rows(
    check_overlap(sub_donors, "V10 subcutaneous"),
    check_overlap(vis_donors, "V10 visceral_omentum")
  ),
  "output/3.heat_signature/5.gtex/data/gtex_v10_expression_v8_metadata_overlap_check.csv"
)
