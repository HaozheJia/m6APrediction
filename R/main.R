# Functions reused from Practical 5. Add documentation in your new package.
#' Encode DNA 5‑mer sequence into factor matrix
#'
#' Convert DNA character sequences (A/T/C/G) into factor data frame for model input.
#' Internal helper function, not exported for end‑users.
#'
#' @param dna_strings Character vector of DNA 5‑mer sequences
#' @return A data frame, each column is factor for one nucleotide position



dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A sites for multiple records
#'
#' Use pre‑trained random forest model to predict m6A modification probability
#' and modification status for many input feature rows.
#' Required columns in feature_df: gc_content, RNA_type, RNA_region, exon_length,
#' distance_to_junction, evolutionary_conservation, DNA_5mer.
#'
#' @param ml_fit Pre‑trained random forest model object
#' @param feature_df Input feature data frame
#' @param positive_threshold Probability cutoff for "Positive" call; default = 0.5
#' @return Input data frame with two appended columns: predicted_m6A_prob, predicted_m6A_status
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' \dontrun{
#' example_csv <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' example_data <- read.csv(example_csv)
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' result <- prediction_multiple(rf_model, example_data)
#' head(result)
#' }



prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A site for one single sample
#'
#' Single‑sample wrapper for prediction_multiple. Takes individual feature values
#' and returns a named vector with predicted probability and modification status.
#'
#' @param ml_fit Pre‑trained random forest model object
#' @param gc_content Numeric GC‑content value
#' @param RNA_type Character string RNA type: mRNA / lincRNA / lncRNA / pseudogene
#' @param RNA_region Character string: CDS / intron / 3'UTR / 5'UTR
#' @param exon_length Numeric exon length
#' @param distance_to_junction Numeric distance to splice junction
#' @param evolutionary_conservation Numeric conservation score
#' @param DNA_5mer Character 5‑mer DNA sequence
#' @param positive_threshold Probability cutoff, default = 0.5
#' @return Named vector: predicted_m6A_prob, predicted_m6A_status
#' @export
#' @examples
#' \dontrun{
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(
#'   ml_fit = rf_model,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA"
#' )
#' }



prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
