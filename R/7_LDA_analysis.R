##Project: Ascaris - Pig Microbiome
##Aim: Liner discriminant analysis (LDA) effect size
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")
##Data
##Question specific
PS.pig<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Rds") ## Data just merged pigs not normalized for alpha diversity plots  
PS.pig.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Norm.Rds") ## Data just merged pigs normalized for beta diversity plots 
PS.Asc<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Rds") ## Data all Ascaris not normalized for alpha diversity plots 
PS.Asc.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Norm.Rds") ## Data all Ascaris normalized for beta diversity plots 
PS.PA<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Rds") ## Data merged pigs and Ascaris (not SH) not normalized for alpha diversity plots 
PS.PA.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Norm.Rds") ## Data merged pigs and Ascaris (not SH) normalized for beta diversity plots 
PS.pig.diff<- readRDS("Data/PS.pig.diff.rds")

##Load libraries 

###Using function from "microbiomeMarker" without the package 
run_lefse <- function(ps,
                      group,
                      subgroup = NULL,
                      taxa_rank = "all",
                      transform = c("identity", "log10", "log10p"),
                      norm = "CPM",
                      norm_para = list(),
                      kw_cutoff = 0.05,
                      lda_cutoff = 2,
                      bootstrap_n = 30,
                      bootstrap_fraction = 2 / 3,
                      wilcoxon_cutoff = 0.05,
                      multigrp_strat = FALSE,
                      strict = c("0", "1", "2"),
                      sample_min = 10,
                      only_same_subgrp = FALSE,
                      curv = FALSE) {
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be phyloseq object", call. = FALSE)
  }
  
  if (!check_rank_names(ps)) {
    stop(
      "ranks of `ps` must be one of ",
      paste(available_ranks, collapse = ", ")
    )
  }
  
  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  strict <- match.arg(strict, c("0", "1", "2"))
  strict <- as.numeric(strict)
  
  # import input from the original lefse python script or galaxy,
  # will be dropped in the next release version
  summarized <- check_tax_summarize(ps)
  if (summarized && norm != "CPM") {
    stop(
      "`norm` must be a 'CPM' or 'none' while `ps` has been summarized",
      call. = FALSE
    )
  }
  
  sample_meta <- sample_data(PS.pig.diff)
  grp_info <- lefse_format_grp(sample_meta, "InfectionStatus")
  grp <- grp_info$group
  subgrp <- grp_info$subgroup
  grp_hie <- grp_info$group_hie
  
  otus <- abundances(PS.pig.diff)
  # transform it for test
  otus_test <- as.data.frame(t(otus), stringsAsFactors = FALSE)
  feature <- tax_table(PS.pig.diff)@.Data[, 6]
  #names(otus_test) <- feature
  
  # kw rank sum test among classes
  kw_p <- purrr::map_dbl(otus_test, ~ kruskal.test(.x, PS.pig.diff@sam_data$InfectionStatus)$p.value)
  
  # remove the taxa, while pvalue is na
  na_ind <- is.na(kw_p)
  if (sum(na_ind) >= 1) {
    otus_test <- otus_test[!na_ind]
    kw_p <- kw_p[!na_ind]
  }
  
  sig_ind <- kw_p <= 0.01
  sig_otus <- otus_test[, sig_ind]
  
  # wilcox test is preformed for each class, if there is no subclass
  features_nms <- names(sig_otus)
  wilcoxon_p <- purrr::map2_lgl(
    sig_otus, features_nms,
    ~ test_rep_wilcoxon(grp_hie,
      .x, .y,
      wilcoxon_cutoff = wilcoxon_cutoff,
      multicls_strat = multigrp_strat,
      strict = 0,
      sample_min = sample_min,
      curv = curv
    )
  )
  sig_otus <- sig_otus[, wilcoxon_p]
  
  # mean abundance in each group
  otus_enriched_group <- get_feature_enrich_group(grp, sig_otus)
  
  # bootsrap iteration of lda
  ldas <- bootstap_lda(
    sig_otus,
    boot_n = 100,
    class = grp,
    sample_fract = 10
  )
  
  lefse_res <- data.frame(
    feature = names(sig_otus),
    enrich_group = otus_enriched_group$group,
    ef_lda = ldas,
    pvalue = kw_p[sig_ind][wilcoxon_p],
    stringsAsFactors = FALSE
  )
  
  lefse_sig <- filter(lefse_res, .data$ef_lda >= lda_cutoff) %>%
    arrange(.data$enrich_group, desc(.data$ef_lda))
  
  lefse_out <- return_marker(lefse_sig, lefse_res)
  lefse_out$padj <- lefse_out$pvalue
  row.names(lefse_out) <- paste0("marker", seq_len(nrow(lefse_out)))
  
  tax <- matrix(feature) %>%
    tax_table()
  row.names(tax) <- row.names(otus)
  
  mm <- microbiomeMarker(
    marker_table = lefse_out,
    norm_method = get_norm_method(norm),
    diff_method = "lefse",
    otu_table = otu_table(otus, taxa_are_rows = TRUE), # normalized
    # new var norm_factor (if it is calculated in normalize)
    sam_data = sample_data(ps_normed),
    tax_table = tax
  )
  
  mm
}

check_tax_summarize <- function(ps) {
  taxa <- row.names(otu_table(ps))
  # whether taxa is separated by `|`,
  # may be required to add extra separate strings in the future
  has_separate <- any(grepl("[|]", taxa))
  
  has_separate
}

check_rank_names <- function(ps) {
  summarized <- check_tax_summarize(ps)
  if (summarized) {
    return(TRUE)
  }
  
  ps_ranks <- rank_names(ps)
  if (!all(ps_ranks %in% available_ranks)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

test<-  run_lefse(
  ps = PS.pig.diff,
  wilcoxon_cutoff = 0.01,
  
  group = "InfectionStatus",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)

lefse_format_grp <- function(sample_meta, group, subgroup = NULL) {
  groups <- sample_meta[[group]]
  group_nms <- unique(groups)
  
  if (is.null(subgroup)) {
    subgroup <- paste0(groups, "_subgrp")
  } else {
    subgroup <- sample_meta[[subgroup]]
  }
  
  group_hie <- split(subgroup, groups) %>%
    purrr::map(unique)
  
  return(list(group = groups, subgroup = subgroup, group_hie = group_hie))
}


bootstap_lda <- function(feature_abundance,
                         boot_n,
                         class,
                         sample_fract,
                         seed = 2020) {
  # Bioconductor not allows set.seed
  ldas <- purrr::rerun(
    boot_n,
    bootstap_lda_one(
      feature_abundance,
      class,
      sample_fract
    )
  ) %>%
    purrr::transpose() %>%
    purrr::map(~ do.call(bind_rows, .x)) %>%
    bind_rows()
  
  mean_lds <- colMeans(ldas)
  mean_lds <- sign(mean_lds) * log10(1 + abs(mean_lds))
  
  mean_lds
}

bootstap_lda_one <- function(feature_abundance,
                             class,
                             sample_fract) {
  sample_groups <- unique(class)
  class_count <- table(class)
  feature_abundance$class <- class
  feature_abundance <- preprocess_feature_all(feature_abundance, class)
  
  sample_n <- nrow(feature_abundance)
  random_n <- floor(sample_n * sample_fract)
  class_n <- length(sample_groups)
  sample_min <- floor(
    min(class_count) * sample_fract * sample_fract * 0.5) %>%
    max(1)
  
  # class vs class
  pairs <- utils::combn(sample_groups, 2, simplify = FALSE) %>%
    purrr::map(sort, decreasing = TRUE)
  
  for (i in seq_len(1000)) {
    # random select samples using bootstrap method
    sample_indx <- sample(sample_n, random_n, replace = TRUE)
    
    is_checked <- check_bootstrap_sample(
      feature_abundance,
      sample_indx,
      sample_min,
      class
    )
    if (is_checked) {
      break
    }
  }
  
  if (!is_checked) {
    stop(
      "Too small samples in each class",
      " or the variance of feature abundances within a",
      " class too small (zero or near zero)",
      call. = FALSE
    )
  }
  
  lda <- purrr::map(
    pairs,
    ~ cal_pair_lda(feature_abundance, sample_indx, .x)
  )
  names(lda) <- purrr::map(pairs, paste, collapse = " -VS- ")
  
  lda
}
