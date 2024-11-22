require(tidyverse)
require(data.table)


# has-map to annotate abbreviation
# https://en.wikipedia.org/wiki/Nucleic_acid_notation
base_abb = list(
  "U" = "T",
  
  "W" = c("A", "T"),
  "S" = c("C", "G"),
  "M" = c("A", "C"),
  "K" = c("G", "T"),
  "R" = c("A", "G"),
  "Y" = c("C", "T"),
  
  "B" = c("C", "G", "T"),
  "D" = c("A", "G", "T"),
  "H" = c("A", "C", "T"),
  "V" = c("A", "C", "G"),
  
  "N" = c("A", "C", "T", "G")
)

binding_dt = rbind(
  # MBNL:  
  # binding motif consists of GC dincleotides embedded in pyrimidines; YGCY (Y = C or U).  
  # The location of motifs in downstream introns are associated with increased exon exclusion 
  # in DM1 and motifs in upstream introns or exons are associated with increased exon inclusion. 
  # Of all pentamers, only one non-MBNL1 biding motif (CUAAU and its variant UAAUC) was a 
  # significant predictor of splicing outcomes, which may indicate a role for Quaking 
  # protein in splicing dysregulation (Tanner, 2021) (YGCU and GCUU; GCUUGCU or GCUUUGC; Lambert, 2014).
  data.table(
    binding_gene = "MBNL",
    binding_site = c("YGCU", "GCUU", "GCUUGCU", "GCUUUGC"),
    source = "Thurman's file",
    note = "",
    group = "intron/exon"
  ),
  data.table(
    binding_gene = "MBNL_3UTR",
    binding_site = c("GCTT", "CGCT", "TGCT", "GCGC", "CCGC", "CTGC", 
                     "GCTA", "ACGC", "CGCA", "AGCT", "TTGC", "CAGC"),
    source = "Wang, et al., 2015",
    note = "MBNL1 motifs for 3' UTR binding",
    group = "3' UTR"
  ),
  # CELF1:   
  # 3 RNA recognition motifs that bind consensus UGU or UGU-rich sequences; UGU, 
  # UGUUGU, or UGUUUGU (Marquis, 2006; Lambert, 2014)
  data.table(
    binding_gene = "CELF1",
    binding_site = c("UGU", "UGUUGU", "UGUUUGU"),
    source = "Marquis, 2006; Lambert, 2014",
    note = "",
    group = "intron/exon"
  ), 
  data.table(
    binding_gene = "CELF1_3UTR",
    binding_site = c("TGTT", "ATGT", "TTGT", "TGTC", "GTGT", "TGTA", "GTTT", "TGTG", "GTCT", "TTTT"),
    source = "Marquis, 2006; Lambert, 2014",
    note = "CELF1 motifs for 3' UTR binding",
    group = "3' UTR"
  ), 
  # RBFOX2:   
  # consensus binding motifs are UGCAUG or simply GCAUG (Lambert, 2014)
  data.table(
    binding_gene = "RBFOX2",
    binding_site = c("UGCAUG", "GCAUG"),
    source = "Thurman's file",
    note = "",
    group = "intron/exon"
  ),
  # add concatenate motif
  data.table(
    binding_gene = "CELF1_concatenate_in_GSTA1",
    binding_site = c("TGT", "TGTTGT", "TGTTTGT", "TGTGGTTT", "ATGTGGTTT", "GTTTGTTT"),
    source = "3UTR_analysis_02032002_TW_ZL.xlsx",
    note = "",
    group = "3' UTR"
  ),
  data.table(
    binding_gene = "CELF1_concatenate_in_PIGT",
    binding_site = c("CAGCAGCT", "TGTGTTGTC"),
    source = "3UTR_analysis_02032002_TW_ZL.xlsx",
    note = "",
    group = "3' UTR"
  )
)
binding_dt[, k := binding_site %>% nchar()]




#' get information from the exon sequence
#'
#' @param seq string; a string of sequence; exon in uppercase
#' @param gene_name string; name of gene
#' @param exon_idxes vector; a vector of exon idexes
#' @param start_with_enxon T/F; if the seqeunce starts with exon or not
#' @param end_with_enxon T/F; if the seqeunce ends with exon or not
#' 
#' @return a data.table annotate the sequence by exon or introns
annot_seq = function(seq, gene_name, exon_idxes, 
                     start_with_exon = T,
                     end_with_exon = T,
                     seqs_r_stopc = NULL){
  # exon -> intron
  exon_ends = str_locate_all(seq, "[A-Z][a-z]")[[1]][, "start"]
  intron_starts = str_locate_all(seq, "[A-Z][a-z]")[[1]][, "end"]
  
  # intron -> exon
  exon_starts = str_locate_all(seq, "[a-z][A-Z]")[[1]][, "end"]
  intron_ends = str_locate_all(seq, "[a-z][A-Z]")[[1]][, "start"]
  
  if(end_with_exon){
    exon_ends = c(exon_ends, nchar(seq))
  }else{
    intron_ends = c(intron_ends, nchar(seq))
  }
  
  if(start_with_exon){
    exon_starts = c(1, exon_starts)
  }else{
    intron_starts = c(1, intron_starts)
  }
  
  # exon coordinates
  exons = c(exon_starts, exon_ends) %>% sort() %>%
    matrix(ncol = 2, byrow = T) %>%
    as.data.frame()
  introns = c(intron_starts, intron_ends) %>% sort() %>%
    matrix(ncol = 2, byrow = T) %>%
    as.data.frame()
  exons[,3] = "enxon"
  exons[, 4] = paste0("exon", exon_idxes)
  introns[,3] = "intron"
  n_exons = length(exon_idxes)
  introns_ids = paste0("intron", exon_idxes[1:(n_exons - 1)], "-", exon_idxes[2:n_exons])
  if(!start_with_exon){
    introns_ids = c(paste0("intron", exon_idxes[1] - 1, "-", exon_idxes[1]), introns_ids)
  }
  if(!end_with_exon){
    introns_ids = c(introns_ids, paste0("intron", exon_idxes[n_exons], "-", exon_idxes[n_exons] + 1))
  }
  introns[, 4] = introns_ids
  
  elements = rbind(exons, introns) %>% setDT()
  colnames(elements) = c("start", "end", "category", "id")
  
  if(!is.null(seqs_r_stopc)){
    split_seq = str_split(seq, seqs_r_stopc) %>% unlist()
    if(length(split_seq)!=2){
      stop("please provide longer sequence around stop codon; 
           it has to end at stop codon, including stop codon\n")
    }
    utr_start = str_locate_all(seq, seqs_r_stopc)[[1]][, "end"] + 1
    elements = rbind(elements, 
      data.table(
        start = utr_start,
        end = utr_start + 999,
        category = "3_UTR",
        id = "3_UTR"
      )
    )
    split_seq = substr(split_seq[2], 1, 1000)
  }
  
  elements[, seq := str_sub(string = seq, start, end), by = 1:nrow(elements)]
  
  # check 3'UTR
  # if(!is.null(seqs_r_stopc)){
  #   if(elements[category == "3_UTR", seq] != split_seq){
  #     stop("3UTR seqeunce error!")
  #   }
  # }
  # 
  
  elements[, list(seq(from = start, to = end, by = 1)), by = 1:nrow(elements)]
  
  # elements[category == "exon", check := grepl("^[[:upper:]]+$", seq)]
  # elements[category == "intron", check := grepl("^[[:lower:]]+$", seq)]
  
  elements[, target := gene_name]
  
  return(elements)
  
}



#' locate all match from genome sequence based on given binding protein
#' 
#' @param binding_sites vector; a vector of string representing the binding site
#' of the geen
#' @param seq string; a string of target genome sequence,
#' @param base_abb list; a list indicating the base abbreviations
#' @param gene_name string; name of the binding protein
#' @param seq_annot data.table; a data.table annotate the seqeunce of target protein
#' 
#' @return a data.table with the coordinates and the id of the binding
hightlight_motif_ = function(binding_sites, seq, base_abb, gene_name, seq_annot){
  all_sites = process_bindingsite(binding_sites = binding_sites,
                                  base_abb = base_abb,
                                  gene_name = gene_name)
  all_sites = all_sites %>% unlist()
  matches = lapply(all_sites, function(site){
    x = str_locate_all(string = seq %>% toupper(), 
                       pattern = site %>% toupper()) %>%
      Reduce(f = rbind) %>%
      as.data.frame()
    setDT(x)
    x[, binding_site := site]
    x[, n_matches := nrow(x)]
    if(nrow(x) != 0){
      x[, double_check := substr(x = seq, start = start, stop = end)]
    }else{
      x = data.table(start = NA,
                     end = NA,
                     binding_site = site,
                     n_matches = 0,
                     double_check = NA)
    }
    return(x)
  })
  matches = do.call(rbind, matches)
  
  matches[, binding_gene := gene_name]
  matches[, match_element := seq_annot[, id] %>% unique()]
  matches[, target_gene := seq_annot[, target] %>% unique()]
  return(matches)
}


#' locate all match from genome sequence based on a list of binding proteins
#' 
#' @param binding_dt data.table; a data.table of binding gene (name) and 
#' their binding sites (value)
#' @param seq string; a string of target genome sequence,
#' @param base_abb list; a list indicating the base abbreviations
#' @param seq_annot data.table; a data.table annotate the seqeunce of target protein
#' 
#' @return a data.table with the coordinates and the id of the binding
hightlight_motif = function(binding_dt, seq, base_abb, seq_annot){
  binding_genes = binding_dt[, binding_gene] %>% unique()
  matches = lapply(setNames(binding_genes, binding_genes), function(gene){
    cat("processing", gene, "..\n")
    hightlight_motif_(
      binding_sites = binding_dt[binding_gene == gene, binding_site],
      base_abb = base_abb,
      seq = seq,
      seq_annot = seq_annot,
      gene_name = gene
    )
  })
  return(do.call(rbind, matches))
}


#' if a match is within intron/exon
#' 
#' @param seq_annot data.table; a data.table annotate the start & end of 
#' @param start_pos integer; start coordinates of match
#' @param end_pos integer; end coordinates of match 
#' 
#' @return id(a string) of intron/exon
locate_match = function(start_pos, end_pos, seq_annot){
  return(seq_annot[start <= start_pos & end >= end_pos, id])
}


#' list all possible binding site
process_bindingsite = function(binding_sites, base_abb, gene_name){
  sites_bases = strsplit(binding_sites, "")
  all_sites = lapply(sites_bases, function(site_bases){
    site_base_elements = lapply(site_bases, function(site_base){
      if(site_base %in% names(base_abb)){
        return(base_abb[[site_base]])
      }else{
        return(site_base)
      }
    })
    all_sites_dt = do.call(CJ, site_base_elements)
    all_sites_dt[, sites := paste(.SD, collapse = ""), by = 1:nrow(all_sites_dt)]
    all_sites_dt[, gene := gene_name]
    return(all_sites_dt[, sites])
  })
  return(all_sites)
}


#' load sequence from .rtf file
#' 
#' @param f string; path to the .rtf file
#' @param load_all T/F; return all or not
#' 
#' @return a string of genome sequence
load_seq = function(f, load_all = F){
  seq = textreadr::read_rtf(f)
  if(load_all){
    # seq = strsplit(seq, "in all CAPS):|Genomic|seq_start|seq_end") %>% unlist() 
    seq = strsplit(seq, "seq_start|seq_end") %>% unlist() 
  }else{
    seq = strsplit(seq, "seq_start|seq_end") %>% unlist() %>%
      .[2]
  }
  return(seq)
}


#' calculate number of k-mers in given length
n_kmer = function(k, length){
  return(sapply(length, function(x) x - k + 1))
}

#' seq_annot for full length seqeunce
#'
#' @param seq string; a string of sequence
#' @param gene string
#'
#' @return a data.table that in format of seq_annot
make_seq_annot = function(seq, gene, id, category){
  dt = data.table(start = 1, 
                  end = seq %>% nchar(),
                  category = category,
                  id = id,
                  seq = seq,
                  target = gene)
  return(dt)
}

merge_bindings = function(bindings, element_id, seq_length){
  merged = bindings[match_element == element_id, ] %>%
    .[, .(n_matches = unique(n_matches)), by = .(binding_site, binding_gene)]
  merged[, k := binding_site %>% nchar()]
  merged[, n_kmers := n_kmer(k, seq_length)]
  return(merged)
}

#' chisqure test based on given binding site and binding data.table
#'
#' @param target_b_merged data.table; a data.table by merge_bindings()
#' @param other_b_merged data.table; a data.table by merge_bindings_l()
#' @param site string; sequence of binding site
#' @param gene string; name of the gene binding to the target
#' binded
chi_test = function(target_b_merged, full_b_merged, 
                    site, gene){
  # contingency table
  ctable = matrix(0, nrow = 2, ncol = 2)
  rownames(ctable) = c("match", "no_match")
  colnames(ctable) = c("target_region", "other_region")
  
  ctable["match", "target_region"] = target_b_merged[binding_site == site & binding_gene == gene, n_matches]
  ctable["no_match", "target_region"] = target_b_merged[binding_site == site & binding_gene == gene, n_kmers] - ctable["match", "target_region"]
  
  ctable["match", "other_region"] = full_b_merged[binding_site == site & binding_gene == gene, n_matches] - ctable["match", "target_region"]
  ### todo
  if(ctable["match", "other_region"] < 0) ctable["match", "other_region"] = 0
  ctable["no_match", "other_region"] = full_b_merged[binding_site == site & binding_gene == gene, n_kmers] - target_b_merged[binding_site == site & binding_gene == gene, n_kmers] - ctable["match", "other_region"]
  
  print(ctable)
  
  # chi-square test
  chisq_test = chisq.test(ctable)
  
  return(list(
    p.value = chisq_test$p.value,
    xsquared = chisq_test$statistic,
    ctable = ctable
  ))
}


#' seqeunce element analysis
element_analysis = function(target_element_id, full_seq,
                            seq_e, target_gene, 
                            binding_dt, base_abb){
  # seqeunce of the target element
  target_seq = seq_e[id == target_element_id, seq]
  
  # binding sites in target region
  cat("analysis target region ..\n")
  target_b = hightlight_motif(
    binding_dt = binding_dt, 
    base_abb = base_abb,
    seq = target_seq,
    seq_annot = make_seq_annot(
      seq = target_seq,
      gene = target_element_id,
      category = target_element_id,
      id = target_element_id
    )
  )
  target_b_merged = merge_bindings(target_b, target_element_id,
                                   seq_length = target_seq %>% nchar())
  
  # binding sites in all region
  cat("analysis full region ..\n")
  full_b = hightlight_motif(
    binding_dt = binding_dt, 
    base_abb = base_abb,
    seq = full_seq,
    seq_annot = make_seq_annot(
      seq = full_seq,
      gene = target_gene,
      category = target_gene,
      id = target_gene
    )
  )
  full_b_merged = merge_bindings(full_b, target_gene,
                                 seq_length = full_seq %>% nchar())
  
  r = target_b_merged[, {site = binding_site; bgene = binding_gene; 
  .(chi_results = list(chi_test(target_b_merged[binding_site == site & binding_gene == bgene, ],
                                full_b_merged, 
                                site = site, 
                                gene = bgene)))},
  by = .(binding_site, binding_gene)]
  r[, ctable := map(chi_results, ~(.x$ctable))]
  r[, chi_p := map(chi_results, ~(.x$p.value)) %>% unlist()]
  r[, adj_chi_p := p.adjust(chi_p, method = "BH")]
  r[, xsquared := map(chi_results, ~(.x$xsquared)) %>% unlist()]
  r[, chi_results := NULL]
  
  r[, target_gene := target_gene]
  r[, target_element_id := target_element_id]
  
  
  return(r)
}

#' seqeunce element analysis
elements_analysis = function(full_seq,
                             seq_e, target_gene, 
                             binding_dt, base_abb){
  target_element_ids = seq_e$id %>% unique()
  l = pblapply(target_element_ids, function(target_element_id){
    cat("\nanalyze", target_element_id, ".\n")
    element_analysis(
      target_element_id = target_element_id,
      full_seq = full_seq,
      seq_e = seq_e,
      target_gene = target_gene,
      binding_dt = binding_dt,
      base_abb = base_abb
    )
  })
  dt = do.call(rbind, l)
  return(dt)
}


substr_tail = function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

seq_analysis = function(full_f, seq_f, gene_name,
                        exon_idxes, start_with_exon = T,
                        end_with_exon = T, 
                        seqs_r_stopc = NULL, 
                        binding_dt, base_abb,
                        exon_e = NULL,
                        data_dir = file.path("..", "results", "Gene_sequences_for_binding_motifs")){
  cat("load full length sequence.\n")
  seq_full = read.fasta(file.path(data_dir, full_f), as.string = T)
  cat("load target sequence.\n")
  exon_seq = load_seq(file.path(data_dir, seq_f))
  
  cat("load sequence elements.\n")
  if(is.null(exon_e)){
    exon_e = annot_seq(seq = exon_seq,
                       gene_name = gene_name,
                       exon_idxes = exon_idxes,
                       start_with_exon = start_with_exon,
                       end_with_exon = end_with_exon,
                       seqs_r_stopc = seqs_r_stopc)
  }
  
  cat("perform analysis.\n")
  tests = elements_analysis(
    full_seq = seq_full,
    seq_e = exon_e,
    target_gene = gene_name,
    binding_dt = binding_dt,
    base_abb = base_abb
  )
  
  return(tests)
}
