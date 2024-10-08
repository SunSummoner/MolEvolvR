# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(rentrez))
# suppressPackageStartupMessages(library(future))
# suppressPackageStartupMessages(library(furrr))
# suppressPackageStartupMessages(library(data.table))

#####################################
## Download Assembly Summary Files ##
#####################################
#' Download the combined assembly summaries of genbank and refseq
#'
#' @author Samuel Chen, Janani Ravi
#'
#' @param outpath String of path where the assembly summary file should be written
#' @param keep Character vector containing which columns should be retained and downloaded
#'
#' @importFrom data.table fwrite setnames
#' @importFrom dplyr bind_rows select
#' @importFrom biomartr getKingdomAssemblySummary
#'
#' @return
#' @export
#'
#' @examples
DownloadAssemblySummary <- function(outpath,
    keep = c(
        "assembly_accession", "taxid",
        "species_taxid", "organism_name"
    )) {
    assembly_kingdom_genbank <- getKingdomAssemblySummary("genbank")
    assembly_kingdom_refseq <- getKingdomAssemblySummary("refseq")

    if (keep == "all") {
        assembly_all <- bind_rows(assembly_kingdom_genbank, assembly_kingdom_refseq)
    } else {
        assembly_all <- bind_rows(assembly_kingdom_genbank, assembly_kingdom_refseq) %>%
            select(all_of(keep))
    }

    assembly_all <- assembly_all %>% data.table::setnames(
        old = c(
            "taxid", "refseq_category", "species_taxid", "organism_name",
            "infraspecific_name", "genome_rep"
        ),
        new = c(
            "TaxID", "RefseqCategory", "Parent.TaxID", "Species",
            "Spp.Strain", "GenomeStatus"
        ),
        skip_absent = T
    )

    # dplyr::rename("AssemblyID"="assembly_accession",
    #               "TaxID"="taxid",
    #               "RefseqCategory"="refseq_category",
    #               "Parent.TaxID"="species_taxid",
    #               "Species"="organism_name",
    #               "Spp.Strain"="infraspecific_name",
    #               "GenomeStatus"="genome_rep")

    fwrite(assembly_all, outpath, sep = "\t")
}



###################################
## Map GCA_ID to TaxID & Lineage ##
###################################
#' Function to map GCA_ID to TaxID, and TaxID to Lineage
#'
#' @author Samuel Chen, Janani Ravi
#' @note
#' Currently configured to have at most kingdom and phylum
#'
#'
#' @param prot_data Dataframe containing a column `GCA_ID`
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the "DownloadAssemblySummary()" function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "create_lineage_lookup()" function
#' @param acc_col
#'
#' @importFrom dplyr pull
#' @importFrom data.table fread setnames
#'
#' @return
#' @export
#'
#' @examples
GCA2lin <- function(prot_data,
    assembly_path = "/data/research/jravilab/common_data/assembly_summary_genbank.txt",
    lineagelookup_path = "/data/research/jravilab/common_data/lineage_lookup.tsv",
    acc_col = "AccNum") {
    assembly_summary <- fread(assembly_path, sep = "\t")
    assembly_summary <- setnames(assembly_summary, "AssemblyID", "GCA_ID")

    mergedTax <- merge(
        x = prot_data, y = assembly_summary,
        by = "GCA_ID", all.x = T
    )
    accessions <- prot_data %>%
        pull(acc_col) %>%
        unique()
    # Prioritize Complete Genome
    best_rows <- integer(length(accessions))
    for (i in 1:length(accessions))
    {
        # browser()
        acc <- accessions[i]
        acc_inds <- which(mergedTax$Protein == acc)
        if (length(acc_inds) > 1) {
            complete <- acc_inds[which(mergedTax[acc_inds, ]$assembly_level == "Complete Genome")]
            if (length(complete) != 0) {
                best_rows[i] <- complete[1]
            } else {
                best_rows[i] <- acc_inds[1]
            }
        }
    }
    mergedTax <- mergedTax[best_rows, ]

    lineage_map <- fread(lineagelookup_path, sep = "\t")
    lineage_map <- lineage_map[, !"Species"]

    mergedLins <- merge(mergedTax, lineage_map,
        by.x = "TaxID", by.y = "TaxID",
        all.x = T
    )

    return(mergedLins)
}

###################################
## !! @SAM why is this called lins?
###################################
#' add_lins
#'
#' @param df
#' @param acc_col
#' @param assembly_path
#' @param lineagelookup_path
#' @param ipgout_path
#' @param plan
#'
#' @importFrom dplyr pull
#' @importFrom rlang sym
#'
#' @return
#' @export
#'
#' @examples
add_lins <- function(df, acc_col = "AccNum", assembly_path,
    lineagelookup_path, ipgout_path = NULL, plan = "multicore") {
    acc_sym <- sym(acc_col)
    accessions <- df %>% pull(acc_sym)
    lins <- acc2lin(accessions, assembly_path,
        lineagelookup_path, ipgout_path,
        plan = plan
    )

    # Drop a lot of the unimportant columns for now? will make merging much easier
    lins <- lins[, c(
        "Strand", "Start", "Stop", "Nucleotide Accession", "Source",
        "Id", "Strain"
    ) := NULL]
    lins <- unique(lins)

    # dup <- lins %>% group_by(Protein) %>% summarize(count=n()) %>% filter(count > 1) %>%
    #   pull(Protein)

    ## !! @SAM: there is no "Protein" column anymore !!
    merged <- merge(df, lins,
        by.x = acc_sym, by.y = "Protein", all.x = TRUE
    )
    return(merged)
}

#######################################
## Map Protein Accessions to Lineage ##
#######################################
#' acc2lin
#'
#' @description
#' Function to map protein accession numbers to lineage
#'
#' @author Samuel Chen, Janani Ravi
#' @description This function combines 'efetch_ipg()' and 'ipg2lin()' to map a set
#' of protein accessions to their assembly (GCA_ID), tax ID, and lineage.
#'
#' @param accessions Character vector of protein accessions
#' @param assembly_path String of the path to the assembly_summary path
#' This file can be generated using the "DownloadAssemblySummary()" function
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' @param ipgout_path Path to write the results of the efetch run of the accessions
#' on the ipg database. If NULL, the file will not be written. Defaults to NULL
#' @param plan
#'
#' @return
#' @export
#'
#' @examples
acc2lin <- function(accessions, assembly_path, lineagelookup_path,
    ipgout_path = NULL, plan = "multicore") {
    tmp_ipg <- F

    if (is.null(ipgout_path)) {
        tmp_ipg <- T
        ipgout_path <- tempfile("ipg", fileext = ".txt")
    }
    efetch_ipg(accessions, out_path = ipgout_path, plan = plan)

    lins <- ipg2lin(accessions, ipgout_path, assembly_path, lineagelookup_path)

    # if(tmp_ipg)
    # {
    #   unlink(tempdir(), recursive=T)
    # }

    # cols <- c("TaxID","GCA_ID", "Protein", "Protein Name", "Species", "Lineage")
    # lins <- unique(lins[,..cols])

    return(lins)
}


#########################################
## Download IPG results for Accessions ##
#########################################
#' efetch_ipg
#'
#' @author Samuel Chen, Janani Ravi
#' @description Perform efetch on the ipg database and write the results to out_path
#'
#' @param accessions Character vector containing the accession numbers to query on
#' the ipg database
#' @param out_path Path to write the efetch results to
#' @param plan
#'
#' @importFrom future future plan
#' @importFrom purrr map
#' @importFrom rentrez entrez_fetch
#'
#' @return
#' @export
#'
#' @examples
efetch_ipg <- function(accessions, out_path, plan = "multicore") {
    if (length(accessions) > 0) {
        partition <- function(v, groups) {
            # Partition data to limit number of queries per second for rentrez fetch:
            # limit of 10/second w/ key
            l <- length(v)

            partitioned <- list()
            for (i in 1:groups)
            {
                partitioned[[i]] <- v[seq.int(i, l, groups)]
            }

            return(partitioned)
        }

        plan(strategy = plan, .skip = T)

        ## ! Note: LS changed it to 600 because she has 5K results and wanted x ≤ 9
        min_groups <- length(accessions) / 600
        groups <- min(max(min_groups, 15), length(accessions))
        partitioned_acc <- partition(accessions, groups)
        sink(out_path)

        a <- map(1:length(partitioned_acc), function(x) {
            # Avoid hitting the rate API limit
            if (plan != "sequential" & x %% 9 == 0) {
                Sys.sleep(1)
            }
            f <- future({
                entrez_fetch(
                    id = partitioned_acc[[x]],
                    db = "ipg",
                    rettype = "xml", # parsed=T,
                    api_key = "YOUR_KEY_HERE"
                )
            })
        })

        for (f in a)
        {
            cat(value(f))
        }
        sink(NULL)
    }
}

#########################################
## Maps IPG results to TaxID + Lineage ##
#########################################
#' ipg2lin
#'
#' @author Samuel Chen, Janani Ravi
#' @description Takes the resulting file of an efetch run on the ipg database and
#' append lineage, and taxid columns
#'
#' @param accessions Character vector of protein accessions
#' @param ipg_file Path to the file containing results of an efetch run on the
#' ipg database. The protein accession in 'accessions' should be contained in this
#' file
#' @param refseq_assembly_path
#' @param genbank_assembly_path
#' @param lineagelookup_path String of the path to the lineage lookup file
#' (taxid to lineage mapping). This file can be generated using the
#' "create_lineage_lookup()" function
#'
#' @importFrom data.table fread setnames
#'
#' @return
#' @export
#'
#' @examples
ipg2lin <- function(accessions, ipg_file,
    refseq_assembly_path, genbank_assembly_path,
    lineagelookup_path) {
    ipg_dt <- fread(ipg_file, sep = "\t", fill = T)

    accessions <- unique(accessions)
    ipg_dt <- ipg_dt[.data$Protein %in% accessions]

    ipg_dt <- setnames(ipg_dt, "Assembly", "GCA_ID")

    # Call GCA2Lins with different assembly_paths depending on refseq or not
    # Select for Refseq rows over other DB rows
    refseq_rows <- integer(length(accessions))
    genbank_rows <- integer(length(accessions))
    for (i in 1:length(accessions))
    {
        # browser()
        acc <- accessions[i]
        acc_inds <- which(.data$mergedTax$.data$Protein == acc)
        if (length(acc_inds) != 0) {
            # refseq inds take precedence
            refseq_inds <- acc_inds[which(.data$mergedTax[acc_inds, ]$Source == "RefSeq")]
            if (length(refseq_inds) != 0) {
                # Take the first first row of the refseq (smallest index)
                refseq_rows[i] <- refseq_inds[1]
            } else {
                # take the first row of whatever is left?
                genbank_rows[i] <- acc_inds[1]
            }
        }
    }

    # Empty values be gone
    refseq_rows <- refseq_rows[which(refseq_rows != 0)]
    genbank_rows <- genbank_rows[which(genbank_rows != 0)]

    # Call GCA2lins using refseq
    ### Possible to run these in parallel if it takes a while
    if (length(refseq_rows) != 0) {
        refseq_ipg_dt <- ipg_dt[refseq_rows, ]
        refseq_lins <- GCA2lin(refseq_ipg_dt,
            .data$assembly_path = refseq_assembly_path,
            lineagelookup_path
        )
    }
    if (length(genbank_rows) != 0) {
        genbank_ipg_dt <- ipg_dt[genbank_rows, ]
        genbank_lins <- GCA2lin(.data$gca_ipg_dt,
            .data$assembly_path = genbank_assembly_path,
            lineagelookup_path
        )
    }


    lins <- GCA2lin(prot_data = ipg_dt, .data$assembly_path, lineagelookup_path)
    lins <- lins[!is.na(.data$Lineage)] %>% unique()

    return(lins)
}


#########################################
## !! @SAM: Add TaxID based on AccNum? ##
#########################################
#' add_tax
#'
#' @param data
#' @param acc_col
#' @param version
#'
#' @importFrom data.table as.data.table
#'
#' @return
#' @export
#'
#' @examples
add_tax <- function(data, acc_col = "AccNum", version = T) {
    if (!is.data.table(data)) {
        data <- as.data.table(data)
    }

    accessions <- data[[acc_col]]

    if (version) {
        data <- data[, AccNum.noV := substr(data[[acc_col]],
            start = 0, stop = nchar(data[[acc_col]]) - 2
        )]
        acc_col <- "AccNum.noV"
    }

    out_path <- tempdir()
    tax <- prot2tax(accessions, "TEMPTAX", out_path, return_dt = TRUE)

    data <- merge.data.table(data, tax,
        by.x = acc_col, by.y = "AccNum.noV", all.x = T
    )
    return(data)
}

##################################
## Maps Protein AccNum to TaxID ##
##################################
#' prot2tax
#'
#' @param accnums
#' @param suffix
#' @param out_path
#' @param return_dt
#'
#' @importFrom data.table fread
#'
#' @return
#' @export
#'
#' @examples
prot2tax <- function(accnums, suffix, out_path, return_dt = FALSE) {
    # Write accnums to a file
    acc_file <- tempfile()
    write(paste(accnums, collapse = "\n"), acc_file)
    script <- "/data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh"
    call <- paste(script, acc_file, suffix, out_path)
    system(call, wait = TRUE)
    if (return_dt) {
        out_file <- paste0(out_path, "/", suffix, ".acc2info.tsv")
        dt <- fread(out_file, sep = "\t", fill = T)
        return(dt)
    }
}


#######################################
## OLD: Maps Protein AccNum to TaxID ##
#######################################
#' prot2tax_old
#'
#' @author Samuel Chen, Janani Ravi
#' @description Perform elink to go from protein database to taxonomy database
#' and write the resulting file of taxid and lineage to out_path
#'
#' @param accessions Character vector containing the accession numbers to query on
#' the ipg database
#' @param out_path Path to write the efetch results to
#' @param plan
#'
#' @importFrom future plan
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
prot2tax_old <- function(accessions, out_path, plan = "multicore") {
    if (length(accessions) > 0) {
        partition <- function(v, groups) {
            # Partition data to limit number of queries per second for rentrez fetch:
            # limit of 10/second w/ key
            l <- length(v)

            partitioned <- list()
            for (i in 1:groups)
            {
                partitioned[[i]] <- v[seq.int(i, l, groups)]
            }

            return(partitioned)
        }

        plan(strategy = plan, .skip = T)

        ## ! Note: LS changed it to 600 because she has 5K results and wanted x to be ≤ 9
        min_groups <- length(accessions) / 600
        groups <- min(max(min_groups, 15), length(accessions))
        partitioned_acc <- partition(accessions, groups)

        out_path <- tempdir()

        a <- map(1:length(partitioned_acc), function(x) {
            # Avoid hitting the rate API limit
            if (plan != "sequential" & x %% 9 == 0) {
                Sys.sleep(1)
            }
            print(x)
            script <- "/data/research/jravilab/molevol_scripts/upstream_scripts/acc2info.sh"
            # script <- "/data/research/jravilab/molevol_scripts/upstream_scripts/prot2tax.sh"

            # accnum_in <- paste(partitioned_acc[[x]], collapse=",")
            accnum_in <- tempfile()
            write(paste(partitioned_acc[[x]], collapse = ","), accnum_in)

            system(call, wait = F)
            # system(paste(script, accnum_in), wait=TRUE)


            # f <- future({
            #   el=entrez_link(dbfrom="protein", id=partitioned_acc[[x]],
            #                    db="taxonomy",
            #                    by_id=FALSE,
            #                    api_key="YOUR_KEY_HERE")
            #   entrez_fetch(db="taxonomy",
            #                id=el, rettype="taxid",
            #                api_key="YOUR_KEY_HERE")
            #   # Calling Janani's shell script would be easier
            # })
        })

        # for( f in a)
        # {
        # cat(value(f))
        # }
        # sink(NULL)
    }
}
