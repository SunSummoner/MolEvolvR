## To convert IPRScan files to a gggenes viz!
## Janani Ravi, Lauren Sosinski, Samuel Chen
## Created: Apr 9, 2020

# suppressPackageStartupMessages(library(here))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(gggenes))
# suppressPackageStartupMessages(library(ggplot2))
# source("../the-approach/R/pre-msa-tree.R") # for "to_titlecase()"
# source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R")

#################################
## Modified gggenes::theme_genes
#################################
## theme_genes2 adapted from theme_genes (w/o strip.text())
## https://github.com/wilkox/gggenes/blob/master/R/theme_genes.R
#' Theme Genes2
#'
#' @importFrom ggplot2 element_blank element_line theme theme_grey
#'
#' @return
#' @export
#'
#' @examples
theme_genes2 <- function() {
    ggplot2::theme_grey() + ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey80", size = 0.2),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(colour = "grey20", size = 0.5),
        axis.ticks.x = ggplot2::element_line(colour = "grey20", size = 0.5),
        strip.background = ggplot2::element_blank()
        # strip.text = ggplot2::element_blank()
    )
}

##################################
## Get Top N AccNum by Lin+DomArch
##################################
#' Group by lineage + DA then take top 20
#'
#' @param infile_full
#' @param DA_col
#' @param lin_col
#' @param n
#' @param query
#'
#' @importFrom dplyr arrange filter group_by select summarise
#' @importFrom shiny showNotification
#' @importFrom stats na.omit
#' @importFrom rlang sym
#'
#' @return
#' @export
#'
#' @examples
find_top_acc <- function(infile_full,
    DA_col = "DomArch.Pfam",
    lin_col = "Lineage_short",
    n = 20,
    query) {
    lin_sym <- sym(lin_col)
    # cln = fread(infile_full, sep ="\t", fill = T)
    cln <- infile_full
    if (query != "All") {
        cln <- cln %>% filter(cln$QueryName == query)
    }
    cols <- colnames(cln)
    domarch_cols <- cols[which(grepl("^DomArch", cols) & !grepl("repeats$", cols))]
    cln_domarch <- cln %>% select(domarch_cols)
    col_counts <- colSums(is.na(cln_domarch))
    DA_sym <- sym(names(which.min(col_counts)))
    showNotification(paste0("Selecting representatives by unique ", DA_sym, " and lineage combinations"))
    ## Group by Lineage, DomArch and reverse sort by group counts
    grouped <- cln %>%
        group_by({{ DA_sym }}, {{ lin_sym }}) %>%
        arrange(desc(PcPositive)) %>%
        summarise(count = n(), AccNum = dplyr::first(AccNum)) %>%
        arrange(-count) %>%
        filter({{ lin_sym }} != "" && {{ DA_sym }} != "")
    top_acc <- grouped$AccNum[1:n]
    top_acc <- na.omit(top_acc)
    return(top_acc)
}


#############################################
## IPR + FULL files --> DomArch Visualization
#############################################
#' IPR2Viz
#'
#' @param infile_ipr
#' @param infile_full
#' @param accessions
#' @param analysis
#' @param group_by
#' @param topn
#' @param name
#' @param text_size
#' @param query
#'
#' @importFrom dplyr distinct filter select
#' @importFrom gggenes geom_gene_arrow geom_subgene_arrow
#' @importFrom ggplot2 aes aes_string as_labeller element_text facet_wrap ggplot guides margin scale_fill_manual theme theme_minimal ylab
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_wider
#'
#' @return
#' @export
#'
#' @examples
ipr2viz <- function(infile_ipr = NULL, infile_full = NULL, accessions = c(),
    analysis = c("Pfam", "Phobius", "TMHMM", "Gene3D"),
    group_by = "Analysis", # "Analysis"
    topn = 20, name = "Name", text_size = 15, query = "All") {
    CPCOLS <- c(
        "#AFEEEE", "#DDA0DD", "#EE2C2C", "#CDBE70", "#B0B099",
        "#8B2323", "#EE7600", "#EEC900", "chartreuse3", "#0000FF",
        "#FFD900", "#32CD32", "maroon4", "cornflowerblue", "darkslateblue",
        "#AB82FF", "#CD6889", "#FFA07A", "#FFFF00", "#228B22",
        "#FFFFE0", "#FFEC8B", "peru", "#668B8B", "honeydew",
        "#A020F0", "grey", "#8B4513", "#191970", "#00FF7F",
        "lemonchiffon", "#66CDAA", "#5F9EA0", "#A2CD5A", "#556B2F",
        "#EEAEEE", "thistle4", "#473C8B", "#FFB6C1", "#8B1C62",
        "#FFE4B5", "black", "#FF7F50", "#FFB90F", "#FF69B4", "#836FFF",
        "#757575", "#CD3333", "#EE7600", "#CDAD00", "#556B2F", "#7AC5CD"
    )
    # in the case of many, many unique domains, we will start re-using colors
    # to prevent errors
    ADDITIONAL_COLORS <- sample(CPCOLS, 1000, replace = TRUE)
    CPCOLS <- append(x = CPCOLS, values = ADDITIONAL_COLORS)
    ## Read IPR file
    ipr_out <- read_tsv(infile_ipr, col_names = T, col_types = .data$iprscan_cols)
    ipr_out <- ipr_out %>% filter(.data$Name %in% accessions)
    analysis_cols <- paste0("DomArch.", analysis)
    infile_full <- infile_full %>% select(analysis_cols, .data$Lineage_short, .data$QueryName, .data$PcPositive, .data$AccNum)
    ## To filter by Analysis
    analysis <- paste(analysis, collapse = "|")
    ## @SAM: This can't be set in stone since the analysis may change!
    ## Getting top n accession numbers using find_top_acc()
    top_acc <- find_top_acc(
        infile_full = infile_full,
        DA_col = "DomArch.Pfam",
        ## @SAM, you could pick by the Analysis w/ max rows!
        lin_col = "Lineage_short",
        n = topn, query = query
    )
    # Filter by Top Accessions per Accession per DomArch and Lineage
    ipr_out <- subset(
        ipr_out,
        ipr_out$.data$AccNum %in% top_acc
    )
    ## Need to fix this eventually based on the 'real' gene orientation! :)
    ipr_out$Strand <- rep("forward", nrow(ipr_out))

    ipr_out <- ipr_out %>% arrange(.data$AccNum, .data$StartLoc, .data$StopLoc)
    ipr_out_sub <- filter(
        ipr_out,
        grepl(pattern = analysis, x = .data$Analysis)
    )
    # dynamic analysis labeller
    analyses <- ipr_out_sub %>%
        select(.data$Analysis) %>%
        distinct()
    analysis_labeler <- analyses %>%
        pivot_wider(names_from = .data$Analysis, values_from = .data$Analysis)

    lookup_tbl_path <- "/data/research/jravilab/common_data/cln_lookup_tbl.tsv"
    lookup_tbl <- read_tsv(lookup_tbl_path, col_names = T, col_types = .data$lookup_table_cols)

    lookup_tbl <- lookup_tbl %>% select(-.data$ShortName) # Already has ShortName -- Just needs SignDesc
    # ipr_out_sub = ipr_out_sub %>% select(-.data$ShortName)
    # TODO: Fix lookup table and uncomment below
    # ipr_out_sub <- merge(ipr_out_sub, lookup_tbl, by.x = "DB.ID", by.y = "DB.ID")

    ## PLOTTING
    ## domains as separate arrows
    # For odering with tree
    # ipr_out_sub$.data$Name <- paste0(" ", ipr_out_sub$.data$Name)
    if (group_by == "Analysis") {
        plot <- ggplot(ipr_out_sub,
            aes_string(
                xmin = 1, xmax = "SLength",
                y = name, label = "ShortName"
            ),
            color = NA, fill = NA
        ) +
            geom_subgene_arrow(data = ipr_out_sub, aes_string(
                xmin = 1, xmax = "SLength", y = name, fill = "SignDesc",
                xsubmin = "StartLoc", xsubmax = "StopLoc"
            ), color = "white") +
            geom_gene_arrow(fill = NA, color = "grey") +
            # geom_blank(data = dummies) +
            facet_wrap(~.data$Analysis,
                strip.position = "top", ncol = 5,
                labeller = as_labeller(analysis_labeler)
            ) +
            # , ncol = 1 + #scales = "free",
            scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
            theme_minimal() +
            theme_genes2() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.key.size = .data$unit(0.02, "npc"),
                legend.box.margin = margin(),
                text = element_text(size = text_size)
            ) +
            ylab("") +
            guides(fill = guide_legend(nrow = 10))
    } else if (group_by == "Query") {
        plot <- ggplot(
            ipr_out_sub,
            aes(
                xmin = 1, xmax = .data$SLength,
                y = .data$Analysis, # y = .data$AccNum
                label = .data$ShortName
            )
        ) +
            geom_subgene_arrow(data = ipr_out_sub, aes_string(
                xmin = 1, xmax = "SLength", y = "Analysis", fill = "SignDesc",
                xsubmin = "StartLoc", xsubmax = "StopLoc"
            ), color = "white") +
            geom_gene_arrow(fill = NA, color = "grey") +
            facet_wrap(as.formula(paste("~", name)),
                strip.position = "top", ncol = 5,
                labeller = as_labeller(analysis_labeler)
            ) +
            scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
            theme_minimal() +
            theme_genes2() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.key.size = .data$unit(0.02, "npc"),
                legend.box.margin = margin(),
                text = element_text(size = text_size)
            ) +
            ylab("") +
            guides(fill = guide_legend(nrow = 10))
    }
    return(plot)
}

#' IPR2Viz Web
#'
#' @param infile_ipr
#' @param accessions
#' @param analysis
#' @param group_by
#' @param name
#' @param text_size
#' @param legend_name
#' @param cols
#' @param rows
#'
#' @importFrom dplyr arrange distinct filter select
#' @importFrom gggenes geom_gene_arrow geom_subgene_arrow
#' @importFrom ggplot2 aes aes_string as_labeller facet_wrap ggplot guides scale_fill_manual theme theme_minimal ylab
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_wider
#'
#' @return
#' @export
#'
#' @examples
ipr2viz_web <- function(infile_ipr,
    accessions,
    analysis = c("Pfam", "Phobius", "TMHMM", "Gene3D"),
    group_by = "Analysis", name = "Name",
    text_size = 15, legend_name = "ShortName", cols = 5, rows = 10) {
    CPCOLS <- c(
        "#AFEEEE", "#DDA0DD", "#EE2C2C", "#CDBE70", "#B0B099",
        "#8B2323", "#EE7600", "#EEC900", "chartreuse3", "#0000FF",
        "#FFD900", "#32CD32", "maroon4", "cornflowerblue", "darkslateblue",
        "#AB82FF", "#CD6889", "#FFA07A", "#FFFF00", "#228B22",
        "#FFFFE0", "#FFEC8B", "peru", "#668B8B", "honeydew",
        "#A020F0", "grey", "#8B4513", "#191970", "#00FF7F",
        "lemonchiffon", "#66CDAA", "#5F9EA0", "#A2CD5A", "#556B2F",
        "#EEAEEE", "thistle4", "#473C8B", "#FFB6C1", "#8B1C62",
        "#FFE4B5", "black", "#FF7F50", "#FFB90F", "#FF69B4", "#836FFF",
        "#757575", "#CD3333", "#EE7600", "#CDAD00", "#556B2F", "#7AC5CD"
    )
    # in the case of many, many unique domains, we will start re-using colors
    # to prevent errors
    ADDITIONAL_COLORS <- sample(CPCOLS, 1000, replace = TRUE)
    CPCOLS <- append(x = CPCOLS, values = ADDITIONAL_COLORS)
    ## To filter by Analysis
    analysis <- paste0(analysis, collapse = "|")

    ## @SAM, colnames, merges, everything neeeds to be done now based on the
    ## combined lookup table from "common_data"
    lookup_tbl_path <- "/data/research/jravilab/common_data/cln_lookup_tbl.tsv"
    lookup_tbl <- read_tsv(lookup_tbl_path, col_names = T, col_types = lookup_table_cols)

    ## Read IPR file and subset by Accessions
    ipr_out <- read_tsv(infile_ipr, col_names = T)
    ipr_out <- ipr_out %>% filter(Name %in% accessions)
    ## Need to fix eventually based on 'real' gene orientation!
    ipr_out$Strand <- rep("forward", nrow(ipr_out))

    ipr_out <- ipr_out %>% arrange(AccNum, StartLoc, StopLoc)
    ipr_out_sub <- filter(
        ipr_out,
        grepl(pattern = analysis, x = Analysis)
    )
    # dynamic analysis labeller
    analyses <- ipr_out_sub %>%
        select(Analysis) %>%
        distinct()
    analysis_labeler <- analyses %>%
        pivot_wider(names_from = Analysis, values_from = Analysis)
    # analysis_labeler[1,] = colnames(analysis_labeler)

    # ipr_out_sub$label <- paste0(" ", ipr_out_sub$Name)
    lookup_tbl <- lookup_tbl %>% select(-ShortName)
    ## @SAM, make sure the following two work with the Lookup Tables!!
    # ipr_out_sub <- merge(ipr_out_sub, lookup_tbl, by = "DB.ID")
    ## PLOTTING
    ## domains as separate arrows
    # ipr_out_sub$label <- paste0(" ", ipr_out_sub$Name)
    if (group_by == "Analysis") {
        plot <- ggplot(
            ipr_out_sub,
            aes_string(
                xmin = 1, xmax = "SLength",
                y = name, label = "ShortName"
            )
        ) +
            geom_gene_arrow(fill = "white", color = "grey") +
            geom_subgene_arrow(data = ipr_out_sub, aes_string(
                xmin = 1, xmax = "SLength", y = name, fill = "SignDesc",
                xsubmin = "StartLoc", xsubmax = "StopLoc"
            ), color = "NA") +
            # geom_blank(data = dummies) +
            facet_wrap(~Analysis,
                strip.position = "top", ncol = cols,
                labeller = as_labeller(analysis_labeler)
            ) +
            # , ncol = 1 + #scales = "free",
            scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
            theme_minimal() +
            theme_genes2() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.key.size = unit(0.02, "npc"),
                legend.box.margin = margin(),
                text = element_text(size = text_size)
            ) +
            ylab("") +
            guides(fill = guide_legend(nrow = rows))
    } else if (group_by == "Query") {
        plot <- ggplot(
            ipr_out_sub,
            aes(
                xmin = 1, xmax = SLength,
                y = Analysis, # y = AccNum
                label = ShortName
            )
        ) +
            geom_subgene_arrow(data = ipr_out_sub, aes_string(
                xmin = 1, xmax = "SLength", y = "Analysis", fill = "SignDesc",
                xsubmin = "StartLoc", xsubmax = "StopLoc"
            ), color = "white") +
            geom_gene_arrow(fill = NA, color = "grey") +
            facet_wrap(as.formula(paste("~", name)),
                strip.position = "top", ncol = cols,
                labeller = as_labeller(analysis_labeler)
            ) +
            scale_fill_manual(values = CPCOLS, na.value = "#A9A9A9") +
            theme_minimal() +
            theme_genes2() +
            theme(
                legend.position = "bottom",
                legend.box = "horizontal",
                legend.key.size = unit(0.02, "npc"),
                legend.box.margin = margin(),
                text = element_text(size = text_size)
            ) +
            ylab("") +
            guides(fill = guide_legend(nrow = rows))
    }
    return(plot)
}
