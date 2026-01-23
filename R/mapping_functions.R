#' Map UniProt IDs to NCBI accessions
#'
#' This function queries UniProt and returns the corresponding NCBI protein accessions.
#' If the UniProt API is unreachable, it returns NA values with a warning.
#'
#' @param uniprot_ids Character vector of UniProt accession IDs.
#'
#' @return A tibble with UniProt IDs and mapped RefSeq accessions.
#'
#' @examples
#' up2Ncbi(c("P12345", "Q9Y263"))
#'
#' @importFrom httr2 request req_url_query req_perform resp_body_json
#' @importFrom tibble tibble
#' @importFrom dplyr right_join
#'
#' @export
up2Ncbi <- function(uniprot_ids) {

    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("accession:(", paste(uniprot_ids, collapse = " OR "), ")")

    response <- tryCatch(
        {
            httr2::request(base_url) |>
                httr2::req_url_query(
                    query = query,
                    format = "json",
                    fields = "accession,xref_refseq"
                ) |>
                httr2::req_perform()
        },
        error = function(e) {
            warning("UniProt API request failed: ", conditionMessage(e))
            NULL
        }
    )

    if (is.null(response)) {
        return(
            tibble::tibble(
                uniprot_id = uniprot_ids,
                ncbi_accession = NA_character_
            )
        )
    }

    results <- httr2::resp_body_json(response)$results

    if (length(results) == 0) {
        return(
            tibble::tibble(
                uniprot_id = uniprot_ids,
                ncbi_accession = NA_character_
            )
        )
    }

    mapped_df <- tibble::tibble(
        uniprot_id = vapply(
            results,
            `[[`,
            character(1),
            "primaryAccession"
        ),
        ncbi_accession = vapply(
            results,
            function(x) {
                refs <- x$uniProtKBCrossReferences
                if (is.null(refs)) return(NA_character_)

                ids <- vapply(
                    refs,
                    function(r) if (r$database == "RefSeq") r$id else NA_character_,
                    character(1)
                )

                ids <- unique(stats::na.omit(ids))
                if (length(ids) == 0) NA_character_ else paste(ids, collapse = ";")
            },
            character(1)
        )
    )

    dplyr::right_join(
        mapped_df,
        tibble::tibble(uniprot_id = uniprot_ids),
        by = "uniprot_id"
    )
}


#' Map NCBI RefSeq Protein IDs to UniProt Accessions
#'
#' This function maps RefSeq protein accessions (from NCBI) to corresponding
#' UniProt accession IDs using the UniProt ID mapping REST API.
#'
#' @param ncbi_ids Character vector of RefSeq protein accessions.
#' @param max_retries Maximum retries while polling UniProt.
#' @param wait_time Seconds between retries.
#'
#' @return A tibble with ncbi_id and uniprot_id.
#'
#' @examples
#' ncbi2Up(c("NP_001026859.1", "XP_002711597.1"))
#'
#' @importFrom httr2 request req_body_form req_perform resp_body_json
#' @importFrom tibble tibble
#' @importFrom dplyr left_join
#' @importFrom purrr map_df
#'
#' @export
ncbi2Up <- function(ncbi_ids, max_retries = 5, wait_time = 3) {

    run_url    <- "https://rest.uniprot.org/idmapping/run"
    status_url <- "https://rest.uniprot.org/idmapping/status/"
    result_url <- "https://rest.uniprot.org/idmapping/results/"

    response <- httr2::request(run_url) |>
        httr2::req_body_form(
            from = "RefSeq_Protein",
            to   = "UniProtKB",
            ids  = paste(ncbi_ids, collapse = ",")
        ) |>
        httr2::req_perform()

    job_id <- httr2::resp_body_json(response)$jobId
    message("Job submitted. Job ID: ", job_id)

    attempt <- 1
    repeat {
        Sys.sleep(wait_time)

        status <- httr2::request(paste0(status_url, job_id)) |>
            httr2::req_perform() |>
            httr2::resp_body_json()

        # IMPORTANT: mirror the working behavior
        if (!is.null(status$results)) break

        if (attempt >= max_retries) {
            warning("Timeout: Results not ready. Returning NA.")
            return(
                tibble::tibble(
                    ncbi_id = ncbi_ids,
                    uniprot_id = NA_character_
                )
            )
        }

        message("Waiting for UniProt mapping results... (Attempt ", attempt, ")")
        attempt <- attempt + 1
    }

    result <- httr2::request(paste0(result_url, job_id)) |>
        httr2::req_perform() |>
        httr2::resp_body_json()

    if (is.null(result$results)) {
        warning("No results returned by UniProt API.")
        return(
            tibble::tibble(
                ncbi_id = ncbi_ids,
                uniprot_id = NA_character_
            )
        )
    }

    mappings <- purrr::map_df(
        result$results,
        ~tibble::tibble(
            ncbi_id = .x$from,
            uniprot_id = .x$to
        )
    )

    tibble::tibble(ncbi_id = ncbi_ids) |>
        dplyr::left_join(mappings, by = "ncbi_id")
}


#' Map NCBI Protein Accessions to IPG (Identical Protein Group) IDs
#'
#' This function maps NCBI protein accession numbers to their corresponding
#' Identical Protein Group (IPG) IDs using the NCBI Entrez API.
#' It first attempts to retrieve the link via `rentrez::entrez_link()` and,
#' if that fails, falls back to a direct E-utilities REST query.
#'
#' @param acc Character vector of NCBI protein accessions.
#' @param api_key Optional NCBI API key.
#'
#' @return A data.frame with accession and IPG ID.
#'
#' @examples
#' \dontrun{
#' acc2Ipg(c("WP_000003915.1", "NP_414543.1"))
#' }
#'
#' @importFrom rentrez entrez_search entrez_link
#' @importFrom httr2 request req_perform resp_body_string
#' @importFrom xml2 read_xml xml_find_first xml_text
#'
#' @export
acc2Ipg <- function(acc, api_key = NULL) {

    ipgFromRest <- function(accession) {

        url <- paste0(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            "?db=ipg&term=", accession, "[Accession]"
        )

        response <- httr2::request(url) |>
            httr2::req_perform()

        xml <- xml2::read_xml(httr2::resp_body_string(response))
        id  <- xml2::xml_text(xml2::xml_find_first(xml, ".//IdList/Id"))

        if (length(id) == 0 || id == "") NA_character_ else id
    }

    results <- lapply(acc, function(accession) {

        tryCatch({
            search <- rentrez::entrez_search(
                db = "protein",
                term = paste0(accession, "[Accession]"),
                api_key = api_key
            )

            if (length(search$ids) == 0) {
                return(
                    data.frame(
                        accession = accession,
                        ipg_id = ipgFromRest(accession)
                    )
                )
            }

            link <- rentrez::entrez_link(
                dbfrom = "protein",
                db = "ipg",
                id = search$ids[1],
                api_key = api_key
            )

            ipg_id <- if (!is.null(link$links$protein_ipg)) {
                link$links$protein_ipg[1]
            } else {
                ipgFromRest(accession)
            }

            data.frame(accession = accession, ipg_id = ipg_id)

        }, error = function(e) {
            data.frame(accession = accession, ipg_id = NA_character_)
        })
    })

    do.call(rbind, results)
}
