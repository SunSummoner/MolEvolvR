% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blastWrappers.R
\name{run_rpsblast}
\alias{run_rpsblast}
\title{Run RPSBLAST to generate domain architectures for proteins of interest}
\usage{
run_rpsblast(
  rpsblast_path,
  db_search_path,
  db = "refseq",
  query,
  evalue = "1e-5",
  out,
  num_threads = 1
)
}
\arguments{
\item{db_search_path}{Path to the BLAST databases}

\item{num_threads}{}
}
\description{
Run RPSBLAST to generate domain architectures for proteins of interest
}
