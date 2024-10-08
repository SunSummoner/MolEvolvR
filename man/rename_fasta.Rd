% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CHANGED-pre-msa-tree.R, R/pre-msa-tree.R
\name{rename_fasta}
\alias{rename_fasta}
\title{Rename the labels of fasta files}
\usage{
rename_fasta(fa_path, outpath, replacement_function = map_acc2name, ...)

rename_fasta(fa_path, outpath, replacement_function = map_acc2name, ...)
}
\arguments{
\item{fa_path}{Path to fasta file}

\item{outpath}{Path to write altered fasta file to}

\item{replacement_function}{Function to apply to lines starting with '>'}

\item{...}{Additional arguments to pass to replacement_function}
}
\description{
Rename the labels of fasta files

Rename the labels of fasta files
}
