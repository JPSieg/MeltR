% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state.drawer.R
\name{state.drawer}
\alias{state.drawer}
\title{Draws all possible states for a specified double helix}
\usage{
state.drawer(states = enumerate.states())
}
\arguments{
\item{States}{A MeltR enumerate.states output for a helix.}
}
\value{
Writes a directory of R2R stockholm files and a master R2R meta file you can use to draw all of the states at once.
}
\description{
Generates R2R inputs for all states in an enumerate.states output object. Run "r2r --disable-usage-warning" on the master R2R meta file in order
to generate secondary structural depictions of all states.
}
