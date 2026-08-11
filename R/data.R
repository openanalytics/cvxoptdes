#' Fullerenes dataset
#'
#' This dataset reports the production of o-xylenyl adducts of
#' Buckminster fullerenes. Three process conditions (temperature, reaction time and
#' ratio of sultine to C60) are varied to maximize the mole fraction of the desired
#' ' product.
#' @format A data.frame containing 120 experiments with columns:
#' \describe{
#'   \item{reaction_time}{numeric, reaction time in flow reactor (min)}
#'   \item{sultine_amount}{numeric, relative concentration of sultine to C60}
#'   \item{temperature}{numeric, temperature of the reaction (Celsius)}
#'   \item{product_amount}{numeric, mole fraction of the desired product}
#' }
#'
#' @source \href{https://github.com/the-matter-lab/olympus/tree/main/src/olympus/datasets/dataset_fullerenes}{Olympus: Fullerenes dataset}
"fullerenes"

#' Imidazoles dataset
#'
#' This dataset reports the yield of the desired product and the reaction cost
#' in a direct arylation reaction of imidazoles. Two continuous process conditions
#' (temperature, initial concentration) and three reaction molecules
#' (base, ligand, solvent) are varied to maximize the yield of the desired
#' product while minimizing the cost of the reaction.
#' @format A data.frame containing 1728 experiments with columns:
#' \describe{
#'   \item{yield}{numeric, experimental yield of the reaction (\%)}
#'   \item{cost}{numeric, cost of the reaction scaled to the unit interval}
#'   \item{concentration}{numeric, initial substrate concentration (M)}
#'   \item{temperature}{numeric, temperature of the reaction (Celsius)}
#'   \item{base}{character, common name of the reaction base molecule}
#'   \item{base_smiles}{character, SMILES string of the reaction base molecule}
#'   \item{solvent}{character, common name of the reaction solvent molecule}
#'   \item{solvent_smiles}{character, SMILES string of the reaction solvent molecule}
#'   \item{ligand}{character, common name of the reaction ligand molecule}
#'   \item{ligand_smiles}{character, SMILES string of the reaction ligand molecule}
#' }
#'
#' @source \doi{doi:10.1038/s41586-021-03213-y}
#' @source \doi{doi:10.1021/jacs.2c08592}
#' @references Shields B.J. et al. (2021) \emph{“Bayesian reaction optimization as a tool for chemical synthesis”}, Nature, Vol. 590 (7844).
#' @references Torres J.A.G. et al. (2022) \emph{“A Multi-Objective Active Learning Platform and Web App for Reaction Optimization”}, J. Am. Chem. Soc., Vol. 144 (43).
"imidazoles"

#' Ethanol dataset
#'
#' This dataset reports 88 measurements from an experiment in which ethanol
#' was burned in a single cylinder automobile test engine. The data originates from the
#' \href{https://CRAN.R-project.org/package=SemiPar}{SemiPar} package, which is now archived.
#'
#' @format A data.frame containing 88 measurements with columns:
#' \describe{
#'    \item{NOx}{numeric, concentration of nitric oxide (NO) and nitrogen dioxide (NO2) in engine exhaust, normalized by the work done by the engine}
#'    \item{C}{numeric, the compression ratio of the engine}
#'    \item{E}{numeric, the equivalence ratio at which the engine was run – a measure of the richness of the air/ethanol mix}
#' }
#'
#' @source Brinkman, N.D. (1981). Ethanol fuel – a single-cylinder engine study of efficiency and exhaust emissions. SAE transactions Vol. 90, No 810345, 1410–1424.
#' @references Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). \emph{Semiparametric Regression}. Cambridge University Press.
"ethanol"


