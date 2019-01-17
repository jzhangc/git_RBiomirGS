.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.
                        To cite in publication: Zhang J, Storey KB. 2018. RBiomirGS: An all-in-one miRNA gene set analysis solution featuring target mRNA mapping and expression profile integration. PeerJ. 6: e4262.
                        For more details, please visit: http://kenstoreylab.com/?page_id=2448")
  return(TRUE)
}
