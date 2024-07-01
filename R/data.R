#' Sample Data
#'
#' @format
#' A list containing
#' - `forcing` as an [xts::xts()] object, and
#' - `soil` as a named double vector.
"ATNeu_example"

#' Bourne SNOTEL station sample data
#'
#' @source <https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=361>
#'
#' @format
#' A list containing
#' - `forcing` as an [xts::xts()] object, and
#' - `max_depth_sm` numeric vector length 1
#' - `soil` as a named double vector.
#' - `md` a metadata [data.frame()] with 1 row and 28 columns
"Bourne"

#' South American sample gridded data
#'
#' Used to construct a [splash.grid()] call. The CRU in the name likely stands for the
#' [Climatic Research
#' Unit](https://research-portal.uea.ac.uk/en/organisations/climatic-research-unit) at
#' the University of East Anglia.
#'
#' @source <https://crudata.uea.ac.uk/cru/data/hrg/>
#'
#' @format
#' A list containing
#' - `tc` a [RasterBrick][raster::Raster-class] of monthly temperature data (mean deg C)
#' - `pn` a [RasterBrick][raster::Raster-class] of monthly precipitation data (mean mm/d or mm/month?)
#' - `sw_in` a [RasterBrick][raster::Raster-class] of monthly precipitation data (mean W m^-2)
#' - `elev` a [RasterLayer][raster::Raster-class] representing a DEM over South America (presumably in average meters above sealevel)
#' - `soil` a [RasterBrick][raster::Raster-class] representing the soil fractions
#'     - ⚠️ sometimes adding up to less, sometimes to more than 100% ⚠️
"SA_cru"
