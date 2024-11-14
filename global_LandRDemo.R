repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
# install.packages(c("SpaDES.project", "Require"), repos = repos)

# install.packages(c("googledrive", "httpuv"), repos = repos)
# googledriveAuthPath <- "~/SpaDES_book/googledrive_auth_cache"
# dir.create(googledriveAuthPath, showWarnings = FALSE)
# googledrive::drive_auth(cache = "~/SpaDES_book/googledrive_auth_cache")

library(SpaDES.project)

out <- setupProject(
  # INPUT OBJECTS -----------------------------------
  # these need to come *before* any formal arguments, as they are needed for params.R
  sppEquivCol = "Boreal",
  vegLeadingProportion = 0,
  successionTimestep = 10L,
  eventCaching = c(".inputObjects", "init"),
  useParallel = FALSE,
  # overwrite = TRUE, # Just run this once. Then remove it.
  useGit = "eliotmcintire",
  Restart = TRUE,
  # standAlone = TRUE,
  paths = list("packagePath" = "packages/",
               "projectPath" = "~/SpaDES_book/LandRDemo_coreVeg"),
  modules = c(
    "PredictiveEcology/Biomass_speciesData@main",
    "PredictiveEcology/Biomass_borealDataPrep@main",
    "PredictiveEcology/Biomass_speciesParameters@main",
    "PredictiveEcology/Biomass_core@main"
  ),
  packages = c(
    "googledrive",
    "httr",
    #   # these are needed but don't load
    #   "DiagrammeR",
    #   "lattice (>= 0.22.5)",
    #   "SpaDES.project (>= 0.1.0.9003)",
    #   "LandR (>=1.1.5.9000)",
    #   "quickPlot (>= 1.0.2.9003)",
    #   "reproducible (>= 2.1.1.9002)",
    #   "eliotmcintire/PSPclean@Eliot (>= 0.1.4.9006)",
    #   "PredictiveEcology/SpaDES.experiment@development (HEAD)",
    #   "SpaDES.core (>= 2.1.5)",
    "terra"
  ),
  options = list(
    "LandR.assertions" = TRUE,
    "reproducible.destinationPath" = paths$inputPath,
    "spades.inputPath" = paths$inputPath,
    "spades.moduleCodeChecks" = FALSE,
    "Require.cloneFrom" = Sys.getenv("R_LIBS_USER"),
    "repos" = unique(repos)
  ),
  sideEffects = {
    googledriveAuthPath <- "~/SpaDES_book/googledrive_auth_cache"
    dir.create(googledriveAuthPath, showWarnings = FALSE)
    googledrive::drive_auth(cache = "~/SpaDES_book/googledrive_auth_cache")
  },
  # SIMULATION SETUP ------------------------------------
  times = list(start = 2001, end = 2031),
  params = "PredictiveEcology/PredictiveEcology.org@training-book/tutos/LandRDemo_coreVeg/params.R",
  # (more) INPUT OBJECTS -----------------------------------
  # these come after, so that we don't need to pre-install/load LandR
  # species lists/traits
  speciesParams = {
    list(
      "shadetolerance" = list(
        # Betu_Pap = 1
        , Lari_Lar = 1
        , Pice_Gla = 2
        , Pice_Mar = 3
        , Pinu_Ban = 1.5
        , Popu_Spp = 1
      )
    )
  },
  studyArea = {
    smallExtent <- cbind(xmin = -104.757, ymin = 55.68663)
    originalcrs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
    Biomass_corecrs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    # a <- rbind(smallExtent)
    a <- terra::vect(smallExtent)
    terra::crs(a) <- originalcrs
    a <- terra::project(a, Biomass_corecrs)
    studyArea <- LandR::randomStudyArea(a, size = 2e8, seed = 1234)

  },
  studyAreaLarge = terra::buffer(studyArea, width = 3e4),
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyAreaLarge)
    species <- LandR::equivalentName(speciesInStudy$speciesList,
                                     df = LandR::sppEquivalencies_CA, sppEquivCol)
    LandR::sppEquivalencies_CA[Boreal %in% species]
  },
  # OUTPUTS TO SAVE -----------------------
  outputs = {
    rbind(
      data.frame(
        objectName = "cohortData",
        saveTime = seq(times$start, times$end)
      ),
      data.frame(
        objectName = "pixelGroupMap",
        saveTime = seq(times$start, times$end)
      )
    )
  }
)

# initialise then run simulation
# simInitOut <- SpaDES.core::simInit2(out)
simInitOut <- do.call(SpaDES.core::simInit, out)
simOut <- SpaDES.core::spades(simInitOut)

simOut <- SpaDES.core::simInitAndSpades2(out)


SpaDES.core::moduleDiagram(simInitOut)
SpaDES.core::objectDiagram(simInitOut)


SpaDES.core::events(simInitOut)
SpaDES.core::events(simOut)


SpaDES.core::completed(simInitOut)
SpaDES.core::completed(simOut)


SpaDES.core::inputs(simOut)
SpaDES.core::outputs(simOut)
SpaDES.core::parameters(simOut)

# spatial inputs from list above
terra::plot(simOut$studyAreaLarge, col = "navyblue", main = "studyArea & studyAreaLarge")
terra::plot(simOut$studyArea, col = "lightblue", add = TRUE)

# spatial outputs from list above
terra::plot(simOut$vegTypeMap,
            col = hcl.colors(palette = "Dynamic", n = length(unique(simOut$vegTypeMap[]))),
            main = "")
terra::plot(simOut$speciesLayers)

# model used to estimate species establishment probabilities
summary(simOut$modelBiomass$mod)
plot(simOut$modelBiomass$mod)

# model used to calibrate Picea glauca's growth parameters
summary(simOut$speciesGrowthCurves$Pice_Gla$NonLinearModel$Pice_Gla)


simInitOut <- SpaDES.core::simInit2(out)
SpaDES.core::P(simInitOut, param = ".plots", module = "Biomass_core") <- "screen"
SpaDES.core::P(simInitOut, param = ".plotMaps", module = "Biomass_core") <- TRUE
simOut <- SpaDES.core::spades(simInitOut)

SpaDES.core::end(simOut) <- 2061
simOut <- SpaDES.core::spades(simOut)







out2 <- out
out2$modules <- out2$modules[out2$modules != "Biomass_speciesParameters"]
out2$paths$outputPath <- normPath(file.path("~/SpaDES_book/LandRDemo_coreVeg", "outputsSim2"))

simOut2 <- SpaDES.core::simInitAndSpades2(out2)


simOut <- SpaDES.core::spades(simInitOut, debug = "plotSummaryBySpecies")


sim <- plotSummaryBySpecies(sim)

# studyArea could be
studyArea = {
  set.seed(123)
  SpaDES.tools::randomStudyArea(size = 200000000)
}

# studyAreaLarge
studyAreaLarge = {
  terra::buffer(studyArea, width = 10000)
}

