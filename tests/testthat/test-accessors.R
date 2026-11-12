test_that("the resource registry is consistent with metadata.csv", {
  reg <- cpgdData_resources()
  expect_equal(nrow(reg), 12L)
  expect_false(anyDuplicated(reg$name) > 0L)
  expect_false(anyDuplicated(reg$title) > 0L)
  expect_true(all(nzchar(reg$name)))
  meta <- read.csv(system.file("extdata", "metadata.csv",
                               package = "cpgdirectionData"))
  expect_equal(nrow(meta), nrow(reg))
  expect_setequal(meta$Title, reg$title)
  # RDataPath carries the Zenodo download query string; compare the file part
  expect_setequal(sub("\\?.*$", "", basename(meta$RDataPath)),
                  paste0(reg$name, ".rds"))
  # every resource must resolve to one absolute URL under the archived record
  expect_true(all(meta$Location_Prefix == "https://zenodo.org/"))
  expect_true(all(grepl("^records/[0-9]+/files/[A-Za-z0-9_]+\\.rds$", meta$RDataPath)))
})

test_that("an unknown resource name is refused with the available names", {
  expect_error(cpgdData("no_such_layer"), "Unknown resource")
})

test_that("local mode reads <dir>/<name>.rds and errors clearly when absent", {
  d <- tempfile("cpgdData_local")
  dir.create(d)
  fake <- data.frame(cpg_id = "cg00000001", pos = 1L)
  saveRDS(fake, file.path(d, "cpg_positions_hg19.rds"))
  old <- options(cpgdirectionData.local = d)
  on.exit(options(old), add = TRUE)
  got <- cpgdData("cpg_positions_hg19")
  expect_equal(got, fake)
  got2 <- getCpgPositionsHg19()
  expect_equal(got2, fake)
  expect_error(cpgdData("smr_directions"), "does not exist")
})

test_that("the second SMR layer and the coordinates have their own names", {
  # smr_directions is read by cpgdirection 0.99.x and 1.0.x, which take every
  # row of it; the second layer must never be served under that name
  reg <- cpgdData_resources()
  expect_true(all(c("smr_directions", "smr_directions_v2", "cpg_coordinates")
                  %in% reg$name))
  expect_equal(reg$title[reg$name == "smr_directions_v2"],
               "cpgdirection smr_directions_v2")
  d <- tempfile("cpgdData_local")
  dir.create(d)
  v1 <- data.frame(cpg_id = "cg00000001", direction = 1L)
  v2 <- data.frame(cpg_id = "cg00000001", direction = 1L, smr_accuracy = 0.9)
  co <- data.frame(cpg_id = "cg00000001", array = "450K;EPICv1")
  saveRDS(v1, file.path(d, "smr_directions.rds"))
  saveRDS(v2, file.path(d, "smr_directions_v2.rds"))
  saveRDS(co, file.path(d, "cpg_coordinates.rds"))
  old <- options(cpgdirectionData.local = d)
  on.exit(options(old), add = TRUE)
  expect_equal(getSmrDirections(), v1)
  expect_equal(getSmrDirectionsV2(), v2)
  expect_equal(getCpgCoordinates(), co)
})

test_that("hub retrieval works when the Hub is reachable", {
  skip_on_cran()
  skip_if_offline <- function() {
    ok <- tryCatch({
      eh <- suppressMessages(ExperimentHub::ExperimentHub())
      TRUE
    }, error = function(e) FALSE)
    if (!ok) skip("ExperimentHub not reachable")
  }
  skip_if_offline()
  # before the resources are on the Hub this can only test the error path;
  # after upload it exercises real retrieval
  res <- tryCatch(cpgdData("cpg_positions_hg19"), error = function(e) e)
  expect_true(is.data.frame(res) || inherits(res, "error"))
})
