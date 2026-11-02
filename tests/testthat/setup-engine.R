# EVOTEMPO_TEST_ENGINE = "R" runs the suite on the pure-R reference engine
# (CI runs it on both engines); tests that compare the engines choose their own.
# The option is restored when the run ends
testEngine <- Sys.getenv("EVOTEMPO_TEST_ENGINE")
if (nzchar(testEngine)) {
  withr::local_options(EvoTempo.engine = TempoEngine(testEngine),
                       .local_envir = testthat::teardown_env())
}
