# EVOTEMPO_TEST_ENGINE = "R" runs the suite on the pure-R reference engine
# (CI runs it on both engines); tests that compare the engines choose their own
testEngine <- Sys.getenv("EVOTEMPO_TEST_ENGINE")
if (nzchar(testEngine)) {
  options(EvoTempo.engine = TempoEngine(testEngine))
}
