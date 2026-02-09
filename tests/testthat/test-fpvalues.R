test_that("Global Null (All lfdr = 1) yields Uniform P-values", {
  n <- 100
  lfdr <- rep(1, n)

  # Run function
  res <- fpvalues(lfdr)

  # 1. Check Range
  expect_equal(min(res$fp), 1 / n)
  expect_equal(max(res$fp), 1)

  # 2. Check Uniformity (Exact Spacing)
  # The p-values should be exactly 1/n, 2/n, ..., n/n
  # This confirms your tie-breaking (seq_along) is working correctly
  expected_p <- seq_len(n) / n
  expect_equal(sort(res$fp), expected_p)

  # 3. Check q-values
  # Since all lfdr=1, the cumulative mean (q-value) should be exactly 1 everywhere.
  expect_equal(res$fq, rep(1, n))
})

test_that("Tie-Breaking using secondary 'p' input works correctly", {
  # Scenario: Two features have identical lfdr (0.5).
  # Feature A (index 1) has p=0.05
  # Feature B (index 2) has p=0.01 (More significant)

  lfdr <- c(0.5, 0.5)
  p_raw <- c(0.05, 0.01)

  res <- fpvalues(lfdr, p = p_raw)

  # Feature B (index 2) must have a smaller functional p-value than Feature A
  expect_lt(res$fp[2], res$fp[1])

  # Verify exact calculations:
  # Rank 1 (Item 2): lfdr=0.5 -> FDR = 0.5/1 = 0.5.   p = (0.5 * 0.5) / 0.75 = 0.333... Wait!
  # Let's trace the math for this specific case:
  # Sorted: Item 2, Item 1
  # CumSum LFDR: 0.5, 1.0
  # FDR: 0.5, 0.5
  # Max(FDR): 0.5
  # Ranks: 1/2, 2/2
  # P-values: (0.5 * 0.5)/0.5 = 0.5, (0.5 * 1.0)/0.5 = 1.0

  expect_equal(res$fp[2], 0.5)
  expect_equal(res$fp[1], 1.0)
})

test_that("Monotonicity is preserved with ordered inputs", {
  lfdr <- c(0.01, 0.05, 0.2, 0.8, 1.0)
  res <- fpvalues(lfdr)

  # P-values should be strictly increasing
  expect_true(all(diff(res$fp) > 0))

  # Q-values should be non-decreasing
  expect_true(all(diff(res$fq) >= 0))
})

test_that("Robustness to NAs in inputs", {
  lfdr <- c(0.1, NA, 0.9)
  p <- c(0.01, 0.5, 0.8)

  res <- fpvalues(lfdr, p)

  # Output length should match input
  expect_equal(length(res$fp), 3)

  # NA in input should result in NA in output
  expect_true(is.na(res$fp[2]))
  expect_true(is.na(res$fq[2]))

  # Valid values should be calculated correctly (ignoring the NA)
  # The remaining set is effectively {0.1, 0.9}
  expect_false(is.na(res$fp[1]))
  expect_false(is.na(res$fp[3]))
})

test_that("Boundary Conditions (lfdr > 1)", {
  # lfdr > 1 should be clamped to 1
  lfdr <- c(0.5, 1.5)
  res <- fpvalues(lfdr)

  # The 1.5 should be treated exactly like 1.0
  # This ensures the clamp `lfdr[lfdr > 1] <- 1` is working
  expect_equal(res$fp[2], 1.0)
})

test_that("Input Validation matches lengths", {
  lfdr <- c(0.1, 0.2)
  p <- c(0.01, 0.02, 0.03) # Mismatch length

  expect_error(fpvalues(lfdr, p), "Length of 'lfdr' and 'p' must match")
})

test_that("Signal Case: Correctly identifies and ranks mixed signal/noise", {
  # Construct a dataset:
  # 1. Strong Signal: lfdr=0.01, p=1e-6
  # 2. Moderate Signal (Tied lfdr): lfdr=0.20, p=0.01
  # 3. Moderate Signal (Tied lfdr): lfdr=0.20, p=0.05 (Weaker raw p)
  # 4. Noise: lfdr=0.90, p=0.60
  # 5. Pure Noise: lfdr=1.00, p=0.90

  lfdr_in <- c(0.01, 0.20, 0.20, 0.90, 1.00)
  p_in <- c(1e-6, 0.01, 0.05, 0.60, 0.90)

  # Expected Order based on (lfdr, p): 1, 2, 3, 4, 5
  # (Already sorted for manual verification ease)

  res <- fpvalues(lfdr_in, p = p_in)

  # --- CHECK 1: Rank Preservation ---
  # The functional p-values must be strictly increasing because our input
  # was constructed to be strictly ordered by (lfdr, p).
  expect_true(all(diff(res$fp) > 0))

  # --- CHECK 2: Signal Identification ---
  # The first item (Strong Signal) should have a very low functional p-value.
  # Manual Calc:
  #   lfdr=0.01 -> FDR_1 = 0.01.
  #   Rank_Prop = 1/5 = 0.2.
  #   Max_FDR approx (0.01+0.2+0.2+0.9+1.0)/5 = 0.462
  #   p = (0.01 * 0.2) / 0.462 = 0.0043...
  expect_lt(res$fp[1], 0.01)

  # --- CHECK 3: Tie Breaking on Secondary P-value ---
  # Items 2 and 3 have identical lfdr (0.20).
  # Item 2 has p=0.01, Item 3 has p=0.05.
  # Item 2 must have a smaller functional p-value than Item 3.
  expect_lt(res$fp[2], res$fp[3])

  # --- CHECK 4: Noise Clamping ---
  # The last item (Pure Noise) should be close to or exactly 1.0
  expect_equal(tail(res$fp, 1), 1.0)
})
