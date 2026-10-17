test_that(".as_dir joins parts, creates the directory, and adds one trailing slash", {
  root <- file.path(tempdir(), "test_as_dir")
  on.exit(unlink(root, recursive = TRUE))

  expect_equal(.as_dir(root), paste0(root, "/"))
  expect_equal(.as_dir(paste0(root, "/")), paste0(root, "/"))
  expect_equal(.as_dir(paste0(root, "//"), "a/", "/b"), paste0(root, "/a/b/"))
  expect_true(dir.exists(file.path(root, "a", "b")))
})
