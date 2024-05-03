test_that("createSeuratObj works not correct", {
  seu <- createSeuratObj(system.file('extdata', 'test.qptiff.tsv.zip',
                                     package='cdQuPath'))
  expect_is(seu, 'Seurat')
})
