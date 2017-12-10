# ComparativeMethodsInR
A book for teaching people how to do comparative methods in R

You will need to 

```
devtools::install_github("rstudio/bookdown")
devtools::install_github("bomeara/phybase")
```

To render it:

`bookdown::render_book("index.Rmd", "bookdown::pdf_book")`

Or for an html version,

`bookdown::render_book("index.Rmd", "bookdown::gitbook")`

And to deploy on bookdown,

`bookdown::publish_book(render = "local", account="bomeara")`

Or just use the RenderAll.R script.

The book.bib references are exported from Mendeley, from the PhyloMeth group.

use the `citr` package with the rstudio addin to cite easily
