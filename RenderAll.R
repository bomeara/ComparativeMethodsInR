bookdown::render_book("index.Rmd", "bookdown::pdf_book")
bookdown::render_book("index.Rmd", "bookdown::gitbook")
bookdown::publish_book(render = "local", account="bomeara")
