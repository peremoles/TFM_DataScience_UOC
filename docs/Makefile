## from https://groups.wfu.edu/ir/tools/data-scientist-handbook/makefiles.html

SRC=.
DOCS=./docs

## fetch dependencies
DOCS_RFILES := $(wildcard $(DOCS)/*.html)

DOCS_OUT_FILES := $(DOCS_RFILES:.Rmd=.html)

all : $(DOCS_OUT_FILES)
	@echo "$(DOCS_OUT_FILES)"

$(DOCS)/%.html : $(SRC)/%.Rmd _site.yml
	@echo compiling R Markdown
	Rscript -e 'rmarkdown::render_site("$<")' || exit
