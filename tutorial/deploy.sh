#!/bin/sh
# Render the tutorial and publish it to the gh-pages branch
# (served at https://yiqunchen.github.io/selectiveHotelling/).
# Run from the package root after installing the package: R CMD INSTALL .
set -e
cd "$(dirname "$0")"
Rscript -e 'bookdown::render_book("index.Rmd", quiet = TRUE)'
tmp=$(mktemp -d)
cp -R _book/. "$tmp"/ && touch "$tmp/.nojekyll"
cd "$tmp" && git init -q -b gh-pages && git add -A && git commit -q -m "Deploy tutorial"
git push -f https://github.com/yiqunchen/selectiveHotelling.git gh-pages
