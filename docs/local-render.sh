#!/bin/sh -x
## Do a Jekyll render of the pages in this directory to see how it will render
## Be sure to run knit-it.sh if any of the Rmd files have changed.


## go where this script lives:
cd $(dirname $(readlink -e $BASH_SOURCE))

cd docs/
rm -fr _site ## clear the cache; usually not needed as jekyll is smart about it

port=4343

echo "Now (re)open 

   http://127.0.0.1:$port/index.html

(NOTE: public version, no figures)

"

set -x
bundle exec jekyll serve --config local_config.yml --incremental --port $port --trace 
