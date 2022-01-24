<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Data](#data)
- [Overview (needs update)](#overview-needs-update)
- [Development](#development)
  - [Local rendering](#local-rendering)
  - [github pages](#github-pages)
- [Structure](#structure)
  - [Configuration](#configuration)
  - [Content](#content)
    - [Header and footer:](#header-and-footer)
  - [Rmarkdown pages](#rmarkdown-pages)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Data

**NOTE: the Overview etc. below is out of date, I will fix this soon**
**This current bit is just about the data that the participants need,
and on which day1.Rmd depends**

The 'lessons' can be found below the `docs/` directory.

The are in rmarkdown format (`.Rmd`) and automatically 'knitted' to
markdown (`.md`) by the `knit-it.sh` script.

The lessons for 'day1' (1-loadingdata .. 8-diffexp) are structured such 
that they build on the previous lessons, but can knitted independently
(they read the previous session and save their own session afterwards). 
Knitting a lesson takes ~ 2 minutes. Total of all `session.rda`s is around
4 GB.

All data is put below `CFG$data_dir`. CFG + contents is a `list` defined
high up in the `day1.R` file. If all is well, the complete `day1.R` file
should work if just `CFG$data_dir` is changed to what it needs to be.

# Overview (needs update)

This repo contains stuff for our course(s). This README.md is aimed at
the tutors of the course. 

**Note to contributers: this website is visible to the entire world. Don't
put any sensitive data here, be it patient data or locations, let alone passwords**

The course practical information, schedule and materials are maintained
in this repository, which is automatically converted to a a website 
using github's `gh-pages` mechanism, explained below.

The README.md's in this directory and below are only there for tutors,
they are not part of the course contents. For convenience when reading
these file in github, they have been given a table of contents using
[doctoc](https://github.com/thlorenz/doctoc). If there are big changes
to any README.md file, rerun doctoc as

```
$    find . -name README.md -exec doctoc --github {} \;
```

(and commit and push).

# Development

The content of the course sits below `docs/`, wich `course.md`
and `schedule.md` being the 'entry points and [0-9]* the 'lessons'.

Content is managed using plain markdown and R markdown (to be rendered
by `knit`) and Jekyll, which is a templating and rendering engine that
produces static `*.html` pages.

Throughout the code I use `@`'s (usually at least 2) to indicate where
things aren't quite ready yet. Often but not always they are in #-comments 
or in \<  \!  -- HTML comments \-- \> so you don't seem them unless editing.

Plain *.md pages are rendered directly by Jekyl (see below). Rmarkdown
pages need to be knitted first usisng the knit-it.sh script, then
rendered by Jekyll, either local or by github.

## Local rendering

To render things locally run the `local-render.sh` script, which does something like

```
cd docs/
rm -fr _site ## clear the cache; usually not needed as jekyll is smart about it
bundle exec jekyll serve --incremental  --port 4321 
```

(note: this is on Mac, and requires `gem`, `yekyll` and `pandoc`).

and then subsequently open http://127.0.0.1:4321/course.html

(without `--port` it defaults to port 4000). 

All valid [^1] `*.md` files under `docs/` are rendered to `*html`, and the 
all hyperlinks to `*.md` files are also translated to the right `.html`
URLs. All this is written to the `_site`  directory.


[^1] Valid markdown pages are those that start with 
'^---\nlayout: default\n---\n'

## github pages

Once you're satisfied, git commit and push it to github. Github will
then render it automatically and then put it up at

https://scgenomics.github.io/course.html

To check if and when it was last rendered succesfully, go to 
https://github.com/scgenomics/scgenomics.github.io/deployments

Rendering may look slightly differently locally and on github.io.

# Site snapshots

(this used to be in schedule.md)

are copied to  SurfDrive `sc-facility/course/site-snapshots`
and are named by timestamp. You should be able to just open e.g. 
`site-2022-01-18_17_04_16/index.html` and see everything.
The link sc-facility/course/site-snapshots/site-LATEST/schedule.html
should automatically take you to the latest and greatest

uploading is done using the site-uploade.sh script.

------------------------------------------------------------------------



# Structure

## Configuration

The main 'recipe' for building the site is `docs/_config.yml`, which sets
title, description, theme, some stuff needed by Jekyll and lastly
file(s) to be rendered (in our case `course.md`).

It specifies the title, layout, and toplevel (here: `course.md`). 

## Content

`course.md` is the organisational / practical homepage of the course.

The content of the course itself reached through `schedule.md`.

Non-R content such as the main page, schedule page etc. are simply rendered by Jekyll. They need a  header and footer:

### Header and footer:

Make sure that it starts with 

```
---
layout: default
title: The Actual Title
---
```

otherwise Jekyll won't recognize it and won't render it to html (it is still served as 'raw' `*.md` though).

Please include the following line at the end of each page:

```
### [Back to main](course.md)
```

or 

```
### [Back to main](../course.md)
```

or similar.

## Rmarkdown pages

Lessons are specfied in Rmarkdown (`*.Rmd`), and the live in in
`docs/NUMBER-THETOPIC/NUMBER-THETOPIC.Rmd` A lesson template file is
given by `999-lesson.Rmd.template`. It contains the necessary header and
footer. Alter the 'title:' and specify the content.  The `knit-it.sh`
script renders the Rmarkdown to plain markdown.

All lessons (apart from 0-preamble) have a symlink to the 0-preamble.R file
(which is purl'ed from `0-preamble.Rmd` by `knit-it.sh`. This is needed
to load the libraries and set globals. To fake having one session for 
all the different knitting, make sure that where applicable a lesson
reads the state of the previous lesson, and save it for the next one.

When called without arguments, `knit-it.sh` knits all the [0-9]*.Rmd
files into markdown. When called with an argument, it will only 
knit that directory.

The header and footer material in each lesson (as based on the template)
take care of loading the state upon start and saving the sate upon finish.
This way all the lessons can be knitted as if they were part of
one big session. Make sure the saving and subsequent loading of state
(to file `session.rda`) in the next lesson is done properply.

