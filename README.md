# RNA-seq lesson

This repository generates the RNA-seq lesson materials based on the website template from [The Carpentries](https://carpentries.org/) Foundation. 

### Preview changes to the lesson locally
The lesson website is built through Github and Jekyll. 

__Option 1:__ follow the Carpentries setup: http://carpentries.github.io/lesson-example/setup.html 
The detailed instructions are listed in the "Jekyll Setup for Lesson Development" section.   

__Option 2:__ use a Docker container
1. Open a Shell window. 
2. Navigate to the `rna-seq-lesson/` folder using the `cd` command.
3. Since the lesson relies Jekyll 3.8.5, type within the Shell `export JEKYLL_VERSION=3.8.5`.
4. Make sure you have Docker for Windows or Mac installed: https://docs.docker.com/install/
5. With the Docker Desktop application running (you should see a little whale with containers at the top of your screen), type `docker run --rm --volume="$PWD:/srv/jekyll" -p 4000:4000 -it jekyll/jekyll:$JEKYLL_VERSION jekyll serve`  
6. Open a web browser and type `http://0.0.0.0:4000/` in the navigation bar. You should see the lesson website. Your changes should be automatically reflected online.  

## Maintainer(s)

Current maintainers of this lesson are 

* Marc Galland, Data analyst and manager (University of Amsterdam, SILS, Plant Physiology Department).
* Tijs Bliek, research technician (University of Amsterdam, SILS, Plant Development and Epigenetics).

## Authors

A list of contributors to the lesson can be found in [AUTHORS](AUTHORS)

## Citation

To cite this lesson, please consult with [CITATION](CITATION)

## Credits
This lesson is heavily based on teaching materials from the [Harvard Chan Bioinformatics Core (HBC) in-depth NGS data analysis course](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/). Materials have been adapted and some exercises created to comply with the [Carpentries Foundation teaching requirements](https://carpentries.github.io/instructor-training/).

## Contributing

We welcome all contributions to improve the lesson! Maintainers will do their best to help you if you have any
questions, concerns, or experience any difficulties along the way.

We'd like to ask you to familiarize yourself with our [Contribution Guide](CONTRIBUTING.md) and have a look at
the [more detailed guidelines][lesson-example] on proper formatting, ways to render the lesson locally, and even
how to write new episodes.

Please see the current list of [issues][FIXME] for ideas for contributing to this
repository. For making your contribution, we use the GitHub flow, which is
nicely explained in the chapter [Contributing to a Project](http://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) in Pro Git by Scott Chacon.
Look for the tag ![good_first_issue](https://img.shields.io/badge/-good%20first%20issue-gold.svg). This indicates that the mantainers will welcome a pull request fixing this issue.  
