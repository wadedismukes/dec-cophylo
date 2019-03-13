# dec-cophylo
Repo for working on using DEC-epoch model for cophylogeny analyses. Currently using jupyter to make workflow.

To use docker to get the ipython notebook type the following:

```
docker run --rm -p 10000:8888 -v "$PWD":/home/jovyan/work jupyter/scipy-notebook:e5c5a7d3e52d
```

Then go to [http://localhost:10000/edit/work](http://localhost:10000/edit/work/README.md), making sure it is going to port 10000 and not 8888 (I have tested this only using Docker for Mac OSX so proceed with caution!

This jupyter notebook will walk you through the workflow for the simulations (as I have done them thus far). 