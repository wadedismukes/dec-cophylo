---
title: "Prototyping simulation study using DEC model for Cophylogeny"
output:
  html_document:
    df_print: paged
---

This notebook is a basic run through of the workflow that will be used to perform simulations testing the utilization of the epoch-DEC model for examining cophylogeny. Here the dispersal part of this model corresponds to host-switching, the extirpation part to extinction of a symbiont within a host, and the cladogenesis part to the various sorting events occurring when the host speciates (including cospeciation). The hosts here act as the biogeographic areas upon which the symbionts live. For now I am testing this using only true trees simulated using [`treeducken`](https://github.com/wadedismukes/treeducken). While these simulations will reflect a model of gene and spceies tree evolution and are thus not the most biologically realistic way to test this model, I think this is a good starting place especially with regards to comparing these models to other methods for examining cophylogeny. 


### Making `treeducken` settings files

Now we need to make our settings files for `treeducken`. To do this we take in some user inputs for the various parameters and then output that in the correct format.


```r
make_settings <- function(sbr, sdr, trr, gbr, gdr, nt, ipp, ne, reps, nl, ng){
    settings_list <- c(paste("-sbr",sbr),
                paste("-sdr",sdr),
                paste("-gbr", gbr),
                paste("-gdr", gdr),
                paste("-lgtr",trr),
                paste("-ipp",ipp),
                paste("-nt",nt),
                paste("-r",reps),
                paste("-ne",ne),
                paste("-nl",nl),
                paste("-ng",ng),
                paste("-o", paste0("hsrate_",trr)),
                paste("-sout", 0))
    return(settings_list)
}

settings_writer <- function(settings_list, file_begin){
    file_name_end <- "_settings.txt"
    out_fn <- paste0(file_begin, file_name_end)
    write(settings_list, file = out_fn, ncolumns = 1)
}
```

### Setting simulation parameters and printing files
Here I have set the model parameters to be printed into settings files at the top.
The for loop generates the number of regimes I want to simulate under. The number
of tips and timing (13 here corresponding to 13 million years ago (mya)) were 
chosen to mimic the Geomyidae gophers and their chewing lice parasites.



```r
nt <- 10
turnover <- 0.5
br <- 0.05
dr <- 0.025
trr <- as.vector(c(0.0, 0.025, 0.05, 0.1))
gbr <- 0.0
gdr <- 0.0
reps <- 10
nl <- 1
ng <- 1
ipp <- 1
ne <- 1
num_setting_regimes <- length(trr)

for (i in 1:num_setting_regimes){
    settings_regime_name <- paste0(trr[i],"_hsrate")
    settings <- make_settings(br, dr, trr[i], gbr, gdr, nt, 1, 1, reps, nl, ng)
    settings_writer(settings, settings_regime_name)
}
```

## Running simulations
Now that I have printed out my settings files, I can use `treeducken` to generate my trees. First, I have to install `treeducken` here.



```bash
if [ ! -d treeducken/ ]; then
    git clone https://github.com/wadedismukes/treeducken.git
    cd treeducken/src && make install
    cd .. && chmod +x treeducken
else
    echo 'already installed!'
fi
```

```
## already installed!
```


```bash
#make directory structure
mkdir -p settings/
mv *_settings.txt settings/
find ./settings/ -name "*_settings.txt" | xargs basename -s "_settings.txt" | xargs -I {} mkdir -p ./data/{}
```

Now run the simulations! 



```bash

for f in $(basename -s "_settings.txt" settings/*)
do
    prefix=$f
    suffix="_settings.txt"
    input_fn=settings/$prefix$suffix
    echo $input_fn
    treeducken/treeducken -i $input_fn
    mv *.tre ./data/$prefix/
done

```

```
## settings/0.025_hsrate_settings.txt
## ############################################################
## ####	treeducken, version 0.1 			####
## ####	e845a82c08ba308f75f94a270b81a36870299b94	####
## ############################################################
## Gene birth rate is 0.0, locus trees will match species trees.
## 		output file name prefix         = hsrate_0.025
## 		Number of extant taxa           = 10
## 		Number of replicates            = 10
## 		Number of loci to simulate      = 1
## 		Number of genes to simulate     = 1
## 		Species birth rate              = 0.05
## 		Species death rate              = 0.025
## 		Gene birth rate                 = 0
## 		Gene death rate                 = 0
## 		Gene transfer rate              = 0.025
## 		Individuals to sample per locus = 1
## 		Effective pop size per locus    = 1
## 		Tree fraction to set outgroup   = 0
## 		Species tree input as newick    = 
## 		Tree scale                      = -1
## 
## Seeds = {25167, 23727}
## settings/0.05_hsrate_settings.txt
## ############################################################
## ####	treeducken, version 0.1 			####
## ####	e845a82c08ba308f75f94a270b81a36870299b94	####
## ############################################################
## Gene birth rate is 0.0, locus trees will match species trees.
## 		output file name prefix         = hsrate_0.05
## 		Number of extant taxa           = 10
## 		Number of replicates            = 10
## 		Number of loci to simulate      = 1
## 		Number of genes to simulate     = 1
## 		Species birth rate              = 0.05
## 		Species death rate              = 0.025
## 		Gene birth rate                 = 0
## 		Gene death rate                 = 0
## 		Gene transfer rate              = 0.05
## 		Individuals to sample per locus = 1
## 		Effective pop size per locus    = 1
## 		Tree fraction to set outgroup   = 0
## 		Species tree input as newick    = 
## 		Tree scale                      = -1
## 
## Seeds = {25184, 23727}
## settings/0.1_hsrate_settings.txt
## ############################################################
## ####	treeducken, version 0.1 			####
## ####	e845a82c08ba308f75f94a270b81a36870299b94	####
## ############################################################
## Gene birth rate is 0.0, locus trees will match species trees.
## 		output file name prefix         = hsrate_0.1
## 		Number of extant taxa           = 10
## 		Number of replicates            = 10
## 		Number of loci to simulate      = 1
## 		Number of genes to simulate     = 1
## 		Species birth rate              = 0.05
## 		Species death rate              = 0.025
## 		Gene birth rate                 = 0
## 		Gene death rate                 = 0
## 		Gene transfer rate              = 0.1
## 		Individuals to sample per locus = 1
## 		Effective pop size per locus    = 1
## 		Tree fraction to set outgroup   = 0
## 		Species tree input as newick    = 
## 		Tree scale                      = -1
## 
## Seeds = {25199, 23727}
## settings/0_hsrate_settings.txt
## ############################################################
## ####	treeducken, version 0.1 			####
## ####	e845a82c08ba308f75f94a270b81a36870299b94	####
## ############################################################
## Gene birth rate is 0.0, locus trees will match species trees.
## 		output file name prefix         = hsrate_0
## 		Number of extant taxa           = 10
## 		Number of replicates            = 10
## 		Number of loci to simulate      = 1
## 		Number of genes to simulate     = 1
## 		Species birth rate              = 0.05
## 		Species death rate              = 0.025
## 		Gene birth rate                 = 0
## 		Gene death rate                 = 0
## 		Gene transfer rate              = 0
## 		Individuals to sample per locus = 1
## 		Effective pop size per locus    = 1
## 		Tree fraction to set outgroup   = 0
## 		Species tree input as newick    = 
## 		Tree scale                      = -1
## 
## Seeds = {25215, 23727}
```

### Creating extra datafiles for use with DEC model
With trees in hand, I need to make a few more files to run the epoch DEC model. 
I need to make a "range" file for the symbionts. To do this we need to parse the
tips of the locus tree which is our stand-in for the symbiont tree. First, we 
need a couple of packages.


```r
library(phytools)
```

```
## Loading required package: ape
```

```
## Error: package or namespace load failed for 'ape' in dyn.load(file, DLLpath = DLLpath, ...):
##  unable to load shared object '/home/waded/R/x86_64-pc-linux-gnu-library/3.5/ape/libs/ape.so':
##   libgfortran.so.4: cannot open shared object file: No such file or directory
```

```
## Error: package 'ape' could not be loaded
```

```r
library(stringr)
```

    also installing the dependencies ‘magick’, ‘igraph’, ‘fastmatch’, ‘ape’, ‘maps’, ‘animation’, ‘clusterGeneration’, ‘coda’, ‘combinat’, ‘expm’, ‘numDeriv’, ‘phangorn’, ‘plotrix’, ‘scatterplot3d’
    
    Warning message in install.packages("phytools", repos = "http://cran.us.r-project.org"):
    “installation of package ‘magick’ had non-zero exit status”Warning message in install.packages("phytools", repos = "http://cran.us.r-project.org"):
    “installation of package ‘animation’ had non-zero exit status”Warning message in install.packages("phytools", repos = "http://cran.us.r-project.org"):
    “installation of package ‘phytools’ had non-zero exit status”Updating HTML index of packages in '.Library'
    Making 'packages.html' ... done


First, I will make a NEXUS data file that takes the locus trees (i.e. our 
symbiont trees as used here) from `treeducken` and decribes which species those 
loci are associated with.


```r
# read in species trees to get the min and max ages of "epochs"
nexus_range_print <- function(fn, host_tree, symb_tree){
    sink(fn)
    
    cat("#NEXUS\n\n")
    cat("Begin data;\n")
    
    nranges <- length(getExtant(host_tree, tol=1e-3))
    nsymb <- length(getExtant(symb_tree, tol=1e-3))
    
    line <- paste0("Dimensions ntax =", nsymb, " nchar =", nranges, ";\n")
    
    cat(line)
    cat("Format datatype=Standard missing=? gap=- labels=\"01\";\n")
    cat("Matrix\n")
    
    hosts <- getExtant(host_tree, tol=1e-3)
    symbs <- getExtant(symb_tree, tol=1e-3)
    for(i in 1:length(symbs)){
        # first need to figure out which species associated with
        curr_tip <- str_split(symbs[i], "_")
        sp_indx <- str_locate(hosts, curr_tip[[1]][1])
        sp_indx <- which(!is.na(sp_indx[,1]))
        # next mkae the bit vector and put into the thing
        range_data <- rep(0, times = nranges)
        range_data[sp_indx[1]] <- 1
        range_data_str <- str_flatten(str_c(range_data))
        
        # now print the data to nexus file
        line <- paste0("\t", symbs[i], "\t", range_data_str, "\n")
        cat(line)
    }
    cat("\t;")
    cat("\nEnd;")
    sink()
}
```

We also need a matrix that describes which tips are able to be dispersed to
at which time points, this is done in this function. This will print out a 
number of files equal to the number of epochs that we are going to be dealing 
with. Each matrix has rows and columns equal to the number of hosts (which in
this case correspond to areas).


```r
connectivity_graph_print <- function(num_epochs, num_hosts, prefix_ofn){
    for(i in sort(1:(num_epochs), decreasing = TRUE)){
        conn_mat <- matrix(0, nrow = num_hosts, ncol = num_hosts)
        conn_mat[1:i, 1:i] <- 1
        ofn <- paste0(prefix_ofn, ".", i, ".txt")
        write.table(conn_mat, file = ofn, row.names=F, col.names=F)
    }
}
```

This next function simply take everything we've done and puts it together.


```r
run_sim_regime <- function(sim_dir, prefix_fn, num_reps){
    host_tree_suffix_fn <- ".sp.tre"
    symb_tree_suffix_fn <- "_0.loc.tre"
    
    brtimes_suffix_fn <- ".times.txt"
    ranges_suffix_fn <- ".range.nex"

    for (i in 0:(num_reps-1)){
        # read in the host and symbiont trees
        host_fn <- paste0(sim_dir, prefix_fn, i, host_tree_suffix_fn)
        symb_fn <- paste0(sim_dir, prefix_fn, i, symb_tree_suffix_fn)
        host_nwck_tr <- read.nexus(host_fn)
        symb_nwck_tr <- read.nexus(symb_fn)
        symb_nwck_tr$node.label <- NULL
        pruned_symb_tree_fn <- paste0(sim_dir, prefix_fn, i, "_pruned", symb_tree_suffix_fn)
        pruned_symb_tree <- drop.tip(symb_nwck_tr, getExtinct(symb_nwck_tr, tol = 1e-3))
        write.nexus(pruned_symb_tree, file = pruned_symb_tree_fn, translate = FALSE)

        # print out epoch times (no uncertainty in age here)
        out_fn <- paste0(sim_dir, prefix_fn, i)
        brtimes_out_fn <- paste0(out_fn, brtimes_suffix_fn)
        br_times <- sort(c(0.0, branching.times(host_nwck_tr)), decreasing = TRUE)
        br_times_mat <- matrix(nrow = length(br_times), ncol = 2)
        br_times_mat[,1] <- br_times
        br_times_mat[,2] <- br_times - (br_times*0.05)
        write.table(br_times_mat, file = brtimes_out_fn, row.names=F, col.names=F)
        
        # print out range files
        nex_range_outfn <- paste0(out_fn, ranges_suffix_fn)
        nranges <- length(getExtant(host_nwck_tr, tol = 1e-3))
        nexus_range_print(nex_range_outfn, host_nwck_tr, symb_nwck_tr)
        
        # print out connectivity graphs
        connectivity_prefix_outfn <- paste0(out_fn, ".connectivity")
        
        connectivity_graph_print(num_epochs = length(br_times),
                                 num_hosts = nranges,
                                 connectivity_prefix_outfn)
        
        # print out "distance" matrix
        dist_mat <- cophenetic(host_nwck_tr)
        distmat_ofn <- paste0(out_fn, ".distances.txt")
        write.table(dist_mat, file = distmat_ofn, row.names=F, col.names=F)

    }
}
```

Output rev scripts for each replicate

```r
print_rev_script <- function(dir_name, prefix_fn, reps){
    symb_tree_suffix_fn <- "_0.loc.tre"
    rev_script_dir <- "rev-scripts/"
    brtimes_suffix_fn <- ".times.txt"
    ranges_suffix_fn <- ".range.nex"
    
    for(i in 0:(reps-1)){
	print("poop")	
        rev_out_fn <- paste0(rev_script_dir, "run_epoch_", prefix_fn, i, ".Rev")
        system2("touch", args = rev_out_fn)
        sink(rev_out_fn)
        
        
        range_fn <- paste0(dir_name, prefix_fn, i, ranges_suffix_fn)
        tree_fn <- paste0(dir_name, prefix_fn, i, "_pruned", symb_tree_suffix_fn)
        out_fn <- paste0(dir_name, "output/", prefix_fn, i)
        geo_fn <- paste0(dir_name, prefix_fn, i)

        cat(paste0("range_fn = \"", range_fn, "\"\n"))
        cat(paste0("tree_fn = \"", tree_fn, "\"\n"))
        cat(paste0("out_fn = \"", out_fn, "\"\n"))
        cat(paste0("geo_fn = \"", geo_fn, "\"\n"))
        cat("times_fn = geo_fn + \".times.txt\"\n")
        cat("dist_fn  = geo_fn + \".distances.txt\"\n")

        cat("moves = VectorMoves()\n")
        cat("monitors = VectorMonitors()\n")
        cat("n_gen = 5000\n")

        cat("dat_range_01 = readDiscreteCharacterData(range_fn)\n")
        cat("n_areas <- dat_range_01.nchar()\n")
        
        cat("max_areas <- 2\n")
        cat("n_states <- 0\n")
        cat("for (k in 0:max_areas) n_states += choose(n_areas, k)\n")
        
        cat("dat_range_n = formatDiscreteCharacterData(dat_range_01, \"DEC\", n_states)\n")
        cat("state_desc = dat_range_n.getStateDescriptions()\n")
        cat("state_desc_str = \"state,range\"\n")
        cat("for (i in 1:state_desc.size())\n")
        cat("{\n")
        cat("\tstate_desc_str += (i-1) + \",\" + state_desc[i] + \"\\n\"\n")
        cat("}\n")
        
        cat("write(state_desc_str, file=out_fn+\".state_labels.txt\")\n")
        
        cat("time_bounds <- readDataDelimitedFile(file=times_fn, delimiter=\" \")\n")
        cat("n_epochs <- time_bounds.nrows()\n")
        
        cat("for (i in 1:n_epochs) {\n")
        cat("\tepoch_fn = geo_fn + \".connectivity.\" + i + \".txt\"\n")
        cat("\tconnectivity[i] <- readDataDelimitedFile(file=epoch_fn, delimiter=\" \")\n")
        cat("}\n")
        
        cat("distances <- readDataDelimitedFile(file=dist_fn, delimiter=\" \")\n")
        cat("tree <- readTrees(tree_fn)[1]\n")
        
        #cat("log10_rate_bg ~ dnUniform(-4,2)\n")
        #cat("log10_rate_bg.setValue(-2)\n")
        #cat("rate_bg := 10^log10_rate_bg\n")
        #cat("moves.append( mvSlide(log10_rate_bg, weight=4) )\n")

	cat("rate_bg <- 1.0")

        cat("log_sd <- 0.5\n")
        cat("log_mean <- ln(1) - 0.5*log_sd^2\n")
        cat("dispersal_rate ~ dnLognormal(mean=log_mean, sd=log_sd)\n")
        cat("moves.append( mvScale(extirpation_rate, weight=5) )\n")

        cat("for (i in 1:n_epochs) {\n")
        cat("\tfor (j in 1:n_areas) {\n")
        cat("\t\tfor (k in 1:n_areas) {\n")
        cat("\t\t\tdr[i][j][k] <- 0.0\n")
        cat("\t\t\tif (connectivity[i][j][k] > 0) {\n")
        cat("\t\t\t\tdr[i][j][k] := dispersal_rate\n")
        cat("\t\t\t}\n")
        cat("\t\t}\n")
        cat("\t}\n")
        cat("}\n")
        
        cat("log_sd <- 0.5\n")
        cat("log_mean <- ln(1) - 0.5*log_sd^2\n")
        cat("extirpation_rate ~ dnLognormal(mean=log_mean, sd=log_sd)\n")
        cat("moves.append( mvScale(extirpation_rate, weight=5) )\n")

        cat("for (i in 1:n_epochs) {\n")
        cat("\tfor (j in 1:n_areas) {\n")
        cat("\t\tfor (k in 1:n_areas) {\n")
        cat("\t\t\ter[i][j][k] <- 0.0\n")
        cat("\t\t}\n")
        cat("\t\ter[i][j][j] := extirpation_rate\n")
        cat("\t}\n")
        cat("}\n")
        
        
        cat("for (i in n_epochs:1) {\n")
        cat("Q_DEC[i] := fnDECRateMatrix(dispersalRates=dr[i],
                                        extirpationRates=er[i],
                                        maxRangeSize=max_areas)\n")
        cat("}\n")

        cat("for (i in 1:n_epochs) {\n")
        cat("\ttime_max[i] <- time_bounds[i][1]\n")
        cat("\ttime_min[i] <- time_bounds[i][2]\n")
        cat("\tif (i != n_epochs) {\n")
        cat("\t\tepoch_times[i] ~ dnUniform(time_min[i], time_max[i])\n")
        cat("\t\tmoves.append( mvSlide(epoch_times[i], delta=(time_max[i]-time_min[i])/2) )\n")
        cat("\t} else {\n")
        cat("\t\tepoch_times[i] <- 0.0\n")
        cat("\t}\n")
        cat("}\n")        
        
        
        cat("Q_DEC_epoch := fnEpoch(Q=Q_DEC, times=epoch_times, rates=rep(1, n_epochs))\n")
        
        cat("clado_event_types <- [ \"s\", \"a\" ]\n")
        cat("p_sympatry ~ dnUniform(0,1)\n")
        cat("p_allopatry := abs(1.0 - p_sympatry)\n")
        cat("clado_type_probs := simplex(p_sympatry, p_allopatry)\n")
        cat("moves.append( mvSlide(p_sympatry, weight=2) )\n")
        cat("P_DEC := fnDECCladoProbs(eventProbs=clado_type_probs,
                                                eventTypes=clado_event_types,
                                                numCharacters=n_areas,
                                                maxRangeSize=max_areas)\n")
        
        cat("rf_DEC <- rep(0, n_states)\n")
        cat("rf_DEC[2] <- 1\n")
        cat("rf_DEC_simp <- simplex(rf_DEC)\n")
        
        cat("m_bg ~ dnPhyloCTMCClado(tree=tree,
                           Q=Q_DEC_epoch,
                           cladoProbs=P_DEC,
                           branchRates=rate_bg,
                           rootFrequencies=rf_DEC_simp,
                           type=\"NaturalNumbers\",
                           nSites=1)\n")
        cat("m_bg.clamp(dat_range_n)\n")
        
        
        cat("monitors.append( mnScreen(printgen=100, rate_bg, extirpation_rate) )\n")
        cat("monitors.append( mnModel(file=out_fn+\".model.log\", printgen=10) )\n")
        cat("monitors.append( mnFile(tree, filename=out_fn+\".tre\", printgen=10) )\n")
        cat("monitors.append( mnJointConditionalAncestralState(tree=tree,
                                                       ctmc=m_bg,
                                                       type=\"NaturalNumbers\",
                                                       withTips=true,
                                                       withStartStates=true,
                                                       filename=out_fn+\".states.log\",
                                                        printgen=10) )\n")
        cat("mymodel = model(m_bg)\n")
        cat("mymcmc = mcmc(mymodel, monitors, moves)\n")
        cat("mymcmc.run(n_gen)\n")
        
        
        sink()
    }
}
```


Now we can run this for each simulation regime. 


```r
dir_name <- c("data/0_hsrate/", "data/0.1_hsrate/", "data/0.05_hsrate/", "data/0.025_hsrate/")
prefix_fn <- c("hsrate_0_", "hsrate_0.1_", "hsrate_0.05_", "hsrate_0.025_")

for(i in 1:length(dir_name)){
    run_sim_regime(dir_name[i], prefix_fn[i], reps)
    print_rev_script(dir_name[i], prefix_fn[i], reps)
}
```

```
## Error in read.nexus(host_fn): could not find function "read.nexus"
```

Now we just need to write up our Rev scripts and run them.


## Rev file bit

Also a test of Rev code chunks in Rmarkdown/knitr.



```bash
mkdir -p tdcken_stats/
mv *.stats.txt tdcken_stats/
```
