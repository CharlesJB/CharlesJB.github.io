---
layout: post
published: true
tags: R, cookbook, GenomicAlignments, TxDb
status: publish
title: "How to extract promoters positions"
author: "Charles Joly Beauparlant"
output: html_document
---
 
# Introduction
 
In this post, I will show how easy it is to extract the genomic positions of every promoters of a specific genome build.
 
For this demo, you will need the `TxDb.Hsapiens.UCSC.hg19.knownGene` package:

{% highlight r %}
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
# To avoid have to type the whole package name every time, we use the variable name txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
{% endhighlight %}
 
# TL;DR
 

{% highlight r %}
promoters(genes(txdb), upstream = 1500, downstream = 500)
{% endhighlight %}



{% highlight text %}
## GRanges object with 23056 ranges and 1 metadata column:
##         seqnames                 ranges strand   |     gene_id
##            <Rle>              <IRanges>  <Rle>   | <character>
##       1    chr19 [ 58873715,  58875714]      -   |           1
##      10     chr8 [ 18247255,  18249254]      +   |          10
##     100    chr20 [ 43279877,  43281876]      -   |         100
##    1000    chr18 [ 25756946,  25758945]      -   |        1000
##   10000     chr1 [244006387, 244008386]      -   |       10000
##     ...      ...                    ...    ... ...         ...
##    9991     chr9 [115095445, 115097444]      -   |        9991
##    9992    chr21 [ 35734823,  35736822]      +   |        9992
##    9993    chr22 [ 19109468,  19111467]      -   |        9993
##    9994     chr6 [ 90538119,  90540118]      +   |        9994
##    9997    chr22 [ 50964406,  50966405]      -   |        9997
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
{% endhighlight %}
 
# `TxDb`
 
The trick is to use a type of packages that are known as [`TxDb`](http://bioconductor.org/packages/release/BiocViews.html#___TxDb). `TxDb` stands for Transcripts Database, and as the name implies, it contains information about the transcripts for a specific genome build of a given specie.
 
If you are lucky, a `TxDb` package is already available on Bioconductor and extracting the promoter information will be very straighforward. Otherwise, it is possible to create your own `TxDb` object, but this is beyond the scope of the current post.
 
The goal of this document is not to describe in details the inner workings of `TxDb` objects. We will only show two helper functions that allow to easily extract relevant information from this type of object.
 

{% highlight r %}
txdb
{% endhighlight %}



{% highlight text %}
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: UCSC
## # Genome: hg19
## # Organism: Homo sapiens
## # UCSC Table: knownGene
## # Resource URL: http://genome.ucsc.edu/
## # Type of Gene ID: Entrez Gene ID
## # Full dataset: yes
## # miRBase build ID: GRCh37
## # transcript_nrow: 82960
## # exon_nrow: 289969
## # cds_nrow: 237533
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2015-03-19 13:55:51 -0700 (Thu, 19 Mar 2015)
## # GenomicFeatures version at creation time: 1.19.32
## # RSQLite version at creation time: 1.0.0
## # DBSCHEMAVERSION: 1.1
{% endhighlight %}
 
# `promoters` and `genes`
 
The `promoters` function can be used to extract the information for the promoters of every transcripts from a `TxDb` object:
 

{% highlight r %}
promoters_txdb <- promoters(txdb)
promoters_txdb
{% endhighlight %}



{% highlight text %}
## GRanges object with 82960 ranges and 2 metadata columns:
##           seqnames               ranges strand   |     tx_id     tx_name
##              <Rle>            <IRanges>  <Rle>   | <integer> <character>
##       [1]     chr1     [  9874,  12073]      +   |         1  uc001aaa.3
##       [2]     chr1     [  9874,  12073]      +   |         2  uc010nxq.1
##       [3]     chr1     [  9874,  12073]      +   |         3  uc010nxr.1
##       [4]     chr1     [ 67091,  69290]      +   |         4  uc001aal.1
##       [5]     chr1     [319084, 321283]      +   |         5  uc001aaq.2
##       ...      ...                  ...    ... ...       ...         ...
##   [82956]     chrY [27605479, 27607678]      -   |     78803  uc004fwx.1
##   [82957]     chrY [27606222, 27608421]      -   |     78804  uc022cpc.1
##   [82958]     chrY [27607233, 27609432]      -   |     78805  uc004fwz.3
##   [82959]     chrY [27635755, 27637954]      -   |     78806  uc022cpd.1
##   [82960]     chrY [59360655, 59362854]      -   |     78807  uc011ncc.1
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
{% endhighlight %}
 
The promoters function returns a `GRanges` object corresponding to the positions of the promoters of every transcripts in the `Txdb` object.
 
This returns 82960 promoter regions. But often we are only interested in the promoters of the genes and not of all the transcripts. This is where the `genes` function becomes handy:
 

{% highlight r %}
genes_txdb <- genes(txdb)
promoters_txdb <- promoters(genes_txdb)
promoters_txdb
{% endhighlight %}



{% highlight text %}
## GRanges object with 23056 ranges and 1 metadata column:
##         seqnames                 ranges strand   |     gene_id
##            <Rle>              <IRanges>  <Rle>   | <character>
##       1    chr19 [ 58874015,  58876214]      -   |           1
##      10     chr8 [ 18246755,  18248954]      +   |          10
##     100    chr20 [ 43280177,  43282376]      -   |         100
##    1000    chr18 [ 25757246,  25759445]      -   |        1000
##   10000     chr1 [244006687, 244008886]      -   |       10000
##     ...      ...                    ...    ... ...         ...
##    9991     chr9 [115095745, 115097944]      -   |        9991
##    9992    chr21 [ 35734323,  35736522]      +   |        9992
##    9993    chr22 [ 19109768,  19111967]      -   |        9993
##    9994     chr6 [ 90537619,  90539818]      +   |        9994
##    9997    chr22 [ 50964706,  50966905]      -   |        9997
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
{% endhighlight %}
 
Both function can also be nested to avoid the intermediate `genes_txdb` object:
 

{% highlight r %}
promoters_txdb <- promoters(genes(txdb))
promoters_txdb
{% endhighlight %}



{% highlight text %}
## GRanges object with 23056 ranges and 1 metadata column:
##         seqnames                 ranges strand   |     gene_id
##            <Rle>              <IRanges>  <Rle>   | <character>
##       1    chr19 [ 58874015,  58876214]      -   |           1
##      10     chr8 [ 18246755,  18248954]      +   |          10
##     100    chr20 [ 43280177,  43282376]      -   |         100
##    1000    chr18 [ 25757246,  25759445]      -   |        1000
##   10000     chr1 [244006687, 244008886]      -   |       10000
##     ...      ...                    ...    ... ...         ...
##    9991     chr9 [115095745, 115097944]      -   |        9991
##    9992    chr21 [ 35734323,  35736522]      +   |        9992
##    9993    chr22 [ 19109768,  19111967]      -   |        9993
##    9994     chr6 [ 90537619,  90539818]      +   |        9994
##    9997    chr22 [ 50964706,  50966905]      -   |        9997
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
{% endhighlight %}
 
By default, the `promoters` function will fetch the 2000 nucleotides before the transcription start site (TSS) and the 200 nucleotides after the TSS. This can be controled with the `upstream` and `downstream` parameters:
 

{% highlight r %}
unique(width(promoters_txdb))
{% endhighlight %}



{% highlight text %}
## [1] 2200
{% endhighlight %}



{% highlight r %}
promoters_txdb <- promoters(genes(txdb), upstream = 1500, downstream = 500)
unique(width(promoters_txdb))
{% endhighlight %}



{% highlight text %}
## [1] 2000
{% endhighlight %}
