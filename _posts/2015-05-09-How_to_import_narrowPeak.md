---
layout: post
published: true
draft: true
tags: R, cookbook, rtracklayer, GenomicRanges 
status: publish
title: "How to import narrowPeak files"
author: "Charles Joly Beauparlant"
output: html_document
---
 
The goal of this post is to show how to import [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) and [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) files into R in a valid `GRanges` format.
 
To do so, you need to have installed the [`rtracklayer`](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package. To replicate the examples in this document, you will also need to install the [`GenomicFormatExamples`](https://github.com/CharlesJB/GenomicFormatExamples) package.
 

{% highlight r %}
require(rtracklayer)
require(GenomicFormatExamples)
{% endhighlight %}
 
# TL;DR
 

{% highlight r %}
# To import narrowPeak files
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
gr_narrowPeak <- import(file_narrowPeak, format = "BED",
                        extraCols = extraCols_narrowPeak)
 
# To import broadPeak files
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric")
gr_broadPeak <- import.bed(file_broadPeak, format = "BED",
                           extraCols = extraCols_narrowPeak)
{% endhighlight %}
 
# The `rtracklayer` package
 
The `rtracklayer` package offers multiple ways to easily import various genomic formats such as BED, WIG or GFF/GTF. For instance, if we want to import a BED file we can use the `import` function:
 

{% highlight r %}
bed_file <- get_demo_file(format = "bed")
gr_bed <- import(bed_file)
gr_bed
{% endhighlight %}



{% highlight text %}
## GRanges object with 1000 ranges and 2 metadata columns:
##          seqnames                 ranges strand   |        name     score
##             <Rle>              <IRanges>  <Rle>   | <character> <numeric>
##      [1]     chr1 [ 62902074,  62903462]      *   |      Rank_1      1859
##      [2]     chr1 [167188356, 167189834]      *   |      Rank_2      1720
##      [3]     chr3 [ 52321760,  52322314]      *   |      Rank_3      1632
##      [4]    chr12 [  2112636,   2113682]      *   |      Rank_4      1606
##      [5]    chr17 [ 37843582,  37845970]      *   |      Rank_5      1597
##      ...      ...                    ...    ... ...         ...       ...
##    [996]    chr12 [111180217, 111180675]      *   |    Rank_996       770
##    [997]    chr12 [ 79941170,  79942669]      *   |    Rank_997       770
##    [998]     chr1 [ 41445363,  41445648]      *   |    Rank_998       770
##    [999]    chr17 [ 55161909,  55162356]      *   |    Rank_999       770
##   [1000]     chr3 [181425978, 181427339]      *   |   Rank_1000       770
##   -------
##   seqinfo: 23 sequences from an unspecified genome; no seqlengths
{% endhighlight %}
 
The `import` function can be used to import the following file formats:
 
* GFF
* BED
* Bed15
* bedGraph
* WIG
* BigWig
 
As shown in the previous example, the file format is derived from the file extension which is why it generally works correctly without have to specify the format.
 
# Importing narrowPeak and broadPeak
 
Unfortunately, the [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) and [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) are not directly supported by the `import` function:
 

{% highlight r %}
narrowPeak_file <- get_demo_file(format = "narrowPeak")
import(narrowPeak_file)
{% endhighlight %}



{% highlight text %}
## Error in import(FileForFormat(con), ...): error in evaluating the argument 'con' in selecting a method for function 'import': Error in FileForFormat(con) : Format 'narrowPeak' unsupported
{% endhighlight %}



{% highlight r %}
broadPeak_file <- get_demo_file(format = "broadPeak")
import(broadPeak_file)
{% endhighlight %}



{% highlight text %}
## Error in import(FileForFormat(con), ...): error in evaluating the argument 'con' in selecting a method for function 'import': Error in FileForFormat(con) : Format 'broadPeak' unsupported
{% endhighlight %}
 
Even if we specify the format to be BED, the `import` function will fail:
 

{% highlight r %}
import(narrowPeak_file, format = "BED")
{% endhighlight %}



{% highlight text %}
## Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings, : scan() expected 'an integer', got '49.55787'
{% endhighlight %}



{% highlight r %}
import(broadPeak_file, format = "BED")
{% endhighlight %}



{% highlight text %}
## Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings, : scan() expected 'an integer', got '49.55787'
{% endhighlight %}
 
The reason is that the `import` function checks the format of the content of every columns to make sure the file is in the good format and that the columns in BED files are not completely identical to those in narrowPeak/broadPeak files.
 
BED files:
 
1. chrom 
2. chromStart
3. chromEnd
4. name
5. score
6. strand
7. thickStart
8. thickEnd
9. itemRgb
10. blockCount
11. blockSizes
12. blockStarts
 
narrowPeak/broadPeak files:
 
1. chrom
2. chromStart
3. chromEnd
4. name
5. score
6. strand
7. signalValue
8. pValue
9. qValue
10. peak (narrowPeak only)
 
The first 6 columns are the same, but the seventh column is different. In the BED format, the `thickStart` is an integer while in the narrowPeak/broadPeak format the `signalValue` is a numeric.
 
In order to solve this problem, we need to use the `extraCols` parameter:
 
```
extraCols: A character vector in the same form as ‘colClasses’ from
	‘read.table’.  It should indicate the name and class of each
	extra/special column to read from the BED file. As BED does
	not encode column names, these are assumed to be the last
	columns in the file. This enables parsing of the various
	BEDX+Y formats.
```
 
In other words, we need to give the name and type of every columns starting at the one that is different from the standard BED format. In our case, we need to give the name and type of the 7th, 8th, 9th and 10th (in the case of narrowPeak) columns. To to so, we have to create a named vector:
 

{% highlight r %}
extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
gr_narrowPeak <- import(narrowPeak_file, format = "BED",
                        extraCols = extraCols_narrowPeak)
gr_narrowPeak
{% endhighlight %}



{% highlight text %}
## GRanges object with 1000 ranges and 6 metadata columns:
##          seqnames                 ranges strand   |        name     score
##             <Rle>              <IRanges>  <Rle>   | <character> <numeric>
##      [1]     chr1 [ 62902074,  62903462]      *   |      Rank_1      1859
##      [2]     chr1 [167188356, 167189834]      *   |      Rank_2      1720
##      [3]     chr3 [ 52321760,  52322314]      *   |      Rank_3      1632
##      [4]    chr12 [  2112636,   2113682]      *   |      Rank_4      1606
##      [5]    chr17 [ 37843582,  37845970]      *   |      Rank_5      1597
##      ...      ...                    ...    ... ...         ...       ...
##    [996]    chr12 [111180217, 111180675]      *   |    Rank_996       770
##    [997]    chr12 [ 79941170,  79942669]      *   |    Rank_997       770
##    [998]     chr1 [ 41445363,  41445648]      *   |    Rank_998       770
##    [999]    chr17 [ 55161909,  55162356]      *   |    Rank_999       770
##   [1000]     chr3 [181425978, 181427339]      *   |   Rank_1000       770
##          singnalValue    pValue    qValue      peak
##             <numeric> <numeric> <numeric> <integer>
##      [1]     49.55787  185.9134  176.4229       497
##      [2]     40.23876  172.0305  163.9549      1285
##      [3]     38.65456  163.2142  155.6267       178
##      [4]     43.39093  160.6900  153.2080       857
##      [5]     46.80020  159.7838  152.3103       636
##      ...          ...       ...       ...       ...
##    [996]     24.79482  77.04454  72.39902       350
##    [997]     24.79482  77.04454  72.39902       470
##    [998]     24.79482  77.04454  72.39902        97
##    [999]     24.79482  77.04454  72.39902       305
##   [1000]     24.79482  77.04454  72.39902       370
##   -------
##   seqinfo: 23 sequences from an unspecified genome; no seqlengths
{% endhighlight %}
 

{% highlight r %}
extraCols_broadPeak <- c(singnalValue = "numeric", pValue = "numeric",
                         qValue = "numeric", peak = "integer")
gr_broadPeak <- import(broadPeak_file, format = "BED",
                       extraCols = extraCols_broadPeak)
gr_broadPeak
{% endhighlight %}



{% highlight text %}
## GRanges object with 1000 ranges and 6 metadata columns:
##          seqnames                 ranges strand   |        name     score
##             <Rle>              <IRanges>  <Rle>   | <character> <numeric>
##      [1]     chr1 [ 62902074,  62903462]      *   |      Rank_1      1859
##      [2]     chr1 [167188356, 167189834]      *   |      Rank_2      1720
##      [3]     chr3 [ 52321760,  52322314]      *   |      Rank_3      1632
##      [4]    chr12 [  2112636,   2113682]      *   |      Rank_4      1606
##      [5]    chr17 [ 37843582,  37845970]      *   |      Rank_5      1597
##      ...      ...                    ...    ... ...         ...       ...
##    [996]    chr12 [111180217, 111180675]      *   |    Rank_996       770
##    [997]    chr12 [ 79941170,  79942669]      *   |    Rank_997       770
##    [998]     chr1 [ 41445363,  41445648]      *   |    Rank_998       770
##    [999]    chr17 [ 55161909,  55162356]      *   |    Rank_999       770
##   [1000]     chr3 [181425978, 181427339]      *   |   Rank_1000       770
##          singnalValue    pValue    qValue      peak
##             <numeric> <numeric> <numeric> <integer>
##      [1]     49.55787  185.9134  176.4229       497
##      [2]     40.23876  172.0305  163.9549      1285
##      [3]     38.65456  163.2142  155.6267       178
##      [4]     43.39093  160.6900  153.2080       857
##      [5]     46.80020  159.7838  152.3103       636
##      ...          ...       ...       ...       ...
##    [996]     24.79482  77.04454  72.39902       350
##    [997]     24.79482  77.04454  72.39902       470
##    [998]     24.79482  77.04454  72.39902        97
##    [999]     24.79482  77.04454  72.39902       305
##   [1000]     24.79482  77.04454  72.39902       370
##   -------
##   seqinfo: 23 sequences from an unspecified genome; no seqlengths
{% endhighlight %}
 
# Conclusion
 
It is very important to avoid to re-invent the wheel. Since the narrowPeak and the broadPeak are not directly supported by the `import` function in the `rtracklayer` package, we could be tempted to import them manually with the `read.table` function.
 
Not only would this be more complicated because we need to add the names of all the columns before converting in the `GRanges` format, but more importantly we need to make sure to convert the 0-based coordinate system of the BED file to the 1-based coordinate system of the `GRanges`:
 

{% highlight r %}
tbl_bed <- read.table(bed_file, header = FALSE)
head(tbl_bed)
{% endhighlight %}



{% highlight text %}
##      V1        V2        V3     V4   V5 V6
## 1  chr1  62902073  62903462 Rank_1 1859  .
## 2  chr1 167188355 167189834 Rank_2 1720  .
## 3  chr3  52321759  52322314 Rank_3 1632  .
## 4 chr12   2112635   2113682 Rank_4 1606  .
## 5 chr17  37843581  37845970 Rank_5 1597  .
## 6 chr22  41487757  41489245 Rank_6 1563  .
{% endhighlight %}



{% highlight r %}
gr_bed
{% endhighlight %}



{% highlight text %}
## GRanges object with 1000 ranges and 2 metadata columns:
##          seqnames                 ranges strand   |        name     score
##             <Rle>              <IRanges>  <Rle>   | <character> <numeric>
##      [1]     chr1 [ 62902074,  62903462]      *   |      Rank_1      1859
##      [2]     chr1 [167188356, 167189834]      *   |      Rank_2      1720
##      [3]     chr3 [ 52321760,  52322314]      *   |      Rank_3      1632
##      [4]    chr12 [  2112636,   2113682]      *   |      Rank_4      1606
##      [5]    chr17 [ 37843582,  37845970]      *   |      Rank_5      1597
##      ...      ...                    ...    ... ...         ...       ...
##    [996]    chr12 [111180217, 111180675]      *   |    Rank_996       770
##    [997]    chr12 [ 79941170,  79942669]      *   |    Rank_997       770
##    [998]     chr1 [ 41445363,  41445648]      *   |    Rank_998       770
##    [999]    chr17 [ 55161909,  55162356]      *   |    Rank_999       770
##   [1000]     chr3 [181425978, 181427339]      *   |   Rank_1000       770
##   -------
##   seqinfo: 23 sequences from an unspecified genome; no seqlengths
{% endhighlight %}
 
Notice how the start value is different depending of the strategy used to import the file. The correct one is the one obtained with the `import` function.
