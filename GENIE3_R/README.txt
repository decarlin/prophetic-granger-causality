R implementation of GEne Network Inference with Ensemble of Trees (GENIE3)

Reference:
Inferring regulatory networks from expression data using tree-based methods
Van Anh HUYNH-THU, Alexandre IRRTHUM, Louis WEHENKEL, Pierre GEURTS
PLoS ONE vol. 5(9): e12776

The source code is provided as a single file "genie3.R".

In addition, we provide an expression data file "dream_multi_100_1.tsv"
from the DREAM4 challenge.

DREAM4 reference:
The DREAM4 In Silico network challenge.
http://wiki.c2b2.columbia.edu/dream09/index.php/D4c2)

The implementation is a research prototype and is provided AS IS. No
warranties or guarantees of any kind are given. Do not distribute this
code or use it other than for your own research without permission of
the author.

Prerequisites
-------------

You will have to install R (http://www.r-project.org).

You must also install the R package randomForest.
Random Forests references:
Breiman L (2001) Random forests. Machine Learning 45: 5-32
http://stat-www.berkeley.edu/users/breiman/RandomForests/

This installation can be done from R with this command (as root under Linux):
> install.packages("randomForest")

Alternatively, under Windows, you can install the package randomForest
from the RGui menu with "Packages > Install package(s)..."

Usage
-----

Example usage (with the provided file "dream_multi_100_1.tsv"):

> source("genie3.R")

This command will load the source file. You must be in the directory where "genie3.R"
is when you issue this command. You can change your current directory with the R
command setwd, or, under Windows, from the RGui menu with "File > Change dir...".
Alternatively, you can give the full path to the file.

> expr.matrix <- read.expr.matrix("dream_multi_100_1.tsv", form="rows.are.samples")

This command will read the expression matrix from file. Again, you must be in the
directory where the expression file is when you issue this command, or give the full
path. The command will automatically recognize if gene names and/or sample names are
provided in the first row and/or column. The form parameter is important to set
correctly. It tells if every row in the expression file corresponds to a gene
("rows.are.genes"), or if every row corresponds to a sample ("rows.are.samples")

> weight.matrix <- get.weight.matrix(expr.matrix)

This command computes the weighted adjacency matrix of the gene network with the
GENIE3 algorithm (Random Forests). In the weight matrix, element w_ij (row i,
column j) gives the "importance" of the link from regulatory gene i to target gene i,
with high scores corresponding to more likely regulatory links. You can specify that
only a subset of the genes are to be used as input genes (regulatory genes) by
passing an array of gene indices through the input.idx parameter:

> weight.matrix <- get.weight.matrix(expr.matrix, input.idx=30:50)

Here, for example, only genes from 30 to 50 (included) will be used as input genes.
This can be useful when you know which genes are transcription factors. If you have
a list of gene names to be used as input genes, you can also use it like this:

> input.genes = c("G3", "G4", "G7", "G8", "G9", "G10", "G11", "G12", "G17", "G33")
> weight.matrix <- get.weight.matrix(expr.matrix, input.idx=input.genes)

You can obtain the list of regulatory links (from more likely to less likely) with
this command:

> link.list <- get.link.list(weight.matrix, report.max=1000)

Usually, one is only interested in extracting the more likely regulatory links. The
optional parameter report.max sets the number of top scoring links to report.

Additional function parameters are explained in the source file.

Contact
-------

Alexandre Irrthum
alexandre.irrthum@ulg.ac.be

Van Anh Huynh-Thu
vhuynht@inf.ed.ac.uk

Department of Electrical Engineering and Computer Science,
Systems and Modeling, and
GIGA-Research, Bioinformatics and Modeling
University of Liege, Belgium
