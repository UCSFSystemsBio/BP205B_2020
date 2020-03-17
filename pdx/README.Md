### VISION analysis of PDX data

VISION genrates a web interface with a representation of single cell data, as well as a score for the activity of a given signature (gene-set) by aggregating the expression, for each single cell.

VISION also conducts a test for locality - that is, asking whether a specific signature enriched in a specific cluster of cells. 

There is VISION result object here (new and improved with 100% less mouse). To use it, first untar / gzip them with tar -xvf file.rda.tar.gz.

Then, follow the file launch_vision.R to install VISION and load the interface.

```
install.packages('devtools')
devtools::install_github("YosefLab/VISION")
library(VISION)

## Loads pdx_vis_result into your environment
load('/path/to/object.rda')

VISION::viewResults(object = pdx_vis_result,port=8787,browser= T)
``` 

Based on Hani + Matt's suggestion, we should use **PDX_agg_med_norm_all.rda**

Note: The same convention for determining signature activity is used here - UP and DOWN signatures are scored by aggregating expression, but for signatures that have both UP and DOWN components, activity is determined as UP + -1*DOWN, to remove any anticorrelations. 

Additionally, cells were labeled using the KNN approach with 10 neighbors, to yield a subtype:
- LM-H (True Mets)
- LM-L (Mets surrounded by parentals)
- p-L (True parentals)
- p-H (parentals surrounded by Mets).