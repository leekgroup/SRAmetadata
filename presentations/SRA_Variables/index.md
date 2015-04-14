---
title: SRA variables
author: José Alquicira Hernández
highlighter: highlight.js
job: null
mode: selfcontained
hitheme: tomorrow
subtitle: 
framework: io2012
widgets: []
---

# Relevant variables in sample table



|Variables   |
|:-----------|
|cell_type   |
|tissue      |
|cell_line   |
|age         |
|disease     |
|population  |
|sex         |
|source_name |

---

# Cell type
[1] "Total number of samples: 6580"
<iframe src="cell_type.html"> </iframe>

---

# Tissue
[1] "Total number of samples: 6927"
<iframe src="tissue.html"> </iframe>

---

# Cell line
[1] "Total number of samples: 5614"
<iframe src="cell_line.html"> </iframe>

---

# Age
[1] "Total number of samples: 4052"
<iframe src="age.html"> </iframe>

---

# Disease
[1] "Total number of samples: 1208"
<iframe src="disease.html"> </iframe>

---


# Population
[1] "Total number of samples: 1670"
<iframe src="population.html"> </iframe>

---

# Sex
[1] "Total number of samples: 2928"
<iframe src="sex.html"> </iframe>

---

# Source name
[1] "Total number of samples: 20513"
<iframe src="source_name.html"> </iframe>

---

```r
options(width = 120)
dbDisconnect(sra_con)
```


```r
devtools::session_info()$packages
```

 package     * version   date       source                           
 codetools   * 0.2-11    2015-03-10 CRAN (R 3.1.3)                   
 DBI           0.3.1     2014-09-24 CRAN (R 3.1.2)                   
 devtools    * 1.7.0     2015-01-17 CRAN (R 3.1.2)                   
 digest      * 0.6.8     2014-12-31 CRAN (R 3.1.2)                   
 DT          * 0.0.34    2015-04-03 Github (rstudio/DT@d446dde)      
 evaluate    * 0.5.5     2014-04-29 CRAN (R 3.1.2)                   
 formatR     * 1.0       2014-08-25 CRAN (R 3.1.2)                   
 htmltools   * 0.2.6     2014-09-08 CRAN (R 3.1.2)                   
 htmlwidgets * 0.3.2     2014-12-09 CRAN (R 3.1.3)                   
 knitr       * 1.9       2015-01-20 CRAN (R 3.1.2)                   
 magrittr      1.5       2015-04-13 Github (smbache/magrittr@89f143d)
 markdown    * 0.7.4     2014-08-24 CRAN (R 3.1.2)                   
 RJSONIO     * 1.3-0     2014-07-28 CRAN (R 3.1.2)                   
 RSQLite       1.0.0     2014-10-25 CRAN (R 3.1.2)                   
 rstudio     * 0.98.1091 2014-12-01 local                            
 rstudioapi  * 0.2       2014-12-31 CRAN (R 3.1.2)                   
 slidify       0.3.52    2015-02-12 Github (ramnathv/slidify@1e11ccf)
 stringr     * 0.6.2     2012-12-06 CRAN (R 3.1.2)                   
 whisker     * 0.3-2     2013-04-28 CRAN (R 3.1.2)                   
 yaml        * 2.1.13    2014-06-12 CRAN (R 3.1.2)                   



