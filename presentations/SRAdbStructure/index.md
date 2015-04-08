---
title: SRA database
author: José Alquicira Hernández
highlighter: highlight.js
job: null
mode: selfcontained
hitheme: tomorrow
subtitle: Tables and fields
framework: io2012
widgets: []
---

## SRA database structure
<center>
<img src="DDBJ_SRA.jpg" style="width:600px;height:550px"/>
</center>

--- 

## Exploring SRA database


```r
# Set path
setwd("/home/joseah/Documents/jeff_leek_lab/SRA")

# Load library
library(RSQLite)

# Path of sqlite file
sqlfile <- "/home/joseah/SRAmetadb.sqlite"

# Create connection
sra_con <- dbConnect(SQLite(),sqlfile)
```

--- 


```r
tables <- data.frame(dbListTables(sra_con)); knitr::kable(tables, col.names ="Tables")
```



|Tables          |
|:---------------|
|col_desc        |
|experiment      |
|fastq           |
|metaInfo        |
|run             |
|sample          |
|sra             |
|sra_ft          |
|sra_ft_content  |
|sra_ft_segdir   |
|sra_ft_segments |
|study           |
|submission      |

---


```r
fields_col_des <- dbListFields(sra_con, "col_desc")
knitr::kable(fields_col_des, col.names ="col_desc table")
```



|col_desc table |
|:--------------|
|col_desc_ID    |
|table_name     |
|field_name     |
|type           |
|description    |
|value_list     |
|sradb_updated  |

```r
query <- 
  "SELECT table_name, field_name, description 
  FROM col_desc 
  WHERE table_name = 'study' OR table_name = 'experiment' OR 
  table_name = 'sample' OR table_name = 'run';"
x <- dbGetQuery(sra_con, query)
```

---

### Table descriptions

<iframe src="table_descriptions.html"> </iframe>

---

### Examples with selected important fields

<iframe src="query.html"> </iframe>

---



## Conclusions (1)

- 57 fields were selected from the following tables:
  - Study 
  - Experiment
  - Sample
  - Run
- The **sample_attribute** field contains important information about ethnicity, sex, cell line, etc. Most of the time, this field is homogenous.
- We can look for the information about ethnicity in other fields.

---

## Conclusions (2)
- We could potentially analyze data divided into some categories as:
  - Sex
  - Ethnicity (differences in expression due to ancestry)
  - Cell line (diferential expression among different cell lines)
  - Center (Looking for batch effects)
  - Platform models
  - Date of submission (differences between older and new platforms)
  - Library selection (differences between the methods used to select and/or enrich the material being sequenced)


---

```r
options(width = 120)
devtools::session_info()
```

 setting  value                       
 version  R version 3.1.3 (2015-03-09)
 system   x86_64, linux-gnu           
 ui       RStudio (0.98.1091)         
 language en_US                       
 collate  en_US.UTF-8                 
 tz       <NA>                        

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
 magrittr    * 1.5       2014-11-22 CRAN (R 3.1.3)                   
 markdown    * 0.7.4     2014-08-24 CRAN (R 3.1.2)                   
 RJSONIO     * 1.3-0     2014-07-28 CRAN (R 3.1.2)                   
 rmarkdown   * 0.5.1     2015-01-26 CRAN (R 3.1.2)                   
 RSQLite       1.0.0     2014-10-25 CRAN (R 3.1.2)                   
 rstudio     * 0.98.1091 2014-12-01 local                            
 rstudioapi  * 0.2       2014-12-31 CRAN (R 3.1.2)                   
 slidify       0.3.52    2015-02-12 Github (ramnathv/slidify@1e11ccf)
 stringr     * 0.6.2     2012-12-06 CRAN (R 3.1.2)                   
 whisker     * 0.3-2     2013-04-28 CRAN (R 3.1.2)                   
 yaml        * 2.1.13    2014-06-12 CRAN (R 3.1.2)                   


