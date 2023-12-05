
<!-- rnb-text-begin -->

---
title: "Plot Supplementary Figure"
output: html_notebook
---

## 01. Load libraries


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShkcGx5cilcbmxpYnJhcnkoZ2dwbG90MilcbmxpYnJhcnkoZ2dwdWJyKVxuYGBgIn0= -->

```r
library(dplyr)
library(ggplot2)
library(ggpubr)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## 02. Define constants


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuT1VUUFVULkRJUiA8LSBcIi9Vc2Vycy9sZWV3b3JrL0RvY3VtZW50cy9SZXNlYXJjaC9wcm9qZWN0cy9wcm9qZWN0X2FjZS9kYXRhL3Byb2Nlc3NlZC9maWd1cmVzXCJcbkNPTU1PTi5USEVNRSA8LSB0aGVtZShheGlzLnRpdGxlID0gZWxlbWVudF90ZXh0KHNpemUgPSAxNCksXG4gICAgICAgICAgICAgICAgICAgICAgYXhpcy50ZXh0ID0gZWxlbWVudF90ZXh0KHNpemUgPSAxNCksXG4gICAgICAgICAgICAgICAgICAgICAgbGVnZW5kLnRpdGxlID0gZWxlbWVudF90ZXh0KHNpemUgPSAxNCksXG4gICAgICAgICAgICAgICAgICAgICAgbGVnZW5kLnRleHQgPSBlbGVtZW50X3RleHQoc2l6ZSA9IDE0KSlcbkdST1VQLkFMUEhBIDwtIDAuNjE4XG5ET0RHRS5XSURUSCA8LSAwLjVcbkVSUk9SQkFSLldJRFRIIDwtIDAuMzg0XG5CT1hQTE9ULldJRFRIIDwtIDAuNjE4XG5MSU5FLkNJUkNMRS5TSVpFIDwtIDIuNjJcbkdST1VQLkNPTE9SUyA8LSBjKFwiaW5ncm91cFwiID0gXCIjN0JBRUQxXCIsXG4gICAgICAgICAgICAgICAgICBcIm91dGdyb3VwXCIgPSBcIiNBNEVFNTdcIilcbk1FVEhPRFMuQ09MT1JTIDwtIGMoXCJhY2Utc1wiID0gXCIjODQ1RUMyXCIsXG4gICAgICAgICAgICAgICAgICAgIFwiYWNlXCIgPSBcIiNENjVEQjFcIixcbiAgICAgICAgICAgICAgICAgICAgXCJyYW5kb21cIiA9IFwiI0ZGOTY3MVwiLFxuICAgICAgICAgICAgICAgICAgICBcInJlcGVhdGVkXCIgPSBcIiNBQkFCQUFcIilcbk1FVEhPRFMuT1JERVIgPC0gYyhcImFjZS1zXCIsXG4gICAgICAgICAgICAgICAgICAgXCJhY2VcIixcbiAgICAgICAgICAgICAgICAgICBcInJhbmRvbVwiLFxuICAgICAgICAgICAgICAgICAgIFwicmVwZWF0ZWRcIilcbmBgYCJ9 -->

```r
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 14),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 14))
GROUP.ALPHA <- 0.618
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.384
BOXPLOT.WIDTH <- 0.618
LINE.CIRCLE.SIZE <- 2.62
GROUP.COLORS <- c("ingroup" = "#7BAED1",
                  "outgroup" = "#A4EE57")
METHODS.COLORS <- c("ace-s" = "#845EC2",
                    "ace" = "#D65DB1",
                    "random" = "#FF9671",
                    "repeated" = "#ABABAA")
METHODS.ORDER <- c("ace-s",
                   "ace",
                   "random",
                   "repeated")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## 03. Plot Sars-Cov2 Dataset


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGxvdC5wcmVjaXNpb24gPC0gUGxvdFByZWNpc2lvbihcbiAgZGYuZXZhbHVhdGlvbi5tZXRyaWNzID0gZGYuZXZhbHVhdGlvbi5tZXRyaWNzXG4pXG5cbmBgYCJ9 -->

```r
plot.precision <- PlotPrecision(
  df.evaluation.metrics = df.evaluation.metrics
)

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiYHN1bW1hcmlzZSgpYCBoYXMgZ3JvdXBlZCBvdXRwdXQgYnkgJ2dyb3VwJy4gWW91IGNhbiBvdmVycmlkZSB1c2luZyB0aGUgYC5ncm91cHNgIGFyZ3VtZW50LlxuIn0= -->

```
`summarise()` has grouped output by 'group'. You can override using the `.groups` argument.
```



<!-- rnb-output-end -->

<!-- rnb-frame-begin eyJtZXRhZGF0YSI6eyJjbGFzc2VzIjpbInRibF9kZiIsInRibCIsImRhdGEuZnJhbWUiXSwibnJvdyI6NCwibmNvbCI6Miwic3VtbWFyeSI6eyJBIHRpYmJsZSI6WyI0IMOXIDIiXX19LCJyZGYiOiJINHNJQUFBQUFBQUFBMTFSelVyRU1CQk9mNkp1d1ZVVWZJdjJJc0xlZW5EWkJ4Q0Z2WWpFTkxzVzI2UWtxVDhYOFN6NExDSWVQT3RqaU5DSDhDcGJwejhqYkE5dDV2c3k4MzJUbVpQcC9EQ1lCNFFRai9pT1F6d0tJYUZucDdOd1FvanZBbkNJVDBiTmVRZEoreEFBNlczRDZiZVhMWWJxRnE4VmJHVGlSbVFHb3QzMnRtTXA0eUkwUGZBQVlMWm1NbEY1ajdhMEtBU3pJaGxvVXA0eGc1TC9SZ3ZHcmRJUXJlQWJOMlp4OWY3eU5QTys0K3JuK09EaHk0bXIxL3ZBVFBiaXovUDQrZTFqT3BTVkxCY282eUs1MUtvc2VqRE9CWk1YaFJZOE5hbVNnL3FSVnJjUmFqU3pjUi9oVjlmMWI1Zm8xWU9XZzRSWkZpMDBsSFJ0cjhsdHFzS0NDWWk1emNEcG9OalJBMktubEkxNUV2S3JVbDZIUjQxQnZ4elNOK1QwQzhQWTdTeDliSXppTElWY3BoS1hRak4yS1RJY0FUeXlmV05VNkZSYWZBbXdKckxLTXN3THVNcVE2VmF5K2dNT2crd1Nad0lBQUE9PSJ9 -->

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["group"],"name":[1],"type":["fctr"],"align":["left"]},{"label":["mean_precision"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"ace-s","2":"0.4955863"},{"1":"ace","2":"0.4992225"},{"1":"random","2":"0.4952071"},{"1":"repeated","2":"0.1434707"}],"options":{"columns":{"min":{},"max":[10],"total":[2]},"rows":{"min":[10],"max":[10],"total":[4]},"pages":{}}}
  </script>
</div>

<!-- rnb-frame-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## 04. Plot Rojas et al., Nature 2023


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuICBkZi5wbG90IDwtIGRmLmV2YWx1YXRpb24ubWV0cmljcyAlPiVcbiAgICBkcGx5cjo6Z3JvdXBfYnkoZ3JvdXAsIHBlcmNfaW1tdW5vZ2VuaWMpICU+JVxuICAgIGRwbHlyOjpzdW1tYXJpc2UoXG4gICAgICBzZCA9IHNkKHByZWNpc2lvbiksXG4gICAgICBwcmVjaXNpb24gPSBtZWFuKHByZWNpc2lvbiksXG4gICAgICByZWNhbGwgPSBtZWFuKHNlbnNpdGl2aXR5KVxuICAgIClcblxuYGBgIn0= -->

```r
  df.plot <- df.evaluation.metrics %>%
    dplyr::group_by(group, perc_immunogenic) %>%
    dplyr::summarise(
      sd = sd(precision),
      precision = mean(precision),
      recall = mean(sensitivity)
    )

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiYHN1bW1hcmlzZSgpYCBoYXMgZ3JvdXBlZCBvdXRwdXQgYnkgJ2dyb3VwJy4gWW91IGNhbiBvdmVycmlkZSB1c2luZyB0aGUgYC5ncm91cHNgIGFyZ3VtZW50LlxuIn0= -->

```
`summarise()` has grouped output by 'group'. You can override using the `.groups` argument.
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


## 05. Merge


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucCA8LSBnZ2FycmFuZ2UocGxvdGxpc3QgPSBsaXN0KHBhbmVsLmEsIHBhbmVsLmIpLCBucm93ID0gMSwgYWxpZ24gPSBcImh2XCIsIGNvbW1vbi5sZWdlbmQgPSBUKVxuYGBgIn0= -->

```r
p <- ggarrange(plotlist = list(panel.a, panel.b), nrow = 1, align = "hv", common.legend = T)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->

