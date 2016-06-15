## Shiny Module Drug Screen
###Input: raw data and/or processed data

####Requirements:

1. Columns

  * Raw Data

    | column name | description | optional/required? |
    |---|---|---|
    | sample | name of the cell line | **required** |
    | replicate | replicate of the sample | optional |
    | conc | concentrantion of the drug | **required**|
    | drug |  name of the drug | **required** |
    | normViability |  normalized viability | **required** |

  * Processed Data

    | column name | description | optional/required? |
    |---|---|---|
    | sample | name of the cell line | **required** |
    | AC50 | value of AC50 in micromolar(uM) | optional |
    | IC50 | value of IC59 in micromolar(uM) | optional |
    | maxResp | value of maximum response in scale 1-100 | **required** |
    | curveClass | defined curve class | optional |
    | drugID | drug ID | optional |
    | drug | name of the drug| **required** |
    | target | target class of the drug | optional |
    | AUC | area under the curve using the Simpson’s rule or the trapezoid method| **required** |

2. Unit/Scale

  * Raw Data
    
      **conc**: molar(M)

  * Processed Data

      **AC50** and **IC50** : micromolar (uM) 

      **maxResp** : scale of 1 - 100



### Output: A data-exploring shiny app	with 3 components: a data selecting tab, a filter tab, and a plot-display window. 

#### Example 1 - Using proccessed data
![alt text](../img/ds_example1_1.png "Example 1")



