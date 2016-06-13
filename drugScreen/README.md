## Shiny Module Drug Screen
###Input: raw data and/or processed data

####Requirements:

1. Columns

  * Raw Data

    | column name | description | optional/required? |
    |---|---|---|
    | sample | name of the cell line | **required** |
    | replicate | blah | optional |
    | conc | concentrantion of the drug | **required**|
    | drug |  name of the drug | **required** |
    | normViability |  normalized viability | **required** |

  * Processed Data

    | column name | description | optional/required? |
    |---|---|---|
    | sample | name of the cell line | **required** |
    | AC50 | blah | optional |
    | IC50 | blah | optional |
    | maxResp | blah | **required** |
    | curveClass | blah | optional |
    | drugID | blah | optional |
    | drug | name of the drug| **required** |
    | target | target class of the drug | optional |
    | AUC | blah | **required** |

2. Unit/Scale

  * Raw Data
    
      **conc**: molar(M)

  * Processed Data

      **AC50** and **IC50** : micromolar (uM) 

      **maxResp** : scale of 1 - 100




### Output: A data-exploring shiny app	

#### Example 1 - Using proccessed data
![alt text](../img/ds_example1_1.png "Example 1")



