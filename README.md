# What is is and what does it do?

Creates Array of 2d Coordiantes from multidimensional data.

# Installation

'npm i Kibusch_Pointcloud --save'

# How to use it


```
import {pca_simple} from 'Kibusch_Pointcloud';


var data = [[1,2,3,4],[1,2,3,4],[1,2,3,4]];

var labels = ["Name1","Name2","Name3"];


var pca = PCA_simple(data, labels, scale=true, maxValue=500);
```



# Options 

* *scale* - _boolean_ (Defaults to true)
* *scale* - _float_ (Defaults to 500)