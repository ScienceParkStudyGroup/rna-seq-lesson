---
title: "Assignments"
---


{% include links.md %}


## Description
Assignments are mainly designed for student curriculums but can be also used to test your working knowledge upon self-study. 


# Assignments
Each assigment comprises two comparisons (comp #1 and comp #2). 

## Assignment 1

| seed                	| infection            	| dpi 	| comparison 	| assignment 	|
|---------------------	|----------------------	|-----	|------------	|------------	|
| MgCl2                	| P. syringae DC3000   	| 2   	| comp #1    	| A1          |
| S. melonis Fr1     	  | P. syringae DC3000 	  | 2   	| comp #1     | A1         	|
| MgCl2               	| P. syringae DC3000   	| 7   	| comp #2    	| A1         	|
| S. melonis Fr1 	      | P. syringae DC3000 	  | 7   	| comp #2  	  | A1         	|

* Comparison 1: 2978 genes (p < 0.01)
* Comparison 2: 3031 genes (p < 0.01)

## Assignment 2

| seed                	| infection            	| dpi 	| comparison 	| assignment 	|
|---------------------	|----------------------	|-----	|------------	|------------	|
| MgCl2                	| Mock                	| 2     | comp #1     | A2     	  |
| MgCl2             	  | P. syringae DC3000    | 2     | comp #1     | A2         	|
| MgCl2                	| Mock               	  | 7     | comp #2  	  | A2         	|
| MgCl2         	      | P. syringae DC3000    | 7     | comp #2     | A2         	|

* Comparison 1: 795 genes (p < 0.01)
* Comparison 2: 3249 genes (__p < 0.001__) and around 4979 at p < 0.01 (comparison from the tutorial).

## Assignment 3

| seed                	| infection            	| dpi 	| comparison 	| assignment 	|
|---------------------	|----------------------	|-----	|------------	|------------	|
| M. extorquens PA1     | Mock               	  | 2   	| comp #1    	| A3        	|
| S. melonis Fr1        | Mock                  | 2   	| comp #1   	| A3          	|
| M. extorquens PA1     | Mock               	  | 7   	| comp #2    	| A3          	|
| S. melonis Fr1        | Mock              	  | 7   	| comp #2    	| A3          	|


* Comparison 1: 1008 genes (p < 0.01)
* Comparison 2: 680 genes (p < 0.01) 

## Assigment 4

| seed                	| infection            	| dpi 	| comparison 	| assignment 	|
|---------------------	|----------------------	|-----	|------------	|------------	|
| M. extorquens PA1     | P. syringae DC3000    | 2   	| comp #1    	| A4    	|
| S. melonis Fr1        | P. syringae DC3000    | 2   	| comp #1    	| A4          	|
| M. extorquens PA1     | P. syringae DC3000   	| 7   	| comp #2    	| A4          	|
| S. melonis Fr1        | P. syringae DC3000   	| 7   	| comp #2    	| A4          	|

* Comparison 1: 3189 genes (p < 0.01)
* Comparison 2: 696 genes (p < 0.01) 

## Tasks
1. Rephrase the scientific question (0.5 point).
2. Compute the tables of differential expressed genes for both comparisons. Create a histogram of corrected p-values for each comparison (QC). Annotate the genes with their gene symbols and annotation using biomartr (2 points). 
3. Create a volcano plot for each comparison to display the differential genes. Justify your choice for the fold change cutoff based on the distribution of your log2 fold change values (e.g. median, 75th percentile, etc). 
