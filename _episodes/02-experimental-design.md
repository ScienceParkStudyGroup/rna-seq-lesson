---
title: "Experimental design"
teaching: 15
exercises: 15
questions:
- "What should I consider when designing an RNA-seq experiment?"
- "What is being measured in a RNA-seq experiment?"
- "What are common sources of variation in RNA-seq?"
- "What is statistical power in the context of RNA-seq?"
objectives:
- "Be able to understand and describe a RNA-seq experimental design."
- "Learn how to generate a reproducible code notebook with Rmarkdown."
keypoints:
- "A."

---

## Table of Contents
1. [Introduction](#introduction)
2. [Statistical power](#statistical-power)


## Introduction

Before diving into the specifics of RNA-seq, it is important to remember that an RNA-seq experiment remains _an experiment_.  


We will exemplify some concepts in experimental design by using the example on Arabidopsis ([see the introduction episode](http://0.0.0.0:4000/01-introduction/index.html).

## Asking a precise question

It starts by asking a question in plain English: ""

## Statistical power

First a small funny way to remember what are type I and type II errors (Taken from the [Essential Guide to Effect Sizes](#references)):

<img src="../img/type-i-and-type-ii-errors.jpg" width="400px" alt="type-i-and-type-ii-errors">

Also taken from the [Essential Guide to Effect Sizes](#references): 
> The power of any test of statistical significance is defined as the probability that it will reject a false null hypothesis. Statistical power is inversely related to beta or the probability of making a Type II error. 
__In short, power = 1 – β__.   
  
> __In plain English, statistical power is the likelihood that a study will detect an effect when there is an effect there to be detected. If statistical power is high, the probability of making a Type II error, or concluding there is no effect when, in fact, there is one, goes down__.


Type I error

Type II error: 

The power of any test of statistical significance is defined as the probability that it will reject a false null hypothesis. Statistical power is inversely related to beta or the probability of making a Type II error. In short, power = 1 – β.

In plain English, statistical power is the likelihood that a study will detect an effect when there is an effect there to be detected. If statistical power is high, the probability of making a Type II error, or concluding there is no effect when, in fact, there is one, goes down.

## References
- The Essential Guide to Effect Sizes: Statistical Power, Meta-Analysis, and the Interpretation of Research Results. 2015. Paul D. Ellis. Cambridge University Press.  
- Replication Concepts for Bioenergy Research Experiments. 
- Concepts of Experimental Design. SAS White Paper. 