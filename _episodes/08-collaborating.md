---
title: "Collaborating with Github"
teaching: 30
exercises: 60 
questions:
- "How can I co-develop and collaborate on code with another scientist?"
- "How can I stay in sync on code"
- "How can I solve conflicts that arise from that collaboration?"
- "What are Github "
objectives:
- "Be able to work together on a R script through RStudio and Github integration."
- ""
keypoints:
- "Github allows you to synchronise work efforts and collaborate with other scientists on (R) code."
- ""
- ""
---

## Introduction

The collaborative power of GitHub and RStudio is really game changing. So far we've been collaborating with our most important collaborator: ourselves. But, we are lucky that in science we have so many other collaborators, so let's learn how to accelerate our collaborations with them through GitHub! 

We are going to teach you the simplest way to collaborate with someone, which is for both of you to have privileges to edit and add files to a repository. GitHub is built for software developer teams, and there is a lot of features that limit who can directly edit files, but we don't need to start there. 

We will do this all with a partner, and we'll walk through some things all together, and then give you a chance to work with your collaborator on your own. 


## Objectives & Resources


### Objectives

We are going to create a website with a collaborator!

- create a new repo and give permission to a collaborator
- open as a new RStudio project!
- collaborate with a partner 
- explore github.com blame, history, issues

<!---**Resources**--->



## Create repo (Partner 1)

Team up with a partner sitting next to you. Partner 1 will create a new repository. We will do this in the same way that we did in Chapter \@ref(github): [Create a repository on Github.com]. 

## Create a gh-pages branch (Partner 1)

We aren't going to talk about branches very much, but they are a powerful feature of git/GitHub. I think of it as creating a copy of your work that becomes a parallel universe that you can modify safely because it's not affecting your original work. And then you can choose to merge the universes back together if and when you want. By default, when you create a new repo you begin with one branch, and it is named `master`. When you create new branches, you can name them whatever you want. However, if you name one `gh-pages` (all lowercase, with a `-` and no spaces), this will let you create a website. And that's our plan. So, Partner 1, do this to create a `gh-pages` branch: 

On the homepage for your repo on GitHub.com, click the button that says "Branch:master". Here, you can switch to another branch (right now there aren't any others besides `master`), or create one by typing a new name. 

<img src="../img/github-branch.png" width="350px">


Let's type `gh-pages`. 

<img src="../img/github_create-branch_gh-pages.png" width="350px"> 

Let's also change `gh-pages` to the default branch and delete the master branch: this will be a one-time-only thing that we do here: 

First click to control branches:

<img src="../img/github-branch2.png" width="350px"> 

And then click to change the default branch to `gh-pages`. I like to then delete the `master` branch when it has the little red trash can next to it. It will make you confirm that you really want to delete it, which I do!

<img src="../img/github-change-branch.png" width="450px"> 


## Give your collaborator administration privileges (Partner 1 and 2)

Now, Partner 1, go into Settings > Collaborators > enter Partner 2's (your collaborator's) username. 

Partner 2 then needs to check their email and accept as a collaborator. Notice that your collaborator has "Push access to the repository" (highlighted below):

![](../img/github_collab.png) 

## Clone to a new Rproject  (Partner 1)

Now let's have Partner 1 clone the repository to their local computer. We'll do this through RStudio like we did before (see Chapter \@ref(github): [Clone your repository using RStudio]), but with a final additional step before hitting "Create Project": select "Open in a new Session".

<img src="../img/github_clone_newproject.png" width="450px"> 


Opening this Project in a new Session opens up a new world of awesomeness from RStudio. Having different RStudio project sessions allows you to keep your work separate and organized. So you can collaborate with this collaborator on this repository while also working on your other repository from this morning. I tend to have a lot of projects going at one time:



![](../img/Rproj_screenshot.jpg)



Have a look in your git tab. 

Like we saw this morning, when you first clone a repo through RStudio, RStudio will add an `.Rproj` file to your repo. And if you didn't add a `.gitignore` file when you originally created the repo on GitHub.com, RStudio will also add this for you. So, Partner 1, let's go ahead and sync this back to GitHub.com. 

Remember: 

![](../img/commit_overview.png)

 

Let's confirm that this was synced by looking at GitHub.com again. You may have to refresh the page, but you should see this commit where you added the `.Rproj` file.

## Clone to a new Rproject  (Partner 2)

Now it's Partner 2's turn! Partner 2, clone this repository following the same steps that Partner 1 just did. When you clone it, RStudio should not create any new files — why? Partner 1 already created and pushed eht the `.Rproj` and `.gitignore` files so they already exist in the repo.  

## Edit a file and sync (Partner 2)

Let's have Partner 2 add some information to the README.md. Let's have them write: 
```
Collaborators: 

- Partner 2's name

```

When we save the README.md, And now let's sync back to GitHub. 

![](../img/commit_overview.png)




When we inspect on GitHub.com, click to view all the commits, you'll see commits logged from both Partner 1 and 2!

> Question: Would you still be able clone a repository that you are not a collaborator on? What do you think would happen? Try it! Can you sync back? 

## State of the Repository

OK, so where do things stand right now? GitHub.com has the most recent versions of all the repository's files. Partner 2 also has these most recent versions locally. How about Partner 1? 

Partner 1 does not have the most recent versions of everything on their computer. 

Question: How can we change that? Or how could we even check? 

Answer: PULL. 

Let's have Partner 1 go back to RStudio and Pull. If their files aren't up-to-date, this will pull the most recent versions to their local computer. And if they already did have the most recent versions? Well, pulling doesn't cost anything (other than an internet connection), so if everything is up-to-date, pulling is fine too. 

I recommend pulling every time you come back to a collaborative repository. Whether you haven't opened RStudio in a month or you've just been away for a lunch break, pull. It might not be necessary, but it can save a lot of heartache later.

## Merge conflicts

What kind of heartache are we talking about? Let's explore. **Stop and watch me create and solve a merge conflict with my Partner 2, and then you will have time to recreate this with your partner.** Here's what I am going to do:

Within a file, GitHub tracks changes line-by-line. So you can also have collaborators working on different lines within the same file and GitHub will be able to weave those changes into each other -- that's it's job! It's when you have collaborators working on *the same lines within the same file* that you can have **merge conflicts**. Merge conflicts can be frustrating, but they are actually trying to help you (kind of like R's error messages). They occur when GitHub can't make a decision about what should be on a particular line and needs a human (you) to decide. And this is good -- you don't want GitHub to decide for you, it's important that you make that decision. 

So let's test this. Let's have both Partners 1 and 2 go to RStudio and pull so you have the most recent versions of all your files. Now, Partners 1 and 2, both go to the README, and on Line 7, write something, anything. I'm not going to give any examples because I want both Partners to write something different. And be sure to save the README. 

OK. Now, let's have Partner 2 sync: pull, stage, commit, push. Great. 

Now, when Partner 2 is done, let's have Partner 1 (me) try. 

Partner 1: pull ---- Error! Merge conflict!

![](../img/github_mergeconflict.png)

So Partner 1 is not allowed to pull, it failed. GitHub is protecting Partner 1 because if they did successfully pull, their work would be overwritten by whatever Partner 2 had written. So GitHub is going to make a human (Partner 1 in this case) decide. GitHub says, either commit this work first, or "stash it" (I interpret that as saving a copy of the README in another folder somewhere outside of this GitHub repository). 

Let's follow their advice and have Partner 1 commit. Great. Now let's pull again. 

Still not happy!

![](../img/github_mergeconflict2.png)



OK, actually, we're just moving along this same problem that we know that we've created: Both Partner 1 and 2 have both added new information to the same line. You can see that the pop-up box is saying that there is a CONFLICT and the merge has not happened. OK. We can close that window and inspect. 

Notice that in the git tab, there are orange `U`s; this means that there is an unresolved conflict, and it is not staged with a check anymore because modifications have occurred to the file since it has been staged. 

Let's look at the README file itself. We got a preview in the diff pane that there is some new text going on in our README file: 

```r
<<<<<<< HEAD
Julie is collaborating on this README.
=======
**Jamie is adding lines here.**
>>>>>>> 05a189b23372f0bdb5b42630f8cb318003cee19b
```

In this example, Partner 1 is Jamie and Partner 2 is Julie. GitHub is displaying the line that Julie wrote and the line Jamie wrote separated by `=======`. So these are the two choices that Partner 2 has to decide between, which one do you want to keep? Where where does this decision start and end? The lines are bounded by `<<<<<<<HEAD` and `>>>>>>>long commit identifier`. 

So, to resolve this merge conflict, Partner 2 has to chose, and delete everything except the line they want. So, they will delete the `<<<<<<HEAD`, `=====`, `>>>>long commit identifier` and one of the lines that they don't want to keep. 

Do that, and let's try again. In this example, we've kept Jamie's line: 

![](../img/github_mergeconflict3.png)



Then be sure to stage, and write a commit message. I often write "resolving merge conflict" or something so I know what I was up to. When I stage the file, notice how now my edits look like a simple line replacement (compare with the image above before it was re-staged): 

![](../img/github_mergeconflict4.png)

### Your turn

Create a merge conflict with your partner, like we did in the example above. And try other ways to get and solve merge conflicts. For example, when you get the following error message, try both ways (commit or stash. Stash means copy/move it somewhere else, for example, on your Desktop temporarily).

![](../img/github_mergeconflict.png)

## How do you avoid merge conflicts?

I'd say pull often, commit and sync often. 

Also, talk with your collaborators. Although our Ocean Health Index project is highly collaborative, we are actually rarely working on the exact same file at any given time. And if we are, we are also on Slack, Gchat, or sitting next to the person. 

But merge conflicts will occur and some of them will be heartbreaking and demoralizing. They happen to me when I collaborate with myself between my work computer and laptop. So protect yourself by pulling and syncing often! 

## Create your collaborative website

OK. Let's have Partner 2 create a new RMarkdown file. Here's what they will do: 

1. Pull!
1. Create a new RMarkdown file **and name it `index.Rmd`**. Make sure it's all lowercase, and named `index.Rmd`. This will be the homepage for our website! 
1. Maybe change the title inside the Rmd, call it "Our website"
1. Knit!
1. Save and sync your .Rmd and your .html files 
    - (pull, stage, commit, pull, push)
1. Go to GitHub.com and go to your rendered website! Where is it? Figure out your website's url from your github repo's url. For example: 
    - my github repo: <https://github.com/jules32/collab-research>
    - my website url: <https://jules32.github.io/collab-research/>
    - note that the url starts with my **username.github.io**
    
So cool! On websites, if something is called `index.html`, that defaults to the home page. So <https://jules32.github.io/collab-research/> is the same as <https://jules32.github.io/collab-research/index.html>. If you name your RMarkdown file `my_research.Rmd`, the url will become <https://jules32.github.io/collab-research/my_research.html>.

## Your turn

Here is some collaborative analysis you can do on your own. We'll be playing around with airline flights data, so let's get setup a bit. 

1. Person 1: clean up the README to say something about you two, the authors.
1. Person 2: edit the `index.Rmd` or create a new RMarkdown file: maybe add something about the authors, and knit it. 
1. Both of you: sync to GitHub.com (pull, stage, commit, push). 
1. Both of you: once you've both synced (talk to each other about it!), pull again. You should see each others' work on your computer.
1. Person 1: in the RMarkdown file, add a bit of the plan. We'll be exploring the `nycflights13` dataset. This is data on flights departing New York City in 2013.
1. Person 2: in the README, add a bit of the plan. 
1. Both of you: sync

## Explore on GitHub.com

Now, let's look at the repo again on GitHub.com. You'll see those new files appear, and the commit history has increased.

### Commit History

You'll see that the number of commits for the repo has increased, let's have a look. You can see the history of both of you. 

### Blame

Now let's look at a single file, starting with the README file. We've explored the "Raw" and "History" options in the top-right of the file, but we haven't really explored the "Blame" option. Let's look now. Blame shows you line-by-line who authored the most recent version of the file you see. This is super useful if you're trying to understand logic; you know who to ask for questions or attribute credit.

### Issues

Now let's have a look at issues. This is a way you can communicate to others about plans for the repo, questions, etc. Note that issues are public if the repository is public.

![](../img/github-issues.png)

Let's create a new issue with the title "NYC flights". 

In the text box, let's write a note to our collaborator. You can use Markdown in this text box, which means all of your header and bullet formatting will come through. You can also select these options by clicking them just above the text box. 

Let's have one of you write something here. I'm going to write: 

```
Hi @jafflerbach! 

# first priority

- explore NYC flights
- plot interesting things
```

Note that I have my collaborator's GitHub name with a `@` symbol. This is going to email her directly so that she sees this issue. I can click the "Preview" button at the top left of the text box to see how this will look rendered in Markdown. It looks good! 

Now let's click submit new issue. 

On the right side, there are a bunch of options for categorizing and organizing your issues. You and your collaborator may want to make some labels and timelines, depending on the project. 

Another feature about issues is whether you want any notifications to this repository. Click where it says "Unwatch" up at the top. You'll see three options: "Not watching", "Watching", and "Ignoring". By default, you are watching these issues because you are a collaborator to the repository. But if you stop being a big contributor to this project, you may want to switch to "Not watching". Or, you may want to ask an outside person to watch the issues. Or you may want to watch another repo yourself!

![](../img/github-collab.png)

Let's have Person 2 respond to the issue affirming the plan.

## NYC flights exploration

Let's continue this workflow with your collaborator, syncing to GitHub often and practicing what we've learned so far. We will get started together and then you and your collaborator will work on your own.

Here's what we'll be doing (from [R for Data Science's Transform Chapter](http://r4ds.had.co.nz/transform.html)):

**Data**: You will be exploring a dataset on flights departing New York City in 2013. These data are actually in a package called `nycflights13`, so we can load them the way we would any other package. 

Let's have Person 1 write this in the RMarkdown document (Person 2 just listen for a moment; we will sync this to you in a moment). 

```{r, message = FALSE}
library(nycflights13) # install.packages('nycflights13')
library(tidyverse)
```

This data frame contains all `r format(nrow(nycflights13::flights), big.mark = ",")` flights that departed from New York City in 2013. The data comes from the US [Bureau of Transportation Statistics](http://www.transtats.bts.gov/DatabaseInfo.asp?DB_ID=120&Link=0), and is documented in `?flights`.

```{r}
flights
```

Let's select all flights on January 1st with:

```{r}
filter(flights, month == 1, day == 1)
```

To use filtering effectively, you have to know how to select the observations that you want using the comparison operators. R provides the standard suite: `>`, `>=`, `<`, `<=`, `!=` (not equal), and `==` (equal). We learned these operations yesterday. But there are a few others to learn as well. 

#### Sync

Sync this RMarkdown back to GitHub so that your collaborator has access to all these notes. Person 2 should then pull and will continue with the following notes: 

### Logical operators

Multiple arguments to `filter()` are combined with "and": every expression must be true in order for a row to be included in the output. For other types of combinations, you'll need to use Boolean operators yourself: 

- `&` is "and" 
- `|` is "or" 
- `!` is "not"

Let's have a look:

The following code finds all flights that departed in November or December:

```{r, eval = FALSE}
filter(flights, month == 11 | month == 12)
```

The order of operations doesn't work like English. You can't write `filter(flights, month == 11 | 12)`, which you might literally translate into  "finds all flights that departed in November or December". Instead it finds all months that equal `11 | 12`, an expression that evaluates to `TRUE`. In a numeric context (like here), `TRUE` becomes one, so this finds all flights in January, not November or December. This is quite confusing!

A useful short-hand for this problem is `x %in% y`. This will select every row where `x` is one of the values in `y`. We could use it to rewrite the code above:

```{r, eval = FALSE}
nov_dec <- filter(flights, month %in% c(11, 12))
```

Sometimes you can simplify complicated subsetting by remembering De Morgan's law: `!(x & y)` is the same as `!x | !y`, and `!(x | y)` is the same as `!x & !y`. For example, if you wanted to find flights that weren't delayed (on arrival or departure) by more than two hours, you could use either of the following two filters:

```{r, eval = FALSE}
filter(flights, !(arr_delay > 120 | dep_delay > 120))
filter(flights, arr_delay <= 120, dep_delay <= 120)
```

Whenever you start using complicated, multipart expressions in `filter()`, consider making them explicit variables instead. That makes it much easier to check your work. 

## Your turn

OK: Person 2, sync this to GitHub, and Person 1 will pull so that we all have the most current information. 

With your partner, do the following tasks. Each of you should work on one task at a time. Since we're working closely on the same document, talk to each other and have one person create a heading and a R chunk, and then sync; the other person can then create a heading and R chunk and sync, and then you can both work safely. 

Remember to make your commit messages useful!

As you work, you may get merge conflicts. This is part of collaborating in GitHub; we will walk through and help you with these and also teach the whole group. 

### Use logicals

1.  Find all flights that:

    1. Had an arrival delay of two or more hours
    1. Flew to Houston (`IAH` or `HOU`)
    1. Were operated by United, American, or Delta
    1. Departed in summer (July, August, and September)
    1. Arrived more than two hours late, but didn't leave late
    1. Were delayed by at least an hour, but made up over 30 minutes in flight
    1. Departed between midnight and 6am (inclusive)

1.  Another useful dplyr filtering helper is `between()`. What does it do?
    Can you use it to simplify the code needed to answer the previous 
    challenges?

### Missing values

To answer these questions: read some background below.

1.  How many flights have a missing `dep_time`? What other variables are 
    missing? What might these rows represent?

1.  Why is `NA ^ 0` not missing? Why is `NA | TRUE` not missing?
    Why is `FALSE & NA` not missing? Can you figure out the general
    rule?  (`NA * 0` is a tricky counterexample!)


One important feature of R that can make comparison tricky are missing values, or `NA`s ("not availables"). `NA` represents an unknown value so missing values are "contagious": almost any operation involving an unknown value will also be unknown.

```{r}
NA > 5
10 == NA
NA + 10
NA / 2
```

The most confusing result is this one:

```{r}
NA == NA
```

It's easiest to understand why this is true with a bit more context:

```{r}
# Let x be Mary's age. We don't know how old she is.
x <- NA

# Let y be John's age. We don't know how old he is.
y <- NA

# Are John and Mary the same age?
x == y
# We don't know!
```

If you want to determine if a value is missing, use `is.na()`:

```{r}
is.na(x)
```



