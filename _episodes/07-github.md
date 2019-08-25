---
title: "Version control with git and Github"
teaching: 45
exercises: 45 
questions:
- "What is version control? How do I use it?"
- "What is the difference between `git`and Github?"
- "What benefits does a version control system brings in for my research?"
objectives:
- "Understand the benefits of using a version control system such as `git`."
- "Understand the basics of `git` and its usage in RStudio."    
keypoints:
- "`git` and Github allow you to version control files and go back in time if needed."
- "In a version control system, file names do not reflect their versions."
- "An RStudio project folder can be fully version controlled and synchronized online with Github."
- "Working locally in RStudio with a synchronised online folder will make your work more stable and understandable for you and others."
---

# Introduction

## Why should scientists use git and Github?

* Ends (or, nearly ends) the horror of keeping track of versions.
  Basically, we get away from this: 

![](../img/MessySaves.png)

![](../img/phdcomics_version_control.png)

When you open your repository, you only see the most recent version.  But, it easy to compare versions, and you can easily revert to previous versions. 

* Improves collaborative efforts.  Different researchers can work on the same files at the same time!
* It is easy to share and distribute files through the Github website.
* Your files are available anywhere, you just need internet connection!  

## git and Github

We will learn about version control using `git` and [GitHub](https://en.wikipedia.org/wiki/GitHub), and we will interface with this through RStudio. This will change your scientific life (for the better!). Github was developed for social coding (i.e., sort of like an open source Wikipedia for programmers). Consequently, much of the functionality and terminology of Github (e.g., branches and pull requests) will not be relevant for most scientists. Therefore, we will skip over all this stuff!    

Github will facilitate your daily coding life when working with your most important collaborator: **you** ! A famous quote that we like to emphasize:
> Your past self from 6 months ago is gone and won't answer emails from your present self!  

**git:**   
The version control program `git` will track and version your files locally on your machine. `git` is locally executed and works on your local machine. It was named created and named by *Linus Torvalds*, the creator of Linux. Torvalds sarcastically quipped about the name _git_ (which means unpleasant person in British English slang): "I'm an egotistical bastard, and I name all my projects after myself. First 'Linux', now 'git'."

`git` is a version control system that lets you track changes to files over time. These files can be any kind of file (e.g. .doc, .pdf, .xls), but free text differences are visible and can be read by humans (eg txt, csv, md). 

**Github:**  
[GitHub](https://github.com/) is a website for storing your git versioned files remotely. It has many nice features to be able visualize differences between [images](https://help.github.com/articles/rendering-and-diffing-images/), [rendering](https://help.github.com/articles/mapping-geojson-files-on-github/) & [diffing](https://github.com/blog/1772-diffable-more-customizable-maps) map data files, [render text data files](https://help.github.com/articles/rendering-csv-and-tsv-data/), and [track changes in text](https://help.github.com/articles/rendering-differences-in-prose-documents/).

> If you are a student you can get the micro account which includes 5 private repositories for free (normally a $7/month value).  You can sign up for the student account [here](https://education.github.com/pack).  Instructors can also request a free organization [account, "Request a discount"](https://education.github.com/). These concepts are more important for coders who want the entire coding community (and not just people working on the same project) to be able to suggest changes to their code.  This isn't how most scientists will use Github. To get the full functionality of Github, you will eventually want to learn other concepts. But, this can wait.  

**git and Github**:  
Although `git` and GitHub are two different things, distinct from each other, I think of them as a bundle since I always use them together. It also helped me to think of GitHub like Dropbox: you make folders that are 'tracked' and can be synced to the cloud. GitHub does this too, but you have to be more deliberate about when syncs are made. This is because GitHub saves these as different versions, with information about who contributed when, line-by-line. This makes collaboration easier, and it allows you to roll-back to different versions or contribute to others' work.

<figure>
    <img src="../img/octocat_GitHub_mascot.png" alt='Github Mascot' width="250" />
    <figcaption><center>The octocat: the official mascot of Github</center></figcaption>
</figure>

## Resources

These materials borrow from: 

- [The `git` version control system wikipedia page](https://en.wikipedia.org/wiki/Git)
- [Github wikipedia page](https://en.wikipedia.org/wiki/GitHub)
- Jenny Bryan's lectures from STAT545 at UBC: [The Shell](http://stat545.com/git09_shell.html)
- Jenny Bryan's [Happy git with R](http://happygitwithr.com) tutorial
- Melanie Frazier's [GitHub Quickstart](https://rawgit.com/nazrug/Quickstart/master/GithubQuickstart.html)
- Ben Best's [Software Carpentry at UCSB](http://remi-daigle.github.io/2016-04-15-UCSB/git/)

Today, we'll only introduce the features and terminology that scientists need to learn to begin managing their projects.

## Outline

Today, we'll interface with GitHub from our local computers using RStudio. There are many other ways to interact with GitHub, including GitHub's Desktop App or the command line ([here is Jenny Bryan's list of git clients](http://stat545.com/git02_git-clients.html)), but today we are going to work from RStudio. You have the largest suite of options if you interface through the command line, but the most common things you'll do can be done through one of these other applications (i.e. RStudio and the GitHub Desktop App).

Here's what we'll do after we set up `git` on your computers: 

1. create a repository on Github.com
2. clone locally using RStudio. 
3. learn the RStudio-GitHub workflow by syncing to Github.com: `pull`, `stage`, `commit`, `push`.
4. explore github.com: files, commit history, file history.
5. practice the RStudio-GitHub workflow by editing and adding files. 
6. practice R Markdown.


### Some Github terminology

* **User**: A Github account for you (e.g., jules32).
* **Organization**: The Github account for one or more user (e.g., datacarpentry).
* **Repository**: A folder within the organization that includes files dedicated to a project.
* **Local Github**: Copies of Github files located your computer.
* **Remote Github**: Github files located on the https://github.com website.
* **Clone**: Process of making a local copy of a remote Github repository.  This only needs to be done once (unless you mess up your local copy).
* **Pull**: Copy changes on the remote Github repository to your local Github repository.  This is useful if multiple people are making changes to a repository.
* **Push**: Save local changes to remote Github
<br />
<br />

![](../img/push_pull_clone.png)
<br />
<br />

## Setup Git & GitHub

We're going to switch gears from R for a moment and set up Git and GitHub, which we will be using along with R and RStudio for the rest of the workshop. This set up is a one-time thing! You will only have to do this once per computer. We'll walk through this together. 

1. Create **Github** account at <http://github.com>, if you don't already have one. For username, I recommend all lower-case letters, short as you can. I recommend using your academic email (e.g. *.uva.nl*), since you can request free private repositories via [GitHub Education](https://education.github.com/) discount.

2. Configure **git** with global commands, which means it will apply 'globally' to all files on your computer, rather than to a specific folder. Open the Git Bash program (Windows) or the Terminal (Mac) and type the following:

~~~
        # display your version of git
        git --version
        
        # replace USER with your Github user account
        git config --global user.name YOUR_USER_NAME
        
        # replace NAME@EMAIL.EDU with the email you used to register with Github
        git config --global user.email YOUR_NAME@@EMAIL.UNIVERSITY
        
        # list your config to confirm user.* variables set
        git config --list
~~~
{: .language-bash}

Not only have you just set up git as a one-time-only thing, you have just used the command line. We don't have time to learn much of the command line today, but you just successfully used it following explicit instructions, which is huge! There are great resources for learning the command line from the [Software Carpentry Shell novice lesson](swcarpentry.github.io/shell-novice). 

### Troubleshooting

If you have problems setting up git, please see the [Troubleshooting section](http://happygitwithr.com/troubleshooting.html) in Jenny Bryan's amazing [HappyGitWithR](http://happygitwithr.com). 

#### New(ish) Error on a Mac
We've also seen the following errors from RStudio: 

~~~
error key does not contain a section --global terminal
~~~
{: .language-bash}

and
~~~
fatal: not in a git directory
~~~
{: .language-bash}

To solve this, go to the Terminal and type:
~~~
which git
~~~
{: .language-bash}

<img src="../img/git_whichgit.png" width="250px">

  
Look at the filepath that is returned. Does it say anything to do with Apple?

-> If yes, then the [Git you downloaded](https://git-scm.com/downloads) isn't installed, please redownload if necessary, and follow instructions to install.  

-> If no, (in the example image, the filepath does not say anything with Apple) then proceed below:

In RStudio, navigate to: Tools > Global Options > Git/SVN. 

<img src="../img/git_options.png" width="250px">


<br>

Does the **“Git executable”** filepath match what the url in Terminal says? 

<br>

<img src="../img/git_options_filepath.png" width="500px">


If not, click the browse button and navigate there.   

>*Note*: on my laptop, even though I navigated to /usr/local/bin/git, it then automatically redirect because /usr/local/bin/git was an alias on my computer. That is fine. Click OK.

Quit RStudio.   

Then relaunch RStudio.  

Try syncing or cloning, and if that works and then you don’t need to worry about typing into the Terminal, you’re all done!

# Synchronising changes from RStudio to Github

## Create a repository on Github.com

First, go to your account on github.com and click "New repository".

<img src="../img/create_repository.png" width="900px"> 


Name it `my-repo`, short for "my-repository" (any short self-describing name would be good).   

Also, add a description, make it public, create a README file, and create your repo!

<img src="../img/create_repository_2.png" width="900px"> 

The *Add gitignore* option adds a document where you can identify files or file-types you want Github to ignore. These files will stay in on the local `git` folder (the one on your computer), but will not be uploaded onto the web version of Github.

The *Add a license* option adds a license that describes how other people can use your Github files (e.g., open source, but no one can profit from them, etc.).  We won't worry about this today.

Check out our new repository!  

Notice how the README.md file we created is automatically displayed at the bottom. The .md means that it is Markdown (remember how .Rmd was RMarkdown?) so the formatting we learned in the last lesson apply here.

<img src="../img/new_repository.png" width="900px"> 

**From here, you will work locally (on your computer).**


## Clone your repository using RStudio

We'll start of by cloning to our local computer using RStudio. We are going to be cloning a copy of our Remote repository on Github.com to our local computers. Unlike downloading, cloning keeps all the version control and user information bundled with the files. 

**Step 0**: Create your `github` folder 

This is really important! We need to be organized and deliberate about where we want to keep all of our GitHub repositories (since this is the first of many in your career). 

Let's all make a folder called `github` (all lowercase!) in our home directories. So it will look like this: 

- Windows: `Users\[User]\Documents\github\`
- Mac: `Users/[User]/github/`

This will let us take advantage of something that is really key about GitHub.com: you can easily navigate through folders within repositories and the urls reflect this navigation. The greatness of this will be evident soon. So let's set ourselves up for easily translating (and remembering) those navigation paths by having a folder called `github` that will serve as our local mirror of the repositories on 'github.com'.

So really. Make sure that you have an all-lowercase folder called `github` in your home directory!!

**Step 1**: Copy the web address of the repository you want to clone.

<img src="../img/clone_step1.png" width="900px"> 

**Step 2**: from RStudio, go to New Project (also in the File menu).

<img src="../img/new_project_1.png" width="600px"> 

**Step 3**: Select Version Control

<img src="../img/new_project_2.png" width="600px"> 

**Step 4**: Select Git
    
<img src="../img/new_project_3.png" width="600px"> 


**Step 5**: Paste it in the Repository URL field, and type **tab** to autofill the Project Directory name. Make sure you keep the Project Directory Name **THE SAME** as the repository name from the URL.

Save it in your github folder (click on Browse) to do this. 

<img src="../img/new_project_4.png" width="600px"> 

If everything went well, the repository will be added to the list located here:

<img src="../img/select_project.png" width="700px"> 

And the repository will be saved to the Github folder on your computer:

<img src="../img/cloned_repository.png" width="900px"> 

<font size="+2">Ta da!!!</font> The folder doesn't contain much of interest, but we are going to change that.

## Inspect your repository

Notice a few things in our repo here: 

1. Our working directory is set to `~/github/my-repo`. This means that I can start working with the files I have in here without setting the filepath. This is that when we cloned this from RStudio, it created an RStudio project, which you can tell because: 
    - `.RProj` file, which you can see in the Files pane. 
    - The project is named in the top right hand corner
1. We have a git tab! This is how we will interface directly to Github.com

<img src="../img/RStudio_IDE_git.png" width="900px"> 

When you first clone a repo through RStudio, RStudio will add an `.Rproj` file to your repo. And if you didn't add a `.gitignore` file when you originally created the repo on GitHub.com, RStudio will also add this for you. These will show up with little yellow `?` icons in your git tab. This is `git` way of saying: "I am responsible for tracking everything that happens in this repo, but I haven't seen these files yet. Do you want me to track them too?" (We'll see that when you click the box to stage them, they will turn into `A`s because they have been added to the repo. 

# A typical workflow: add, commit and push 

## Add files to our local repo

The repository will contain:

* one `.gitignore` file.
* one `README.md` file.
* one `Rproj`file.

And, I typically create the following:

* folders for `data/` and `figures/`.  
* R scripts.
* etc.

I'm going to go to the Finder (Windows Explorer on a PC) and copy a file into my repository from there. And then I'm going to go back to RStudio -- it shows up in the git tab! So the repository is being tracked, no matter how you make changes to it (changes do not have to be done only through RStudio). 

To make changes to the repository, you will work from your computer in the local `git` folder (mirror of the online Github `my-repo`). 

When files are changed in the local repository, these changes will be reflected in the git tab of RStudio:

<img src="../img/modify_files.png" width="900px"> 

## Inspect what has changed

These are the codes RStudio uses to describe how the files are changed, (from the RStudio [cheatsheet](http://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf)):

![](../img/modified.png)

## Sync from RStudio to GitHub


When you are ready to commit your changes, you follow these steps:

![](../img/commit_overview.png)

We walk through this process below:

### Pull 
From the Git tab, "Pull" the repository.  This makes sure your local repository is synced with the remote repository.  This is very important if other people are making changes to the repository or if you are working from multiple computers.

<img src="../img/pull.png" width="700px"> 

### Stage
Stage the files you want to commit.  In RStudio, this involves checking the "Staged" boxes:

<img src="../img/staged.png" width="700px"> 

### Commit

<img src="../img/commit.png" width="700px"> 

### Push

<img src="../img/push.png" width="700px"> 

## Explore remote Github
The files you added should be on github.com:

<img src="../img/Github_remote.png" width="700px"> 

Let's also explore commit history, file history.

# Your turn!

> ## Exercise: update the README file
>
>  1. Open your README file by clicking on it in the `Files` panel (lower right corner).
>  2. Write a few line of text and save it.
>  3. See what happens in your Git tab.
>  4. Sync it to your remote repository (`my-repo` on Github.com).
{: .challenge}

> ## Exercise: add a new file
>
>  1. Open your file explorer (Finder on Mac, Explorer on Windows).
>  2. Copy-paste a file into your local `my-repo` git folder.
>  3. Go back to RStudio.
>  4. Confirm that `git` can see this file.
>  5. Add/stage this file.
>  6. Commit this new file with a commit message. 
>  7. Sync it to your remote repository (`my-repo` on Github.com).
>
> > ## Solution
> > * 4) `git` sees the file as untracked and displays a question mark icon: **?** 
> > * 5) In the Git panel of RStudio, click on the **Staged** tick box. The status of the file should change. 
> > * 6) Click on the **Commit** button and write a small descriptive message. It should be the end of the sentence: "with this commit, ...". 
> > 
> {: .solution}
{: .challenge}  
 
Remember, `git` will track anything within that folder (the way Dropbox does), it's not specific to RStudio!

## Committing - how often? Tracking changes in your files

Whenever you make changes to the files in Github, you will walk through the Pull -> Stage -> Commit -> Push steps.

I tend to do this every time I finish a task (basically when I start getting nervous that I will lose my work).  Once something is committed, it is very difficult to lose it.

One thing that I love about about Github is that it is easy to see how files have changed over time.  Usually I compare commits through github.com:

![](../img/commit_history.png)
<br />
<br />
![](../img/commit_compare_2.png)
<br />
<br />

You can click on the commits to see how the files changed from the previous commit:
<br />
<br />

![](../img/commit_compare_3.png)
<br />
<br />


## Happy Git with R 

Jenny Bryan's [HappyGitWithR](http://happygitwithr.com) is very useful for troubleshooting, particularly the sections on [Detect Git from RStudio](http://happygitwithr.com/rstudio-see-git.html) and [RStudio, Git, GitHub Hell (troubleshooting)](http://happygitwithr.com/troubleshooting.html).


