[//]: ## (grip -b git_workflow.md)

# Git Workflow for COMPAS software developers

---

## Contents of this document

 Introduction

 Getting Set Up

 Day to Day Commands

 Typical Workflow

 Terminology

---

## 1. Introduction: Git & Github for COMPAS developers

### Git and Github
For those who are unfamiliar, git and github are popular tools in the software development community for sharing and collaborating on software projects. 

Git is a light-weight command line tool for maintaining different versions of software locally, and sharing those versions to remote servers. Github is a website that stores git-managed projects and enables developers to collaborate centrally on many projects. 

It is a bigger topic than we can get into here, but if you are curious you should [read more here.](https://www.atlassian.com/git/tutorials/what-is-version-control) 

Learning git is somewhat similar to learning a new language, and it can be difficult to fully grasp the vocabulary when starting out (which makes searching the internet for help significantly more challenging!). Some of the most fundamental terms are [described below](#terminology) to assist new users.

### Purpose of this document
The purpose of this document is to:
- Help COMPAS users new to git to get setup,
- Outline a consistent workflow for COMPAS developers in their day-to-day use of git, and
- Provide some of the commands that are required for this workflow 

Git is very powerful, so this is only a very small subset of the available git commands. 

This is, in some sense, a living document, meaning we are always open to [suggestions and criticism](mailto:reinhold.willcox@monash.edu) with the workflow, and seek only to find the best option for everybody. 

With that said, all developers should commit to learning the agreed upon workflow, to ensure consistency and protect against conflicts which may derail development.

### Outline of the COMPAS code repository

*Note:* If anything below doesn't make sense, try looking at the end of this document for relevant [Terminology.](#terminology)

COMPAS Users who are not developers can download the source code from the Main Repository, found at [github.com/TeamCOMPAS/COMPAS](github.com/TeamCOMPAS/COMPAS) (details can be found below). You will only need the default `master` branch and do not need to worry about what branches are. 

For developers, this repository (or 'repo') is considered "pristine", meaning that any work done here should be in a mature stage. 

The repository contains 2 permanent branches, `master` and `dev`. All other branches are either feature or hotfix branches, whose purpose is to either introduce some new functionality or fix a bug, respectively, and then be deleted. 

Feature branches on the Main Repository (also called the Main Fork or simply Main) should be ready to be tested by others. The Main Fork is not a "sandbox" for new, experimental ideas. You should [create your own fork] off of the Main Repository if you want to have public-facing experimental work.

This approach to the repository and workflow below are based on the [Feature Branch Workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) which is in common use in industry.

---

## 2. Getting Set Up: Step-by-step directions for how to configure your local and remote git repositories

### *COMPAS Users and Developers*

### Setup a Github account and git

If you have not already, go to [github.com](https://github.com/) and setup an account.

Check that you have a [working install of git.](https://www.atlassian.com/git/tutorials/install-git)

It is recommended, though not necessary, that you configure [Github with ssh](https://help.github.com/en/articles/connecting-to-github-with-ssh) as well.

### Clone the COMPAS repository to your personal computer

The COMPAS repo has a green button to clone your repository to your local computer. Alternatively, you can clone the repo with ssh (if configured, see above) or https by typing either of the following comands into a terminal window. First, go into the directory where you want to install COMPAS, e.g `cd ~/codes`, then type:

- SSH: `git clone git@github.com:TeamCOMPAS/COMPAS.git`

- HTTPS: `git clone https://github.com/TeamCOMPAS/COMPAS.git`

### Confirm that it worked

To test that cloning worked, run the following two commands:

`cd COMPAS`

`git branch`

If the clone finished without error, you should see as output: 

`* master`

At this point, if you do not plan to do any COMPAS development, you're all set. See [getting_started.md](getting_started.md) to see how to compile and run COMPAS. If you run into issues or would like to see new features implemented, you can [contact us here.](compas-user@googlegroups.com). You should read on if you are curious, but if you are not invited to be a collaborator, you will only have read-access to the repository.


### *COMPAS Developers Only*

*Note:* This section is very technical. Take a look at the section below on [Terminology.](#terminology) if you get stuck!

### Join as a collaborator 

In order to contribute to COMPAS, you will need to be added as a collaborator. Non-collaborators have read-only access to all of the branches.

[Contact us here](compas-dev@googlegroups.com) to inquire about collaborating, or reach out to one of us directly (see the [COMPAS homepage](https://compas.science/) for an up-to-date list).

### Fork the main repo.

As a COMPAS developer, you are highly encouraged to create your own personal fork of the Main repo. This serves as a public-facing 'sandbox' of your current work, where you can share partially-developed ideas and projects with others who might be interested in assisting. 

On Github, go to the TeamCOMPAS/COMPAS repo and click on `Fork` in the upper-right corner. This will create a copy of the current state of the TeamCOMPAS/COMPAS repo, including all branches and all commit histories, and place it in your profile as your personal Fork, identified as <your-username>/COMPAS.

Since this is your personal repo, you can be as organized or scatter-brained as you wish here. If you work best with 50 branches, obscure names, and code scraps everywhere, have at it. You can also give or take away access to any other collaborators who you might wish to contribute. Note that for public repositories, your code will still be read-only for everyone who is not a collaborator. 

### Clone from your remote fork to your local repo

Once your fork is created, you'll want to connect it to your local repository. In the terminal, navigate to your COMPAS git repo and type:

`git remote add <fork-nickname> <remote-fork-url>`

The <remote-fork-url> can be found on your remote repo under the same green 'Clone or Download' button as before. If you have ssh configured, it will be similar to `git@github.com:reinhold-willcox/COMPAS.git`. The <fork-nickname> is your choice, but should be informative, e.g `reinhold_fork`.

---

## 3. Day to Day commands

### Basic commands for navigating local git 

Branches allow a developer to experiment with multiple new features simultaneously on the same code-base. In git, branches are very lightweight and easy to manage, making them incredibly useful.

To view, create, and switch branches, use: (similar to `ls`, `mkdir`, and `cd`)

`git branch` 

`git checkout -b <newbranch>`

`git checkout <branchname>`

*Note:* Many git commands require that you are on the correct branch before executing the command - using these 3 commands regularly before running more complicated commands will save you headaches down the road!

### Committing changes

What git does best is to record all the small changes and edits that accumulate as we modify code. After many small changes, you might have a feature that you decide isn't actually what you want, and you want to get rid of it. Or you might have introduced a bug at some point that spans many files, and you need to remove it without undoing all the work you've accomplished since then. Git makes this incredibly easy by storing small edits as "commits". Commits, like branches, are incredibly versatile and powerful, but can be conceptually tricky to grasp at first. 

Committing is the process of adding a selection of changes to the history of your branch. It is effectively saving your work, and should be done often (every time any small fix has been made). To perform a commit, you first need to add the relevant files to your "index", then submit the commit with a commit message. The message should describe every change you made in some detail, so that in the event that we decide to undo (or "revert") a previous commit, we can identify exactly when the mistake occured.

`git branch` (check that you're on the correct branch)

`git add <file1> <file2> <...>` (whatever files you've just edited)

`git commit -m "really clear message indicating all the changes you made in this commit."`

*Note:* A single commit should capture an entire "fix" of one kind. If, for example, you want to add a function to a C file and it's header, and you also want to update the internal contents of a completely different function in the same C file, you should do 2 commits. First, make the edits to the first function and header, then `git add file.C file.h`, `git commit -m "created function myFunction to do someStuff and added it to the header file"`. Then make the second edits for the contents of the existing function, and run `git add file.C`, `git commit -m "updated internal contents of thisOtherFunction to allow for specificUseCase"`. 

You can undo a `git add` before you have done a `git commit` with:

`git reset <file>` 

- You can also use `git commit --amend` to substitute the previous commit with a new one:

```
git commit -m "first commit"
git add <fogotten_file_name>
git commit --amend
```

will open an editor and make it possible to modify the commit message. 

### Status of the Index

We can check the status of our git repository with:

`git status`

The printout is self explanatory and tells you which files have been added and which ones we need to add before committing


### F. Deleting branches 

You should be comfortable deleting branches, or else your repos might pile up with old branches that are no longer active. Branches are also very easy to manage in git (relative to other version control systems), so you should practice creating new branches, making quick edits, testing and updating, and deleting again without worry. To delete a branch, first navigate to any other branch, then:

`git branch -d <branch-name>`

- If that throws an error, likely there were some uncommited changes (work that would be completely lost if the branch gets deleted). Either commit the branch before deleting, or if you decided against all the changes, you can force the delete with:

`git branch -D <branch-name>`

### E. Fetch other branches from a remote

If you followed the above workflow, you can verify that the COMPAS repo is a designated remote fork in your local repo, nicknamed `origin`. 

`git remote -v` should output

```
origin	git@github.com:TeamCOMPAS/COMPAS.git (fetch)
origin	git@github.com:TeamCOMPAS/COMPAS.git (push)
```

To see all of the other branches on this fork:

`git branch -a` should output something similar to

```
* master
  remotes/origin/HEAD -> origin/master
  remotes/origin/dev
  remotes/origin/hotfix-input-file-bug  
  remotes/origin/master
  remotes/origin/release
```

All of the branches found under `remotes/origin/` are available to be copied locally with:

`git checkout -b <local-branch-name> origin/<remote-branch-name>`





To view which remote repositories you are tracking, and all the branches on those remotes:

`git remote -v` 

`git branch -r` 



""" Probably not useful now
### E. Create new branches locally 
- To create a new branch off of the current one (Note: you need to switch to the correct base branch first for this to work):
`git checkout -b new-branch`
"""



- To set the upstream parent (tells git from what branch to git pull, git fetch etc.)

`git push --set-upstream-to <remote-repo> <name-of-branch-on-remote>` 

an example: git push --set-upstream-to dev new-branch










---

## 3. Typical Workflow

### New project idea 

When beginning a new project, you will typically want to branch off of the most updated version of the `dev` branch

```git checkout dev 
git pull
git checkout -b <new-project>
``` 


- Pushing to your remote repository is a way to save all of your commits (i.e the history of edits) somewhere off your local computer. This is good practice because it acts as a backup in the event something happens to your local machine, and it also allows other collaborators to see your work (without having to give them access to your personal device). This should also be done often, but does not need to be as frequent as commits. A good rule of thumb is to push any updated branches at the end of every day. 

`git branch -vv` (check that you're on the correct branch, and that you are pointed to the correct remote branch)

`git push`

- Another good rule of thumb is to pull every morning (at least on branches that might have multiple contributors)!

`git branch -vv` (check that you're on the correct branch, and that you are pointed to the correct remote branch)

`git pull`

### Work collaboratively across forks 

- The purpose of forks is to give you a place to work on projects privately and without stepping on the toes of other people (or letting anyone else step on your toes!). However, there will likely still be times when you want to check out a branch that someone else worked on. This may happen before the branch is polished enough to be sent to the Main repo, but they would still like feedback or edits. In this case, the other developer needs to add you as a collaborator on their personal repo (this is done on Github, in the settings menu of their repo). Once you are added, you can set up a local branch to point upstream to their remote branch.

`git remote add <collaborator_name-repo> <url>` where \<url\> is the https or ssh url you can copy from Github, and \<name\> is your nickname for that fork, e.g `reinhold_fork`. 

`git checkout dev`

`git checkout -b <collaborator_name-project1>`

`git branch --set-upstream-to=<collaborator_name-repo>`

-- Note: now you can push and pull to that remote branch (assuming the other developer has given you permission) to see what work they have done as well.

### Mature projects - adding to the Main Repo

- When a project is nearing completion (e.g when the code is nearly ready to be joined into the Main Repository), the author of the branch should do a final push to their personal remote repo, then submit a pull-request onto the Main Repo. This branch should cover only the scope of the named project and should not include work on any other project or bug fixes. 

- In order to merge your remote branch into a branch on COMPAS, you will first need to create a new branch on the Main Repo (the reason for this is that we decided to keep all branches off of the Main Repo until they were mature, but you need to have an existing branch to receive a pull request). 

    - 1. On the Main Repo, click the Code tab, click the Branch dropdown and select Dev, click the Branch dropdown again and into the Find or Create a Branch textbox, and type in the name of the new branch. In most situations, this will be the name of your mature branch (but you should ensure that it is clear what edits were made in this branch).

    - 2. If you did the initial fork properly, there should be a clear yellow box right in the center of your repo's main page with recently pushed branches, and a button to "Compare and pull request". Click the button to see the Open a Pull Request GUI. 

    - 3. The title of the pull request should be the name of the branch you are merging in, **only if this is the first pull-request**. If a pull-request is rejected and resubmitted, or if it is accepted but additions are made to it later on, the newer titles should broadly reflect the changes.

    - 4. You should leave a detailed description of the branch in the comment section, and if you would like to request any Reviewers, feel free to do so. You can leave Assignees, Labels, Projects, and Milestones blank for now.

    - 5. Once you have created the pull request, it is up to the other team members to review it (see below). 

---

## 3. Typical Workflow


### Overview 
First, any and all changes to the main repo are done through pull-requests, not pushes (a pull-request really represents a push, with the added step of confirmation from a third party that the push meets certain standards). This is to ensure that the main repo is "clean" and contains code that has at least been somewhat tested. Branches will first be added into the main repo from the remote forks of a given collaborator, then (after further testing) added to the Dev branch. Finally, when sufficient sub-branches have been added to the Dev branch to warrant a new version release, the Dev branch will be added to the Master branch and published to the public. Ultimately, we should always be thinking of the next Version Release and which projects/concepts we would like to include in it, and define our timelines around those expectations.  

### Q&A Protocol
A mature branch will have to go through 3 rounds of Q&A review before it is ready for deployment in a new COMPAS version. The Primary Review occurs when a pull-request is sent to have the branch included as a sub-branch in the main repo. The Secondary Review occurs when a pull-request is sent to have that sub-branch merged into the Dev branch on the main repo. The Final Review occurs when a pull-request is sent to have the Dev branch (which might contain several sub-branches) merged into the Master branch. **Testing and validation of a branch should never be done by someone who worked on it extensively** (we can decide on a case-by-case basis who falls into that category). 

1. Primary Review: To be included as a sub-branch in the main repo, the branch must be reviewed by one person in the broader COMPAS collaboration. It must compile and run, producing somewhat sensible output (e.g output files should not be empty, but at this stage they may contain data which is "wrong", or physically inconsistent) 

2. Secondary Review: To be merged into the Dev branch on the main repo, the branch must then be reviewed by two people in the broader COMPAS collaboration, one of whom is a core developer. Additionally, one of the reviewers should focus on the consistency of the code layout and structure with the current code, while the other reviewer should focus primarily on the accuracy of the physics and its implementation, although ultimately both reviewers are responsible for ensuring that the branch is clearly and correctly written. This is the primary review step, with the expectation that code that is accepted to be merged into Dev is ready for the public. 

3. Final Review: When Dev is updated with a new branch, all collaborators should be testing it, and in particular should ensure that all of the different projects which have been added to Dev work with each other. When Dev is nearly ready to be merged into the Master branch, the core collaborators will host a COMPAS Powwow to discuss and assign any final revisions/tests, and set a firm release date. Other collaborators are invited to join (over Zoom or in person), details decided on a case-by-case basis.

---

## Quick Fixes

- This workflow should work well for larger projects and features which need oversight and testing before merging into the Main Repo. However, there will inevitibly be times when only a simple bux-fix or last-minute addition is warranted. In these situations, we will use a stream-lined workflow for merging with the Main Repo. 

    - 1. Determine the seriousness of the fix: is this a bug in the published Master branch that we need to fix urgently, or is there simply a typo in one of the sub-branches? Once this is decided, go onto the relevant branch on the Main Repo, and create a new branch (see Step 1 under Mature Projects, above). Give the branch a descriptive name that starts with "hotfix-", e.g "hotfix-typo_in_speed_of_light". 

    - 2. Create a branch locally (from an up-to-date branch, e.g the local master), and set the new hotfix branch as the upstream of the local branch (see Work Collaboratively Across Forks, above). 

	- 3. Make the edits on the local branch

	- 4. `git push`, to merge your edits into the remote branch

	- 5. Submit a pull request on the Main repo from the hotfix branch onto it's parent

	- 6. Notify another collaborator to have them quickly review your hotfix and accept the pull request

-- Note: This will only work for branches that start with "hotfix-"

---

## Terminology

- **Branch**: Branches in git separate work-streams for different features (e.g front-end developers might have a branch for a fancy new button for their website, while back-end developers might have a branch to make database-entry easier). In our case, branches will distinguish different projects or concepts (e.g Supernova-Kicks, White-Dwarf-Accretion, etc.). Branches should _not_ be used to distinguish developers. As mentioned previously, only the Master and Dev branches are permanent, and any new branches should be created with the intention of contributing some new feature or physics, and being [deleted once that is done.](https://rickandmorty.fandom.com/wiki/Mr._Meeseeks)

- **Repository**: A Repository (or Repo) is a collection of all the different branches of a given project which are kept in the same location. To be specific, a location might be your local computer or part of a remote server. A Repo can be public (often called Open Source) or private, with a selected list of collaborators who have read and possibly write access. A single github user may have many Repos for all of their different projects, all of which might have any number of branches. Note: COMPAS development will be done in a private repo, but there will be a second public repo which will hold only a copy of the Master branch for the public to download.

- **Fork**: A Fork is a full copy of a Repo, including all its branches, to another location. Most of the time, "another location" will mean elsewhere on the github servers, since we will be Forking from the Main Repo to our Personal Repo when we are setting up. In our case, Forks will distinguish different users, or perhaps groups of users (e.g Copenhagen/COMPAS). All core developers should have a personal fork. If you are familiar with the `git clone` command, this is identical to Forking from a remote server onto your own personal computer. 

- **Local/Remote**: Local refers to what's on your personal computer, while Remote refers to anything that isn't. Github will be remote for everyone. My local computer is only local to me, and would be considered remote to anyone else. This should not come up often, because it would be very foolish to give anyone access to my local computer, even if it's a developer I trust. The purpose of your personal remote fork is to be a proxy for your local fork, where you can add things you've worked on that you wish to share around.

- **Tracking**: The word tracking refers to either remote repositories, or local files in a single branch, but they have slightly different meanings. 
    - A tracked repository is one which contains a branch which is "upstream" of a branch in your local repository. By default, all the branches on a cloned or forked repository track the branches they were cloned from, and new branches track their (local) parent branches. You can modify the upstream branch of a given branch to point at any other branch you like, whether local or remote, and can even track multiple remote repositories (see below). This is useful if you want to check out a branch that sits on a colleague's Fork. You can view all tracked repositories with `git remote -v`
    - A tracked file is one in your current working branch which has been either added or previously commited. If you make an edit to a file, it will be untracked until you add it. If a file (or collection of files) is in your .gitignore, it is considered neither tracked nor untracked, but ignored. If you try to change branches while you have untracked files, you will run into errors.

- **Origin**: Origin is a shorthand for the most relevant remote repository. It is configured by default whenever you clone from a repository, so yours will (and usually should) point to your personal remote fork. If you track multiple remote forks, you should give them all helpful, distinguishing names (e.g jeff_fork, main_compas_repo)

- **Commit**: A commit is a complete set of edits to one or more files that you want added to the history of the branch. The branch itself is just made up of the sequence of commits throughout it's history, so committing to the branch just adds one more commit onto the end of it. Commits can be reverted if it is decided later on that they are not desired. Commits also include a message which should be very detailed. 

- **Index**: The Index (aka Staging Area) exists only in the intermediate step between editing local files and committing those files. Historically, other Version Control systems only allowed editting files, and then committing those files one by one. The issue with that is that sometimes a collection of edits of different files logically make up one full "commit-worthy-edit". The classic example of this is adding a function to a .C file and it's header .h file. If you need to revert this commit back for any reason, it makes sense to remove both of those edits at once - you would virtually never need to remove the function from the C file but leave it in the header. Adding files to the index is the way to collect all of the files that were involved in a given series of edits that you want to treat as one big Edit. 

- **Push, Pull, and Pull Request**: These commands form the backbone of file-sharing across repositories. They all pretty much cover the same action, to get a given branch on _this_ repo over to _that_ repo. Whether you use `push` or `pull` really depends on where you are in relation to the branch. If the branch is on your local repo, then you will want to `push` it somewhere else. If it is on a remote repository, you will want to `pull` it down. The confusing one, `pull request` is used when there is an extra step of validation on the receiver side. The donor is not allowed to simply `push` some possibly broken code, so they have to request that the receiver pull it, hence `pull request`. For our setup, we will almost exclusively require pull-requests onto Main, to maintain a high-quality of branches. 

- **Revert**: A revert is used when the chain of commits that make up a branch has gone too far - you have decided that you don't like the latest edits and we want to remove them from the branch. In this case, you revert the HEAD of the branch (the latest commit) to an earlier commit, identified by it's unique SHA hash. This can get quite complicated though, so make sure to use this one with caution, and do lots of testing before you try anything. 

- **Working Directory**
