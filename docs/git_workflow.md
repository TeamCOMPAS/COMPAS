[//]: ## (grip -b git_workflow.md)

# Git Workflow for COMPAS software developers

## Introduction 

### 1. Git and Github
For those who are unfamiliar, git and github are popular tools in the software development community for sharing and collaborating on software projects. Git is a light-weight command line tool for maintaining different versions of software locally, and distributing those versions to remote servers. Github is a website that centralizes for storage of git-managed projects. It is a big topic, and not worth getting into here, but if you are curious you can [read more](https://www.atlassian.com/git/tutorials/what-is-version-control).

### 2. Purpose
The purpose of this document is to outline a consistent workflow for COMPAS developers in their day-to-day use of git, protocols for whenever new projects are started or completed, and the commands that are required for this workflow (git is very powerful, so this is only a very small subset of the availalbe git commands). This is, in some sense, a living document, meaning we are always open to [suggestions and criticism](mailto:reinhold.willcox@monash.edu) with the workflow, and seek only to find the best option for everybody. With that said, everyone should commit to learning the agreed upon workflow, to ensure consistency between developers and protect against user-error which may derail development.

### 3. Outline

- The workflow here is based on the [Feature Branch Workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) in common use in industry, in which there are 2 permanent branches, Master and Dev. All other branches are considered "sub-branches", based on projects and/or concepts, whose purpose is to add a feature and then be deleted.

- The Main Repository (also called the Main Fork, or simply Main) is considered "pristine", and should only contain the Master and Dev branches, and any sub-branches that are nearly ready to be made public (more on this later). 

- Developers should each have their own personal Forks of the Main Repository, and any work done locally (on their personal computers) should be pushed up through their personal Repository before being added to Main (more on this later). 

### 4. Sections of this Document

- Getting Set Up: Step-by-step directions for how to configure your personal local and remote git repos

- Lifetime of a Project: Walkthrough of creating, committing, pushing/pulling, and setting parents of branches

- Test and Release New Versions: **Important** The specific COMPAS workflow around merges, pull-requests, and updates to the Main Git Repo

- Terminology: Important git keywords that may be unfamiliar

---

## Getting Set Up 

### 1. Join as a collaborator 

- If you have not already, go to [github.com](https://github.com/) and setup an account.

- The Repo is found at [TeamCOMPAS/COMPAS](https://github.com/TeamCOMPAS/COMPAS), however this will be hidden to you if you are not yet a Collaborator

- Reach out to a member of the core developer team (see [COMPAS homepage](https://compas.science/) for an up-to-date list) to request Collaborator access.

- You will receive an email notifying you that you have been added, at which point you can browse the Repo.

### 2. Fork the Main Repo into your Personal Remote Repo

- While logged in, go to the TeamCOMPAS/COMPAS repo and click on `Fork` in the upper-right corner to create your personal Fork (found at user-name/COMPAS)

- This is your personal, private Repo. You can be as organized or scatter-brained as you wish here. If you work best with 50 branches, all nested within each other, have at it. You can also give or take away access to any other collaborators who you might wish to share your work with, though you should not make your Repo public (only the Master branch should be public).

### 3. Clone from your remote repo locally

- If you have not yet configured [Github with ssh](https://help.github.com/en/articles/connecting-to-github-with-ssh), you can clone over http: 

`git clone https://github.com/user-name/COMPAS.git`

- With ssh configured, you can clone with: 

`git clone git@github.com:user-name/COMPAS.git`

### 4. Basic commands for navigating local git 

- To view and switch to available local branches, (equivalent of `ls` and `cd`) (Note: many git commands require that you are on the correct branch before executing the command - using these 2 commands regularly will save you headaches down the road): 

`git branch` 

`git checkout <branchname>`

- To view which remote repositories you are tracking, and all the branches on those remotes:

`git remote -v` 

`git branch -r` 

### 5. Create new branches locally 

- To create a new branch off of the current one (Note: you need to switch to the correct base branch first for this to work):

`git checkout -b new-branch`

- To set the upstream parent 

`git push --set-upstream <remote-repo> <name-of-branch-on-remote>` 

### 6. Delete branches locally

- You should be comfortable deleting branches, or else your repos might pile up with old branches that are no longer active. Branches are also very easy to manage in git (relative to other version control systems), so you should practice creating new branches, making quick edits, testing and updating, and deleting again without worry. To delete a branch, first navigate to any other branch, then:

`git branch -d <branch-name>`

- If that throws an error, likely there were some uncommited changes (work that would be completely lost if the branch gets deleted). Either commit the branch before deleting, or if you decided against all the changes, you can force the delete with:

`git branch -D <branch-name>`


---

## Lifetime of a project 

### New project idea 

- When beginning a new project, you should start your first branch locally, branching off of your local Dev branch (which should ideally be identical to the upstream Dev branch on the Main Repo)

`git checkout dev`

`git pull` 

`git checkout -b <new-project>`

### Commit often and push/pull 

- Committing is the process of adding a selection of changes to the history of your branch. It is effectively saving your work, and should be done often (every time any small fix has been made). To perform a commit, you first need to add the relevant files to your "index", then commit with a message. The message should describe every change you made in some detail, so that in the event that we decide we want to revert back to a previous commit, we know exactly what happened at each step.

`git branch` (check that you're on the correct branch)

`git add <file1> <file2> <...>`

`git commit -m "really clear message indicating all the changes you made in this commit."`

-- Note: A single commit should capture an entire "fix" of one kind. If, for example, you want to add a function to a C file and it's header, and you also want to update the internal contents of a completely different function in the same C file, you should do 2 commits. First, make the edits to the first function and header, then `git add file.C file.h`, `git commit -m "created function myFunction to do someStuff and added it to the header file"`. Then make the second edits for the contents of the existing function, and run `git add file.C`, `git commit -m "updated internal contents of thisOtherFunction to allow for specificUseCase"`. 

- Pushing to your remote repository is a way to save all of your commits (i.e the history of edits) somewhere off your local computer. This is good practice because it acts as a backup in the event something happens to your local machine, and it also allows other collaborators to see your work (without having to give them access to your personal device). This should also be done often, but does not need to be as frequent as commits. A good rule of thumb is to push any updated branches at the end of every day. 

`git branch -vv` (check that you're on the correct branch, and that you are pointed to the correct remote branch)

`git push`

- Another good rule of thumb is to pull every morning (at least on branches that might have multiple contributors)!

`git branch -vv` (check that you're on the correct branch, and that you are pointed to the correct remote branch)

`git pull`

### Work collaboratively across forks 

- The purpose of forks is to give you a place to work on projects privately and without stepping on the toes of other people (or letting anyone else step on your toes!). However, there will likely still be times when you want to check out a branch that someone else worked on. This may happen before the branch is polished enough to be sent to the Main repo, but they would still like feedback or edits. In this case, the other developer needs to add you as a collaborator on their personal repo (this is done on Github, in the settings menu of their repo). Once you are added, you can set up a local branch to point upstream to their remote branch.

`git remote add <collaborator_name-repo> <url>` where \<url\> is the https or ssh url you can copy from Github, and \<name\> is the shorthand for what 

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

## Test and Release New Versions 

### Overview 
First, any and all changes to the main repo are done through pull-requests, not pushes (a pull-request really represents a push, with the added step of confirmation from a third party that the push meets certain standards). This is to ensure that the main repo is "clean" and contains code that has at least been somewhat tested. Branches will first be added into the main repo from the remote forks of a given collaborator, then (after further testing) added to the Dev branch. Finally, when sufficient sub-branches have been added to the Dev branch to warrant a new version release, the Dev branch will be added to the Master branch and published to the public. Ultimately, we should always be thinking of the next Version Release and which projects/concepts we would like to include in it, and define our timelines around those expectations.  

### Q&A Protocol
A mature branch will have to through 3 rounds of Q&A review before it is ready for deployment in a new COMPAS version. The Primary Review occurs when a pull-request is sent to have the branch included as a sub-branch in the main repo. The Secondary Review occurs when a pull-request is sent to have that sub-branch merged into the Dev branch on the main repo. The Final Review occurs when a pull-request is sent to have the Dev branch (which might contain several sub-branches) merged into the Master branch. **Testing and validation of a branch should never be done by someone who worked on it extensively** (we can decide on a case-by-case basis who falls into that category). 

1. Primary Review: To be included as a sub-branch in the main repo, the branch must be reviewed by one person in the broader COMPAS collaboration. It must compile and run, producing somewhat sensible output (e.g output files should not be empty, but at this stage they may contain data which is "wrong", or physically inconsistent) 

2. Secondary Review: To be merged into the Dev branch on the main repo, the branch must then be reviewed by two people in the broader COMPAS collaboration, one of whom is a core developer. Additionally, one of the reviewers should focus on the consistency of the code layout and structure with the current code, while the other reviewer should focus primarily on the accuracy of the physics and its implementation, although ultimately both reviewers are responsible for ensuring that the branch is clearly and correctly written. This is the primary review step, with the expectation that code that is accepted to be merged into Dev is ready for the public. 

3. Final Review: When Dev is updated with a new branch, all collaborators should be testing it, and in particular should ensure that all of the different projects which have been added to Dev work with each other. When Dev is nearly ready to be merged into the Master branch, the core collaborators will host a COMPAS Powwow to discuss and assign any final revisions/tests, and set a firm release date. Other collaborators are invited to join (over Zoom or in person), details decided on a case-by-case basis.

---

## Terminology

- **Branch**: Branches in git separate work-streams for different features (e.g front-end developers might have a branch for a fancy new button for their website, while back-end developers might have a branch to make database-entry easier). In our case, branches will distinguish different projects or concepts (e.g Supernova-Kicks, White-Dwarf-Accretion, etc.). Branches should _not_ be used to distinguish developers. As mentioned previously, only the Master and Dev branches are permanent, and any new branches should be created with the intention of contributing some new feature or physics, and being [deleted once that is done.](https://rickandmorty.fandom.com/wiki/Mr._Meeseeks)

- **Repository**: A Repository (or Repo) is a collection of all the different branches of a given project which are kept in the same location. To be specific, a location might be your local computer or part of a remote server. A Repo can be public (often called Open Source) or private, with a select list of collaborators who have read and possibly write access. A single github user may have many Repos for all of their different projects, all of which might have any number of branches. Note: COMPAS development will be done in a private repo, but there will be a second public repo which will hold only a copy of the Master branch for the public to download.

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
