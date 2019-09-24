[//]: ## (grip -b git_workflow.md)
	
# Git Workflow for COMPAS software developers

## Introduction 

### 1. Git and Github
For those who are unfamiliar, git and github are popular tools in the software development community for sharing and collaborating on software projects. Git is a light-weight command line tool for maintaining different versions of software locally, and distributing those versions to remote servers. Github is a website that centralizes for storage of git-managed projects. It is a big topic, and not worth getting into here, but if you are curious you can [read more](https://www.atlassian.com/git/tutorials/what-is-version-control).

### 2. Purpose
The purpose of this document is to outline a consistent workflow for COMPAS developers in their day-to-day use of git, protocols for whenever new projects are started or completed, and the commands that are required for this workflow (git is very powerful, so this is only a very small subset of the availalbe git commands). This is, in some sense, a living document, meaning we are always open to [suggestions and criticism](mailto:reinhold.willcox@monash.edu) with the workflow, and seek only to find the best option for everybody. With that said, everyone should commit to learning the agreed upon workflow, to ensure consistency between developers and protect against user-error which may derail development.

### 3. Outline
Broadly speaking, the setup is outlined below:

- The workflow here is based on the [Feature Branch Workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow) in common use in industry, in which there are 2 permanent branches, Master and Dev. All other branches are considered "sub-branches", based on projects and/or concepts, whose purpose is to add a feature and then be deleted.

- The Main Repository (also called the Main Fork, or simply Main) is considered "pristine", and should only contain the Master and Dev branches, and any sub-branches that are nearly ready to be made public (more on this later). 

- Developers should each have their own personal Forks of the Main Repository, and any work done locally (on their personal computers) should be pushed up through their personal Repository before being added to Main (more on this later). 

### 4. Terminology

- *Branch*: Branches in git separate work-streams for different features (e.g front-end developers might have a branch for a fancy new button for their website, while back-end developers might have a branch to make database-entry easier). In our case, branches will distinguish different projects or concepts (e.g Supernova-Kicks, White-Dwarf-Accretion, etc.). Branches should _not_ be used to distinguish developers. As mentioned previously, only the Master and Dev branches are permanent, and any new branches should be created with the intention of contributing some new feature or physics, and being [deleted once that is done.](https://rickandmorty.fandom.com/wiki/Mr._Meeseeks)

- *Repository*: A Repository (or Repo) is a collection of all the different branches of a given project which are kept in the same location. To be specific, a location might be your local computer or part of a remote server. A Repo can be public (often called Open Source) or private, with a select list of collaborators who have read and possibly write access. A single github user may have many Repos for all of their different projects, all of which might have any number of branches. Note: COMPAS development will be done in a private repo, but there will be a second public repo which will hold only a copy of the Master branch for the public to download.

- *Fork*: A Fork is a full copy of a Repo, including all its branches, to another location. Most of the time, "another location" will mean elsewhere on the github servers, since we will be Forking from the Main Repo to our Personal Repo when we are setting up. In our case, Forks will distinguish different users, or perhaps groups of users (e.g Copenhagen/COMPAS). All core developers should have a personal fork. If you are familiar with the `git clone` command, this is identical to Forking from a remote server onto your own personal computer. 

- *Local/Remote*: Local refers to what's on your personal computer, while Remote refers to anything that isn't. Github will be remote for everyone. My local computer is only local to me, and would be considered remote to anyone else. This should not come up often, because it would be very foolish to give anyone access to my local computer, even if it's a developer I trust. The purpose of your personal remote fork is to be a proxy for your local fork, where you can add things you've worked on that you wish to share around.

- *Tracking*: A tracked repository is one which has a branch which is "upstream" of a branch in your local repository. By default, all the branches on a cloned or forked repository track the branches they were cloned from, and new branches track their parent branches. You can modify the upstream branch of a given branch to point at any other branch you like, whether local or remote, and can even track multiple remote repositories (see below). This is useful if you want to check out a branch that sits on a colleague's Fork.

- *Origin*: Origin is a shorthand for the most relevant remote repository. It is configured by default whenever you clone from a repository, so yours will point to your personal remote fork. If you track multiple remote forks, you should give them all helpful, distinguishing names (e.g jeff_fork, main_compas_repo)

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

- If you have not yet configured Github with ssh, you can clone over http: 

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

- When beginning a new project, you should start your first branch locally


### Commit often and push to your repo 
### Work collaboratively across forks 

---

## Testing and Release of New Versions 

### Overview 
First, any and all changes to the main repo are done through pull-requests, not pushes (a pull-request really represents a push, with the added step of confirmation from a third party that the push meets certain standards). This is to ensure that the main repo is "clean" and contains code that has at least been somewhat tested. Branches will first be added into the main repo from the remote forks of a given collaborator, then (after further testing) added to the Dev branch. Finally, when sufficient sub-branches have been added to the Dev branch to warrant a new version release, the Dev branch will be added to the Master branch and published to the public. Ultimately, we should always be thinking of the next Version Release and which projects/concepts we would like to include in it, and define our timelines around those expectations.  

A mature branch will have to through 3 rounds of Q&A before it is ready for deployment in a new COMPAS version. The first occurs when a pull-request is sent to have the branch included as a sub-branch in the main repo. The second occurs when a pull-request is sent to have that sub-branch merged into the Dev branch on the main repo. The final occurs when a pull-request is sent to have the Dev branch (which might contain several sub-branches) merged into the Master branch. Testing of a branch can never be done by someone who worked on it extensively (we can decide on a case-by-case basis who falls into that category). 

1. To be included as a sub-branch in the main repo, the branch must be reviewed by one person in the broader COMPAS collaboration. It must compile and run, producing somewhat sensible output (e.g output files should not be empty, but at this stage they may contain data which is "wrong", or physically inconsistent) 

2. To be merged into the Dev branch on the main repo, the branch must then be reviewed by two people in the broader COMPAS collaboration, one of whom is a core developer. Additionally, one of the reviewers should focus on the consistency of the code layout and structure with the current code, while the other reviewer should focus primarily on the accuracy of the physics and its implementation, although ultimately both reviewers are responsible for ensuring that the branch is clearly and correctly written. This is the primary review step, with the expectation that code that is accepted to be merged into Dev is ready for the public. 

3. When Dev is updated with a new branch, all collaborators should be testing it, and in particular should ensure that all of the different projects which have been added to Dev work with each other. When Dev is nearly ready to be merged into the Master branch, the core collaborators will host a COMPAS powwow to discuss and assign any final revisions/tests, and set a firm release date. Other collaborators are invited to join (over Zoom or in person), details decided on a case-by-case basis.

---

## Example 
