Git Workflow for COMPAS software developers
===========================================



Contents of this document


`Introduction <#introduction>`__


`Getting Set Up <#getting-set-up>`__


`Day to Day Commands <#day-to-day-commands>`__


`Lifetime of a Project <#lifetime-of-a-project>`__


`COMPAS Git Workflow <#the-compas-git-workflow>`__


`Terminology <#terminology>`__ 




Introduction
============

Git & Github for COMPAS developers


For those who are unfamiliar, git and github are popular tools in the
software development community for sharing and collaborating on software
projects.

Git is a light-weight command line tool for maintaining different
versions of software locally, and sharing those versions to remote
servers.
Github is a website that stores git-managed projects and enables
developers to collaborate centrally on many projects.

It is a bigger topic than we can get into here, but if you are curious
you should `read more
here. <https://www.atlassian.com/git/tutorials/what-is-version-control>`__

Learning git is somewhat similar to learning a new language, and it
can be difficult to fully grasp the vocabulary when starting out (which
makes searching the internet for help significantly more challenging!).
Some of the most fundamental terms are `described
below <#terminology>`__ to assist new users.



Purpose of this document


The purpose of this document is to:

-  Help COMPAS users who are new to git to get set up,
-  Outline a consistent workflow for COMPAS developers in their
   day-to-day use of git, and
-  Provide some of the commands that are required for this workflow

Git is very powerful, so this is only a very small subset of the
available git commands.

This is, in some sense, a living document, meaning we are always open
to suggestions and criticism with the workflow, and seek only to find
the best option for everybody.
Please send any feedback to the `COMPAS User Google
Group <mailto:compas-user@googlegroups.com>`__.

With that said, all developers should commit to learning the agreed upon
workflow, to ensure consistency and protect against conflicts which may
derail development.



Outline of the COMPAS code repository


*Note:* If anything below doesn't make sense, try looking at the end of this document for relevant `Terminology. <#terminology>`__

COMPAS users who are not developers can download the source code from
the Main Repository, found at
`github.com/TeamCOMPAS/COMPAS <github.com/TeamCOMPAS/COMPAS>`__ (details
can be found below).
You will only need the default ``production`` branch and do not need
to worry about what branches are.

For developers, this repository (or 'repo') is considered "pristine",
meaning that any work done here should be in a mature stage.

The repository contains 2 permanent branches, ``production`` and
``dev``.
All other branches are either feature or hotfix branches, whose
purpose is to either introduce some new functionality or fix a bug,
respectively, and then be deleted.

Feature branches on the Main Repository (also called the Main Fork or
simply Main) should be ready to be tested by others.
The Main Fork is not a "sandbox" for new, experimental ideas.
You should `create your own fork <#fork-the-main-repo>`__ off of the
Main Repository if you want to have public-facing experimental work.

This approach to the repository and workflow below are based on the
`Feature Branch
Workflow <https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`__
which is in common use in industry.



Getting Set Up
==============

**Step-by-step directions for how to configure your local and remote git
repositories**

*COMPAS Users and Developers*


Set up git and create a Github account


#. If you have not already, go to `github.com <https://github.com/>`__
   and setup an account.

#. Check that you have a `working install of
   git. <https://www.atlassian.com/git/tutorials/install-git>`__

#. It is recommended, though not necessary, that you configure `Github
   with
   ssh <https://help.github.com/en/articles/connecting-to-github-with-ssh>`__
   as well.
   This will allow you to bypass frequent login verification requests.



Clone the COMPAS repository to your personal computer


Cloning a git repository is slightly different to downloading a code
database.
Git includes a history of changes made to a repo, and these changes
are included with a clone.
Cloning also includes pointers to the original (remote) repo, so that
changes to the remote can be easily imported to your local machine.

#. To clone the COMPAS repo, go to the public COMPAS Github page and
   click on the green ``code`` button (see below).
   You can copy the repo address for ssh if you have it configured,
   otherwise click ``Use HTTPS`` and copy that address.

.. raw:: html

    <p align="center">
    <img src="./media/git_clone_button.png" width="600" />
    </p>

#. In a terminal window, change directory into the location where you
   plan to install COMPAS, e.g ``cd ~/git_repos/``

#. Type the appropriate of the following two commands into a terminal
   window (pasting the repo name copied above).

   -  SSH: ``git clone git@github.com:TeamCOMPAS/COMPAS.git``

   -  HTTPS: ``git clone https://github.com/TeamCOMPAS/COMPAS.git``

#. To confirm that it worked, run the following two commands:

::

    cd COMPAS
    git branch

#. If the clone finished without error, you should see as output:

``* production``

At this point, if you do not plan to do any COMPAS development, you're
all set.
See `getting\_started.md <getting_started.md>`__ to see how to compile
and run COMPAS.
If you run into issues or would like to see new features implemented,
you can `contact us here. <mailto:compas-user@googlegroups.com>`__.
You can read on if you are curious, but if you are not invited to be a
collaborator, you will only have read-access to the repository.

*COMPAS Developers Only*


*Note:* This section is very technical.  Take a look at the section below on `Terminology, <#terminology>`__ if you get stuck!



Join as a collaborator


In order to contribute to COMPAS, you will need to be added as a
collaborator.
Non-collaborators have read-only access to all of the branches.

`Contact us here <mailto:compas-dev@googlegroups.com>`__ to inquire
about collaborating, or reach out to one of us directly (see the `COMPAS
homepage <https://compas.science/>`__ for an up-to-date list).



Fork the main repo


As a COMPAS developer, you are highly encouraged to create your own
personal github fork of the Main repo.
This is a second, remote repository (distinct from your local repo),
but is managed by your github account.
It serves as a public-facing 'sandbox' of your current work, where you
can share partially-developed ideas and projects with others who might
be interested in assisting, without interferring with or clogging up the
Main repo.

On Github, go to the TeamCOMPAS/COMPAS repo and click on ``Fork`` in
the upper-right corner.
This will create a copy of the current state of the TeamCOMPAS/COMPAS
repo, including all branches and all commit histories, and place it on
your profile, identified as ``<your-username>/COMPAS``.

Since this is your personal repo, you can be as organized or
scatter-brained as you wish here.
If you work best with 50 branches, obscure names, and code scraps
everywhere, have at it.
You can also give or take away access to any other collaborators who
you might wish to contribute.
Note that for public repositories, your code will still be read-only
for everyone who is not a collaborator.



Connect to your remote fork from your local repo


Once your fork is created, you'll want to connect it to your local
repository.
In the terminal, navigate to your COMPAS git repo and type:

``git remote add <fork-nickname> <remote-fork-url>``

The ``<remote-fork-url>`` can be found on your remote repo under the
same green 'Clone or Download' button as before.
If you have ssh configured, it will be similar to
``git@github.com:reinhold-willcox/COMPAS.git``.
The ``<fork-nickname>`` is your choice, but should be informative, e.g
``reinhold_fork``.



Day to Day commands
===================

Basic commands for navigating local git


Branches allow a developer to experiment with multiple new features
simultaneously on the same code-base.
In git, branches are very lightweight and easy to manage, making them
incredibly useful.

To view, switch, and create branches (akin to ``ls``, ``cd``, and
``mkdir``), use:

::

    git branch
    git checkout <branch-name>
    git checkout -b <new-branch>

*Note:* Many git commands require that you are on the correct branch before executing the command - using these 3 commands regularly before running more complicated commands will save you headaches down the road!

**Important:** A new branch is already created as a copy of current
branch, so you always need to double check that you're on the branch you
want to copy (typically, ``dev``).



Committing changes


What git does best is to record all the small changes and edits that
accumulate as we modify code.
After many small changes, you might have a feature that you decide
isn't actually what you want, and you want to get rid of it.
Or you might have introduced a bug at some point that spans many
files, and you need to remove it without undoing all the work you've
accomplished since then.
Git makes this incredibly easy by storing small edits as "commits".
Commits, like branches, are incredibly versatile and powerful, but can
be conceptually tricky to grasp at first.

Committing is the process of adding a selection of changes to the
history of your branch.
It is effectively saving your work, and should be done often (every
time any small fix has been made).
To perform a commit:

#. Check that you're on the correct branch!

``git branch``

#. Add the relevant files to your "index" (whatever files you've just
   edited)

``git add <file1> <file2> <...>``

#. Then submit the commit with a commit message. The message should be
   clear and concise, to help identify exactly when certain changes were
   made and undo them if necessary.

``git commit -m "really clear message indicating all the changes you made in this commit."``

*Note:* A single commit should capture an entire "fix" of one kind.

*Example:* Say you want to add a function to
a C file and its header, and you also want to update the internal
contents of a completely different function in the same C file, you
should do 2 commits.

#. First, make the edits to the first function and header, then add and
   commit

::

    git add file.C file.h
    git commit -m "created function myFunction to do someStuff and added it to the header file"

#. Then do the same for the second edits

::

    git add file.C
    git commit -m "updated internal contents of thisOtherFunction to allow for specificUseCase"

You can (and should) check the status of the current index regularly
with:

``git status``

The printout is pretty self explanatory: it tells you which files have
been added, and which have been changed that you might consider adding
before committing.

If you accidentally staged a file to your index, you can undo a
``git add`` before you have done a ``git commit`` with:

``git reset <file>``

You can also use ``git commit --amend`` to fix the most recent,
erroneous commit.

::

    git commit -m "previous commit which had the wrong files staged"
    git add <fogotten-file>
    git reset <file-that-does-not-belong>
    git commit --amend

which will open an editor where you can modify the commit message.

The takeaway message is that branches are made up of many commits
strung together, one after another, forming a history of minor edits to
a given branch.
You can view the commit history of a branch with any of:

::

    git log
    git log --pretty=oneline
    git log --pretty=medium --graph
    git log --all --decorate --oneline --graph

(If you have some spare time/ interest, there are actually quite a few
elaborate git log setups online you can look at for inspiration).

Looking through your git log, you may begin to appreciate the value of
clear, concise, commit messages.



Merging branches


Creating a branch for every new idea is great, but at some point
you'll have two somewhat-finalized, distinct features on different
branches that you will want to combine into one.
To do that, you need to merge the branches.
Merging a separate branch onto your current branch adds a 'merge
commit' to the tip of your current branch, and leaves the other branch
as it was.
The two original branches are called parent branches, and the result,
appropriately, the child.
Typically, once you successfully merge, it is desirable to delete the
separate branch to keep things tidy.

::

    git checkout branch1
    git merge branch2
    git branch -d branch2

Merging can be difficult at first because, unless you are very good at
thinking ahead or very lucky, you probably have some overlap in the two
branches that you were working on.
Git has some pretty clever tools to figure out which changes to pull
into the merge commit, but if it is ambiguous (e.g if you've edited the
same part of a file in both parents), you will get a merge conflict.
You will have to manually edit the files to choose how to resolve the
conflict.
You are encouraged to make backup copies of both parent branches until
you are more comfortable.
Git has several `ways to deal with merge
conflicts; <https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts>`__
the best option for you may depend on the particular IDE you are using.



Comparing branches


Often it is useful to see differences between branches and workspaces
without actually making any changes to either.
To accomplish this, we use the ``git diff`` command.
The arguments (or lack thereof) determine which objects are compared.

To see all the recent changes in your working directory:

::

    git diff            # compare working directory to index
    git diff HEAD       # compare working directory to tip of the current branch
    git diff --cached   # compare index to tip of the current branch

To compare two branches ``<b1>`` and ``<b2>`` (or even a single file on
separate branches):

::

    git diff <b1> <b2>                         # compare the tips of two branches
    git diff <b1> <remote-fork>/<b2>           # compare local branch to a remote branch
    git diff <b1>:./file/path <b2>:./file/path # compare the same file on different branches

For even more flexibility and control over branch/file comparisons, you
should checkout ``git difftool`` and its customizations for your
preferred text editor.



Deleting branches


You should become comfortable deleting branches, or else your repos
might pile up with old branches that are no longer active.
Branches are also very easy to manage in git (relative to other
version control systems), so you should practice creating new branches,
making quick edits, committing, and deleting again without worry.
To delete a branch,

#. Navigate to any other branch

``git checkout <unrelated-branch>``

#. Try deleting the branch

``git branch -d <branch-name>``

#. If that throws an error, likely there were some uncommited/unmerged
   changes (work that would be completely lost if the branch gets
   deleted).
   Either commit/merge the branch before deleting, or if you don't want
   to keep the changes, you can force the delete with:

``git branch -D <branch-name>``



Fetch other branches from a remote


If you followed the above workflow, you can verify that the COMPAS
repo is a designated remote fork in your local repo, nicknamed
``origin``.
You can also see any other remote forks that you have linked from your
local repo:

``git remote -v``

should output something like:

::

    origin  git@github.com:TeamCOMPAS/COMPAS.git (fetch)
    origin  git@github.com:TeamCOMPAS/COMPAS.git (push)
    reinhold_fork   git@github.com:reinhold-willcox/COMPAS.git (fetch)
    reinhold_fork   git@github.com:reinhold-willcox/COMPAS.git (push)
    another_fork    git@github.com:another-user/COMPAS.git (fetch)
    another_fork    git@github.com:another-user/COMPAS.git (push)

To see all of the available branches across all your linked forks:

``git branch -a``

should output something similar to

::

    * production
    local_feature_branch
    remotes/another_fork/dev
    remotes/another_fork/production
    remotes/another_fork/pythonSubmit
    remotes/origin/HEAD -> origin/production
    remotes/origin/dev
    remotes/origin/production
    remotes/origin/release
    remotes/reinhold_fork/dev
    remotes/reinhold_fork/git_workflow
    remotes/reinhold_fork/production

where anything not starting with "remotes/" is a local branch, and the
\* indicates your current branch.

*Note:* The remote branch named ``origin/HEAD`` is a pointer to the ``origin/production`` branch.  HEAD, when used locally, is a pointer to the most recent commit, or "tip", of the current branch.  `Read more. <https://stackoverflow.com/questions/2529971/what-is-the-head-in-git>`__

All of the remote branches are available to be copied locally with:

``git checkout -b <new-local-branch-name> <remote-name>/<remote-branch-name>``

*Example:*

``git checkout -b myPySubmit another_fork/pythonSubmit``



Configuring remote tracking branches - pushing & pulling


**Important:** This section is crucially important, but it contains
some of the more confusing subtleties of git.
I tried to make these explicit throughout, but as a result this
section is a bit dense (sorry about that).
I highly recommend trying the commands yourself as you read through.

It's often useful, though not required, to point local branches to a
branch on a remote repo, from which it will inherit changes.
For example, when changes occur on the ``dev`` branch of the Main
repo, you will probably want to pull them into your local ``dev`` branch
to keep up to date.

If changes occur on the remote, your local git repo will not
automatically know about it (git does not regularly ping the remote
server with update requests like, e.g, most phone apps).
You can check for remote changes on a fork with:

``git fetch <remote-fork>``

*Warning:* This is a bit subtle - ``git fetch`` only updates
git's "local knowledge" of the remote branches, it does not affect your
local branches.
That makes it very "safe" - you can't overwrite any of your own work
with ``fetch``.
This is not true of ``git pull`` `(see below). <#git-pull>`__

To see which local branches are tracking remote branches, use:

``git branch -vv``

which will have an output that looks similar to:

::

    * compas_hpc_updates eea656f [origin/compas_hpc_updates: behind 14] Removed references to dead files:
    dev                  a110d38 [origin/dev: ahead 2, behind 12] Remove unwanted demo files (#150)
    production           d379be5 [origin/production] Jeff's defect repairs from previous commits that had to be readded (#82)
    new_branch           b6aee96 generic branch to test git branch -vv, don't keep this

#. The first column lists your local branches (the \* indicates your
   current branch).
#. The second column is the unique hash that identifies the commit of
   the tip of that branch (technically, it's only the beginning of the
   hash, but it suffices to identify the commit).
#. If the local branch is tracking a remote branch, this will be
   specified in brackets in the third column as
   ``[<remote_repo>/<remote_tracking_branch>]``.

   -  If there is a colon after the branch name with either "ahead N" or
      "behind M" (or both), this describes whether the tip of the local
      branch has additional commits that the remote does not, and vice
      versa.

#. If there are no brackets, the branch is not tracking anything.



git pull


If you have a branch which is "behind" the remote branch it is tracking
by some number of commits, then yours is out of date and you should
update it with:

::

    git checkout <outdated-branch>
    git pull

The ``git pull`` command defaults to the remote tracking branch of the
current branch (whatever was in the brackets above).
If the current branch is not tracking anything, or if you want to pull
from a different remote branch
(e.g, if ``origin/dev`` was updated and you want your
``<local-feature-branch>`` to pull in those updates), you can set it
explicitly:

::

    git checkout <local-feature-branch>
    git pull <remote-fork> <remote-branch>

*Note:* You should regularly check that your branches are updated. If not, you should pull to avoid larger conflicts later on.



git push


To share your local work with the other collaborators, you need to
"push" your changes to a remote repository.
Similar to ``git pull``, ``git push`` defaults to the designated
remote tracking branch, if it exists.
If not, or if you want to push to a different remote branch, you can
set it manually:

::

    git checkout <local-branch-to-push>
    git push <remote-fork> <remote-branch>

Pushing to your personal remote repository is a way to save all of
your commits (i.e the history of edits) somewhere off your local
computer.
This is good practice because it acts as a backup in the event
something happens to your local machine, and it also allows other
collaborators to see your work
(without having to explicitly send them your work all the time).
This should also be done often, but not necessarily for every commit.
A good rule of thumb is to push any updated branches at the end of the
day.



pull requests


We will briefly introduce here the concept of pull requests. If
working on a remote repo, especially a shared one, it is often desirable
to block direct push access, as this could
potentially lead to bad code being introduced without proper vetting.
The solution is pull requests: the user who wrote the new code will
submit the changes as a pull request,
for another developer to review. If they pass inspection, the reviewer
can then approve the pull request and merge the changes into the remote
repo.

Clarification of the difference between push, pull, and pull requests
can be found in the `Terminology <#terminology>`__ section below.



set remote tracking branch


You can add or update a branch's remote tracking branch (sometimes
called the "upstream" branch) with:

::

    git checkout <branch-to-update>`
    git branch --set-upstream-to=<remote-fork>/<remote-branch-to-track>

*Note:* The syntax may vary slightly depending on your version of git.  ``man git branch`` should be able to shed some light.



Lifetime of a New Feature
=========================

New feature branches


When beginning a new feature, you will typically want to branch off of
the most updated version of the ``dev`` branch.
Ultimately, the feature will be merged back into ``dev`` (or else
abandoned), and this will facilitate the merge later on.

::

    git checkout dev 
    git pull
    git checkout -b <new-project>

The name of your branch should *clearly* describe the feature you plan
to implemented.
This will help you to keep track of where different bits of code live
once the number of branches gets large.



Ongoing feature branches


Commit regularly as you make changes.

::

    git status
    git add <file1> <file2> <...>
    git commit -m "useful message"

When you have made many commits and want to push your work up to the
remote, first check that you have the correct current and target
branches

::

    git branch -vv
    git push

If you are working on a shared remote branch, you should also pull
regularly to keep up with any changes that are made there. A safe way to
check if there are any changes, without risking overwriting your local
work, is to fetch and diff.

::

    git fetch <remote-fork>
    git diff HEAD <remote-fork>/<remote-branch>



Finalized features


When a feature branch is nearing completion (e.g when the code is
nearly ready to be submitted into the Main Repository and tested), you
will want to ensure that it is fully up-to-date with the Main repo.
Then, push your branches up to your personal remote repo before
submitting a Pull Request.

#. Ensure that your branch has the latest updates from ``dev``.

::

    git checkout dev
    git pull
    git checkout <mature-branch>
    git merge dev

#. Push to your personal remote repo

::

    git checkout <mature-branch>
    git push --set-upstream <your-remote-repo>

#. Submit a Pull Request to the Main repo

   -  Login to github and go to your personal remote repo
      ``<your-username>/COMPAS``.

   -  Click ``Pull request`` (If you recently pushed your branch, you
      could also click on ``Compare & pull request``)

   -  Double check that you have selected the correct feature and target
      branches. In almost all cases, the base should be
      ``TeamCOMPAS/COMPAS`` with branch ``dev``, which will probably not
      be the default. Then click ``Create pull request``

   -  Add a comment describing your feature and what changes you made.
      If you have any particular reviewers in mind, or your feature
      solves one of the Git Issues, you should link those here. Then
      click ``Create pull request``, and you're all set!

.. raw:: html

    <p align="center">
    <img src="./media/git_pr_button.png" width="600" />
    </p>

Once you have created the pull request, it is up to the other team
members to review it (see below). They may ask you to fix some parts
before accepting it, so keep an eye on the pull request conversation.



The COMPAS Git Workflow
=======================

The above sections go over many of the available git commands that you
might find useful.
This section delves into how we apply these specifically to the COMPAS
workflow.

Overview


There should always be only 2 branches on the Main Repo:
``production`` and ``dev``.
They are both permanent, and both can only be modified with pull
requests which must be approved by another COMPAS developer.

The ``production`` branch is the current "long term service release",
meaning that it should be well-tested.
Of course, code is never truly bug-free, but this branch is the one
that the public will use, so updates should be extensively tested.

The ``dev`` branch is where new features are joined together in
preparation for the next release.
Pull requests to ``dev`` should be made from feature branches sitting
on other remote repos (e.g the personal repo of the author).
Presumably, these new features have each been tested in isolation and
correctly do what they propose to do.
But ``dev`` is a place to confirm that all the new features combined
together still produce sensible output.



Reviewing Pull Requests


Typically, a new feature branch will be formally reviewed when it is
submitted as a pull request into ``dev``.
Reviewers have a responsibility to check the following:

-  The code compiles without error on the usual assortment of Operating
   Systems.
-  The code runs without error using all default values (``./COMPAS``).
-  The code runs without error on a medium-sized population of binaries.
-  The new feature(s) do what they propose to do.
-  All new features are explicitly mentioned (i.e nothing is fixed
   quietly).
-  Documentation has been updated appropriately.
-  Formatting conforms to the rest of COMPAS.

This does not all have to be done by one reviewer, but there should be a
consensus among all reviewers that all tests have been passed.

A new release is defined by a pull request from ``dev`` to
``production`` and should involve most of the active developers.
The ``dev`` branch should be tested heavily for a variety of potential
bugs, including speed tests, different package and OS versions, and
comparisons of key plots from different papers.



Terminology
===========

-  **Commit**: A single commit records a collection of edits to one or
   more files, with an associated commit message.
   You can make and undo many changes before making a commit, and you
   can similarly revert commits which are later deemed unnecessary.
   As a verb, committing changes means to create a commit of the
   changes and append that commit onto a sequence of previous commits (a
   "branch", see below).

-  **Branch**: A single branch is an ordered sequence of commits.
   A new commit is always appended onto the tip of a branch, and the
   name of the branch is really just a pointer to this most updated
   branch tip.
   When a new branch is created from an old one, they initially still
   point at the same commit, the tip of both (currently identical)
   branches.
   New commits can be applied to one branch or the other, leading to a
   divergent history (which is not a bad thing).
   The imagery of the shared history of commits, followed by the split
   into two separate histories, readily leads to the name "branches".
   A branch will often represent a place to experiment with changes in
   a way that doesn't risk destroying the existing code.
   Major branches will add some new functionality or some new physical
   prescription, while sub-branches may pop-up to quickly test some
   variation to the new functionality.
   These sub-branches might be merged in to the major feature branch,
   destroyed, or possibly continue on their own to be expanded into a
   more major feature (and then merged in later on).
   Whether the branch is merged or scrapped, it should always
   `ultimately be deleted <#deleting-branches>`__
   `[1] <https://rickandmorty.fandom.com/wiki/Mr._Meeseeks>`__ (aside
   from the permanent ``production`` and ``dev`` branches).

.. raw:: html

    <p align="center">
    <a href="https://nvie.com/posts/a-successful-git-branching-model/">
    <img src="./media/git_branches.png" width="600" />
    </a>
    </p>

-  **Repository**: A Repository (or Repo) is a single storage location
   for a given code base.
   A single github user may have many repos for all of their different
   software projects.
   In our case, we have the Main Repository hosted by on Github at
   `TeamCOMPAS/COMPAS. <github.com/TeamCOMPAS/COMPAS>`__
   There are often many repositories for a given development project -
   these can be local or remote repositories (see below), each (usually)
   hosted by one the developers.
   Each repo can contain different branches each with slight
   variations on the code base, and these branches can be readily shared
   between repos, along with their history of commits.
   A Repo can be public (often called Open Source) or private.
   COMPAS is Open Source, but the general public has only read-access.
   Prospective contributors need to be added as a collaborator in
   order to make changes and submit pull requests.

-  **Local/Remote**: Local refers to the repository on your personal
   computer, while Remote refers to any repo that isn't.
   Github repos (whether Main or someone else's) will be remote for
   everyone.
   My local computer is only local to me; from a purely git
   perspective, it would be considered remote to anyone else, though
   this should not come up often because other users should never have
   even remote access to your personal computer.
   The purpose of your personal remote fork is to be a public proxy
   for your local fork, where you can add things you've worked on that
   you wish to share around.

-  **Fork**: A Fork is a full copy of a repo, including all its
   branches, to another location.
   Most of the time, "another location" will mean elsewhere on the
   github servers, since we will be Forking from the Main Repo to our
   Personal Repo when we are setting up.
   In our case, Forks will distinguish different users, or perhaps
   groups of users (e.g Copenhagen/COMPAS).
   All core developers should have a personal fork.
   If you are familiar with the ``git clone`` command, this is
   identical to Forking from a remote server onto your own personal
   computer.

-  **Origin**: Origin is the name commonly used for the primary remote
   repository.
   It is configured by default whenever you clone from a repository,
   so yours will probably point to the Main Repo.
   If you track multiple remote forks, you should give them all
   helpful, distinguishing names (e.g ``jeff_fork``, ``reinhold_fork``,
   etc.)

-  **Working Directory**: The Working Directory is where a user makes
   edits to files.
   It has meaning in git only in reference to the Index and the most
   recent commit (i.e the tip of the current branch).
   Files are editted in the working directory, before being added to
   the Index (or "staged"), and then finally committed to the current
   branch, or HEAD (see below).

-  **Index**: The Index (aka Staging Area) exists only in the
   intermediate step between editing local files and committing those
   files.
   Historically, other Version Control systems only allowed editting
   files, and then committing those files one by one.
   The issue with that is that sometimes a collection of edits of
   different files logically make up one full "commit-worthy-edit".
   The classic example of this is adding a function to a .C file and
   it's header .h file.
   If you need to revert this commit back for any reason, it makes
   sense to remove both of those edits at once - you would virtually
   never need to remove the function from the C file but leave it in the
   header.
   Adding files to the index is the way to collect all of the files
   that were involved in a given series of edits that you want to treat
   as one big Edit.

-  **Tracking**: The word tracking has two meanings, and could refer to
   either tracked remote repositories, or tracked local files in the
   current branch, and they have slightly different meanings.

   - A tracked repository is one which contains a branch which is
     currently being tracked, or "upstream", of a branch in your local
     repository.
     By default, all the branches on a forked repository track the
     branches they were forked from.
     You can modify the upstream branch of a given branch to point at
     any other branch you like, whether local or remote. You can also
     have multiple tracked remote repositories, though any given branch
     can only track at most one other branch at a time.
     This is useful if you want to check out and keep up-to-date with a
     branch that sits on a colleague's fork.
     You can view all tracked repositories with ``git remote -v``
   - A tracked file is one that git "knows about", meaning it was
     included in the last commit.
     You can have other files in the same folders as your git repo
     which are not tracked (if, e.g, you want to have output files from
     COMPAS runs but do not want to share those around).
     If you make modifications to a tracked file but don't commit it,
     git will not let you leave the branch.

-  **Push, Pull, and Pull Request**: These commands form the backbone
   of file-sharing across repositories.
   They all cover the same conceptual idea of "taking a branch and
   copying it over to a different branch on another repo." The
   difference is where you are relative to the target.
   You ``pull`` from a remote into your local, and you ``push`` from
   your local into a remote.
   For many remotes, there are protections in place to keep arbitrary
   users from pushing changes ad hoc.
   ``Pull-requests`` are the polite version of a ``push`` - instead of
   forcing your changes onto a remote, you are asking the manager of the
   remote to review your changes, and hopefully pull them into the
   remote if they approve.

-  **Revert**: A revert is used when it is decided that a particular
   previous commit (or perhaps several) have introduced bugs or are
   otherwise no longer undesired, and we want to remove them from the
   branch.
   A ``git revert`` will attempt to identify the changes made in those
   commits, and create a new commit to undo them.
   This is a fairly advanced git command and can easily become quite
   complicated, so make sure to use this one with caution, make backups
   of your work, and do lots of testing before you try anything.

-  **HEAD**: HEAD is a pointer to a commit, but the specific commit it
   points to moves around regularly.
   In general, it refers to the tip of whichever is the current
   branch.
   When you make a commit to the branch, HEAD updates to the new tip.

-  **Log**: The log of a branch is the history of that branch in terms
   of its commits.
   The log shows when the commits occured, who authored them, and what
   the commit message stated.


