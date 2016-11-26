# Using the terminal, git and GitHub

If you're not yet familiar with working with the terminal, this is a
critical skill you'll have to learn. Some quick tutorials can be found
by googling some combination of: command line / linux / basics / cheat sheet.

## Text editor

It's good to learn how to use a terminal-based editor such as *emacs*
or *vi*. The difficulty here is to learn a number of keyboard-based
shortcuts, but in the long run it's advantageous to be able to edit
scripts without using the mouse. Emacs can be connected to work with R
via the [ESS](https://ess.r-project.org/) package (although much of
this functionality is available through RStudio).

## Config file

Having a good configuration file helps a lot: you can define *aliases*
for commands which you use often. My `~/.bashrc` file has the
following aliases:

```
alias ls='ls -G'
alias ll='ls -lh'
alias rm='rm -i'
alias mv='mv -i'
alias cp='cp -i'
```

The first two are useful changes to the way the `ls` program lists
files and directories. The last three change the `rm`, `mv`, and `cp`
function such that commands which would overwrite or remove a file are
first verified.

Note that R and emacs also have config files: `~/.Rprofile` and
`~/.emacs`.

My config files are:

* [.bashrc](https://gist.github.com/mikelove/d96fb988db039250fb8d)
* [.emacs](https://gist.github.com/mikelove/b0f4eb15a21387ddb534)
* [.Rprofile](https://gist.github.com/mikelove/c3f7ff05ce18541b8b92)

## git

*git* is a useful program for saving text files or other small files
(e.g. images) and tracking changes to these over time. There are
alternative version control systems (*svn* for example is used by the
R and Bioconductor projects), but *git* is very popular and commonly
used among computational biologists.
 
It's helpful to download a *git cheat sheet* off the internet so you
can remember the basic functions. There are many functions of git,
which can get confusing, but you really only need to use a few basic
functions, which I will discuss here. You can highlight these on the
cheat sheet and just focus on these.

## Install git

* The easiest way to install git (and other Linux software) on a Mac is
using [homebrew](http://brew.sh/). 
* The easiest way to install git on a Windows machine is to download
[this program](https://git-for-windows.github.io/). 
* If you use Linux, git is probably already installed, or if you are
  on a cluster, you may need to load it with `module load git` or
  something similar.

### Using git locally

git can be hooked up to external repositories to backup your code, but
in its most basic form, you can use git just to keep track of code on
your own machine.

Suppose you are working in a directory with a file called `file.R`.
To start a git repository or *repo* from scratch on your machine, and
to add the file you would do:

```
git init
git add file.R
git commit -m "adding the first file"
```

The first step you only need to run once, but then in the future any
change made to `file.R` will be tracked and can be commited to the git
history by using `git commit` with the argument `-a` (all changes) and
`-m "here is a message"`. You always need to add a message to each
commit action (this makes the history much more useful anyway). You
can combine the `-a` and `-m` into a single command like so:

```
git commit -am "changes to variance estimation"
```

You can use functions like `git log` to see previous commits. Using
just these functions, you will always have a backup in case some text
or files are lost or overwritten. The history is stored in binary
files in a hidden `.git` directory in the top directory of the repo.
As long as this `.git` directory is not removed, any previous state
can be recovered. It is also possible to have multiple branches, but
you can explore this functionality on your own.

## GitHub (or other external repos)

