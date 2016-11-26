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

## GitHub (or other external repos)

