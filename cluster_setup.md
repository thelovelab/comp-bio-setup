# Cluster setup (specific to UNC)

Our lab primarily uses the longleaf cluster for computing. While you
may sometimes work locally on your laptop, you will eventually want to
also have the same code on the cluster, so you can run simulations for
example, without having to keep your laptop open and running. You may
also need to work on longleaf if you have data that cannot be
downloaded to a laptop.

The main pieces you need to work on longleaf are:

* Use OnDemand for interactive work with RStudio (new as of 2020), or...
* X11 forwarding for when you SSH (so you can see plot windows)
* Some way of editing and running R code on the cluster,
  either [ESS](https://ess.r-project.org/)
  or [RStudio](https://www.rstudio.com/products/RStudio/) 
* Version control using [git](terminal_git_github.md).

## OnDemand for interactive work

As of 2020, Research Computing at UNC has made a very nice solution
for interactive work on the cluster, which makes the piece below about
X11 forwarding and ESS irrelevant. For various data science applications,
first see if they are supported here, as this will be a much easier interface
for most students. If you are off campus you will need to connect via VPN first.

<https://ondemand.rc.unc.edu>

## X11 forwarding

For Windows machines, you need to download and install an X11
application, such as Xming. (You may need to open Xming manually for
X11 forwarding to work when you connect to longleaf.)

<https://sourceforge.net/projects/xming/>

For Mac, you should install XQuartz. XQuart will open automatically
when you start an X11 forwarding session.

<https://www.xquartz.org/>

Then when you ssh to longleaf, you can use the following shell
command, where the -Y flag enables trusted X11 forwarding

```
ssh -Y username@longleaf.unc.edu
```

To test if this works, you can try the following in an interactive
shell and see if a plot open up.

```
module load r
R
> plot(1:10)
```

## Editing and running R code on the cluster

First, as always, you need to load an interactive shell with:

```
srun --x11=first --pty --mem=5000 --time=360 /bin/bash
```

I have these commands in my `.bashrc` file as aliases, so I don't have
to type this out each time:

```
alias inter='srun --x11=first --pty --mem=5000 --time=360 /bin/bash'
alias interbig='srun --x11=first --pty --mem=20000 --time=360 /bin/bash'
```

For editing scripts and running R code this, I personally use emacs
and [ESS](https://ess.r-project.org/), but you can also load the
RStudio module on the cluster and work within an RStudio window. From
my experience there is too much lag in the window update, so it
restricts the speed of my typing / data analysis.

Note that there are multiple versions of R actually on the
cluster. You can see the available versions:

```
module avail 2>&1 >/dev/null | grep ' r/'
```

To load the latest version you can use, e.g.:

```
module load r/x.y.z
```

You can put the `module load` commands into your `~/.bashrc` file so
that it loads every time.

To use emacs with ESS, put `module load ess` into your `~/.bashrc` file.

Then in your `~/.emacs` file, add:

```
(require 'ess-site)
```

You can see
my [.emacs](https://gist.github.com/mikelove/b0f4eb15a21387ddb534)
file for more ESS customization. 

Assuming you know about emacs keybindings, to load R within emacs,
either type `M-x R` or you can start by running a line of code from an
R script with `C-c C-n`. See this reference card for ESS keybindings:

<http://ess.r-project.org/refcard.pdf>


## Version control using git

For editing data analysis R scripts or working on a new method, you
should be saving your code in git repositories, and typically also
syncing this with a BitBucket or GitHub remote server.

You will have to set up SSH keys on the cluster,
to sync git repositories on the cluster with GitHub or BitBucket.
You can follow the steps described on the [the git page](terminal_git_github.md).

At the end, the ideal setup is to have GitHub repos on your laptop and
the same repo on the cluster, and you will use `git pull` to keep all
code up to date on all locations. You should `commit` and `push` your
code daily, to avoid any lost work.
