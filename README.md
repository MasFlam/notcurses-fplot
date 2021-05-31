# Notcurses function plotter
Function plotting program for the terminal, powered by
[Notcurses](https://github.com/dankamongmen/notcurses). A practical example
use of `ncvisual` with `NCBLIT_PIXEL`. It also implements the Wu's line
drawing algorithm.

# Configuration
The configuration is changed through [`config.h`](/config.h). You can change the function
to plot and the domain to plot it in. The function can be any function that takes a `double`
ans returns a `double`.

# Building
Make sure Notcurses is installed (on Debian it's `libnotcurses-dev`). Then run `make`
to compile and link the `fplot` executable.
