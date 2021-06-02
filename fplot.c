#include <float.h>
#include <math.h>
#include <notcurses/notcurses.h>
#include "functions.h"

#define SWAPDBL(a, b) do { double t_m_p = a; a = b; b = t_m_p; } while (0)

struct g {
	struct notcurses *nc;
	struct ncplane *stdp, *drawp;
	int termw, termh;
	int wpx, hpx;
	void *emptybuf;
};

static inline double fpart(double x); // For Wu's algorithm
static void fplot(double (*f)(double), double argbeg, double argend);
static void init();
static inline uint32_t linepixel(double intensity);
static inline double rfpart(double x); // For Wu's algorithm

struct g g;

#include "config.h"

double
fpart(double x)
{
	return x - floor(x);
}

void
fplot(double (*f)(double), double argbeg, double argend)
{
	double step = (argend - argbeg) / g.wpx;
	double minval = DBL_MAX, maxval = DBL_MIN;
	double *vals = malloc(g.wpx * sizeof(double));
	for (int i = 0; i < g.wpx; ++i) {
		double val = f(argbeg + step * i);
		if (val < minval) minval = val;
		if (val > maxval) maxval = val;
		vals[i] = val;
	}
	double valspan = maxval - minval;
	struct ncvisual *ncv = ncvisual_from_rgba(g.emptybuf, g.hpx, g.wpx * 4, g.wpx);
	// 1. Draw lines in between plot points
	for (int i = 0; i < g.wpx-1; ++i) {
		double y0 = (g.hpx-1) * (1 - (vals[i] - minval) / valspan);
		double y1 = (g.hpx-1) * (1 - (vals[i+1] - minval) / valspan);
		double x0 = i;
		double x1 = i+1;
		// Wu's line drawing algorithm (adapted from Wikipedia pseudocode)
		bool steep = fabs(y1 - y0) > 1; // ... > abs(i+1 - i)
		if (steep) {
			SWAPDBL(x0, y0);
			SWAPDBL(x1, y1);
		}
		if (x0 > x1) {
			SWAPDBL(x0, x1);
			SWAPDBL(y0, y1);
		}
		double dx = x1 - x0;
		double dy = y1 - y0;
		double gradient = dy / dx;
		if (dx == 0) gradient = 1;
		// First endpoint of the line
		double xend = round(x0);
		double yend = y0 + gradient * (xend - x0);
		double xgap = rfpart(x0 + 0.5);
		int xpxl1 = xend;
		int ypxl1 = floor(yend);
		if (steep) {
			ncvisual_set_yx(ncv, xpxl1, ypxl1, linepixel(rfpart(yend) * xgap));
			ncvisual_set_yx(ncv, xpxl1, ypxl1+1, linepixel(fpart(yend) * xgap));
		} else {
			ncvisual_set_yx(ncv, ypxl1, xpxl1, linepixel(rfpart(yend) * xgap));
			ncvisual_set_yx(ncv, ypxl1+1, xpxl1, linepixel(fpart(yend) * xgap));
		}
		double intery = yend + gradient;
		// Second endpoint of the line
		xend = round(x1);
		yend = y1 + gradient * (xend - x1);
		xgap = fpart(x1 + 0.5);
		int xpxl2 = xend;
		int ypxl2 = floor(yend);
		if (steep) {
			ncvisual_set_yx(ncv, xpxl2, ypxl2, linepixel(rfpart(yend) * xgap));
			ncvisual_set_yx(ncv, xpxl2, ypxl2+1, linepixel(fpart(yend) * xgap));
		} else {
			ncvisual_set_yx(ncv, ypxl2, xpxl2, linepixel(rfpart(yend) * xgap));
			ncvisual_set_yx(ncv, ypxl2+1, xpxl2, linepixel(fpart(yend) * xgap));
		}
		// Main loop
		if (steep) {
			for (int x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
				ncvisual_set_yx(ncv, x, floor(intery), linepixel(rfpart(intery)));
				ncvisual_set_yx(ncv, x, floor(intery)+1, linepixel(fpart(intery)));
				intery += gradient;
			}
		} else {
			for (int x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
				ncvisual_set_yx(ncv, floor(intery), x, linepixel(rfpart(intery)));
				ncvisual_set_yx(ncv, floor(intery)+1, x, linepixel(fpart(intery)));
				intery += gradient;
			}
		}
	}
	// 2. Draw dots at plot points (the endings of the lines)
	for (int i = 0; i < g.wpx; ++i) {
		double y = (g.hpx-1) * (1 - (vals[i] - minval) / valspan);
		ncvisual_set_yx(ncv, round(y), i, linepixel(1));
	}
	free(vals);
	ncvisual_render(g.nc, ncv, &(struct ncvisual_options) {
		.n = g.drawp,
		.x = 0, .y = 0,
		.scaling = NCSCALE_NONE,
		.blitter = NCBLIT_PIXEL
	});
	ncvisual_destroy(ncv);
}

void
init()
{
	g.nc = notcurses_core_init(&(struct notcurses_options) {
		.flags = NCOPTION_SUPPRESS_BANNERS,
//		.loglevel = NCLOGLEVEL_DEBUG
	}, stdout);
	g.stdp = notcurses_stddim_yx(g.nc, &g.termh, &g.termw);
	if (notcurses_check_pixel_support(g.nc) < 1) {
		notcurses_stop(g.nc);
		fputs("No pixel support!", stderr);
		exit(2);
	}
	// We need this because Notcurses doesn't allow pixel blitting over the standard plane.
	g.drawp = ncplane_create(g.stdp, &(struct ncplane_options) {
		.x = 0, .y = 1,
		.rows = g.termh - 1,
		.cols = g.termw
	});
	ncplane_pixelgeom(g.drawp, NULL, NULL, NULL, NULL, &g.hpx, &g.wpx);
	uint8_t *buf = malloc(g.wpx * g.hpx * 4);
	for (int i = 0; i < g.wpx * g.hpx; ++i) {
		buf[4*i + 0] = 0;
		buf[4*i + 1] = 0;
		buf[4*i + 2] = 0;
		buf[4*i + 3] = 255;
	}
	g.emptybuf = buf;
}

uint32_t
linepixel(double intensity)
{
	int gray = round(intensity * 255);
	return ncpixel(gray, gray, gray);
}

double
rfpart(double x)
{
	return 1 - fpart(x);
}

int
main()
{
	init();
	fplot(g_function, g_xstart, g_xend);
	notcurses_render(g.nc);
	notcurses_getc_blocking(g.nc, NULL);
	notcurses_stop(g.nc);
	free(g.emptybuf);
	return 0;
}
