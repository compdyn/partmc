#!/usr/bin/env python
# Authors: 
#    Wenhan TANG - 01/2021 (Original version)
#    Wenhan TANG - 05/2021 (Add multicases_subplot_scheme class)
#    ...
if __name__ == "__main__":
    print("This script can only be imported by another python script using \"import Bench\"")
    exit()

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mpl.use("PDF")

class multicases_subplot_scheme(object):
    def __init__(self, Nrow_max, Ncol_max, span_dim = None, auto_reverse = True, remain_adjust = True):#, prev = "col"):
        self.Nrow_max = Nrow_max
        self.Ncol_max = Ncol_max
        if span_dim is None:
            span_dim = "row" if Nrow_max >= Ncol_max else "col"
        span_dim = span_dim.lower()
        assert span_dim in ["row", "col"]
        self.span_dim = span_dim
        self.auto_reverse = auto_reverse
        self.subloc_generator = None
        self.remain_adjust = remain_adjust
        #prev = prev.lower()
        #assert prev in ["row", "col"]
        #self.prev = prev

    def layout(self, Ncase):
        self.Nsub_max = self.Nrow_max * self.Ncol_max
        self.Ncase = Ncase
        if self.Ncase >= self.Nsub_max:
            layout_Nrow = self.Nrow_max
            layout_Ncol = self.Ncol_max
        else:
            if self.span_dim == "row":
                layout_Nrow = int(Ncase / self.Ncol_max) + (0 if np.mod(Ncase, self.Ncol_max) == 0 else 1)
                layout_Ncol = self.Ncol_max
            else:
                layout_Nrow = self.Nrow_max
                layout_Ncol = int(Ncase / self.Nrow_max) + (0 if np.mod(Ncase, self.Nrow_max) == 0 else 1)
        if self.auto_reverse and (self.Nrow_max - self.Ncol_max) * (layout_Nrow - layout_Ncol) < 0:
            layout_Nrow, layout_Ncol = layout_Ncol, layout_Nrow
        self.layout_Nrow = layout_Nrow
        self.layout_Ncol = layout_Ncol
        self.layout_Nsub = layout_Nrow * layout_Ncol
        #self.Ncase = Ncase

    def subloc_gen_unlimited(self):
        ic = 0
        while True:
            ic += 1
            if ic > self.layout_Nsub:
                # Just turn to next page
                ic = 1
                isNextPage = True
            else:
                isNextPage = False
            #isNextPage = True if ic == self.layout_Nsub else False
            yield isNextPage, self.layout_Nrow, self.layout_Ncol, ic

    def subloc_gen_limited(self):
        icase = 0
        ic = 0
        Ncase_orig = self.Ncase
        while(icase < Ncase_orig):
            ic += 1
            if ic > self.layout_Nsub:
                # Just turn to next page
                self.layout(Ncase_orig - icase)
                ic = 1
                isNextPage = True
            else:
                isNextPage = False
            #isNextPage = True if ic == self.layout_Nsub else False
            icase += 1
            yield isNextPage, self.layout_Nrow, self.layout_Ncol, ic
    
    def get_subloc(self):
        if self.subloc_generator is None:
            self.subloc_generator = self.subloc_gen_limited() if self.remain_adjust else self.subloc_gen_unlimited()
        return next(self.subloc_generator)

class BenchPlots(object):

    def __init__(self, figsize = (8, 11), lsc = None):

        if lsc == "l":
            figsize = (11, 8)
        if lsc == "p":
            figsize = (8, 11)

        self.figsize = figsize
        self.cfig = plt.figure(figsize = self.figsize)

    def change_pageSize(self, figsize):
        self.figsize = figsize
        plt.close(self.cfig)
        self.cfig = plt.figure(figsize = self.figsize)

    def page_reverse(self):
        figsize = self.figsize[1], self.figsize[0]
        self.change_pageSize(figsize = figsize)

    def close(self):
        pass

    def __del__(self):
        pass

    def print(self):
        pass
 
    def clear(self):
        self.cfig.clf()

    def next_page(self):

        self.print()
        self.clear() 

    #def vpage(self, xs, ys, xl, yl, kwargs = None):
    def vpage(self, xs, ys, xl, yl, **kwargs):
        #if kwargs == None:
        #    return self.cfig.add_axes([xs,ys,xl,yl])
        #else:
        return self.cfig.add_axes([xs,ys,xl,yl], **kwargs)
    
    def subplot(self, nrows, ncols, index, kwargs = None):
        if kwargs == None:
            return  self.cfig.add_subplot(nrows, ncols, index)
        else:
            return self.cfig.add_subplot(nrows, ncols, index, **kwargs)

    def subplots(self, nrows = 1, ncols = 1, **kwargs):
        return self.cfig.subplots(nrows, ncols, **kwargs)

    def subplot2grid(self, shape, loc, rowspan = 1, colspan = 1, kwargs = None):
        if kwargs == None:
            return plt.subplot2grid(shape, loc, rowspan = rowspan, colspan = colspan, fig = self.cfig)
        else:
            return plt.subplot2grid(shape, loc, rowspan = rowspan, colspan = colspan, fig = self.cfig, **kwargs)


class BenchPlots_pdf(BenchPlots):

    def __init__(self, figsize = (8, 11), lsc = None, save_name = "bench"):

        BenchPlots.__init__(self, figsize, lsc)
        self.pdf = PdfPages(save_name + ".pdf")
        self.pdf_open = True
    
    def close(self):
        self.pdf.close()
        self.pdf_open = False

    def __del__(self):
        #self.pdf.close()
        if self.pdf_open:
            self.close()
    
    def print(self):
        self.pdf.savefig(self.cfig)
    

class BenchPlots_png(BenchPlots):
    
    def __init__(self, figsize = (8, 11), lsc = None, save_name = "bench"):

        BenchPlots.__init__(self, figsize, lsc)
        self.pngdir = save_name
        if not os.path.exists(self.pngdir + "__"):
            os.makedirs(self.pngdir + "__")
        else:
            os.system("rm -f " + self.pngdir + "__/BenchPlots*.png")
        self.png_drawing = True
        self.pngNum = 1

    def close(self):
        pdfName = self.pngdir + ".pdf"
        os.system("convert " + self.pngdir + "__/BenchPlots*.png " + pdfName)
        os.system("rm -f " + self.pngdir + "__/BenchPlots*.png")
        os.system("rmdir "  + self.pngdir + "__")

    #def __del__(self):
    #    if self.png_drawing:
    #        self.close()

    def print(self):
        plt.savefig(self.pngdir + "__/BenchPlots_" + str(self.pngNum).zfill(3) + ".png", dpi = 300)
        self.pngNum += 1



