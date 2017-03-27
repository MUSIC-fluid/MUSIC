#! /usr/bin/env python

from numpy import *


def getPlotElements(idx):
    """
    return line style, marker type, and color for curve. It loops over the
    whole lists of options.
    """
    lineStyleList = ['-', '--', '-.', ':']
    colorList = ['k', 'r', '#009900', 'b', '#990099', '#CC6600', '#009999',
                 '#FF8000', '#00d2ff', '#f12b6e', '#cfc755', '#598254']
    shadowColorList = ['#C0C0C0', '#FFA8A8', '#B3FFB8', '#A4CCFF', '#EF9BFF',
                       '#FFDE80', '#99FFFF', '#FFB266', '#00d2ff', '#fb6eb1',
                       '#cfc755', '#598254']
    MarkerList = ['s', 'o', '^', 'v', 'D', 'p', '*', 'H', '.', ',', '<',
                  '>', '1', '2', '3', '4', 'h', '+', 'x', 'd', '|', '_']
    return (lineStyleList[idx % len(lineStyleList)],
            MarkerList[idx % len(MarkerList)], colorList[idx % len(colorList)],
            shadowColorList[idx % len(shadowColorList)])


def getBinnedAveragedDatawithErrorbars(dataMatrix, nbin, setBinBoundary=False,
                                       binBoudary=(0.0, 1.0), bincol=0):
    """
      Return the binned average of data, together with the count of the number 
      of data and normalized probability distribution. It returns another
      matrix for the statistical error bars for each binned quantities
    """
    if dataMatrix.ndim == 1:
        ncol = 1
        dataMatrix = dataMatrix.reshape(len(dataMatrix), 1)
    else:
        ncol = len(dataMatrix[0, :])
    binColumn = dataMatrix[:, bincol]
    ntotal = len(binColumn)

    if setBinBoundary:
        binMin = binBoudary[0]
        binMax = binBoudary[1]
    else:
        binMin = min(binColumn)
        binMax = max(binColumn) + 1e-8
    binEdges = linspace(binMin, binMax, nbin + 1)

    binnedData = zeros([nbin, ncol + 2])
    binnedData_err = zeros([nbin, ncol + 2])

    for idxbin in range(nbin):
        binWidth = binEdges[idxbin + 1] - binEdges[idxbin]
        binmid = (binEdges[idxbin + 1] + binEdges[idxbin]) / 2.
        idxdata = logical_and(binColumn >= binEdges[idxbin],
                              binColumn < binEdges[idxbin + 1])
        nsamples = len(dataMatrix[idxdata, bincol])
        for icol in range(ncol):
            if nsamples == 0:
                binnedData[idxbin, icol] = 0
                binnedData_err[idxbin, icol] = 0
            elif nsamples == 1:
                binnedData[idxbin, icol] = mean(dataMatrix[idxdata, icol])
                binnedData_err[idxbin, icol] = 0
            else:
                binnedData[idxbin, icol] = mean(dataMatrix[idxdata, icol])
                binnedData_err[idxbin, icol] = std(
                    dataMatrix[idxdata, icol]) / sqrt(nsamples - 1)
        if nsamples == 0:
        # if there is no data in the bin, middle value of the bin is used
            binnedData[idxbin, bincol] = binmid

        binnedData[idxbin, ncol] = float(nsamples) / ntotal
        binnedData[idxbin, ncol + 1] = (float(nsamples) / ntotal) / binWidth
        # calculate statistical errors according to Poisson distribution
        if ntotal > 1:
            binnedData_err[idxbin, ncol] = (sqrt(float(nsamples) / ntotal)
                                            / sqrt(ntotal - 1))
            binnedData_err[idxbin, ncol + 1] = sqrt(
                float(nsamples) / ntotal) / binWidth / sqrt(ntotal - 1)
    return binnedData, binnedData_err

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y
    
    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)
    def ufunclike(xs):
        return array(map(pointwise, array(xs)))
    return ufunclike
