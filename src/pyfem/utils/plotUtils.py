from pylab import plot, show


def plotCurve(output):
    plot([x[0] for x in output], [x[1] for x in output], 'r-o')

    show()
