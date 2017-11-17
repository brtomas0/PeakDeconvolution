import numpy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def getData(filename):
    with open(filename) as file:
        data = [[], []]
        titles = [c.strip() for c in file.readline().split(",")]
        print titles
        for line in file:
            point = [c.strip() for c in line.split(",")]
            data[0].append(float(point[0]))
            data[1].append(float(point[1]))
    return data


def findPeaks(yvalues, window=None):
    max_peak_index = 0
    decrease_counter = 0
    increase_counter = 0
    peak_indices = []
    if window == None:
        window = len(yvalues) / 100

    for i in range(0, len(yvalues)):
        if yvalues[i] > yvalues[max_peak_index]:
            max_peak_index = i
            decrease_counter = 0
            increase_counter += 1
        elif yvalues[i] < yvalues[max_peak_index]:
            decrease_counter += 1

        if decrease_counter > window and increase_counter > window:
            peak_indices.append(max_peak_index)
            max_peak_index = i
            decrease_counter, increase_counter = 0, 0

    if len(peak_indices) == 0:
        return [None]
    return peak_indices


def gauss(x, *p):
        # A is height, mu is center, and sigma is sigma
    A, mu, sigma = p
    return A * numpy.exp(-(x - mu)**2 / (2. * sigma**2))


def fitPeak(xvalues, yvalues, peak_index, fit_window=None):
    # returns a array where the 0th index contain the xvalues and the 1 index contains the yvalues
    if fit_window == None:
        fit_window = len(xvalues) / 50

    # This is where the peak initialization is for each peak
    p0 = [yvalues[peak_index], xvalues[peak_index], 1.]

    xvalues_to_fit = xvalues[peak_index - fit_window:peak_index + fit_window]
    yvalues_to_fit = yvalues[peak_index - fit_window:peak_index + fit_window]

    coeff, var_matrix = curve_fit(gauss, xvalues_to_fit, yvalues_to_fit, p0=p0)
    # print(var_matrix)

    return [xvalues, gauss(xvalues, *coeff)]


def subtractAndZero(curve1, curve2):
    if len(curve1) != len(curve2):
        raise Exception("Curves are not the same length")

    final = [0 for i in range(len(curve1))]
    for i in range(len(curve1)):
        val = curve1[i] - curve2[i]
        if val <= 0:
            continue
        else:
            final[i] = val
    return final


def initialAllPeakFits(data, expected_peaks=10):
    yvalue_copy = [y for y in data[1]]
    first_peak_index = findPeaks(yvalue_copy)[0]
    gaussian_peaks = []
    color_index = 0
    # colors = 'rgbcm'

    while first_peak_index != None and len(gaussian_peaks) < expected_peaks:
        # print(first_peak_index, (data[0][first_peak_index], data[1][first_peak_index]))
        gauss_data = fitPeak(data[0], yvalue_copy, first_peak_index)
        gaussian_peaks.append(gauss_data[1])

        yvalue_copy = subtractAndZero(yvalue_copy, gauss_data[1])
        first_peak_index = findPeaks(yvalue_copy)[0]

        # plt.plot(data[0][index], data[1][index], marker='o', markersize=3, color='red')
        # plt.plot(gauss_data[0], gauss_data[1], '{}--'.format(colors[color_index]), label='fit-with-bounds')
        color_index += 1

    return gaussian_peaks


def refinePeaksOnce(data, peak_fits):
    new_fits = []
    color_index = 0
    # colors = 'rgbcm'

    for i in range(len(peak_fits)):
        yvalue_copy = [y for y in data[1]]
        for j in range(len(peak_fits)):
            if i != j:
                yvalue_copy = subtractAndZero(yvalue_copy, peak_fits[j])
        peak_index = max(findPeaks(yvalue_copy))
        if peak_index == None:
            continue
        gauss_data = fitPeak(data[0], yvalue_copy, peak_index)
        new_fits.append(gauss_data[1])

        # plt.plot(gauss_data[0], gauss_data[1], '{}--'.format(colors[color_index]), label='fit-with-bounds')
        # plt.plot(data[0], yvalue_copy, '{}--'.format(colors[color_index]), label='fit-with-bounds')
        color_index += 1

    return new_fits


def fitFitness(exp_yval, fit_yval):
    squared_residuals = 0
    for i in range(len(exp_yval)):
        squared_residuals += (exp_yval[i] - fit_yval[i])**2
    return squared_residuals


def fitPeaksFitness(data, peak_fits):
    fit_sums = [0 for i in range(len(data[0]))]
    for i in range(len(data[0])):
        fit_sums[i] = sum([fit[i] for fit in peak_fits])
    return fitFitness(data[1], fit_sums)


def peakDeconvolution(data, expected_peaks=10, max_iterations=100):
    peak_fits = initialAllPeakFits(data, expected_peaks)
    best_fit = fitPeaksFitness(data, peak_fits)
    best_peaks = peak_fits

    for i in range(max_iterations):
        peak_fits = refinePeaksOnce(data, peak_fits)
        fitness = fitPeaksFitness(data, peak_fits)
        if fitness < best_fit:
            best_fit = fitness
            best_peaks = peak_fits

    return best_peaks


def main():
    data = getData("test-deconvolve.csv")
    peak_fits = peakDeconvolution(data, 4)

    # The following is for coloring each fit and the total, along with the data set
    color_index = 0
    colors = 'rgbcm'
    fit_sums = [0 for i in range(len(data[0]))]
    for fit in peak_fits:
        plt.plot(data[0], fit, '{}--'.format(colors[color_index]), label='fit-with-bounds')
        fit_sums = [fit_sums[i] + fit[i] for i in range(len(fit))]
        color_index += 1

    plt.plot(data[0], fit_sums, 'black')

    plt.plot(data[0], data[1])
    plt.show()


if __name__ == '__main__':
    main()
