"""
trajectory.plot
===============

A module with statistical functions for trajectories.
"""

__author__ = ("KT")
__version__ = "6.21, September 2013"


# python modules

# other modules
import numpy
import scipy.stats
import pandas
import statsmodels.formula.api
import statsmodels.stats.anova
import statsmodels.stats.multicomp


# custom made modules
import signale.tools

# package modules


def decisionLine(n, h=None, p0=.5, p1=.75, alpha=.05, beta=.01):
    """
    h ... hypothesis type, 0 for H0 and 1 for H1

    Returns:
    decision line
    """
    if not h in [0, 1]:
        print 'Calculate for H0 or H1?'
        return 0

    x = numpy.log((1 - p0) / (1 - p1))
    nenner = numpy.log(p1 / p0) + x

    b = x / nenner
    if h:
        a = numpy.log((1 - beta) / alpha) / nenner
    else:
        a = numpy.log(beta / (1 - alpha)) / nenner

    return b * n + a



def testLearningCurve(arr):
    ''' For running several statistical tests on a learning curve.

        Expects a list of 1D arrays that contain the learning curve of an
        individual subject.
    '''

    tarr = signale.tools.transposeArr(arr)
    sarr = signale.tools.sameLengthArr(arr)

    tests_subjects = {}
    tests_subjects['one-way ANOVA'] = scipy.stats.f_oneway(*arr)
    tests_subjects['one-way Kruskal-Wallis'] = scipy.stats.kruskal(*arr)
    tests_subjects['Friedman'] = scipy.stats.friedmanchisquare(*sarr)


    tests_treatments = {}
    tests_treatments['one-way ANOVA'] = scipy.stats.f_oneway(*tarr)
    tests_treatments['one-way Kruskal-Wallis'] = scipy.stats.kruskal(*tarr)
    tests_treatments['Friedman'] = scipy.stats.friedmanchisquare(*sarr.T)

    print 'Subjects'
    for t in tests_subjects:
        print t, 'stat:', tests_subjects[t][0], ' p:', tests_subjects[t][1]

    print 'Treatments'
    for t in tests_treatments:
        print t, 'stat:', tests_treatments[t][0], ', p:', tests_treatments[t][1]



def anovaLearningCurve(arr, multicomp=False):
    ''' For running anova on a learning curve.

        Expects a list of 1D arrays that contain the learning curve of an
        individual subject.
    '''

    tarr = signale.tools.transposeArr(arr)
    sarr = signale.tools.sameLengthArr(arr)


    # First, check if the variances are equal, with the "Levene"-test
    (W, p) = scipy.stats.levene(*arr)
    if p < 0.05:
        print('Warning: the p-value of the Levene test on subjects is <0.05: p={0}'.format(numpy.round(p, 4)))
    (W, p) = scipy.stats.levene(*tarr)
    if p < 0.05:
        print('Warning: the p-value of the Levene test on sessions is <0.05: p={0}'.format(numpy.round(p, 4)))

    dummy = []
    for i, r in enumerate(arr):
        dummy.append(numpy.array([r, numpy.ones(r.size) * i, numpy.arange(r.size)]).T)
    df = pandas.DataFrame(numpy.concatenate(dummy), columns=['values', 'subjects', 'sessions'])


    formula = 'values ~ C(subjects) + C(sessions)'
    lm = statsmodels.formula.api.ols(formula, df).fit()
    anovaResults = statsmodels.stats.anova.anova_lm(lm)
    print anovaResults
    print 'number of subjects =', numpy.unique(df['subjects'].values).size

    if multicomp:
        # multicomp for the animals
        print "Multiple comparisons for subjects:"
        mod = statsmodels.stats.multicomp.MultiComparison(df['values'], df['subjects'])
        print(mod.tukeyhsd().summary())
        # print(mod.allpairtest(scipy.stats.ttest_rel, method='bonferroni')[0])

        # multicomp for the sessions
        print "Multiple comparisons for sessions:"
        mod = statsmodels.stats.multicomp.MultiComparison(df['values'], df['sessions'])
        print(mod.tukeyhsd().summary())
        # print(mod.allpairtest(scipy.stats.ttest_rel, method='bonferroni')[0])


    return anovaResults, df
