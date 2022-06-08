#!/usr/bin/env python
""" imf.py

Author: Chong-Chong He (che1234@umd.edu)
Written on Sat,  4 Apr 2020.
"""

from random import random
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from .center import MASS_SHIFT
from . import center, utilities
from .ramses import Ramses, NoSinkParticle

# from ramses.py

def plot_step_func(ax, n, bins, **kwargs):
    """
    Accurate and precise. The bins are representing the value between the
    pivots.
    """

    n = [n[0]] + n.tolist()
    ax.step(bins, n, **kwargs)

def plot_imf_from_nbins(ax, n, bins, isLog=True, **kwargs):
    """ Plot a mass spectrum on ax, given n and bins. y: d N / d log(m)
    """

    # bins_log = np.log10(bins)
    bins_log = np.log10(bins) if isLog else bins
    bins_size = bins_log[1:] - bins_log[:-1]
    n = np.double(n) / bins_size
    # n = [n[0]] + n.tolist()
    # ax.step(bins, n, **kwargs)
    plot_step_func(ax, n, bins, **kwargs)
    ax.set_yscale('log')
    if isLog:
        ax.set_xscale('log')
        ax.set_xlabel(r"$M$ $({\rm M}_\odot)$")
        ax.set_xlim([1e-2, 1e2])
    else:
        ax.set_xlabel(r"$\log \, M$ $({\rm M}_\odot)$")
    # ax.set_ylabel(r"${\rm d}N/{\rm d}\log_{10}(M/{\rm M}_{\odot})$")
    # ax.set_ylabel(r"$dN/d\log M$ $({\rm M}_\odot)$")
    return

def plot_imf(ax, mass, bins=None, xlim=[-2, 2], ylim=[0.7, 3.3e2],
             is_over_dlogm=True, is_text=True, **kwargs):
    """ Plot the mass spectrum of mass on ax.
    x: mass bins; y: d N / d log10(m)
    """

    if bins is None or type(bins) is str:
        bins = 10**np.arange(-2.0, 3.1, 0.2)
    n, bins = np.histogram(mass, bins=bins)  # or bins='doane')
    log_bins = np.log10(bins)
    if is_over_dlogm:
        dlogm = log_bins[1:] - log_bins[:-1]
        n = n / dlogm
        ylabel = r"${\rm d}N/{\rm d}\log(M)$"
    else:
        ylabel = "N"
    plot_step_func(ax, n, log_bins, **kwargs)
    ax.set(xlim=xlim, ylim=ylim,
           xlabel=r"$\log \, M$ $({\rm M}_\odot)$",
           ylabel=ylabel,
           yscale='log',
    )
    if not is_over_dlogm:
        if ylim[1] > 100:
            ax.set_yticks([1, 10, 100])
            ax.set_yticklabels(['1', '10', '100'])
        else:
            ax.set_yticks([1, 10])
            ax.set_yticklabels(['1', '10'])
    ax.minorticks_on()

    if is_text:
        # write texts
        text = r"$m_{{\rm tot}} = {:.3g}$".format(np.sum(mass))
        text += "\n"
        text += "$N = {}$".format(len(mass))
        ax.text(0.98, 0.98, text, va='top', ha='right',
                # color='k',
                #fontsize='small',
                transform=ax.transAxes)

def imf(mass, ax=None, **kwargs):
    if ax is None:
        _, ax = plt.subplots()
    plot_imf(ax, mass, **kwargs)
    return

def plot_imf_N_vs_m(ax, m, color='k', istext=True, bin_width=None, bins=None):
    """Given m, plot the mass function where the y axis is the number of
    stars and x axis is log(m), with bin size 0.1 deg.

    """

    ax.set(yscale='log', ylim=[0.7, 3.3e2], xlim=[-1, 3],
           ylabel='Number of stars',
           xlabel=r"$\log \, M$ $({\rm M}_\odot)$")
    ax.set_yticks([1, 10, 100])
    ax.set_yticklabels(['1', '10', '100'])
    ax.minorticks_on()
    # write texts
    text = r"$m_{{\rm tot}} = {:.3g}$".format(np.sum(m))
    text += "\n"
    text += "$N = {}$".format(len(m))
    ax.text(0.98, 0.98, text, va='top', ha='right', color='k', fontsize='small',
            transform=ax.transAxes)
    if len(m)==0:
        return

    # set bins and get n
    if bins is None:
        if bin_width is None:
            bin_width = 0.1
        bins = 10**np.arange(-2.1, 3, bin_width)
    else:
        logbin = np.log10(bins)
        bin_width = logbin[1:] - logbin[:-1]
    n, bins = np.histogram(m, bins=bins)

    # plot mass function
    numb = n * bin_width        # so that dn = numb / dlogm
    bins_log = np.log10(bins)
    plot_imf_from_nbins(ax, numb, bins_log, isLog=False, color=color)

    if not istext:
        return
    return

# sink (from sink.py)
def get_bins(mass, maxBins=None, m_min=0.1):
    """ Get the bins for IMF. Apply on shifted mass """

    logmass = np.log10(mass[mass > m_min])
    # logmass = logmass[logmass > np.log10(m_min)]
    n, bins = np.histogram(logmass, bins='auto')
    if maxBins is None:
        maxBins = 5 + int(len(logmass) / 100.0)
    # else:
    #     maxBins = 16
    if len(n) > maxBins:
        n, bins = np.histogram(logmass, bins=maxBins)
    return 10**bins

def likeli(m1, m2, alpha, mu):
    """ Return the log likelihood. mu = np.log10(mass) is all masses. """

    return len(mu) * (np.log(alpha) - np.log(m1**(-alpha) - m2**(-alpha))) \
           - alpha * np.sum(mu) * np.log(10)

def find_sigma(logl, alpha, deltalogl=-0.5):
    """ Return the most probable alpha and 1 sigam region of alpha,
    given the log likelihood logl.

    Return:
    -------
    alpha_max, alpha_minus, alpha_plus
    """

    assert deltalogl < 0.0
    region = alpha[logl > np.max(logl) + deltalogl]
    alpha_max = alpha[np.argmax(logl)]
    return alpha_max, region[0] - alpha_max, region[-1] - alpha_max

def bayesian_imf_fit(mass, m1, m2):
    """ Do a Bayesian fitting to mass with a power-law slope at -alpha
    in range m1 to m2. Return alpha.
    mass is the shifted sink mass.

    Parameters
    ----------
    mass: [array], all the masses to fit with
    m1, m2: fitting range

    Returns
    -------
    (A, a_peak, a_peak_minus, a_peak_plus)
    """

    mu = np.log10(mass)
    alphas = np.arange(0.5, 1.8, 0.001)
    L = [likeli(m1, m2, a, mu) for a in alphas]
    a_peak, a_minus, a_plus = find_sigma(L, alphas)
    A = np.log(10) * len(mass) * a_peak / (m1**(-a_peak) - m2**(-a_peak))

    is_use_accum = 1
    if is_use_accum:
        # use 16% and 84% point as the 1-sigma region
        L = np.array(L) - np.max(L)  # raise L up so that 10**L is not too small
        like = np.exp(np.array(L))
        like_accum = np.cumsum(like)
        like_sum = like_accum[-1]
        one_sigma = .6827
        left = np.searchsorted(like_accum, like_sum * (1-one_sigma)/2)
        right = np.searchsorted(like_accum, like_sum * (1+one_sigma)/2)
        a_minus = alphas[left] - a_peak
        a_plus = alphas[right] - a_peak

    return A, a_peak, a_minus, a_plus

def bayesian_imf_fit_new(mass, m1, m2=None):
    """Based on and replacing bayesian_imf_fit. Now the mass fed in could
    be all the masses in a simulation. Masses outside of [m1, m2] will
    be cut out.

    """
    if m2 is None:
        pick = mass >= m1
        m2 = np.max(mass)
    else:
        pick = np.logical_and(mass >= m1, mass <= m2)
    mass_cut = mass[pick]
    return bayesian_imf_fit(mass_cut, m1, m2)

def get_bins_fixed():

    # return 10.**np.arange(-1, 3, 0.2)
    # return 10.**np.arange(-1.04, 3, 0.2)
    return 10**np.arange(-1.1, 3, 0.2)  # This bin makes Job3.3.2 look much better

def overplot_kroupa_new(self, ax, m_sum, m_min=0.08, m_max=100, textxy=None,
                        ls='--', isxloged=False, scale=1., **kwargs):
    """ over plot a Kroupa IMF curve. parameters are in log scale
    Normalized to total stellar mass

    Notes
    -----
    self is useless. I put INTERP here for INTERP to be copied in the Sink class
    """

    print("here 23")
    def kroupa_lin(m):
        if m < 0.5:
            return m**(-1.3) * m
        else:
            return 0.5 * m**(-2.3) * m
    kroupaNorm = integrate.quad(kroupa_lin, m_min, m_max)[0]
    # print("kroupaNorm = ", kroupaNorm)
    normFactor = m_sum / kroupaNorm
    def kroupa_log(m):
        """ Generate the Kroupa IMF """
        if m < 0.5:
            return normFactor * np.log(10) * m**(-0.3)
        else:
            return normFactor * 0.5 * np.log(10) * m**(-1.3)
    if m_min < 0.5:
        masses = [m_min, 0.5, m_max]
    else:
        masses = [m_min, m_max]
    dNdlogm = np.array([kroupa_log(m) for m in masses])
    x = masses if not isxloged else np.log10(masses)
    ax.plot(x, scale * dNdlogm, ls=ls, **kwargs)
    print("x =", x, "y =", scale * dNdlogm)
    if textxy is not None:
        ax.text(10**textxy[0], 10**textxy[1], "Kroupa single",
                **kwargs)

    return

def overplot_salpeter(ax, m_sum, m_min=0.08, m_max=100, textxy=None, **kwargs):
    """ Overplot a Salpeter line. textxy=[textx, texty] in
    log scale """

    alpha = 2.35
    m_tot = (m_max**(2 - alpha) - m_min**(2 - alpha)) / (2 - alpha)
    norm = m_sum / m_tot

    def salpeter_log(m):
        return norm * np.log(10) * m**-1.35

    xs = [m_min, m_max]
    ys = [salpeter_log(mi) for mi in xs]
    # xs = [10**(xi) for xi in xs]
    # ys = [10**(yi) for yi in ys]
    ax.plot(xs, ys, '--', **kwargs)
    if textxy is not None:
        ax.text(10**textxy[0], 10**textxy[1], r"$\Gamma = -1.35$", **kwargs)

    return

def mass_vertical_line(ax, masses, thresh, color='k'):
    """ add mass divider. Apply to masses directly (will not shift the mass)
    """

    # masses = self.get_sink_masses(outputid)
    mass_sort = np.sort(masses)
    mass_cum = np.cumsum(mass_sort)
    indice = np.searchsorted(mass_cum, np.sum(masses) * thresh)
    mass_thr = mass_sort[indice]
    ax.axvline(mass_thr, color=color, ls='--')

def bayes_get_param(masses, thresh):
    """Given star masses and the percentile cummanative mass fraction to
    fit with, return the bayes parameters and the mass cut

    """

    mass_sort = np.sort(masses)
    mass_cum = np.cumsum(mass_sort)
    indice = np.searchsorted(mass_cum, np.sum(masses) * thresh)
    mass_thr = mass_sort[indice]
    # choose mass-percentage threshhold as the starting point of Bayes fit
    bayes = bayesian_imf_fit_new(masses, mass_thr)
    return bayes, mass_thr

def plot_bayes_fit_any_width(ax, masses, thresh, yscaleup=1,
                             color='g', is_xlog=False):
    """ plot a bayes fit on top of a IMF"""

    bayes, mass_thr = bayes_get_param(masses, thresh)
    x = np.array([mass_thr, np.max(masses)])
    y = bayes[0] * x**(-bayes[1])
    y *= yscaleup
    if is_xlog:
        ax.plot(np.log10(x), y, '-', color=color)
    else:
        ax.plot(x, y, '-', color=color)
    plus = bayes[3] if bayes[3] >= 0.05 else 0
    minus = -bayes[2] if bayes[2] <= -0.05 else 0
    text_string = r'$ \Gamma ={:.1f}^{{+{:.1f}}}_{{-{:.1f}}}$'.format(
        bayes[1], plus, minus)
    tx = np.sqrt(x[0] * x[-1]) / 1.3
    if is_xlog:
        tx = np.log10(tx)
    # ty = y[0] / 4 if i==0 else y[0] / 4
    ty = np.sqrt(y[0] * y[-1]) * 2
    font = {'family': 'sans-serif',
            'color':  color,
            'weight': 'extra bold',
            'size': 'small',
            }
    ax.text(tx, ty, text_string, fontdict=font,)

    return

def plot_bayes_fit(ax, masses, thresh, color='g'):

    bayes, mass_thr = bayes_get_param(masses, thresh)
    x = np.array([mass_thr, np.max(masses)])
    y = bayes[0] * x**(-bayes[1])
    ax.plot(x, y, '-', color=color)
    plus = bayes[3] if bayes[3] >= 0.05 else 0
    minus = -bayes[2] if bayes[2] <= -0.05 else 0
    text_string = r'$ \Gamma ={:.1f}^{{+{:.1f}}}_{{-{:.1f}}}$'.format(
        bayes[1], plus, minus)
    tx = np.sqrt(x[0] * x[-1]) / 1.3
    # ty = y[0] / 4 if i==0 else y[0] / 4
    ty = np.sqrt(y[0] * y[-1]) * 2
    # ax.text(tx, ty, text_string, )
    font = {'family': 'sans-serif',
            'color':  color,
            'weight': 'extra bold',
            'size': 'small',
            }
    ax.text(tx, ty, text_string, fontdict=font,)


class Sink(Ramses):
    """ Read sink info from outputs, make plots of SFE, sink number, and IMF
    """

    overplot_kroupa = overplot_kroupa_new

    def find_SNs(self, lastoutput=None, mSN=8.0):
        """Return the raw time (in Myr) of all SN explosions
        !!! Apply to the version of RAMSES where all stars never die.
        Return
        ------
        time: {list}, the time of all SNe explosions, not necessary on
        output snapshots, before the lastoutput, when the first SN
        explode
        i:  the outputID when the first SN explode

        Exeptions:
        -1, 0

        """

        if lastoutput is None:
            lastoutput = self.get_last_output()
        # mass, bad = self.get_sink_mass(lastoutput)
        # if bad:
        #     raise SystemExit("This outputID does not have sink particles")
        mass = self.get_sink_masses(lastoutput)
        time = self.get_time(lastoutput)
        lifetime = utilities.mass_to_lifetime(
            MASS_SHIFT * mass)
        age = self.get_age(lastoutput)
        time_ago = age - lifetime
        is_SN = lifetime < age
        time_ago = time_ago[is_SN]
        # optimize this
        time_of_SNs = [time - i_time_ago for i_time_ago in time_ago]
        indices_of_SNs = np.where(is_SN)

        # # No SN found
        # print(self.jobPath, "No SN found!")
        # return -1, 0

        return time_of_SNs, indices_of_SNs

    def find_first_SN(self, mSN=8.0):
        """
        Return
        ------
        i: the output of the first SN. If no SN found, return -1
        time: the time (in Myr) of the first SN (not - tRelax)
        """

        for i in range(1, 100):
            try:
                mass = self.get_sink_masses(i)
            except NoSinkParticle:
                continue
            except FileNotFoundError:
                return -1, -1e10
            lifetime = utilities.mass_to_lifetime(MASS_SHIFT * mass)
            age = self.get_age(i)
            for lf, ag, mas in zip(lifetime, age, mass):
                if mas > mSN and ag > lf:
                    return i, self.get_time(i)

    def metal_yield(self, outputID=None):
        maxOutputID = self.get_last_output()
        if outputID is None:
            outputID = maxOutputID
        else:
            if outputID > maxOutputID:
                outputID = maxOutputID
        masses, bad = self.get_sink_mass(outputID)
        assert not bad, "{}, output_{:05d}".format(self.jobID, outputID)
        massesNotAlive = masses[np.logical_not(self.is_sink_alive(outputID))]
        # print('massesNotAlive =')
        # print(massesNotAlive)
        if not len(massesNotAlive):
            return 0.0
        # scale = 0.3 if para['is_shift_mass'] else 0.0
        metalYield = utilities.mass_to_metal_yield(
            MASS_SHIFT * massesNotAlive)
        return np.sum(metalYield)

    def num_of_SN(self, outputID=None):
        maxOutputID = self.get_last_output()
        if outputID is None:
            outputID = maxOutputID
        else:
            if outputID > maxOutputID:
                outputID = maxOutputID
        return np.sum(np.logical_not(self.is_sink_alive(outputID)))

    def sfe_time_length(self):
        time = self.get_times()
        t_plot = time[time > self.tRelax] - self.tRelax
        total_mass = np.array(self.get_tot_sink_mass())
        assert len(time) == len(total_mass)
        total_mass = total_mass[time > self.tRelax]
        sf_indices = self.find_SF_length(total_mass)
        return t_plot[sf_indices[1]] - t_plot[sf_indices[0]]

    def plot_SFE(self, ax=None, label=None, color=None, yscale='log', ylim=None,
                 is_disp_SN=True, **kwargs):
        """ Plot a single SFE curve. """

        if ax is None:
            ax = plt.gca()

        # Set labels and lims """
        # ax.set_ylabel(r"$f_\star=m_\star/M_{\rm gas}$")
        ax.set_yscale(yscale)
        ax.set_xlabel(r"$t/t_{\rm ff}$")
        ax.set_ylabel(r"$f_*\;(\%)$")
        ax.set_xlim(center.PLOT_RANGE_IN_TFF)
        if ylim is not None:
            ax.set_ylim(ylim)

        t = self.get_t_over_tff()
        total_mass = self.get_tot_sink_mass_list()
        assert len(t) == len(total_mass)
        # plotline, = ax.plot(t_plot / self.tff, total_mass / gas_mass, **kwargs)
        # ax.plot(t_plot/self.tff, total_mass/gas_mass, label=label, color=color,
        #         **kwargs)
        ax.plot(t[t>0], 100 * total_mass[t>0] / self.gas_mass,
                label=label, color=color, **kwargs)
        # plotline, = ax.plot(time, self.total_mass_cut / self.gasMass,
        #                     **kwargs)

        is_disp_sf_indices = False
        if is_disp_sf_indices:
            sf_indices = self.find_SF_length(total_mass)
            if sf_indices[0] < len(total_mass):
                ax.scatter(t_plot[sf_indices]/self.tff,
                        100 * total_mass[sf_indices]/self.gas_mass,
                        marker='*', color=color)

        if not is_disp_SN:
            return
        time_of_SN, _ = self.find_SNs()
        t_raw = self.get_times()
        mass_interp = np.interp(time_of_SN, t_raw, total_mass)
        max_SN = 2
        for i in range(len(time_of_SN)):
            if i > max_SN-1:
                break
            tt = (time_of_SN[i] - self.tRelax) / self.tff
            y = 100 * mass_interp[i] / self.gas_mass
            ax.scatter(tt, y, s=6, c='k', zorder=101)
            height = .05
            pos = ax.transLimits.transform([tt, y])
            # print(tt, y)
            # print(pos)
            # import sys
            # sys.exit()
            pos_x_ax = [pos[0], pos[0]]
            pos_y_ax = [pos[1] - height, pos[1] + height]
            ax.plot(pos_x_ax, pos_y_ax, 'k', zorder=100,
                    transform=ax.transAxes,
                    **kwargs)

            # yy = np.array([mass_interp[i] / height_ratio,
            #                mass_interp[i] * height_ratio])
            # yy = yy / self.gas_mass
            # ax.plot([tt, tt], yy, 'k', zorder=100, **kwargs)

        #TODO: remove this
        # time_of_SN -= self.tRelax
        # # time_of_SN = (time_of_SN - self.tRelax) / self.tff
        # t_plot = t * self.tRelax  # Wrong!
        # logy = np.interp(time_of_SN, t_plot, np.log(total_mass))
        # half_height = np.log(ax.get_ylim()) / 15
        # # print(self.jobid)
        # # print(t_plot)
        # # print(total_mass)
        # # print(time_of_SN)
        # # print(np.exp(logy))
        # max_SN = 2
        # for i in range(len(time_of_SN)):
        #     if i > max_SN-1:
        #         break
        #     t = time_of_SN[i] / self.tff
        #     ax.scatter(t, np.exp(logy[i])/self.gas_mass, s=6, c='k', zorder=101)
        #     # ax.plot([t, t], [sfeSN / scale, sfeSN * scale],
        #     #         **kwargs)
        #     yy = np.exp([logy[i]-half_height, logy[i]+half_height])/self.gas_mass
        #     ax.plot([t, t], yy, 'k', zorder=100, **kwargs)

        # ylim is defined in the outer scope
        # if yscale == 'log':
        #     ax.set_ylim([1e-3, 1])
        # else:
        #     ax.set_ylim(0)
        return

    def find_SF_length(self, sink_masses):

        # time = self.get_times()
        # total_mass = np.array(self.get_tot_sink_mass())
        # assert len(time) == len(total_mass)
        # t_plot = time[time > self.tRelax] - self.tRelax
        # total_mass = total_mass[time > self.tRelax]
        # gas_mass = self.get_gas_mass()

        low = 0.01
        high = 0.5
        # max_m = 0.29 * sink_masses**0.68
        # starting = np.searchsorted(max_m, 20)
        starting = np.searchsorted(sink_masses, low * sink_masses[-1])
        ending = np.searchsorted(sink_masses, high * sink_masses[-1])

        return [starting, ending]

    def display_summary(self, ax=None, **kwargs):
        """ Summary of SFE """

        if ax is None:
            ax = self.ax1
        t = "Summary of {}:\n{: <34}{:16.4f}\n{: <34}{:16.4f}\n" \
            "{:<34}{:16d}\n{: <34}{:16.4f}".format(
            self.jobID, "Final time (Myr):", self.time_cut[-1],
            "Maximum SFE:", self.total_mass_cut[-1] / self.gasMass,
            "Final sink number:", self.num_of_sink_cut[-1],
            "Final mean sink mass(Msun):",
                            self.total_mass_cut[-1] / self.num_of_sink_cut[-1])
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        textx = xlim[0] + 0.02 * (xlim[1] - xlim[0])
        texty = ylim[1] - 0.02 * (ylim[1] - ylim[0])
        ax.text(textx, texty, t, wrap=True,
                verticalalignment='top', size='small', **kwargs)

    def plot_SFE_and_MSM(self, ax=None, hold=False, isDisplayLegend=False,
                         **kwargs):

        self.plotflag = 'SFE'
        self.fetch_all_outputs()
        if ax is None:
            ax = self.ax1

        """  Plot SFE  """
        # self.ax1.set_ylim([10**(-5.2), 10**0.2])
        # if ylim is not None:
        #     self.ax1.set_ylim(10**np.array(ylim))
        ls = '-' # if len(time) < 10 else '-'
        plotline1 = self.plot_SFE(ax=ax) # use self.ax1 as plot axis

        ##### Plot mean sink mass
        if not hold:
            self.ax2 = ax.twinx()
            self.ax2.axis['right'].toggle(all=True)
            self.ax2.set_ylabel(r"$\overline{m}_{\rm s} \;(M_\odot)$")
        y2 = self.total_mass_cut / self.num_of_sink_cut
        legend_y2 = r"$\overline{m}_{\rm s} / M_\odot$"
        ls = ':'  # if len(self.time) < 10 else ':'
        plotline2, = self.ax2.plot(self.time_cut, y2, ls,
                                   color=plotline1.get_color())

        self.plotlines.append([plotline1, plotline2])
        return

        #self.ax2.set_ylim(ylim2)

        # [self.ax2_min, self.ax2_max] = self.ax2.get_ylim()
        # self.ax2.set_ylim([self.ax2_min * 3 / 2, self.ax2_max * 3 / 2])
        # self.ax2.set_ylabel("Number of sink", color='g')
        # self.ax2.tick_params('y', width=1)
        # self.ax2.tick_params('y', colors='g')

        """ Legends """
        if isDisplayLegend:
            legend1 = self.ax1.legend(self.plotlines[0],
                                      [r"$\overline{m}_{\rm s} (M_\odot)$",
                                       legend_y2],
                                      loc='best')
            self.ax1.legend(loc='best')
            # self.ax1.add_artist(legend1)

        """  Set labels  """
        # self.ax1.yaxis.set_ticks_position('both')
        # self.ax1.tick_params('y', width=3)
        # self.ax1.set_ylabel("SFE", color='r')
        return

    def MassFunc(self, ax=None, bins='auto', outputID=None, tff=None,
                 is_x_log_number=False, is_plot_kroupa=False, isFill=False, newsink=False, label=None, plotstyle={}):
        """
        Plot the mass function
        kwargs could be: refPoint = [x_start, x_end], bin_max, bin_min,
        nbins, xlim, ylim, (ALL IN DEX), isDisplayMSM(bool)

        is_x_log_number: If True, plot x as x = log10(masses).
        """

        if ax is None:
            plt.figure()
            ax = plt.gca()
            # ax = self.ax1
        self.plotflag = 'IMF'
        if outputID is None:
            # outputID = self.GetFinalOutputID()
            outputID = self.get_last_output()

        if label=="time":
            # t = self.fetch_single_output(outputID)[0]
            t = self.get_time(outputID)
            t -= self.tRelax
            if self.tff is not None:
                ##### Edit here to change the number of digits in the legend
                tStr = "{:.1f}".format(t) if t < 10.0 else "{:.0f}".format(t)
                label = r"$t = {}\,$Myr $({:.0f} \, t_{{\rm ff}})$". \
                    format(tStr, t / self.tff)
            else:
                label = r"$t = {:.2g}\,$Myr".format(t)

        # mass, bad = self.get_sink_mass(outputID)
        # if bad:
        #     return 0
        try:
            mass = self.get_sink_masses(outputID)
            print(mass.shape)
            if newsink:
#                mass_before = self.get_sink_masses(outputID-5)
                mass_before = self.get_sink_masses(24)
                print(mass_before.shape[0])
                mass = mass[mass_before.shape[0]:]
                print(mass*MASS_SHIFT)
        except NoSinkParticle:
            return 0
        mass *= MASS_SHIFT

        # set the bins if not given
        if type(bins) is str:
            if bins=='auto':
                bins = get_bins(mass)
        n, bins = np.histogram(mass, bins=bins)  # or bins='doane')
        bins = (bins[:-1] + bins[1:]) / 2
        normf = np.log10(bins[1]) - np.log10(bins[0])
        n = np.double(n) / normf
        # self.bins = bins
        # self.n = n

        if isFill:
            ax.fill_between(bins, n, step='pre', label=label, **plotstyle)
        else:
            if is_x_log_number:
                bins = np.log10(bins)
            else:
                ax.set_xscale('log')
            ax.step(bins, n, label=label, **plotstyle)
        ax.set_yscale('log')
        ax.set_xlim([1e-1, 1e3])
        ax.set_ylim([6, 2e3])
        
        if is_plot_kroupa:
            self.overplot_kroupa(ax, self.get_sink_masses(outputID).sum(), m_min=0.08,
                                 m_max=mass.max(),
                                 # m_max=1000,
                                 color='b')


        ##### Labels #####
        ax.set_xlabel(r"$M$ (M$_{\odot}$)")
        ax.set_ylabel(r"${\rm d}N/{\rm d}\log(M/{\rm M}_{\odot})$")
        # if label is not None:
        #     ax.legend(loc='best')

        print("{} output_{:05d} mass range (log): {:.2f} to {:.2f}\n".format
              (self.jobPath, outputID, np.log10(np.min(mass)),
               np.log10(np.max(mass))))

        return 1

    def plot_imf_one_grid_single_line(self, ax, outputid,
                                      # t_max=None,
                                      color='k',
                                      ylim=None,
                                      mass_shift=MASS_SHIFT,
                                      is_plot_kroupa=False,
                                      is_title=True,
                                      is_text=False,
                                      is_met=False,
                                      is_text_cloud_param=True,
                                      **kwargs):
        """Plot one single line of IMF, Kroupa line (optional), set title
        (optional), and add text (optional)
        Parameters

        ----------
        is_text: bool
            Control if or not to add texts at topleft (gas mass and den) and
            topright corner.
        is_met: bool
            Only when is_text = True. Control if or not only add text at
            topright corner.

        """

        # s = sink.Sink(jobid)
        # t = s.get_t_over_tff()
        try:
            mass = self.get_sink_masses(outputid)
        except NoSinkParticle:
            return
        bins = get_bins_fixed()
        shifted_mass = mass * mass_shift
        # if t_max is None:
        #     color = 'k'
        # else:
        #     ti = self.norm_time(self.get_time(outputid))
        #     cmap = matplotlib.cm.get_cmap('plasma')
        #     color = cmap(1 - ti/t_max)
        ylim = [4, 5e3] if ylim is None else ylim
        plot_imf(ax, shifted_mass, bins=bins, ylim=ylim, color=color,
                 **kwargs)

        # add title
        if is_title:
            ax.set_title(center.JOBNAMES[self.jobid])
        # plot a Kroupa IMF
        if is_plot_kroupa:
            self.overplot_kroupa(ax, mass.sum(), m_min=0.08,
                                 m_max=shifted_mass.max(),
                                 # m_max=1000,
                                 color=color)

        if is_text:
            text1, text2 = self.imf_summary_text(outputid)
            bbox = dict(boxstyle="round", fill=False, edgecolor="0.6")
            if not is_met and is_text_cloud_param:
                ax.text(.04, .96, text1, wrap=False, va='top', ha='left',
                        multialignment='left', fontsize=10,
                        transform=ax.transAxes, #bbox=bbox,
                )
            ax.text(.96, .94, text2, wrap=False, va='top', ha='right',
                    multialignment='left', fontsize=9,
                    transform=ax.transAxes, bbox=bbox,
            )

    def imf_fit(self, outputID=None):
        if outputID is None:
            outputID = self.get_last_output()

    def imf_summary_text(self, outputID=None):
        """ Return the texts of the mean sink mass """

        if outputID is None:
            outputID = self.get_last_output()

        # t = self.get_time(outputID) - self.tRelax
        masses = self.get_sink_masses(outputID)
        mass_sh = masses * MASS_SHIFT

        # m_mean = np.mean(mass_sh)
        # m_mean_str = "{:.1f}".format(m_mean) if m_mean < 10 \
        #     else "{:.0f}".format(m_mean)
        # m_max = np.max(mass_sh)
        # m_max_str = "{:.1f}".format(m_max) if m_max < 10 \
        #     else "{:.0f}".format(m_max)

        # text = r"$m_{\rm gas}/$M$_\odot =" + "{}$".format(
        #     utilities.to_latex(self.gas_mass, 0))
        # text += "\n"
        # text += r"$\overline{n}_{\rm gas}/$cm$^{-3} =" + "{}$".format(
        #     utilities.to_latex(self.gas_den * center.den_c_to_m, 0),
        # )
        text = r"${}$ M$_\odot$".format(utilities.to_latex(self.gas_mass, 0))
        text += "\n"
        text += r"${}$ cm$^{{-3}}$".format(
            utilities.to_latex(self.gas_den * center.den_c_to_m, 0))

        # text += r"$M_{{\rm mean}}={}$ M$_\odot$".format(m_mean_str)
        text2 = "$t = {:.1f}$ Myr".format(self.get_time(outputID) - self.tRelax)
        text2 += "\n"
        # m_sum = int(tools.round_sci_digits(np.sum(masses), 2))
        m_sum = int(np.sum(masses))
        text2 += r"$M_*={:d}$ M$_\odot$".format(m_sum)
        text2 += "\n"
        text2 += r"$N_*={:d}$".format(len(mass_sh))
        # text += "\n"
        # text += r"$M_{{\rm max}}={}$ M$_\odot$".format(m_max_str)
        return text, text2

    def imf_summary(self, ax, outputID=None, pos=None, **kwargs):
        """ Write the mean sink mass at the top right corner """

        if outputID is None:
            # outputID = self.GetFinalOutputID()
            outputID = self.get_last_output()

        # t, massSum, num, mass = self.fetch_single_output(outputID)
        t = self.get_time(outputID) - self.tRelax
        # mass, _ = self.get_sink_mass(outputID)
        mass = self.get_sink_masses(outputID)

        # def round3(x):
        #   return round(x, -int(np.floor(np.log10(abs(x))))+1)
        # # "mean sink mass = {:.2f} {}\n".format(meanMass, r"$M_{\odot}$") + \
        # text = r"$t$ = {:.1f} Myr ({:.1f} $t_{{\rm ff}}$)".format(t, t/self.tff) + \
        #     "\nnumber of sinks = {}\n".format(len(mass)) + \
        #     r"$\Sigma m_i = {:.0f} M_{{\odot}}$".format(round3(np.sum(mass)))
        #     # "total sink mass = {:.0f} {}".format(massSum, r"$M_{\odot}$")
        text = "Number of sinks = {}\n".format(len(mass)) + \
            r"Mean sink mass = {:.1g} $M_\odot$".format(np.mean(mass))
        pos = [0.02, 0.98] if pos is None else pos
        ax.text(pos[0], pos[1], text, wrap=True, verticalalignment='top',
                transform=ax.transAxes, **kwargs)

        return

    def quick_plot_IMF(self, ax=None, out=None, masses=None, bin_width=None,
                       bins=None, is_over_dlogm=False, shift=None,
                       is_kroupa=False, kroupa_kwargs={}, **kwargs):

        if out is None:
            out = self.get_last_output()
        if masses is None:
            try:
                m_orig = self.get_sink_masses(out)
            except NoSinkParticle:
                m_orig = []
            except FileNotFoundError:
                return
        else:
            m_orig = masses
        if shift is not None:
            m = np.array(m_orig) * shift
        is_ret = False
        if ax is None:
            is_ret = True
            f, ax = plt.subplots()
        # plot_imf_N_vs_m(ax, m, color='k', bin_width=bin_width, bins=bins)
        if bin_width is not None:
            bins = 10**np.arange(-2, 3.1, bin_width)
        plot_imf(ax, m, bins=bins, is_over_dlogm=is_over_dlogm, **kwargs)

        # write time and outputid tag
        text = "t = {:.3g} Myr".format(self.get_time(out, readinfo=True))
        ax.text(0.02, 0.98, text, va='top', transform=ax.transAxes)
        # ax.set_title('Job' + self.jobid)
        # utilities.add_out_tag(ax, out)

        if is_kroupa: # overplot Kroupa IMF
            overplot_kroupa_new(None, ax, np.sum(m_orig), color='k',
                                isxloged=True, m_max=1000, **kroupa_kwargs,
                                )
        if is_ret:
            return f, ax

    def plot_IMF(self, ax, out=None, bins='auto', is_over_dlogm=False, **kwargs):
        """ Plot MF of an output. Tag with time, out id, m_tot, and N

        Parameters
        ----------
        ax: plot axis
            description
        out: int
            output number. Default: self.get_last_output()
        bins: str or list/array
            description
        is_over_dlogm: bool
            False: y = dN. True: y = dN/dlogm
        kwargs:
            **kwargs for plot_imf()
        """

        if out is None:
            out = self.get_last_output()
        try:
            m = self.get_sink_masses(out)
        except NoSinkParticle:
            m = []
        except FileNotFoundError:
            return

        plot_imf(ax, m, bins, is_over_dlogm=is_over_dlogm, is_text=0, **kwargs)
        text = "t = {:.3g} Myr".format(self.get_time(out, readinfo=True))
        ax.text(0.02, 0.98, text, va='top', transform=ax.transAxes)
        utilities.add_out_tag(ax, out)
        text = r"$m_{{\rm tot}} = {:.3g}$".format(np.sum(m))
        text += "\n"
        text += "$N = {}$".format(len(m))
        ax.text(0.98, 0.98, text, va='top', ha='right', color='k', fontsize='small',
                transform=ax.transAxes)


# IMF plotting
def set_imf_xticks(axes):

    xtick_left = [1e-1, 1e0, 1e1, 1e2]
    xtick_rightmost = [1e-1, 1e0, 1e1, 1e2, 1e3]
    for iaxes in axes[:-1]:
        iaxes.set_xticks(xtick_left)
    axes[-1].set_xticks(xtick_rightmost)

def set_imf_xticks_single(ax, is_remove_rightmost=False):
    if not IsXLogNum:
        xtick = [1e-1, 1e0, 1e1, 1e2, 1e3]
        if is_remove_rightmost:
            xtick = xtick[:-1]
        ax.set_xticks(xtick)
    else:
        if not is_remove_rightmost:
            major = 10**(np.arange(-1, 3.1, 1))
            ax.xaxis.set_ticks(major)
            ax.set_xticklabels(['-1', '0', '1', '2', '3'])
        else:
            major = 10**(np.arange(-1, 2.1, 1))
            ax.xaxis.set_ticks(major)
            ax.set_xticklabels(['-1', '0', '1', '2'])
        ax.xaxis.set_ticks(10**(np.arange(-1, 3.1, 0.2)), minor=True)

def plot_single_imf(ax, jobid, is_set_title=True, is_plot_kroupa=True):
    """ Outdated. Do not use.
    Plot the MFs of one job.
    The indxes must be in reversed order (from big to small) """

    s = sink.Sink(jobid)
    # tff = center.TFFS[jobid[0]]
    # t = (s.get_times() - s.tRelax) / tff
    t = s.get_t_over_tff()
    # indices = pick_indices(t)
    indices = [s.get_last_output()]
    # print(jobid, 'indices =', indices)
    # times = [t[indice] for indice in indices]
    # mass, _ = self.get_sink_mass(indices[-1])
    mass = s.get_sink_masses(indices[-1])
    # bins = sink.get_bins(mass * center.PARAMS['mass_shift'], maxBins=18)
    bins = get_bins_fixed()
    if is_set_title:
        ax.set_title(center.JOBNAMES[jobid])
    t_max = 7.0
    is_x_log_number = False     # Use False. True option is not working properly
    for i, outputID in enumerate(indices):
        # color = colors[i]
        cmap = matplotlib.cm.get_cmap('plasma')
        color = cmap(1 - t[outputID-1]/t_max)
        signal = s.MassFunc(ax=ax, bins=bins, outputID=outputID,
                            label='time', tff=tff, isFill=False,
                            is_x_log_number=is_x_log_number,
                            # plotstyle={'linewidth':3, 'color':color})
                            plotstyle={'color':color})
        if not signal:
            continue

        # plot a Kroupa IMF
        # masses, _ = self.get_sink_mass(outputID)
        masses = s.get_sink_masses(outputID)
        s.overplot_kroupa(m_sum=masses.sum(), ax=ax, color=color)

    text = s.imf_summary_text(indices[-1])
    bbox = dict(boxstyle="round", fill=False, edgecolor="0.6")
    ax.text(.96, .94, text, wrap=False, va='top', ha='right',
            multialignment='left', fontsize='small',
            transform=ax.transAxes, bbox=bbox,
    )
    # ax.legend(fontsize='small', loc='upper right')
    ax.set(xlim=[1e-1, 1e3] if not is_x_log_number else [-1, 3],
           ylim=[1, 8e3])
    if IsXLogNum:
        set_imf_xticks_single(ax)

    # add mass divider
    mass_sort = np.sort(masses)
    mass_cum = np.cumsum(mass_sort)
    # for i in range(5):
    #     thresh = 0.5 - 0.1 * i
    # threshs = [0.1, 0.2]
    threshs = [0.1, 0.3]
    for thresh in threshs:
        indice = np.searchsorted(mass_cum, np.sum(masses) * thresh)
        mass_thr = mass_sort[indice] * MASS_SHIFT
        ax.axvline(mass_thr, color='k', ls='--')
    # choose 0.2 line to use as the starting point of Bayes fitting
    # Replacing all lines below inside this function
    masses_shifted = masses * MASS_SHIFT
    bayes = sink.bayesian_imf_fit_new(masses_shifted, mass_thr)
    # x = np.logspace(np.log10(bayes[0]), np.log10(bayes[1]), 20)
    x = np.array([mass_thr, np.max(masses_shifted)])
    y = bayes[0] * x**(-bayes[1])
    color = 'g'
    ax.plot(x, y, '-', color=color)
    plus = bayes[3] if bayes[3] >= 0.05 else 0
    minus = -bayes[2] if bayes[2] <= -0.05 else 0
    # ts = "$-\Gamma="
    # ts = "-Γ $=$"
    ts = "jiba $="
    ts += "{:.1f}^{{+{:.1f}}}_{{-{:.1f}}}".format(bayes[1], plus, minus)
    ts += "$"
    ts = "%r" % ts
    raise SystemExit(ts)
    tx = np.sqrt(x[0] * x[-1]) / 1.3
    # ty = y[0] / 4 if i==0 else y[0] / 4
    ty = np.sqrt(y[0] * y[-1]) * 1.3
    ax.text(tx, ty, ts, color=color, ) #fontsize='small')

    # # add the Baysian fit
    # bays_fit = True
    # if not bays_fit:
    #     return
    # bayes_list = s.imf_fit(indices[-1])
    # # bayes_list = [mass_thr]

    # # colors = ['magenta', 'orangered']
    # # colors = [CMAP(1), CMAP(0.6)]
    # ls = ['-', '--']
    # for i in range(len(bayes_list)):
    #     bayes = bayes_list[i]
    #     x = np.logspace(np.log10(bayes[0]),
    #                     np.log10(bayes[1]), 20)
    #     y = bayes[2] * x**(-bayes[3])
    #     # color = 'm'
    #     # ax.plot(x, y,
    #     #         # color=colors[i],
    #     #         alpha=0.6,
    #     #         # color="0.6",
    #     #         # color='k',
    #     #         color=color,
    #     #         # ls=ls[i],
    #     #         linewidth=3,)
    #     color = 'g'
    #     ax.plot(x, y, '-', color=color)
    #     num = len(bayes_list)
    #     subscript = '' if num==1 else '_{}'.format(i+1)
    #     plus = bayes[5] if bayes[5] > 0.01 else 0
    #     text_string = r"$\alpha{}={:.1f}^{{+{:.1f}}}_{{{:.1f}}}$".format(
    #         subscript, bayes[3], plus, bayes[4])
    #     # ax.text(10**(x[0]+0.2), 10**(a*x[0]+b), text_string)
    #     # tx = x[0] * 3 if i==0 else x[0] * 3
    #     tx = np.sqrt(x[0] * x[-1]) / 1.3
    #     # ty = y[0] / 4 if i==0 else y[0] / 4
    #     ty = np.sqrt(y[0] * y[-1]) * 1.3
    #     ax.text(tx, ty, text_string,
    #             color=color,
    #             fontsize='small',
    #             # color=colors[i],
    #     )

    return

def plot_imf_from_output(ax, r, out, isbay=False):
    """
    Parameters
    ----------
    ax: plt axis
        the axis to plot on
    r: Sink class object
        Instance of ramses.Ramses
    out: int
        output ID
    isbay: bool
        If True, plot Bayesian fit to the power-law slopw
    """

    Bw = 0.2
    Bay_th = 0.2

    # check if the star is in zoom region
    try:
        r.get_star_pos(out)
    except NoSinkParticle:
        print("No sink particle found. Not plotting")
        return
#     print("Out {}, Are all sinks inside zoom region? {}".format(out,
#         all(np.linalg.norm(poss - zoomc_pc, axis=1) < 0.5)))
    r.quick_plot_IMF(ax, out, bin_width=Bw)
    if not isbay:
        return
    try:
        masses = r.get_sink_masses(out)
    except FileNotFoundError:
        return
    except NoSinkParticle:
        return
    if len(masses) < 5:
        return
    sink.plot_bayes_fit_any_width(ax, masses, Bay_th, color='g', is_xlog=1, yscaleup=Bw)


# IMF Sampling (from imf_sample_pkg.py and imf_sample.py)
def sample_filename(jobid, outputID, version, kind=0):
    """ When version='' and kind=0, returns the original sink filename
    Parameters
    ----------
    jobid: str
        e.g. "3.3.2", or "XL-VC"
    version: str
        description
    kind: int
        0: sink, 1: Kroupa, 2: Salpeter

    """

    root_dir = "../work_outputs/imf_sample/"
    utilities.makedir(root_dir + "Job" + jobid)
    fname = root_dir + "Job{}/output_{:05d}{}{}.txt".format(
        jobid, outputID, version, KIND_DICT[kind])
    return fname

def pow_int(a, b, beta):
    """Return the integral of x^-alpha from a to b. beta = 1 - alpha.

    """

    # beta = 1.0 - alpha
    return (b**beta - a**beta) / beta;

def pow_formula(x, a, b, beta):

    return (a**beta + x * (b**beta - a**beta))**(1. / beta)

def kroupa_int_two_side(a, b, alpha1, alpha2):
    """ Sample Kroupa IMF. The coefficient is set by assuming alpha1 =
    0.3, and alpha2 = 1.3 """

    assert a > 0.0 and a < 0.5 and b > 0.5;
    beta1 = 1. - alpha1
    beta2 = 1. - alpha2
    return pow_int(a, 0.5, beta1) + 0.5 * pow_int(0.5, b, beta2);

def pick_salpeter(a, b, alpha=2.35, rand=None):
    """For salpeter, rand is set for faster sampling based on numpy
arrays.

    """

    assert alpha > 1.0
    assert a > 0.0 and b > a, "a = {:.2g}, b = {:.2g}".format(a, b)
    beta = 1. - alpha
    if rand is None:
        x = random()
    else:
        x = rand
    # return (a**beta * (1. - x) + b**beta * x)**(1. / beta)
    return pow_formula(x, a, b, beta)

def sample_two_powerlaw(rand, a1, a2, m1, m2, m3):
    """Given a random number from 0 to 1, return a sampled mass based the
following distribution: a broken power law from m1 to m2, m2 to m3,
with negative slopes a1 and a2

    """

    # m^a1, if m1 < m < m2
    # k2 m^a2, if m2 < m < m3, where k2 = m2^(a1 - a2)
    int1 = pow_int(m1, m2, 1 - a1)
    int2 = m2**(a2 - a1) * pow_int(m2, m3, 1 - a2)
    intab = int1 + int2
    if rand <= int1 / intab:
        # return (x * beta1 * intab + a**beta1)**(1. / beta1)
        return (rand * intab * (1-a1) + m1**(1 - a1))**(1 / (1 - a1))
    else:
        # return ((x * intab - pow_int(a, 0.5, beta1)) * 2. * beta2 + .5**beta2)**(1. / beta2)
        return ((rand * intab - int1) * (1-a2) / m2**(a2-a1) + m2**(1 - a2))**(1. / (1 - a2))

def sample_imf(pdf, m_tot, m_mean=0.2):
    """ Sample a cluster of stars with total mass m_tot with a given PDF

    Parameters
    ----------
    pdf: function
        The MF
    m_tot: double
        total mass to sample to (e.g. the mass of a sink particle to fragment)
    m_mean: double (not used)
        The estimated mean mass from pdf
    """

    msum = 0
    masses = []
    while True:
        m = pdf(random())
        msum += m
        if msum > m_tot:
            masses.append(m_tot - msum + m)
            break
        masses.append(m)
    return masses

    n_est = int(m_tot / m_mean * 3)
    if n_est == 0:
        raise SystemExit('Fail to sample m_tot = {:.4f}. '
                         'Try decreasing m_mean'.format(m_tot))
    rands = np.random.random(n_est)
    mass_samples = np.vectorize(pdf)(rands)
    mcum = np.cumsum(mass_samples)
    idx = np.argmax(mcum > m_tot)
    if idx == 0:
        if mass_samples[0] > m_tot:
            return [m_tot]
        else:
            raise SystemExit('Fail to sample m_tot = {:.4f}. '
                             'Try decreasing m_mean'.format(m_tot))
    return mass_samples[:idx]

def kroupa(m):
    if m < 0.5:
        return 2. * m**-1.3
    else:
        return m**-2.3

def pick_kroupa(a, b, rand=None):
    """ Pick a Kropam sampling
    rand is set for faster sampling based on numpy arrays.
    """

    assert a > 0.0 and b > a, "a = {:.2g}, b = {:.2g}".format(a, b)
    alpha1 = 1.3
    alpha2 = 2.3
    beta1 = 1. - alpha1
    beta2 = 1. - alpha2
    x = random()
    if b <= 0.5:
        return pow_formula(x, a, b, beta1)
    elif a >= 0.5:
        return pow_formula(x, a, b, beta2)
    else:
        intab = pow_int(a, 0.5, beta1) + 0.5 * pow_int(0.5, b, beta2)
        if x < pow_int(a, 0.5, beta1) / intab:
            return (x * beta1 * intab + a**beta1)**(1. / beta1)
        else:
            return ((x * intab - pow_int(a, 0.5, beta1)) * 2. * beta2
                    + .5**beta2)**(1. / beta2)

def lognormal(logm):
    """
    A lognormal IMF.  The default parameters correspond to the Chabrier IMF
    """

    # return scale * np.exp(-1.*(np.log10(m)-np.log10(offset))**2/(2*width**2))
    return 0.086 * np.exp(-1.*(logm - np.log10(0.22))**2/(2*0.57**2))

def chabrier2003(logm):
    """ The chabrier 2003 system IMF
    """
    if logm < 0:
        return lognormal(logm) / lognormal(0)
    else:
        return 10.**(-1.3 * logm)

def plot_chab03(ax, msum, m_min=0.1, m_max=100, **kwargs):

    assert m_min > 0 and m_min < 1 and m_max > 1
    scale = integrate.quad(lambda logm: 10**logm * chabrier2003(logm),
                           np.log10(m_min), np.log10(m_max))[0]
    norm = msum / scale
    masses = np.logspace(np.log10(m_min), np.log10(m_max), 100)
    logm = np.log10(masses)
    vch03 = np.vectorize(chabrier2003)
    ax.plot(masses, norm * vch03(logm), **kwargs)

def chabrier05(m):
    """ The Chabrier 2015 single-star IMF. Reference: Eq. 3 of  Krumholz:2014.
    The value of mc is missing. Refer to the original paper: Chabrier:2015

    :return float: dndlogm
    """

    alpha = 2.35
    sigma = 0.55
    mb = 1.
    mc = 0.2
    const1 = 0.093
    const2 = 0.041
    if m <= mb:
        dndlogm = const1 * np.exp(-(np.log10(m) - np.log10(mc))**2 /
                                  (2 * sigma**2))
    else:
        dndlogm = const2 * (m / mb)**(1. - alpha)
    return dndlogm

def chabrier05_logm(logm):
    """
    :return: float: dndlogm
    """

    return chabrier05(10**logm)

def plot_chab05(ax, msum, m_min=0.1, m_max=100, isxloged=0, **kwargs):

    assert m_min > 0 and m_min < 1 and m_max > 1
    scale = integrate.quad(lambda logm: 10**logm * chabrier05_logm(logm),
                           np.log10(m_min), np.log10(m_max))[0]
    norm = msum / scale
    masses = np.logspace(np.log10(m_min), np.log10(m_max), 100)
    chabrier05_v = np.vectorize(chabrier05)
    x = masses if not isxloged else np.log10(masses)
    ax.plot(x, norm * chabrier05_v(masses), **kwargs)

def test_mean():
    """ A driver to test the mean mass of sampling. For Salpeter, INTERP's ~0.6.
    """

    N = 100000
    rands = np.random.random(N)
    print('rand mean =', np.mean(rands))
    mass = pick_salpeter(0.08, 100, rand=rands)
    # mass = [pick_salpeter(99, 100) for i in range(N)]
    mmean = np.mean(mass)
    mmean_t = 0.283
    success = np.abs((mmean - mmean_t) / mmean_t) < 1e-1
    print("Salpeter: mean, max, min: ", mmean, np.max(mass), np.min(mass))
    print("Test mean: {}".format({True: "Success", False: "Fail"}[success]))

    mass = [pick_kroupa(0.08, 100) for i in range(N)]
    mmean = np.mean(mass)
    mmean_t = 0.574
    success = np.abs((mmean - mmean_t) / mmean_t) < 1e-1
    # mass = [pick_salpeter(99, 100) for i in range(N)]
    print("Kroupa: mean, max, min: ", np.mean(mass), np.max(mass), np.min(mass))
    print("Test mean: {}".format({True: "Success", False: "Fail"}[success]))

def sample_masses(masses, kind='sal', m_min=0.1, m_lower=0.08,
                  m_max=-1, **kwargs):
    """
    Parameters
    ----------
    masses: array
        The masses to sample from
    kind: str
        'sal' for salpeter and 'kro' for Kroupa sampling recipe
    m_min: double
        The lower bound of IMF. 0.0 < m_min < 0.5
    m_lower: double
        Only useful when m_min is a string. m_min = max(m_lower, x * m_sink)
    m_max: double
        If m_max < 0: use m_sink as m_max
    **kwargs: passed to pick_salpeter (not to pick_kroupa) to set alpha.

    Return
    ------
    mass_samp: array
        The sampled masses
    # n: list of integers
    #     Counts of the histogram of the sampled masses
    # bins: list of double
    #     Bins of the histogram of the sample masses
    """

    mass_samp = []
    # m_min = 0.08
    # m_min = 0.1
    for massi in masses:
        msum = 0.0
        m_max_local = massi if m_max < 0 else m_max
        if type(m_min) is float:
            m_min_local = m_min
        elif type(m_min) is str:
            m_min_local = max(float(m_min) * massi, m_lower)
        if m_max_local < m_min_local:
            continue
        while True:
            if kind == 'sal':
                # salpeter, may pass an alpha
                m = pick_salpeter(m_min_local, m_max_local, **kwargs)
            elif kind == 'kro':
                # kropua, may not pass an alpha
                m = pick_kroupa(m_min_local, m_max_local)
            else:
                raise SystemExit("kind not in ['sal', 'kro']")
            msum += m
            if msum > massi:
                break
            mass_samp.append(m)
    # n, bins = np.histogram(mass_samp, bins=10.0**np.arange(-1, 3, 0.2))
    # return n, bins
    return mass_samp

def plot_sp_imf_one_job(ax, jobid, gammas, m_min=0.1, m_lower=0.08,
                        sample_label_ext=''):
    """
    Parameters
    ----------
    ax: plt axis
        the axis to plot in
    jobid: str
        jobid, e.g. '4.3.1'
    # is_sal: bool
    #     If True, add Salpeter sampling along with Kroupa sampling.

    """

    bins = 10**np.arange(-1, 3, 0.2)
    # bins = imf.get_bins_fixed()

    s = Ramses(jobid)
    last = s.get_last_output()
    mass = s.get_sink_masses(last)
    # fn = sample_filename(jobid, last, '')
    # mass = np.loadtxt(fn)

    # Plot the original SMF (Sink Mass Functin) and sSMF (shifted SMF)
    # color = plotutils.indigo
    color = 'C2'
    plot_imf(ax, mass, bins=bins, label='Sink', color=color,
                  linewidth=1)
    plot_imf(ax, 0.4 * mass, bins=bins, linestyle='--',
                  # dashes=plotutils.get_dashes(1),
                  linewidth=1,
                  label=r'0.4 $\times$ Sink',
                  color=color)

    # if is_sal:
    #     mass = np.loadtxt(fn + "_min{:d}_sal".format(MinType))
    #     plot_imf(ax, mass, m_min=m_min, label='Sal Sample', color='C1')
    #     sink.overplot_salpeter(ax, mass.sum(),
    #                            m_min=m_min,
    #                            # m_max=max(mass),
    #                            m_max=100,
    #                            label='Salpeter IMF',
    #                            # textxy=[1e-1, 5e5],
    #                            color='C1')

    # Do Salpeter sampling for multiple times and take the average
    def sample_alpha(alpha, **kwargs):
        # alpha = 2.35
        kind = 'sal'
        n_sample = 20
        bins_lower = -1 if m_lower > 0.07 else -2
        bins = 10**np.arange(bins_lower, 3, 0.2)
        n_mean = np.zeros(len(bins) - 1)
        for i in range(n_sample):
            masses_samp = sample_masses(mass, m_min=m_min, m_lower=m_lower,
                                        kind=kind, alpha=alpha)
            ni, _ = np.histogram(masses_samp, bins=bins)
            n_mean += ni
        n_mean /= n_sample
        plot_imf_from_nbins(
            ax, n_mean, bins,
            # dashes=plotutils.get_dashes(3),
            # label=r'Sampled with $\Gamma = {:.1f}$'.format(alpha-1))
            label=r'$\Gamma = {:.2f}$ sample{}'.format(alpha-1,
                                                       sample_label_ext),
            # linewidth=1,
            **kwargs,
        )
    colors = ['C1', 'C0']
    for count, gamma in enumerate(gammas):
        sample_alpha(gamma + 1., color=colors[count])

    # Overplot analytic IMF line
    m_max = np.max(mass) * MASS_SHIFT
    # sink.overplot_kroupa_new(1, ax, np.sum(mass), m_min=0.1, m_max=m_max,
    #                          color='0.2',
    #                          # label='Mass-normalized Kroupa',
    #                          label='Kroupa',
    #                          zorder=-10)
    # plot_chab03(ax, np.sum(mass), m_max=m_max, color='C0',
    #             label='Chabrier03', zorder=-9,
    #             dashes=plotutils.get_dashes(3),
    #             )
    plot_chab05(ax, np.sum(mass), m_max=m_max, color='k', #color='C0',
                label='Chabrier05', zorder=-9,
                # dashes=plotutils.get_dashes(3),
    )

    ax.set(xlim=[1e-1, 1e3], ylim=[2e0, 1e5])


def imf_compare(m1, m2, ax, c1=None, c2=None, groups1=None, groups2=None):
    """Compare two MFs in a top-bottom view.

    Args:
        m1: (array )
        m2: (array )
        ax: plt axis
        c1: color for m1 (default None)
        c2: color for m2 (default None)
        groups1: (list of lists) particle groups for m1 (default None)
        groups2: (list of lists) particle groups for m2 (default None)

    Returns:
        f, ax

    """

    plot_imf(ax, m1, is_over_dlogm=False, is_text=1,
             # label='l=14',
             color=c1, ylim=[7e-1, 1e1],)
    plot_imf(ax, m2, is_over_dlogm=False, is_text=1,
             # label='l=18',
             color=c2, ylim=[7e-1, 1e1],)

if __name__ == '__main__':
    pass
