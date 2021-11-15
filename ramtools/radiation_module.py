# from utilities import *
import os
from math import log10
import numpy as np
from scipy import integrate, interpolate
import astropy.constants as C

from .utilities import TLIM
from .utilities import read_quant_from_ramses_info as rinfo
from . import center, center2, units
# import ramses

# TAUH_TO_TAUD = 7.9e-4        # not used anymore.
# TAUH_TO_TAUD = 5e-5

# QSCALE = 1e-44
# logtauLim = [-1.0, 2.0] # TODO fix logtauLim
# LOG_TAO_LIM = [-1.0, 2.0] # TODO fix logtauLim

def QVacca(Mass):
    """ Return VaccaQHI in 1e44 s^-1. Mass range: 1 - 10000. """

    # s**(-1) then normalised in code units in read_stellar
    stf_k=9.634642584812752e48
    # Msun then normalised in code units in read_stellar
    stf_m0=2.728098824280431e1
    stf_a=6.840015602892084e0
    stf_b=4.353614230584390e0
    stf_c=1.142166657042991e0
    return 1. / center.unit_lum * (stf_k * (Mass /stf_m0)**stf_a / (1 + (
            Mass/stf_m0)**stf_b)**stf_c)

def QVaccaRaw(m):
    # return QVacca in s^-1
    return QVacca(m) * center.unit_lum

def QVacca_sp(m):
    """ Modified Vacca, single power law """
    stf_k=9.634642584812752e48
    stf_m0=2.728098824280431e1
    stf_a=6.840015602892084e0
    stf_b=4.353614230584390e0
    stf_c=1.142166657042991e0
    return stf_k * (m / stf_m0)**(stf_a - stf_b * stf_c)

def scluster_fit(mtot):
    return 4.6e46 * mtot**1.0746

def rescale_ratio_do_not_use(mass):
    """Calculate the rescaling ratio -- scluster / wcluster, where
    wcluster = Sigma ssink -- of a given sink mass function. This
    rescaling is used in all the simulations shown in He2019 and
    He2019a.

    """

    scluster = scluster_fit(mass.sum())
    wcluster = np.sum(QVaccaRaw(0.3 * mass))
    return scluster / wcluster

def luminosity(masses, recipe, info_path, masstot=None):
    """
    Args
        masses (list): list of sink mass in Msun
        recipe (str): one of ["frig_he", "sp"]
        info_path (str): of the of info_000##.txt file, which is used to
            get unit_d and unit_l
        masstot (float): total mass. When masstot = None (default), this
            is np.sum(masses)
    """

    assert recipe in["frig_he", "sp"]
    if recipe == "frig_he":
        scale_d = rinfo(info_path, "unit_d")
        scale_l = rinfo(info_path, "unit_l")
        scale_msun = scale_d * scale_l**3 / units.Msun
        print(f"scale_msun = {scale_msun}")
        sraw = QVaccaRaw(0.3 * masses / scale_msun)  # wrong scale/units use in ramses_frig or ramses_frig_he
        wcluster = sraw.sum()  # wrong scale used in _frig
        if masstot is None:
            masstot = masses.sum()
        scluster = scluster_fit(masstot)
        raw_scale = scluster / wcluster
        scales = raw_scale if raw_scale < 1.0 else 1.0
        return sraw * scales
    if recipe == "sp":
        return QVacca_sp(0.4 * masses)

def poly(logm, a):
    ret = 0.0
    for i, ai in enumerate(a):
        ret += ai * logm**i
    return ret

def QHeI(mass):
    """ Return Schaerer QHeI * center.unit_lum. """

    a = [16.05, 48.87, -24.70, 4.29]
    logm = np.log10(mass)
    if mass < 6:
        return 0.
    elif mass < 150:
        return 10 ** poly(logm, a) * center.qscale
    else:
        if mass > 1e4:
            hardness = 0.0
        else:
            hardness1 = poly(np.log10(150), a) - np.log10(QVaccaRaw(150))
            hardness2 = 0.0
            x1 = np.log10(150)
            x2 = 4
            hardness = hardness1 + (hardness2 - hardness1) / \
                       (x2 - x1) * (logm - x1)
        return QVacca(mass) * 10**hardness

def QHeI_raw(mass):
    return QHeI(mass) * center.unit_lum

def QHeII(mass):
    logm = np.log10(mass)
    a = [34.65, 8.99, -1.4]
    x0 = np.log10(150)
    if mass < 6:
        return 0.
    elif mass < 150:
        return 10 ** poly(logm, a) * center.qscale
    elif mass < 7e4:
        return 10 ** (poly(x0, a) + (logm - x0) * (a[1] + 2*a[2]*x0)) * center.qscale
    else:
        return QVacca(mass)

def QHeII_raw(mass):
    return QHeII(mass) * center.unit_lum

def sigmaHI(nu):
    """ [nu] = eV
    Return sigma (cm^2)
    """

    sigma = 0.0
    if nu > 13.6 - 1e-10:
        sigma = center.sigma0 * (nu / 13.6)**-3
    return sigma

def sigmaHeI(nu):
    """ [nu] = eV """

    epsilon = 1e-10
    if nu > 65.4 - epsilon:
        return sigmaHI(nu) * (37.0 - 19.1 * (nu / 65.4)**(-0.76))
    if nu > 24.6 - epsilon:
        return sigmaHI(nu) * (6.53 * (nu / 24.6) - 0.22)
    return 0.0

def sigmaHeII(nu):
    """ [nu] = eV """

    if nu > 54.4 - 1e-10:
        return center.sigma0 / 4 * (nu / 13.6 / 4)**-3
    return 0.0

# def sigmaAll(nu):
#     return sigmaHI(nu) + sigmaHeI(nu) + sigmaHeII(nu)

# sigma1 = sigmaAll(13.6)
# sigma2 = sigmaAll(24.6)
# sigma3 = sigmaAll(54.4)
# print(sigma1, sigma2, sigma3)

# def tau2ef(tau):
#     """ tau is an array of optical depth """
#     return np.exp(-tau)

class DustExtinction():

    def __init__(self, Z=1.0, cloud='SMC'):
        Z0 = 0.2 if cloud=='SMC' else 0.32
        self.ratioZ = Z / Z0
        cur_dir = os.path.dirname(os.path.realpath(__file__))
        param = np.loadtxt('{}/../work/dust_params_{}.txt'.format(cur_dir, cloud), delimiter='\t')
        self.paraml = param[:, 0]
        self.a = param[:, 1]
        self.b = param[:, 2]
        self.p = param[:, 3]
        self.q = param[:, 4]
        if cloud=='SMC':
            self.sigma0 = 1e-22
        else:
            self.sigma0 = 3e-22

    def get_sigma(self, lam):
        """
        Parameters
        ----------
        lam: double or array
             wavelength in um
        """

        sigma = 0.0
        for i in range(7):
            x = lam / self.paraml[i]
            sigma += self.a[i] / (x**self.p[i] + x**(-self.q[i]) + self.b[i])
        sigma *= self.sigma0
        return sigma

    def get_tau(self, lam, N_H):
        """ Formula: N_dust = Z / Z_0 * N_H """
        sigma = self.get_sigma(lam)
        N_dust = self.ratioZ * N_H
        tau = N_dust * sigma
        return tau


def dust_model(model):
    if model==1:
        return DustExtinction(Z=0.1, cloud='SMC')
    elif model==2:
        return DustExtinction(Z=0.2, cloud='SMC')
    elif model==3:
        return DustExtinction(Z=1.0, cloud='SMC')
    raise SystemExit('Unrecognized dust model: {}'.format(model))

def find_bin(nu):
    epsilon = 1e-10
    if nu < 13.6 - epsilon:
        return 0
    if nu < 24.6 - epsilon:
        return 1
    if nu < 54.4 - epsilon:
        return 2
    return 3

def colden_to_tau_ev(eVs, colden_HI=0.0, colden_HeI=0.0,
                     colden_HeII=0.0):
    """ Only apply to gas (not dust)! eVs can be scaler or array
    Remember to multiply tau by self.n_colden_H to get actual tau
    """

    tau = np.zeros((*colden_HI.shape, len(eVs)))
    for i, eV in enumerate(eVs):
        taui = colden_HI * sigmaHI(eV) \
               + center.nH_to_nHe * colden_HeI * sigmaHeI(eV) + \
               + center.nH_to_nHe * colden_HeII * sigmaHeII(eV)
        tau[:, :, i] = taui
    return tau


# def NH_to_taud_scaler(cd_H, hu, dust_model):
#     """ Convert N_H (cm^-2) to tau_d. """
#     NH21 = cd_H * 1e-21
#     const = 5
#     return const * NH21

# def NH_to_taud(cd_H, hus, dust_model):
#     """cd_H is usually a ndim=2 array (stars and directions in each
#     dimension).
#     Parameters
#     ----------
#     cd_H: array
#         description
#     hus: array
#         description
#     dust_model: int (1, 2, ...)
#         Dust models.

#     """

#     dust = DustExtinction(Z=1.0, cloud='SMC')
#     ret = np.zeros((*cd_H.shape, len(hus)))
#     for i, hu in enumerate(hus):
#         # ret[:, :, i] = NH_to_taud_scaler(cd_H, hu, dust_model)
#         ret[:, :, i] = dust.get_tau(hu, cd_H)
#     return ret

class InterpolateTau():
    # def interpolate_tau(logtaulim=[-1, 2.5], Tlim=[5000, 50000]):
    """
    Perform realistic tau-to-get_ef conversion
    :param tau: array, ndim=1. log10 of tau
    :return: get_ef: array, ndim=1. log10 of get_ef
    """

    e = C.e.value
    h = C.h.value
    k = C.k_B.value
    c = C.c.value
    vH = 13.6 * e / h
    # logefV = None

    def planck(self, nu, T):
        """ Return the planck function """
        planck = 2 * self.h * nu**3 / self.c**2
        planck /= np.exp(self.h * nu / self.k / T) - 1
        return planck

    def deno(self, x, T):
        """ x = nu/nuH """
        return self.planck(x * self.vH, T) / (self.h * x * self.vH)

    def nume(self, x, tau, T, field):
        """ x = nu/nuH """
        # return self.deno(x, T) * np.exp(-tau * x**-3)
        # assert field in ['HI', 'HeI', 'HeII']
        if field == 'HI':
            return self.deno(x, T) * np.exp(-sigmaHI(x*center.eHI) / sigmaHI(
                center.eHI) * tau)
        elif field == 'HeI':
            return self.deno(x, T) * np.exp(-sigmaHeI(x*center.eHI) / sigmaHeI(
                center.eHeI) * tau)
        else:
            return self.deno(x, T) * np.exp(-sigmaHeII(x*center.eHI) / sigmaHeII(
                center.eHeII) * tau)

    def logef(self, tau, T, field):
        """ the calculated escape fraction
        Return:
        -------
        log(ef)
        """

        def numeTau(x):
            return self.nume(x, tau, T, field)
        def denoT(x):
            return self.deno(x, T)
        vMax = 5
        if field == 'HI':
            vMin = 1
        elif field == 'HeI':
            vMin = center.eHeI / center.eHI
        else:
            vMin = center.eHeII / center.eHI
        down = integrate.quad(denoT, vMin, vMax)
        top = integrate.quad(numeTau, vMin, vMax)
        return np.log10(top[0] / down[0])

    # def logef(tau, T):
    #     return -tau / np.log(10)

    def interpolate(self, logtauLim=[-1, 2.5], Tlim=[5000, 50000],
                    isNewStyle=True, field='HI'):
        N = 100
        # N = 10
        # yT = np.logspace(*np.log10(np.array(Tlim)), N)
        yT = np.logspace(np.log10(Tlim[0]), np.log10(Tlim[1]), N)
        ylogT = np.log10(yT)
        xtau = np.logspace(logtauLim[0], logtauLim[1], N)
        xlogtau = np.log10(xtau)
        xtau, yT = np.meshgrid(xtau, yT)
        xlogtau, ylogT = np.meshgrid(xlogtau, ylogT)
        # print(self.logef(1, 1000))
        # logefV = np.vectorize(self.logef())
        # zlogef = logefV(xtau, yT)
        nx, ny = xtau.shape
        # zlogef = np.array([[self.logef(xtau[i, j], yT[i, j])
        #                     for i in range(nx)] for j in range(ny)])
        zlogef = self.logefV(xtau, yT, field)
        if isNewStyle:
            # return logtau, logT --> logef / tau
            return interpolate.bisplrep(xlogtau, ylogT, zlogef / xtau)
        else:
            # return logtau, logT --> logef
            return interpolate.bisplrep(xlogtau, ylogT, zlogef)

    def do_interpolation(self):
        self.logefV = np.vectorize(self.logef)
        print("Doing interpolation for T =", TLIM)
        # intTau = InterpolateTau()
        self.logtau2logefovertau = \
            self.interpolate(center2.LOG_TAO_LIM, TLIM)
        # self.logtau2logefovertau_HeI = \
        #     self.interpolate(center2.LOG_TAO_LIM, TLIM, field='HeI')
        # self.logtau2logefovertau_HeII = \
        #     self.interpolate(center2.LOG_TAO_LIM, TLIM, field='HeII')
        print("Done with interpolation.")

    # def ef_inte(self, tau, T, field):
    #     """
    #     TODO: here
    #     """

    #     return np.exp(-tau)     # not used

    #     # try:
    #     #     self.logefV
    #     # except AttributeError:
    #     #     self.do_interpolation()
    #     # assert field in ['HI', 'HeI', 'HeII']

    #     if field == 'HI':
    #         return 10**(interpolate.bisplev(
    #             log10(tau), log10(T), self.logtau2logefovertau) * tau)
    #     elif field == 'HeI':
    #         return 10**(interpolate.bisplev(
    #             log10(tau), log10(T), self.logtau2logefovertau_HeI) * tau)
    #     else:
    #         return 10**(interpolate.bisplev(
    #             log10(tau), log10(T), self.logtau2logefovertau_HeII) * tau)

    def tau_to_ef_H(self, tau, T):
        """ tau and T are doubles """

        return 10**(interpolate.bisplev(log10(tau), log10(T),
                                        self.logtau2logefovertau) * tau)

    def tau_to_ef(self, tau, temps, field):
        """
        Convert tau to log(ef). np.ndim(tau) must be 2
        :param tau:  tau
        :param T:       T in Kelvin
        :return: log10 of ef
        """

        assert np.ndim(temps) == 1
        assert np.ndim(tau) == 2
        assert tau.shape[0] == len(temps)
        ef_tau = np.zeros(np.shape(tau))

        if field != 'HI': # 'HeI' or 'HeII'
            return np.exp(-tau)

        # assert tau.shape[0] == len(temps)
        # ef_tau = np.zeros(np.shape(tau))
        for i in range(tau.shape[0]):
            if temps[i] < TLIM[0] or temps[i] > TLIM[1]:
                ef_tau[i, :] = np.exp(-tau[i, :])
                continue
            for j in range(tau.shape[1]):
                if tau[i, j] < 10**center2.LOG_TAO_LIM[0]:
                    ef_tau[i, j] = np.exp(-tau[i, j])
                elif tau[i, j] > 10**center2.LOG_TAO_LIM[1]:
                    ef_tau[i, j] = 0.0
                else:
                    # ef_tau[i, j] = self.ef_inte(tau[i, j], temps[i], field)
                    ef_tau[i, j] = self.tau_to_ef_H(tau[i, j], temps[i])
        return ef_tau


############## Do the interpolation ##################
# logtauLim = [-1.0, 3.0] # TODO fix logtauLim
# print("Doing interpolation for T =", Tlim)
# intTau = InterpolateTau()
# logtau2logefovertau = intTau.interpolate(logtauLim, Tlim)
# logtau2logefovertau_HeI = intTau.interpolate(logtauLim, Tlim, field='HeI')
# logtau2logefovertau_HeII = intTau.interpolate(logtauLim, Tlim, field='HeII')
# print("Done with interpolation.")

# class Radiation(ramses.Ramses):
#
#     def colden_to_tau0(self, field, colden_HI, colden_HeI=None,
#                        colden_HeII=None, eV=13.6):
#         """ Input colDenHI in pymses code units (unit_density * unit_length).
#         Return tau_0, which is the tau at 13.6, 24.6, or 54.4 eV """
#
#         # if field == 'HI':       # Sep 8, 2019
#         #     return colden_HI * self.n_colden_H * sigmaHI(eV)
#         if field == 'HI':
#             return colden_HI * self.n_colden_H * sigmaHI(13.6)
#         if field == 'H':
#             cd_H = colden_HI
#             cd_H *= self.n_colden_H
#             NH21 = cd_H * 1e-21
#             const = 5
#             return const * NH21
#         # if field in ['HI', 'H']:
#         #     tau0H = colden_HI * self.n_colden_H * sigmaHI(13.6)
#         #     if field == 'HI':
#         #         return tau0H
#         #     else:
#         #         return tau0H * TAUH_TO_TAUD
#         if field == 'HeI':    # 24.6
#             # assert(colden_HeI is not None)
#             tau = colden_HI * self.n_colden_H * sigmaHI(24.6)
#             tau += center2.nH2nHe * colden_HeI * self.n_colden_H * sigmaHeI(24.6)
#             return tau
#         if field == 'HeII':
#             assert(colden_HeI is not None)
#             assert(colden_HeII is not None)
#             tau = colden_HI * self.n_colden_H * sigmaHI(54.4)
#             tau += center2.nH2nHe * colden_HeI * self.n_colden_H * sigmaHeI(54.4)
#             tau += center2.nH2nHe * colden_HeII * self.n_colden_H * sigmaHeII(54.4)
#             return tau
#
#         # if field == 'HI' or field == 'H':       # 13.6
#         #     tau = hydrogenFrac*colDenHI*self.colDenUnit / mH * sigmaHI(13.6)
#         # if field == 'HeI':    # 24.6
#         #     tau = hydrogenFrac*colDenHI*self.colDenUnit / mH * sigmaHI(24.6)
#         #     tau += (1 - hydrogenFrac) * colDenHeI * self.colDenUnit / mHe * \
#         #           sigmaHeI(24.6)
#         # if field == 'HeII':
#         #     tau = hydrogenFrac*colDenHI*self.colDenUnit / mH * sigmaHI(54.4)
#         #     tau += (1 - hydrogenFrac) * colDenHeI * self.colDenUnit / mHe * \
#         #            sigmaHeI(54.4)
#         #     tau += (1 - hydrogenFrac) * colDenHeII * self.colDenUnit / mHe * \
#         #           sigmaHeII(54.4)
#
#     def get_Q_tot(self, outputid):
#         try:
#             masses = self.get_sink_masses(outputid)
#         except ramses.NoSinkParticle:
#             return 0
#         return np.sum(QVacca(masses)) * center.unit_lum
#
#     def get_Q_tot_all(self):
#         Q_tot_all = []
#         for i in range(1, self.get_last_output()+1):
#             Q_tot_all.append(self.get_Q_tot(i))
#         return np.array(Q_tot_all)


