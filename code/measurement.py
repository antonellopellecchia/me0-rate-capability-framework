#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import ROOT as rt
import pandas as pd

import root_style_cms

qe, primaries, errPrimaries = 1.6e-19, 418, 9

class GainCurve:
    def __init__(self, inputFile):
        df = pd.read_excel(inputFile, sheet_name='Data Summary', usecols='E,L', skiprows=29, nrows=15, header=None)
        dividerCurrent, gain = np.array(df[4], dtype='float'), np.array(df[11])
        self.gainPlot = rt.TGraph(len(gain), dividerCurrent, gain)

    def Get(self, dividerCurrent):
        return self.gainPlot.Eval(dividerCurrent)

class DataTaking:
    # parse single data-taking at fixed distance

    def __init__(self, name, dataTakingDf, gainCurve, chamberDividerCurrent, linearizationMethod='saturation2part'):
        # set column names as in root and excel files:
        channelXrayCurrent = 'XrayCurrent'
        channelXrayCurrentError = 'ERRXrayCurrent'
        channelAnodeCurrent = 'Ianode'
        channelAnodeCurrentError = 'ERRIanode'
        
        # sort dataframe by xray currents in ascending order:
        dataTakingDf.sort_values(channelXrayCurrent, inplace=True)

        xrayCurrent, errXrayCurrent = np.array(dataTakingDf[channelXrayCurrent]), np.array(dataTakingDf[channelXrayCurrentError])
        anodeCurrent, errAnodeCurrent = abs(np.array(dataTakingDf[channelAnodeCurrent])), np.array(dataTakingDf[channelAnodeCurrentError])

        self.name = name
        self.measurement = {
            'xray': [xrayCurrent, errXrayCurrent],
            'anode': [anodeCurrent, errAnodeCurrent]
        }
        self.gainCurve = gainCurve
        self.chamberDividerCurrent = chamberDividerCurrent
        self.linearizationMethod = linearizationMethod # saturation, firstpoints, saturation2part
        #self.chamberGain = gainCurve.Get(chamberDividerCurrent)

    def FromExcelFile(name, inputFile, gainCurve, chamberDividerCurrent):
        rootFile = rt.TFile(inputFile, 'READ')
        df = pd.read_excel(inputFile, sheet_name='XRay-off-substracted')

        channelXrayCurrent = 'XrayCurrent'
        channelXrayCurrentError = 'ERR XrayCurrent'
        channelAnodeCurrent = 'Ianode'
        channelAnodeCurrentError = 'ERR Ianode'
        xrayCurrent, errXrayCurrent = np.array(df[channelXrayCurrent]), np.array(df[channelXrayCurrentError])
        anodeCurrent, errAnodeCurrent = abs(np.array(df[channelAnodeCurrent])), np.array(df[channelAnodeCurrentError])

        return DataTaking(name, xrayCurrent, errXrayCurrent, anodeCurrent, errAnodeCurrent, gainCurve, chamberDividerCurrent)

    def SetCollimator(self, radius):
        self.collimatorRadius = radius

    def SetLinearizeMethod(self, method):
        self.linearizationMethod = method

    @property
    def xrayCurrent(self):
        return self.measurement['xray'][0], self.measurement['xray'][1]

    @property
    def anodeCurrent(self):
        return self.measurement['anode'][0], self.measurement['anode'][1]

    @property
    def anodeCurrentLinearized(self):
        try: return self._anodeCurrentLinearized, self._errAnodeCurrentLinearized
        except AttributeError: pass
        
        xray, errXray = self.xrayCurrent
        anodeCurrentPlot = self.anodePlot
        if self.linearizationMethod=='saturation':
            fit = rt.TF1('p', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 0, 200) # fit as (A+Bx)/(1+tau(A+Bx))
        elif self.linearizationMethod=='firstpoints':
            fit = rt.TF1('p', '[0]*x+[1]', 0, 30) # linear fit on first points
        elif self.linearizationMethod=='saturation2part':
            # fit as sum of two saturating functions:
            fit1 = rt.TF1('p1', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 0.1, 100)
            fit2 = rt.TF1('p2', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 101, 200)
        else: raise ValueError('Unrecognized current linearization method')

        if self.linearizationMethod in ['saturation', 'firstpoints']:
            anodeCurrentPlot.Fit(fit, 'R')
            A, B = fit.GetParameter(0), fit.GetParameter(1)
            errA, errB = fit.GetParError(0), fit.GetParError(1)
            self._anodeCurrentLinearized, self._errAnodeCurrentLinearized = A*xray+B, np.sqrt(errB**2 + (A*xray)**2*((errA/A)**2+(errXray/xray)**2))
        else:
            anodeCurrentPlot.Fit(fit1, 'R')
            anodeCurrentPlot.Fit(fit2, 'R')
            A1, B1, A2, B2 = fit1.GetParameter(0), fit1.GetParameter(1), fit2.GetParameter(0), fit2.GetParameter(1)
            errA1, errB1, errA2, errB2 = fit1.GetParError(0), fit1.GetParError(1), fit2.GetParError(0), fit2.GetParError(1)
            xray1, errXray1 = xray[xray<=100], errXray[xray<=100]
            xray2, errXray2 = xray[xray>100], errXray[xray>100]
            linearized1, errLinearized1 = A1*xray1+B1, np.sqrt(errB1**2 + (A1*xray1)**2*((errA1/A1)**2+(errXray1/xray1)**2))
            linearized2, errLinearized2 = A2*xray2+B2, np.sqrt(errB2**2 + (A2*xray2)**2*((errA2/A2)**2+(errXray2/xray2)**2))
            self._anodeCurrentLinearized = np.concatenate([linearized1, linearized2])
            self._errAnodeCurrentLinearized = np.concatenate([errLinearized1, errLinearized2])
        return self._anodeCurrentLinearized, self._errAnodeCurrentLinearized

    @property
    def rate(self):
        try: return self._rate, self._errRate
        except AttributeError: pass

        anodeCurrent, errAnodeCurrent = self.anodeCurrentLinearized
        nominalGain = self.gainCurve.Get(self.chamberDividerCurrent)
        self._rate, self._errRate = anodeCurrent/(qe*primaries*nominalGain), errAnodeCurrent/(qe*primaries*nominalGain)
        return self._rate, self._errRate

    @property
    def flux(self):
        try: return self._flux, self._errFlux
        except AttributeError: pass

        try: spotArea = np.pi*self.collimatorRadius**2
        except AttributeError: spotArea = 100*100
        rate, errRate = self.rate
        self._flux, self._errFlux = rate/spotArea, errRate/spotArea
        return self._flux, self._errFlux

    @property
    def effectiveGain(self):
        try: return self._effectiveGain, self._errEffectiveGain
        except AttributeError: pass

        xray, errXray = self.xrayCurrent
        anodeCurrent, errAnodeCurrent = self.anodeCurrent
        rate, errRate = self.rate
        self._effectiveGain = anodeCurrent/(qe*primaries*rate)
        self._errEffectiveGain = self._effectiveGain * ((errAnodeCurrent/anodeCurrent)**2+(errRate/rate)**2)**0.5
        return self._effectiveGain, self._errEffectiveGain

    @property
    def anodePlot(self):
        try: return self._anodePlot
        except AttributeError: pass

        xray, errXray = self.xrayCurrent
        anode, errAnode = self.anodeCurrent
        self._anodePlot = rt.TGraphErrors(len(xray), xray, abs(anode), errXray, errAnode)
        self._anodePlot.SetName('AnodePlot')
        self._anodePlot.SetTitle(';X-ray current (#muA);Anode current (A)')
        return self._anodePlot

    @property
    def anodePlotLinearized(self):
        try: return self._anodePlotLinearized
        except AttributeError: pass

        xray, errXray = self.xrayCurrent
        linearized, errLinearized = self.anodeCurrentLinearized
        self._anodePlotLinearized = rt.TGraphErrors(len(xray), xray, linearized, errXray, errLinearized)
        self._anodePlotLinearized.SetName('AnodePlotLinearized')
        self._anodePlotLinearized.SetTitle(';X-ray current (#muA);Anode current (A)')
        self._anodePlotLinearized.SetMarkerStyle(5)
        return self._anodePlotLinearized

    @property
    def ratePlot(self):
        try: return self._ratePlot
        except AttributeError: pass

        xray, errXray = self.xrayCurrent
        rate, errRate = self.rate
        self._ratePlot = rt.TGraphErrors(len(xray), xray, rate, errXray, errRate)
        self._ratePlot.SetName('RatePlot')
        self._ratePlot.SetTitle(';X-ray current (#muA);Rate (Hz)')
        return self._ratePlot

    @property
    def effectiveGainPlot(self):
        try: return self._effectiveGainPlot
        except AttributeError: pass

        xray, errXray = self.xrayCurrent
        effectiveGain, errEffectiveGain = self.effectiveGain
        self._effectiveGainPlot = rt.TGraphErrors(len(xray), xray, effectiveGain, errXray, errEffectiveGain)
        self._effectiveGainPlot.SetName('EffectiveGainPlot')
        self._effectiveGainPlot.SetTitle(';X-ray current (#muA);Effective gain')
        return self._effectiveGainPlot

    @property
    def rateCapabilityPlot(self):
        try: return self._rateCapabilityPlot
        except AttributeError: pass

        flux, errFlux = self.flux
        effectiveGain, errEffectiveGain = self.effectiveGain
        self._rateCapabilityPlot = rt.TGraphErrors(len(flux), flux, effectiveGain, errFlux, errEffectiveGain)
        self._rateCapabilityPlot.SetName('RateCapabilityPlot')
        self._rateCapabilityPlot.SetTitle(';Rate (Hz/mm^{2});Effective gain')
        return self._rateCapabilityPlot

    def SaveCurrents(self, epsPath, rootPath):
        rootFile = rt.TFile(rootPath, 'RECREATE')
        self.anodePlot.Write()
        self.anodePlotLinearized.Write()
        rootFile.Close()

        currentCanvas = rt.TCanvas('CurrentCanvas', '', 800, 600)
        currentMultiGraph = rt.TMultiGraph()
        currentMultiGraph.SetTitle(';X-ray current (#muA);Anode current (A)')
        currentMultiGraph.Add(self.anodePlot, 'p')
        currentMultiGraph.Add(self.anodePlotLinearized, 'p')
        currentMultiGraph.Draw('a')
        currentCanvas.SaveAs(epsPath)

    def SaveRate(self, epsPath, rootPath):
        ratePlot = self.ratePlot
        ratePlot.SaveAs(rootPath)
        rateCanvas = rt.TCanvas('RateCanvas', '', 800, 600)
        ratePlot.Draw('ap')
        rateCanvas.SaveAs(epsPath)

    def SaveEffectiveGain(self, epsPath, rootPath):
        effectiveGainPlot = self.effectiveGainPlot
        effectiveGainPlot.SaveAs(rootPath)
        gainCanvas = rt.TCanvas('EffectiveGainCanvas', '', 800, 600)
        effectiveGainPlot.Draw('ap')
        gainCanvas.SaveAs(epsPath)

    def SaveRateCapability(self, epsPath, rootPath):
        rateCapabilityPlot = self.rateCapabilityPlot
        rateCapabilityPlot.SaveAs(rootPath)
        rateCapabilityCanvas = rt.TCanvas('RateCapabilityCanvas', '', 800, 600)
        rateCapabilityCanvas.SetLogx()
        rateCapabilityCanvas.SetGrid()
        rateCapabilityPlot.Draw('ap')
        rateCapabilityCanvas.SaveAs(epsPath)


class Measurement:
    # groups one or more DataTakings

    def __init__(self, dataTakingList):
        self.dataTakingList = dataTakingList

    def FromFile(file, gainCurve, dividerCurrent):
        rootFile = rt.TFile(file, 'READ')
        dataTakingList = list()
        for treeKey in rootFile.GetListOfKeys():
            treeName = treeKey.GetName()
            tree = rootFile.Get(treeName)
            data, cols = tree.AsMatrix(return_labels=True)
            measurementDf = pd.DataFrame(data=data, columns=cols)
            dataTakingList.append(DataTaking(treeName, measurementDf, gainCurve, dividerCurrent))
        return Measurement(dataTakingList)

    def __iter__(self): return self.dataTakingList.__iter__()

    @property
    def rateCapabilityPlot(self):
        try: return self._rateCapabilityPlot
        except AttributeError: pass

        self._rateCapabilityPlot = rt.TMultiGraph()
        self._rateCapabilityPlot.SetName('RateCapabilityPlot')
        self._rateCapabilityPlot.SetTitle('Rate capability;Rate (Hz/mm^{2});Effective gain')
        for dataTaking in self: self._rateCapabilityPlot.Add(dataTaking.rateCapabilityPlot, 'p')
        return self._rateCapabilityPlot

    def SaveRateCapability(self, epsPath, rootPath):
        rateCapabilityPlot = self.rateCapabilityPlot
        rateCapabilityPlot.SaveAs(rootPath)
        rateCapabilityCanvas = rt.TCanvas('RateCapabilityCanvas', '', 800, 600)
        rateCapabilityCanvas.SetLogx()
        rateCapabilityCanvas.SetGrid()
        rateCapabilityPlot.Draw('a')
        rateCapabilityCanvas.SaveAs(epsPath)