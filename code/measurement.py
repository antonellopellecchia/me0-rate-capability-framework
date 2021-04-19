#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import ROOT as rt
import pandas as pd

from physlibs.root import root_style_cms

class GainCurve:
    def __init__(self, inputFile):
        df = pd.read_excel(inputFile, sheet_name='Data Summary', usecols='E,L', skiprows=29, nrows=15, header=None)
        dividerCurrent, gain = np.array(df[4], dtype='float'), np.array(df[11])
        self.gainPlot = rt.TGraph(len(gain), dividerCurrent, gain)

    def Get(self, dividerCurrent):
        return self.gainPlot.Eval(dividerCurrent)

class DataTaking:
    # parse single data-taking at fixed distance

    def __init__(self, inputFile, gainCurve, chamberDividerCurrent):
        self.qe, self.np = 1.6e-19, 404
        df = pd.read_excel(inputFile, sheet_name='XRay-off-substracted')

        channelXrayCurrent = 'XrayCurrent'
        channelXrayCurrentError = 'ERR XrayCurrent'
        channelAnodeCurrent = 'Ianode'
        channelAnodeCurrentError = 'ERR Ianode'

        self.measurement = {
            'xray': [np.array(df[channelXrayCurrent]), np.array(df[channelXrayCurrentError])],
            'anode': [abs(np.array(df[channelAnodeCurrent])), np.array(df[channelAnodeCurrentError])]
        }
        print(self.measurement['xray'])

        self.gainCurve = gainCurve
        self.chamberDividerCurrent = chamberDividerCurrent
        self.linearizationMethod = 'saturation'
        #self.chamberGain = gainCurve.Get(chamberDividerCurrent)

    def SetCollimator(self, radius):
        self.collimatorRadius = radius

    def SetLinearizeMethod(self, method):
        self.linearizationMethod = method

    def GetXrayCurrent(self):
        return self.measurement['xray'][0], self.measurement['xray'][1]

    def GetAnodeCurrent(self):
        return self.measurement['anode'][0], self.measurement['anode'][1]

    def GetAnodeCurrentLinearized(self):
        try: return self.anodeLinearized, self.errAnodeLinearized
        except AttributeError: pass
        
        xray, errXray = self.GetXrayCurrent()
        anodeCurrentPlot = self.GetAnodePlot()
        if self.linearizationMethod=='saturation':
            fit = rt.TF1('p', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 0, 200) # fit as (A+Bx)/(1+tau(A+Bx))
        elif self.linearizationMethod=='firstpoints':
            fit = rt.TF1('p', '[0]*x+[1]', 0, 30) # linear fit on first points
        else: raise ValueError('Unrecognized current linearization method')

        anodeCurrentPlot.Fit(fit, 'R')
        A, B = fit.GetParameter(0), fit.GetParameter(1)
        errA, errB = fit.GetParError(0), fit.GetParError(1)
        self.anodeLinearized, self.errAnodeLinearized = A*xray+B, np.sqrt(errB**2 + (A*xray)**2*((errA/A)**2+(errXray/xray)**2))
        return self.anodeLinearized, self.errAnodeLinearized

    def GetRate(self):
        try: return self.rate, self.errRate
        except AttributeError: pass

        xray, errXray = self.GetXrayCurrent()
        anodeCurrent, errAnodeCurrent = self.GetAnodeCurrentLinearized()
        nominalGain = self.gainCurve.Get(self.chamberDividerCurrent)
        self.rate, self.errRate = anodeCurrent/(self.qe*self.np*nominalGain), errAnodeCurrent/(self.qe*self.np*nominalGain)
        return self.rate, self.errRate

    def GetFlux(self):
        try: return self.flux, self.errFlux
        except AttributeError: pass

        try: spotArea = np.pi*self.collimatorRadius**2
        except AttributeError: spotArea = 100*100
        rate, errRate = self.GetRate()
        self.flux, self.errFlux = rate/spotArea, errRate/spotArea
        return self.flux, self.errFlux

    def GetEffectiveGain(self):
        try: return self.effectiveGain, self.errEffectiveGain
        except AttributeError: pass

        xray, errXray = self.GetXrayCurrent()
        anodeCurrent, errAnodeCurrent = self.GetAnodeCurrent()
        rate, errRate = self.GetRate()
        self.effectiveGain = anodeCurrent/(self.qe*self.np*rate)
        self.errEffectiveGain = self.effectiveGain * ((errAnodeCurrent/anodeCurrent)**2+(errRate/rate)**2)**0.5
        return self.effectiveGain, self.errEffectiveGain

    def GetAnodePlot(self):
        try: return self.anodeCurrentPlot
        except AttributeError: pass

        xray, errXray = self.GetXrayCurrent()
        anode, errAnode = self.GetAnodeCurrent()
        self.anodeCurrentPlot = rt.TGraphErrors(len(xray), xray, abs(anode), errXray, errAnode)
        self.anodeCurrentPlot.SetName('AnodePlot')
        self.anodeCurrentPlot.SetTitle(';X-ray current (#muA);Anode current (A)')
        return self.anodeCurrentPlot

    def GetAnodePlotLinearized(self):
        try: return self.anodeCurrentPlotLinearized
        except AttributeError: pass

        xray, errXray = self.GetXrayCurrent()
        linearized, errLinearized = self.GetAnodeCurrentLinearized()
        self.anodeCurrentPlotLinearized = rt.TGraphErrors(len(xray), xray, linearized, errXray, errLinearized)
        self.anodeCurrentPlotLinearized.SetName('AnodePlotLinearized')
        self.anodeCurrentPlotLinearized.SetTitle(';X-ray current (#muA);Anode current (A)')
        self.anodeCurrentPlotLinearized.SetMarkerStyle(5)
        return self.anodeCurrentPlotLinearized

    def GetRatePlot(self):
        try: return self.ratePlot
        except AttributeError: pass

        xray, errXray = self.GetXrayCurrent()
        rate, errRate = self.GetRate()
        self.ratePlot = rt.TGraphErrors(len(xray), xray, rate, errXray, errRate)
        self.ratePlot.SetName('RatePlot')
        self.ratePlot.SetTitle(';X-ray current (#muA);Rate (Hz)')
        return self.ratePlot

    def GetEffectiveGainPlot(self):
        try: return self.effectiveGainPlot
        except AttributeError: pass

        xray, errXray = self.GetXrayCurrent()
        effectiveGain, errEffectiveGain = self.GetEffectiveGain()
        self.effectiveGainPlot = rt.TGraphErrors(len(xray), xray, effectiveGain, errXray, errEffectiveGain)
        self.effectiveGainPlot.SetName('EffectiveGainPlot')
        self.effectiveGainPlot.SetTitle(';X-ray current (#muA);Effective gain')
        return self.effectiveGainPlot

    def GetRateCapabilityPlot(self):
        try: return self.rateCapabilityPlot
        except AttributeError: pass

        flux, errFlux = self.GetFlux()
        effectiveGain, errEffectiveGain = self.GetEffectiveGain()
        self.rateCapabilityPlot = rt.TGraphErrors(len(flux), flux, effectiveGain, errFlux, errEffectiveGain)
        self.rateCapabilityPlot.SetName('RateCapabilityPlot')
        self.rateCapabilityPlot.SetTitle(';Rate (Hz/mm^{2});Effective gain')
        return self.rateCapabilityPlot

    def SaveCurrents(self, epsPath, rootPath):
        rootFile = rt.TFile(rootPath, 'RECREATE')
        self.GetAnodePlot().Write()
        self.GetAnodePlotLinearized().Write()
        rootFile.Close()

        currentCanvas = rt.TCanvas('CurrentCanvas', '', 800, 600)
        currentMultiGraph = rt.TMultiGraph()
        currentMultiGraph.SetTitle(';X-ray current (#muA);Anode current (A)')
        currentMultiGraph.Add(self.GetAnodePlot(), 'p')
        currentMultiGraph.Add(self.GetAnodePlotLinearized(), 'p')
        currentMultiGraph.Draw('a')
        currentCanvas.SaveAs(epsPath)

    def SaveRate(self, epsPath, rootPath):
        ratePlot = self.GetRatePlot()
        ratePlot.SaveAs(rootPath)
        rateCanvas = rt.TCanvas('RateCanvas', '', 800, 600)
        ratePlot.Draw('ap')
        rateCanvas.SaveAs(epsPath)

    def SaveEffectiveGain(self, epsPath, rootPath):
        effectiveGainPlot = self.GetEffectiveGainPlot()
        effectiveGainPlot.SaveAs(rootPath)
        gainCanvas = rt.TCanvas('EffectiveGainCanvas', '', 800, 600)
        effectiveGainPlot.Draw('ap')
        gainCanvas.SaveAs(epsPath)

    def SaveRateCapability(self, epsPath, rootPath):
        rateCapabilityPlot = self.GetRateCapabilityPlot()
        rateCapabilityPlot.SaveAs(rootPath)
        rateCapabilityCanvas = rt.TCanvas('RateCapabilityCanvas', '', 800, 600)
        rateCapabilityCanvas.SetLogx()
        rateCapabilityCanvas.SetGrid()
        rateCapabilityPlot.Draw('apc')
        rateCapabilityCanvas.SaveAs(epsPath)