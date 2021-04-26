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
        self.gainPlotFit = rt.TF1('e', 'expo(0)')
        self.gainPlot.Fit(self.gainPlotFit)
        self.A, self.B = self.gainPlotFit.GetParameter(0), self.gainPlotFit.GetParameter(1)
        self.errA, self.errB = self.gainPlotFit.GetParError(0), self.gainPlotFit.GetParError(1)
        '''c = rt.TCanvas('c', '', 800, 600)
        self.gainPlot.Draw('ACP')
        c.SaveAs('gain.eps')'''

    def GetGain(self, dividerCurrent):
        return np.exp(self.A+dividerCurrent*self.B)
        #return self.gainPlotFit.Eval(dividerCurrent)
        #return self.gainPlot.Eval(dividerCurrent)

    def GetError(self, dividerCurrent):
        return self.GetGain(dividerCurrent)*np.sqrt(self.errA**2+dividerCurrent**2*self.errB**2)
        #return self.gainPlotFit.Eval(dividerCurrent)
        #return self.gainPlot.Eval(dividerCurrent)

class DataTaking:
    # parse single data-taking at fixed distance

    def __init__(self, name, dataTakingDf, gainCurve, chamberDividerCurrent, linearizationMethod='piecewiseSaturation'):
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
        '''self.measurement = {
            'xray': [xrayCurrent, errXrayCurrent],
            'anode': [anodeCurrent, errAnodeCurrent]
        }'''
        self._xrayCurrent, self._errXrayCurrent = xrayCurrent, errXrayCurrent
        self._anodeCurrent, self._errAnodeCurrent = anodeCurrent, errAnodeCurrent
        self.gainCurve = gainCurve
        self.chamberDividerCurrent = chamberDividerCurrent
        self.linearizationMethod = linearizationMethod # saturation, firstpoints, saturation2part
        #self.chamberGain = gainCurve.GetGain(chamberDividerCurrent)

    def FromExcelFile(name, inputFile, gainCurve, chamberDividerCurrent):
        rootFile = rt.TFile(inputFile, 'READ')
        dataTakingDf = pd.read_excel(inputFile, sheet_name='XRay-off-substracted')
        return DataTaking(name, dataTakingDf, gainCurve, chamberDividerCurrent)

        '''channelXrayCurrent = 'XrayCurrent'
        channelXrayCurrentError = 'ERR XrayCurrent'
        channelAnodeCurrent = 'Ianode'
        channelAnodeCurrentError = 'ERR Ianode'
        xrayCurrent, errXrayCurrent = np.array(df[channelXrayCurrent]), np.array(df[channelXrayCurrentError])
        anodeCurrent, errAnodeCurrent = abs(np.array(df[channelAnodeCurrent])), np.array(df[channelAnodeCurrentError])'''


    def SetCollimator(self, radius):
        self.collimatorRadius = radius

    def SetLinearizeMethod(self, method):
        self.linearizationMethod = method

    @property
    def nominalGain(self):
        try: return self._nominalGain
        except AttributeError: pass
        self._nominalGain = self.gainCurve.GetGain(self.chamberDividerCurrent)
        return self._nominalGain

    @property
    def xrayCurrent(self):
        #return self.measurement['xray'][0], self.measurement['xray'][1]
        return self._xrayCurrent, self._errXrayCurrent

    @property
    def anodeCurrent(self):
        return self._anodeCurrent, self._errAnodeCurrent

    @property
    def anodeCurrentLinearized(self):
        try: return self._anodeCurrentLinearized, self._errAnodeCurrentLinearized
        except AttributeError: pass
        
        xray, errXray = self.xrayCurrent
        anodeCurrent, errAnodeCurrent = self.anodeCurrent
        anodeCurrentPlot = self.anodePlot
        if self.linearizationMethod=='saturation':
            fit = rt.TF1('p', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 0, 200) # fit as (A+Bx)/(1+tau(A+Bx))
            fit.FixParameter(1, 0)
        elif self.linearizationMethod=='firstpoints':
            fit = rt.TF1('p', '[0]*x+[1]', 0, 30) # linear fit on first points
        elif self.linearizationMethod=='saturation2part':
            xrayCurrentMax = 99.2
            # fit as two saturating functions separately:
            separationXrayCurrent = max(xray[xray<=100]) # divide x-ray currents in two sets at about 100
            #separationAnodeCurrent = anodeCurrent[xray==separationXrayCurrent][0]
            fit1 = rt.TF1('p1', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 0.1, separationXrayCurrent)
            fit1.FixParameter(1, 0)
            fit2 = rt.TF1('p2', '([0]*(x-[3])+[1])/(1+([0]*(x-[3])+[1])*[2])+[4]', separationXrayCurrent+1, 200)
            fit2.SetParameter(1, 0)
        elif self.linearizationMethod=='piecewiseSaturation':
            # fit as two saturating functions piecewise:

            f1 = '([0]*x+[1])/(1+([0]*x+[1])*[2])*(x<[3])'
            '''f2 = '+([0]*[3]+[1])/(1+([0]*[3]+[1])*[2])*(x>[3])'
            f3 = '+([4]*(x-[3])+[5])/(1+([4]*(x-[3])+[5])*[6])*(x>[3])'''
            f2 = '+([0]*[3]+[1]+[4]*(x-[3])+[5])/(1+([0]*[3]+[1])*[2]+([4]*(x-[3])+[5])*[6])*(x>[3])'
            fit = rt.TF1('p1', f1+f2, 0, 200)
            fit.SetParNames('A1', 'B1', 't1', 'x0', 'A2', 'B2', 't2')
            
            # get initial guess for parameters A and B from fit on restricted range:
            print('Guessing initial parameters...')
            fitPartial1 = rt.TF1('p1', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 0, 100)
            fitPartial2 = rt.TF1('p2', '([0]*x+[1])/(1+([0]*x+[1])*[2])', 100, 200)
            anodeCurrentPlot.Fit(fitPartial1, 'R')
            anodeCurrentPlot.Fit(fitPartial2, 'R')
            A1, B1, t1 = fitPartial1.GetParameter(0), fitPartial1.GetParameter(1), fitPartial1.GetParameter(2)
            A2, B2, t2 = fitPartial2.GetParameter(0), fitPartial2.GetParameter(1), fitPartial2.GetParameter(2)
            fit.SetParameters(A1, B1, t1, 99, A2, B2, t2)
            #fit.SetParameters(4e-9, 1e-8, 8e4, 99, 1e-9, 1e-8, 1e4)
            fit.FixParameter(3, 99.2)
            #fit.FixParameter(1, 0)
        else: raise ValueError('Unrecognized current linearization method')

        if self.linearizationMethod in ['saturation', 'firstpoints']:
            anodeCurrentPlot.Fit(fit, 'R')
            A, B = fit.GetParameter(0), fit.GetParameter(1)
            errA, errB = fit.GetParError(0), fit.GetParError(1)
            self._anodeCurrentLinearized, self._errAnodeCurrentLinearized = A*xray+B, np.sqrt(errB**2 + (A*xray)**2*((errA/A)**2+(errXray/xray)**2))
        elif self.linearizationMethod=='saturation2part':
            anodeCurrentPlot.Fit(fit1, 'R')
            fit2.FixParameter(3, xrayCurrentMax)
            fit2.FixParameter(4, fit1.Eval(xrayCurrentMax))
            anodeCurrentPlot.Fit(fit2, 'R')
            A1, B1, A2, B2 = fit1.GetParameter(0), fit1.GetParameter(1), fit2.GetParameter(0), fit2.GetParameter(1)
            errA1, errB1, errA2, errB2 = fit1.GetParError(0), fit1.GetParError(1), fit2.GetParError(0), fit2.GetParError(1)
            xray1, errXray1 = xray[xray<=separationXrayCurrent], errXray[xray<=separationXrayCurrent]
            xray2, errXray2 = xray[xray>separationXrayCurrent], errXray[xray>separationXrayCurrent]
            linearized1, errLinearized1 = A1*xray1+B1, np.sqrt(errB1**2 + (A1*xray1)**2*((errA1/A1)**2+(errXray1/xray1)**2))
            #shift2, offset2 = fit2.GetParameter(3), A1*xrayCurrentMax+B1
            shift2, offset2 = fit2.GetParameter(3), fit2.GetParameter(4)
            linearized2, errLinearized2 = offset2 + A2*(xray2-shift2)+B2, np.sqrt(errB2**2 + (A2*xray2)**2*((errA2/A2)**2+(errXray2/xray2)**2))
            self._anodeCurrentLinearized = np.concatenate([linearized1, linearized2])
            self._errAnodeCurrentLinearized = np.concatenate([errLinearized1, errLinearized2])
        elif self.linearizationMethod=='piecewiseSaturation':
            print('Fitting with saturating function...')
            anodeCurrentPlot.Fit(fit, 'R')
            A1, B1, x0, A2, B2,  = fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(3), fit.GetParameter(4), fit.GetParameter(5)
            errA1, errB1, errX0, errA2, errB2 = fit.GetParError(0), fit.GetParError(1), fit.GetParError(3), fit.GetParError(4), fit.GetParError(5)
            self._anodeCurrentLinearized = (A1*xray+B1)*(xray<=x0) + (A1*x0+B1 + A2*(xray-x0)+B2)*(xray>x0)
            self._errAnodeCurrentLinearized = np.sqrt(errA1**2*errXray**2+errB1**2)*(xray<=x0) + np.sqrt(errA1**2*errX0**2+B1**2 + errA2**2*(errXray**2+errX0**2)+errB2**2)*(xray>x0)
            #self._errAnodeCurrentLinearized = np.sqrt(errB1**2 + (A1*xray)**2*((errA1/A1)**2+(errXray/xray)**2))
        return self._anodeCurrentLinearized, self._errAnodeCurrentLinearized

    @property
    def rate(self): # hit rate on chamber in Hz
        try: return self._rate, self._errRate
        except AttributeError: pass

        anodeCurrent, errAnodeCurrent = self.anodeCurrentLinearized
        nominalGain, errNominalGain = self.gainCurve.GetGain(self.chamberDividerCurrent), self.gainCurve.GetError(self.chamberDividerCurrent)
        self._rate = anodeCurrent/(qe*primaries*nominalGain)
        self._errRate = self._rate * np.sqrt((errAnodeCurrent/anodeCurrent)**2 + (errNominalGain/nominalGain)**2)
        return self._rate, self._errRate

    @property
    def flux(self): # flux on chamber in Hz/cm2
        try: return self._flux, self._errFlux
        except AttributeError: pass

        try: spotArea = np.pi*self.collimatorRadius**2
        except AttributeError: spotArea = 10*10
        rate, errRate = self.rate
        self._flux, self._errFlux = rate/spotArea, errRate/spotArea
        return self._flux, self._errFlux

    @property
    def effectiveGain(self):
        try: return self._effectiveGain, self._errEffectiveGain
        except AttributeError: pass

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
        self._ratePlot = rt.TGraphErrors(len(xray), xray, rate/1e3, errXray, errRate/1e3)
        self._ratePlot.SetName('RatePlot')
        self._ratePlot.SetTitle(';X-ray current (#muA);Rate (kHz)')
        return self._ratePlot

    @property
    def effectiveGainPlot(self):
        try: return self._effectiveGainPlot
        except AttributeError: pass

        xray, errXray = self.xrayCurrent
        effectiveGain, errEffectiveGain = self.effectiveGain
        self._effectiveGainPlot = rt.TGraphErrors(len(xray), xray, effectiveGain/1e3, errXray, errEffectiveGain/1e3)
        self._effectiveGainPlot.SetName('EffectiveGainPlot')
        self._effectiveGainPlot.SetTitle(';X-ray current (#muA);Effective gain (#times 10^{4})')
        return self._effectiveGainPlot

    @property
    def rateCapabilityPlot(self):
        try: return self._rateCapabilityPlot
        except AttributeError: pass

        flux, errFlux = self.flux
        effectiveGain, errEffectiveGain = self.effectiveGain
        self._rateCapabilityPlot = rt.TGraphErrors(len(flux), flux/1e3, effectiveGain/1e3, errFlux/1e3, errEffectiveGain/1e3)
        self._rateCapabilityPlot.SetName('RateCapabilityPlot')
        self._rateCapabilityPlot.SetTitle(';Rate (kHz/cm^{2});Effective gain (#times 10^{4})')
        return self._rateCapabilityPlot

    def SaveCurrents(self, epsPath, rootPath):
        rootFile = rt.TFile(rootPath, 'RECREATE')
        self.anodePlot.Write()
        self.anodePlotLinearized.Write()
        rootFile.Close()

        currentCanvas = rt.TCanvas('CurrentCanvas', '', 600, 600)
        currentMultiGraph = rt.TMultiGraph()
        currentMultiGraph.SetTitle(';X-ray current (#muA);Anode current (A)')
        currentMultiGraph.Add(self.anodePlot, 'p')
        currentMultiGraph.Add(self.anodePlotLinearized, 'p')
        currentMultiGraph.Draw('a')
        currentCanvas.SaveAs(epsPath)

    def SaveRate(self, epsPath, rootPath):
        ratePlot = self.ratePlot
        ratePlot.SaveAs(rootPath)
        rateCanvas = rt.TCanvas('RateCanvas', '', 600, 600)
        ratePlot.Draw('ap')
        rateCanvas.SaveAs(epsPath)

    def SaveEffectiveGain(self, epsPath, rootPath):
        effectiveGainPlot = self.effectiveGainPlot
        effectiveGainPlot.SaveAs(rootPath)
        gainCanvas = rt.TCanvas('EffectiveGainCanvas', '', 600, 600)
        effectiveGainPlot.Draw('ap')
        gainCanvas.SaveAs(epsPath)

    def SaveRateCapability(self, epsPath, rootPath):
        rateCapabilityPlot = self.rateCapabilityPlot
        rateCapabilityPlot.SaveAs(rootPath)
        rateCapabilityCanvas = rt.TCanvas('RateCapabilityCanvas', '', 600, 600)
        rateCapabilityCanvas.SetLogx()
        rateCapabilityPlot.Draw('ap')
        rateCapabilityCanvas.SaveAs(epsPath)


class Measurement:
    # groups one or more DataTakings

    def __init__(self, dataTakingList):
        self.dataTakingList = dataTakingList

    def FromFile(file, gainCurve, dividerCurrent, linearizationMethod='saturation2part'):
        rootFile = rt.TFile(file, 'READ')
        dataTakingList = list()
        for treeKey in rootFile.GetListOfKeys():
            treeName = treeKey.GetName()
            tree = rootFile.Get(treeName)
            data, cols = tree.AsMatrix(return_labels=True)
            measurementDf = pd.DataFrame(data=data, columns=cols)
            dataTakingList.append(DataTaking(treeName, measurementDf, gainCurve, dividerCurrent, linearizationMethod))
        return Measurement(dataTakingList)

    def __iter__(self): return self.dataTakingList.__iter__()

    @property
    def rateCapabilityPlot(self):
        try: return self._rateCapabilityPlot
        except AttributeError: pass

        self._rateCapabilityPlot = rt.TMultiGraph()
        self._rateCapabilityPlot.SetName('RateCapabilityPlot')
        self._rateCapabilityPlot.SetTitle('Rate capability;Rate (kHz/cm^{2});Effective gain (#times 10^{4})')
        for dataTaking in self: self._rateCapabilityPlot.Add(dataTaking.rateCapabilityPlot, 'p')
        return self._rateCapabilityPlot

    def SaveRateCapability(self, epsPath, rootPath):
        rateCapabilityPlot = self.rateCapabilityPlot
        rateCapabilityPlot.SaveAs(rootPath)
        rateCapabilityCanvas = rt.TCanvas('RateCapabilityCanvas', '', 600, 600)
        rateCapabilityCanvas.SetLogx()
        rateCapabilityPlot.Draw('a')
        rateCapabilityCanvas.SaveAs(epsPath)