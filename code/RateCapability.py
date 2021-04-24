#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import ROOT as rt
import pandas as pd

import root_style_cms

import measurement

def main():
    ap = argparse.ArgumentParser(add_help=True)
    #ap.add_argument('--input', nargs='+')
    #ap.add_argument('--output')
    #ap.add_argument('--qc5')
    ap.add_argument('--dividerCurrent', type=float)
    ap.add_argument('--verbose', action='store_true')
    options = ap.parse_args(sys.argv[1:])

    measurementFile = os.environ['RATE_CAPABILITY_DATA']+'/RateCapability.root'
    gainFile = os.environ['RATE_CAPABILITY_DATA']+'/EffectiveGain.xlsx'
    resultsDirectory = os.environ['RATE_CAPABILITY_RESULTS']+'/RateCapability'
    outputDirectory = os.environ['RATE_CAPABILITY_OUTPUT']+'/RateCapability'

    for d in [resultsDirectory, outputDirectory]:
        try: os.makedirs(d)
        except FileExistsError: pass

    chamberGainCurve = measurement.GainCurve(gainFile)
    chamberGain = chamberGainCurve.Get(options.dividerCurrent)
    print('%d chamber gain'%(chamberGain))

    meas = measurement.Measurement.FromFile(measurementFile, chamberGainCurve, options.dividerCurrent)
    
    for dataTaking in meas: # fit current plots separately for each xray-to-chamber distance
        path = f'{outputDirectory}/{dataTaking.name}'
        dataTaking.SaveCurrents(f'{path}_AnodeCurrent.eps', f'{path}_AnodeCurrent.root')
        dataTaking.SaveRate(f'{path}_Rate.eps', f'{path}_Rate.root')
        dataTaking.SaveEffectiveGain(f'{path}_EffectiveGain.eps', f'{path}_EffectiveGain.root')
        dataTaking.SaveRateCapability(f'{path}_RateCapability.eps', f'{path}_RateCapability.root')
    meas.SaveRateCapability(f'{resultsDirectory}/RateCapability.eps', f'{resultsDirectory}/RateCapability.root')

    return

    rateCapabilityGraph = rt.TMultiGraph()
    rateCapabilityGraphs = list()
    legend = rt.TLegend(0.2, 0.2, 0.8, 0.35)
    colors = [rt.kBlack, rt.kRed, rt.kGreen, rt.kBlue]
    for iMeasurement,measurementFile in enumerate(options.input):
        title = measurementFile.split('/')[-1][:-5]
        try: os.makedirs('%s/%s/'%(options.output, title))
        except FileExistsError: pass

        meas = measurement.DataTaking(measurementFile, chamberGainCurve, options.dividerCurrent)
        meas.SetCollimator(1)
        #meas.SetLinearizeMethod('firstpoints')
        meas.SaveCurrents('%s/%s/AnodeCurrent.eps'%(options.output, title), '%s/%s/AnodeCurrent.root'%(options.output, title))
        meas.SaveEffectiveGain('%s/%s/EffectiveGain.eps'%(options.output, title), '%s/%s/EffectiveGain.root'%(options.output, title))
        meas.SaveRateCapability('%s/%s/RateCapability.eps'%(options.output, title), '%s/%s/RateCapability.root'%(options.output, title))
        
        rateCapabilityGraphs.append(meas.GetRateCapabilityPlot())
        g = rateCapabilityGraphs[-1]
        g.SetMarkerColor(colors[iMeasurement])
        rateCapabilityGraph.Add(g, 'p')
        legend.AddEntry(g, title, 'p')

    rateCapabilityCanvas = rt.TCanvas('RateCapabilityCanvas', '', 800, 600)
    rateCapabilityGraph.SetTitle(';Rate (Hz/mm^{2});Effective gain')
    rateCapabilityGraph.Draw('a')
    rateCapabilityCanvas.SetGrid()
    rateCapabilityCanvas.SetLogx()
    legend.Draw()
    rateCapabilityCanvas.SaveAs('%s/RateCapability.eps'%(options.output))

if __name__=='__main__': main()
