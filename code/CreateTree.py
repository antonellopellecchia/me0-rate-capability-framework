#!/usr/bin/python3

import os, sys
import argparse
import re

import numpy as np
import ROOT as rt
import pandas as pd

import root_style_cms

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input', nargs='+')
    ap.add_argument('--labels', nargs='+', type=str)
    ap.add_argument('--output')
    options = ap.parse_args(sys.argv[1:])

    outDirectory = '/'.join(options.output.split('/')[:-1])
    try: os.makedirs(outDirectory)
    except FileExistsError: pass

    outFile = rt.TFile(options.output, 'RECREATE')

    # create tree from each input file
    for label,inputFile in zip(options.labels,options.input):
        treeDf = pd.read_excel(inputFile, sheet_name='XRay-off-substracted')

        treeName = f'Tree{label}'
        tree = rt.TTree(treeName, 'Current measurement')
        columnNames = treeDf.columns
        branchNames = [ col.replace(' ', '').replace('-', '') for col in columnNames ]
        branchVariables = dict()
        for branch in branchNames:
            branchVariables[branch] = np.zeros(1)
            tree.Branch(branch, branchVariables[branch], f'{branch}/D')

        print(columnNames, branchNames)
        for irow,row in treeDf.iterrows():
            for branch,col in zip(branchNames,columnNames):
                branchVariables[branch][0] = float(row[col])
            tree.Fill()

        tree.Print()
        outFile.Write(treeName)

    outFile.Write()
    outFile.Close()

if __name__=='__main__': main()