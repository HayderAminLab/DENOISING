# -*- coding: utf-8 -*-
"""
Created on Dec 12 2022
@author:  BIONICS_LAB
@company: DZNE
"""

import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
"""need change these parameters based on different datasets"""
SpikesMax = 15               # (Hz) - maximum acceptable Spikes Rate(Spikes events/min)
SpikesMin = 0.1                # (Hz) - minumum acceptable Spikes Rate(Spikes events/min)
mbrMin = 0.1              # burst/min - minumum acceptable MBR
mbrMax = 10              # burst/min - minumum acceptable MBR
Synchrony_Threshold = 0.05 #Those averageEventRateHists less than Synchrony_Threshold of all the events count in hist are regarded as random Events
Firing_Electrodes_Count_Threshold = 3500  # the maximum acceptable active channel numbers at the same time
class Spikes_DENOISING_Function:
    def __init__(self,srcfilepath):
        self.srcfilepath = srcfilepath  # main path
        self.SpikesMax = SpikesMax
        self.SpikesMin = SpikesMin

    def get_filename(self, filepath, filetype):
        """
        Getting the reference file name and path based on the file type and reference file path.
        Parameters
        ----------
        filepath: string
            the path of the file
            filetype: string
            the type of the file
            Returns
            -------
            filename: string
            the name of the file
        """
        filename = []
        for root, dirs, files in os.walk(filepath):
            for i in files:
                if filetype in i:
                    filename.append(i)
        return filename


    def histogram_filter(self,SpikesTime=None,SpikesID=None):
        """
        The function is used to filter the spikes based on the histogram of the spikes, the main idea is to remove the random spikes, and keep the spikes with high frequency, the threshold is set based on the histogram of the spikes time.
        Parameters
        ----------
        SpikesTime: array
            the time of the spikes
        SpikesID: array
            the channel id of the spikes
        Returns
        -------
        SpikesChIDs_all: array
            the channel id of the spikes after filtering
        SpikesTimes_all: array
            the time of the spikes after filtering
        """
        binsDistr = np.arange(0, max(SpikesTime), 0.1)  # fixed bin size
        averageEventRateHist, averageEventRateBinsEdge = np.histogram(SpikesTime, bins=binsDistr, normed=False,weights=None, density=None)
        Hist_filter = [i for i in range(len(averageEventRateHist)) if
                       averageEventRateHist[i] >= Synchrony_Threshold * np.max(
                           [count for count in averageEventRateHist if
                            count < Firing_Electrodes_Count_Threshold])]
        Hist_filter = np.unique(Hist_filter)
        Remain_Time = []
        for ID in Hist_filter:
            time_bin = [time for time in SpikesTime if
                        time >= averageEventRateBinsEdge[ID] and time <= averageEventRateBinsEdge[ID + 1]]
            Remain_Time.extend(time_bin)
        Remain_ID = [i for i in range(len(SpikesTime)) if SpikesTime[i] in Remain_Time]
        SpikeChIDs_all = [SpikesID[i] for i in Remain_ID]
        SpikeTimes_all = [SpikesTime[i] for i in Remain_ID]
        return np.array(SpikeChIDs_all), np.array(SpikeTimes_all)

    def MBR_filter(self,SamplingRate = 0,SpikeTimes = None,SpikeChIDs = None,ISImax=100):
        '''
        The function is used to filter the spikes based on the MBR(Mean Burst Rate) of the spikes, the main idea is to remove the random spikes, and keep the spikes with high frequency, the threshold is set based on the MBR of the spikes time.
        :param SamplingRate: float
            The sampling rate of the recording reading from the bxr file
        :param SpikeTimes: array
            the time of the spikes in seconds unit
        :param SpikeChIDs: array
            the channel ids of the spikes
        :param ISImax: float
            the maximum acceptable ISI(inter-spike-interval) in ms
        Returns
        -------
        SpikesChIDs: array
            the channel id of the spikes after filtering
        SpikesTimes: array
            the time of the spikes after filtering

        '''
        spikes_unique = np.unique(SpikeChIDs)
        recordingLength = np.ceil(max(SpikeTimes))
        keep_channel_IDs = []
        for indST in spikes_unique:
            indOfST = []
            for i in range(len(SpikeChIDs)):
                if indST == SpikeChIDs[i]:
                    indOfST.append(i)
            ST = SpikeTimes[indOfST] * 1000 / SamplingRate
            ISI = np.diff(np.append(ST, ST[-1] + 1e6))  # add a fake element at the end to close last potential burst
            indISI = np.where(ISI > ISImax)[0]
            LindISI = len(indISI)
            numburst = 0
            for k in range(LindISI):
                numburst += 1
            if 60.0 * numburst / recordingLength >=mbrMin and 60.0 * numburst / recordingLength <= mbrMax: # number of burst / min
                keep_channel_IDs.append(indST)
        new_ids = [i for i in range(len(SpikeChIDs)) if SpikeChIDs[i] in keep_channel_IDs]
        return SpikeChIDs[new_ids],SpikeTimes[new_ids]


    def AnalyzeExp(self,expFile = None):
        """
        The denoising codes main function to filter the spikes based on the histogram and MBR of the spikes time. The main idea is to remove the random spikes, and keep the spikes with high frequency, the threshold is set based on the histogram and MBR of the spikes time. The results are saved in .npy file.
        Parameters
        ----------
        expFile: string
            the name of the file
            Returns
            -------
            SpikesChId_filter: array
                the channel id of the spikes after filtering
                SpikesTimes_filter: array
                the time of the spikes after filtering

        """
        if os.path.exists(self.srcfilepath + expFile[:-4] + '_denoised_SpikesChIDs' + '.npy') and os.path.exists(self.srcfilepath + expFile[:-4] + '_denoised_SpikesTimes' + '.npy'):
            SpikesChId_filter = np.load(self.srcfilepath + expFile[:-4] + '_denoised_SpikesChIDs' + '.npy')
            SpikesTimes_filter = np.load(self.srcfilepath + expFile[:-4] + '_denoised_SpikesTimes' + '.npy')
        else:
            filehdf5_bxr = h5py.File(self.srcfilepath + expFile, 'r')  # read LFPs bxr files
            ##########################################################
            try:
                # Read the data from the file generated from the BrainWave 4.0
                NRecFrames = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["NRecFrames"])[0]
                samplingRate = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["SamplingRate"])[0]
                recordingLength = NRecFrames / samplingRate  # recording duraton in [s]
                SpikesTimes = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["SpikeTimes"]) / samplingRate
                SpikesChId = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["SpikeChIDs"])
                maxVolt = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]['MaxVolt'])[0]
                minVolt = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]['MinVolt'])[0]
                stepVolt = (maxVolt - minVolt) / 4096
                signalInversion = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["SignalInversion"])[0]
            except:
                # Read the data from file generated from the BrainWave 5.0
                import json
                ExperimentSettings = filehdf5_bxr["ExperimentSettings"][:][0]
                settings = json.loads(ExperimentSettings)
                data = filehdf5_bxr[list(filehdf5_bxr.keys())[-1]]
                ##########################################################
                samplingRate = settings['TimeConverter']['FrameRate']

                SpikesChId = np.asarray(data['SpikeChIdxs'][:])
                SpikesTimes = np.asarray(data['SpikeTimes'][:]) / samplingRate
                valueConverter = settings['ValueConverter']
                maxVolt = valueConverter['MaxAnalogValue']
                minVolt = valueConverter['MinAnalogValue']
                QLevel = valueConverter[
                    'MaxDigitalValue']  # quantized levels corresponds to 2^num of bit to encode the signal
                # print('QLevel',QLevel)
                fromQLevelToUVolt = (maxVolt - minVolt) / QLevel
                stepVolt = (maxVolt - minVolt) / maxVolt
                signalInversion = maxVolt / QLevel
                chIdx = data['StoredChIdxs'][:]
                ID = np.asarray([i for i in range(QLevel + 1)])
                Row = np.asarray([int(i / 64) + 1 for i in ID])
                Col = np.asarray([i % 64 + 1 for i in ID])
                MeaChs2ChIDsVector = {'Row': Row, 'Col': Col, 'ID': ID}
                MeaStreams = [[int(j / 64) + 1, int(j % 64) + 1] for j in ID]
                SpikesTimes = [float(i) for i in  SpikesTimes]
                recordingLength =  SpikesTimes[-1]  # recording duraton in [s]
            recordingLength = max(SpikesTimes)
            SpikesChId_uni = np.unique(SpikesChId)

            keep_id = [i for i in SpikesChId_uni if len([j for j in SpikesChId if j == i])/recordingLength >= self.SpikesMin and len([j for j in SpikesChId if j == i])/recordingLength<= self.SpikesMax]
            filter_ID = [i for i in range(len(SpikesChId)) if SpikesChId[i] in keep_id]
            SpikesChId_filter = np.array(SpikesChId)[filter_ID]
            SpikesTimes_filter = np.array(SpikesTimes)[filter_ID]

            ######################################hist filter
            SpikesChId_hist, SpikesTimes_hist = self.histogram_filter(SpikesTime=SpikesTimes_filter,SpikesID=SpikesChId_filter)

            SpikesChId_last, SpikesTimes_last = self.MBR_filter(SamplingRate = samplingRate,SpikeTimes = SpikesTimes_hist,SpikeChIDs = SpikesChId_hist)

            ##############save the new results in .npy file
            np.save(self.srcfilepath + expFile[:-4] + '_denoised_SpikesChIDs', SpikesChId_last)
            np.save(self.srcfilepath + expFile[:-4] + '_denoised_SpikesTimes', SpikesTimes_last)

        print('Denoise is Done')
        return np.array(SpikesChId_filter), np.array(SpikesTimes_filter)