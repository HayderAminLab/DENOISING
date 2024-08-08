# -*- coding: utf-8 -*-
"""
Created on Dec 12 2022
@author:  BIONICS_LAB
@company: DZNE
"""

import h5py
import numpy as np
import os
"""need change these parameters based on different datasets"""
threshold_low = 50   #uV
threshold_high = 4000  #uV
LFPMax= 30               # (Hz) - maximum acceptable LFP Rate(LFP events/min)
LFPMin= 0.1                # (Hz) - minumum acceptable LFP Rate(LFP events/min)
Firing_Electrodes_Count_Threshold = 6000 # the maximum acceptable active channel numbers in Frames
Duration_threshold = 10 # duration thrshold(ms)
Synchrony_Threshold = 0.1 #Those averageEventRateHists less than Synchrony_Threshold of all the events count in hist are regarded as random Events

class LFP_DENOISING_Function:
    def __init__(self,srcfilepath):
        self.srcfilepath = srcfilepath  # main path

    def get_filename(self, filepath, filetype):
        """get the reference file name and path"""
        filename = []
        for root, dirs, files in os.walk(filepath):
            for i in files:
                if filetype in i:
                    filename.append(i)
        return filename

    def denoising(self,lfpChId_raw=None, lfpTimes_raw=None, LfpForms=None,ChsGroups = None,stepVolt=0,signalInversion =0,threshold_low = 210,threshold_high = 2000,samplingRate=None):
        """
        1.denoising based on wavefroms amplitudes threshold
        2.denoising based on duration
        3.denoising based on the active channels at the same time
        Parameters
        ----------
        :param lfpChId_raw: numpy array,
                the channel IDs of the LFPs read from the .bxr file
        :param lfpTimes_raw: numpy array,
                the timepoints of the LFPs read from the .bxr file
        :param LfpForms: numpy array,
                the waveforms of the LFPs read from the .bxr file
        :param ChsGroups: dict,
                the information of the clusters
        :param stepVolt: float,
                the step of the voltage
        :param signalInversion: float,
                the signal inversion
        :param threshold_low: float,
                the low threshold of the waveforms amplitudes
        :param threshold_high: float,
                the high threshold of the waveforms amplitudes
        :param samplingRate: float,
                the sampling rate
        Returns
        ----------------
        :return:
                lfpChId_new: numpy array,
                    the channel IDs of the LFPs after denoising
                lfpTimes_new: numpy array,
                    the timepoints of the LFPs after denoising
                tempWaveForms_new: numpy array,
                    the waveforms of the LFPs after denoising

        """
        numLFPs = len(lfpChId_raw)  # extract the total number of detected LFPs
        LFPLength = len(LfpForms) // numLFPs  # extract the LFP length
        tempWaveForms = LfpForms
        tempWaveForms = np.array(tempWaveForms)
        tempWaveForms = tempWaveForms.reshape(numLFPs, LFPLength)[:]
        cluster_ID = []
        for i in range(len(ChsGroups['Chs'])):
            for j in range(len(ChsGroups['Chs'][i])):
                k = (ChsGroups['Chs'][i][j][0] - 1) * 64 + (ChsGroups['Chs'][i][j][1] - 1)
                cluster_ID.append(k)
        Id = [i for i in range(len(lfpChId_raw)) if lfpChId_raw[i] in cluster_ID]
        lfpChId = [lfpChId_raw[i] for i in Id]
        lfpTimes = [lfpTimes_raw[i] for i in Id]
        tempWaveForms = [tempWaveForms[i] for i in Id]
        ##########################################Duration filter
        #all the duration
        Duration_all = [np.nonzero(el)[0][-1]*1000/samplingRate for el in tempWaveForms]
        #duration threshold(The duration for each event>= mean of all the duration)
        Normal_ID_duration = [j for j in range(len(lfpTimes)) if Duration_all[j] >= Duration_threshold]
        lfpChId_duration_filter = [lfpChId[i] for i in Normal_ID_duration]
        lfpTimes_duration_filter = [lfpTimes[i] for i in Normal_ID_duration]
        tempWaveForms_duration_filter = [tempWaveForms[i] for i in Normal_ID_duration]
        ##############denoise low amplitudes channels
        Normal_ID = [j for j in range(len(lfpChId_duration_filter)) if (max(list(abs((tempWaveForms[Normal_ID_duration[j]][:np.nonzero(tempWaveForms[Normal_ID_duration[j]])[0][-1]] - (4096.0 / 2)) * stepVolt * signalInversion)))>=threshold_low and max(list(abs((tempWaveForms[Normal_ID_duration[j]][:np.nonzero(tempWaveForms[Normal_ID_duration[j]])[0][-1]] - (4096.0 / 2)) * stepVolt * signalInversion)))<=threshold_high)]
        lfpChId_new = [lfpChId_duration_filter[i] for i in Normal_ID]
        lfpTimes_new = [lfpTimes_duration_filter[i] for i in Normal_ID]
        tempWaveForms_new = [tempWaveForms_duration_filter[i] for i in Normal_ID]
        return lfpChId_new,lfpTimes_new,tempWaveForms_new

    def MFR_denoise(self,lfpChId_raw=None, lfpTimes_raw=None,recordingLength = 0):
        """delete the event that actve channels appear at the same time over Firing_Electrodes_Count_Threshold"""
        lfpChId_unique = np.unique(lfpChId_raw)
        filter_ID = [i for i in lfpChId_unique if list(lfpChId_raw).count(i)*60/recordingLength >= LFPMin and list(lfpChId_raw).count(i)*60/recordingLength <= LFPMax]
        lfpChId_filter_index = [i for i in range(len(lfpChId_raw)) if lfpChId_raw[i] in filter_ID]
        lfpChId_filter = [lfpChId_raw[i] for i in lfpChId_filter_index]
        lfpTimes_filter = [lfpTimes_raw[i] for i in lfpChId_filter_index]
        ################deleterd all channels active at the same time
        lfpTimes_filter_uni = np.unique(lfpTimes_filter)
        noise_time = [i for i in lfpTimes_filter_uni if list(lfpTimes_filter).count(i) > Firing_Electrodes_Count_Threshold]
        lfpTimes_filter_remain_index = [i for i in range(len(lfpTimes_filter)) if lfpTimes_filter[i] not in noise_time]
        lfpChId_filter = [lfpChId_filter[i] for i in lfpTimes_filter_remain_index]
        lfpTimes_filter = [lfpTimes_filter[i] for i in lfpTimes_filter_remain_index]
        return lfpChId_filter, lfpTimes_filter
    def duplicate_channels(self,ChsGroups = None):
        """
        Find the wrong culster name and duplicate channels in difference clusters and empty cluster in the ChsGroups dict.
        Parameters
        ----------
        :param ChsGroups: dict,
                the information of the clusters
        Returns
        ----------------
        :return:
                None
        """
        # print(ChsGroups['Name'],self.clusters)
        for clu in range(len(ChsGroups['Name'])):
            if ChsGroups['Name'][clu] in self.clusters:
                channels  = ChsGroups['Chs'][clu]

                if len(channels)>0:
                    continue
                elif len(channels)< 5:
                    print(ChsGroups['Name'][clu] + ' only have '+str(len(channels)) + ' channels!')
                else:
                    print('Error! '+ ChsGroups['Name'][clu] +' is empty!')
            else:
                print('Error! '+ ChsGroups['Name'][clu] +' not match the clusters. Please check the name of cluster.')

        for i in range(len(ChsGroups['Chs'])):
            for j in range(i+1,len(ChsGroups['Chs'])):
                com = [i for i in ChsGroups['Chs'][i] if i in ChsGroups['Chs'][j]]
                if len(com)>0:
                    for id in com:
                        print('Error! '+'Channel: '+str(id) + 'belong to cluster:' + ChsGroups['Name'][i] + ' and '+ChsGroups['Name'][j])
                else:
                    continue

    def histogram_filter(self,recordingLength=0,LFPTime=None,lfpChId =None,LfpForms=None):
        '''
        The function is used to filter the LFPs based on the histogram of the LFPs, the main idea is to remove the random LFPs, and keep the LFPs with high frequency, the threshold is set based on the histogram of the LFPs time.
         Parameters
            ----------
        :param recordingLength: float,
                the length of the recording
        :param LFPTime: numpy array,
                the timepoints of the LFPs read from the .bxr file
        :param lfpChId: numpy array,
                the channel IDs of the LFPs read from the .bxr file
        :param LfpForms: numpy array,
                the waveforms of the LFPs read from the .bxr file
        Returns
        ----------------
            :return:
                    lfpChId_remain: numpy array,
                        the channel IDs of the LFPs after filtering
                    LFPTime_remain: numpy array,
                        the timepoints of the LFPs after filtering
                    LfpForms_remian: numpy array,
                        the waveforms of the LFPs after filtering
        '''
        binsDistr = np.arange(0, recordingLength, 0.15)  # fixed bin size
        averageEventRateHist, averageEventRateBinsEdge = np.histogram(LFPTime, bins=binsDistr,normed=False, weights=None, density=None)
        Hist_filter = [i for i in range(len(averageEventRateHist)) if averageEventRateHist[i] >= Synchrony_Threshold * np.max([count for count in averageEventRateHist if count < Firing_Electrodes_Count_Threshold])]
        Hist_filter = np.unique(Hist_filter)
        Remain_Time = []
        for ID in Hist_filter:
            time_bin = [time for time in LFPTime if time>= averageEventRateBinsEdge[ID] and time <= averageEventRateBinsEdge[ID+1]]
            Remain_Time.extend(time_bin)
        Remain_ID = [i for i in range(len(LFPTime)) if LFPTime[i] in Remain_Time]
        lfpChId_remain = [lfpChId[i] for i in Remain_ID]
        LFPTime_remain = [LFPTime[i] for i in Remain_ID]
        LfpForms_remian = [LfpForms[i] for i in Remain_ID]
        return lfpChId_remain,LFPTime_remain,LfpForms_remian

    def Class_ID_detect(self,LfpChIDs=None,ChsGroups = None,expFile=None,desfilepath =None):
        """Classify the LFPs based on the clusters
        Parameters
        ----------
        :param LfpChIDs: numpy array,
                the channel IDs of the LFPs read from the .bxr file
        :param ChsGroups: dict,
                the information of the clusters
        :param expFile: str,
                the name of the .bxr file
        :param desfilepath: str,
                the path of the .bxr file
        Returns
        ----------------
        :return:
                New_ID: numpy array,
                    the new channel IDs of the LFPs after classification
                New_cluster_count: numpy array,
                    the number of the LFPs in each cluster
                clusters_names: list,
                    the names of the clusters
                number_list: numpy array,
                    the number of the LFPs in each cluster
        """
        colorMapTitle = expFile[:-4] + '_Rater_Plots_Sorted'
        file_path = desfilepath + colorMapTitle + ".txt"
        f = open(file_path, "a")
        f.seek(0)
        f.truncate()
        Class_LFPs_ID = np.unique(LfpChIDs)
        cluster_ids_all = []
        clusters_names = []
        number_list = [-1]
        for i in range(len(ChsGroups['Chs'])):
            clusters_names.append(ChsGroups['Name'][i])
            cluster_ids = [(ChsGroups['Chs'][i][j][0] - 1) * 64 + (ChsGroups['Chs'][i][j][1] - 1) for j in range(len(ChsGroups['Chs'][i]))]
            cluster_ids_all.extend(cluster_ids)
            msg = ChsGroups['Name'][i] + ' goes from ' + str(number_list[-1]+1) + ' to ' + str(len(cluster_ids) + number_list[-1] + 1)
            number_list.append(len(cluster_ids) + number_list[-1] + 1)
            f.write(msg + '\n')
        No_cluter = [i for i in Class_LFPs_ID if i not in cluster_ids_all]
        msg = 'Not in Cluster' + ' goes from ' + str(number_list[-1] + 1) + ' to ' + str(len(No_cluter) + number_list[-1] + 1)
        number_list.append(len(No_cluter) + number_list[-1] + 1)
        clusters_names.append('Not in Cluster')
        f.write(msg + '\n')
        f.close()

        all_ids = cluster_ids_all + No_cluter
        New_ID = [list(all_ids).index(i) for i in LfpChIDs]
        New_cluster_count = [int((number_list[:-1][i]+number_list[1:][i])/2) for i in range(len(number_list[:-1]))]
        number_list[0] = 0
        return New_ID,New_cluster_count,clusters_names,number_list

    def AnalyzeExp(self,expFile = None):
        """ The denoising codes main function based on the .bxr file and the parameters set in the beginning of the script.

        Parameters
        ----------
        :expFile: str,
                the name of the .bxr file
        Returns
        ----------------
        :return:
                lfpChId_last: numpy array,
                    the channel IDs of the LFPs after denoising
                lfpTimes_last: numpy array,
                    the timepoints of the LFPs after denoising
                LfpForms_last: numpy array,
                    the waveforms of the LFPs after denoising
                All results are saved in the .npy files.
        """
        filehdf5_bxr = h5py.File(self.srcfilepath + expFile, 'r')  # read LFPs bxr files
        try:
            NRecFrames = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["NRecFrames"])[0]
            samplingRate = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["SamplingRate"])[0]
            # recordingLength = NRecFrames / samplingRate  # recording duraton in [s]
            ChsGroups = np.asarray(filehdf5_bxr["3BUserInfo"]["ChsGroups"])
            LfpForms = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["LfpForms"])
            lfpChId = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["LfpChIDs"])
            lfpTimes = np.asarray(filehdf5_bxr["3BResults"]["3BChEvents"]["LfpTimes"]) / samplingRate
            maxVolt = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]['MaxVolt'])[0]
            minVolt = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]['MinVolt'])[0]
            stepVolt = (maxVolt - minVolt) / 4096
            signalInversion = np.asarray(filehdf5_bxr["3BRecInfo"]["3BRecVars"]["SignalInversion"])[0]
        except:
            import json
            ExperimentSettings = filehdf5_bxr["ExperimentSettings"][:][0]
            settings = json.loads(ExperimentSettings)
            data = filehdf5_bxr[list(filehdf5_bxr.keys())[-1]]
            ##########################################################
            samplingRate = settings['TimeConverter']['FrameRate']
            NRecFrames = data['FpTOC'][:][-1]
            # recordingLength = NRecFrames / samplingRate  # recording duraton in [s]
            Names = []
            Chs = []
            for i in range(len(settings['MeaPlate']['UnitGroups'])):
                Group_Information = settings['MeaPlate']['UnitGroups'][i]
                Names.append(Group_Information['UserDefinedName'])
                Chs.append([[int(j / 64), int(j % 64)] for j in Group_Information['UnitIndexes']])
            ChsGroups = {'Name': Names, 'Chs': Chs}
            LfpForms = np.asarray(data['FpForms'][:])
            lfpChId = np.asarray(data['FpChIdxs'][:])
            lfpTimes = np.asarray(data['FpTimes'][:]) / samplingRate
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
        if type(ChsGroups['Name'][0]) != np.str:
            ChsGroups['Name'] = [i.decode("utf-8") for i in ChsGroups['Name']]

        recordingLength = max(lfpTimes) # recording duraton in [s]
        self.clusters = list(ChsGroups['Name'])
        if os.path.exists(self.srcfilepath + expFile[:-4] + '_denoised_LfpChIDs' + '.npy') and os.path.exists(
                self.srcfilepath + expFile[:-4] + '_denoised_LfpTimes' + '.npy') and os.path.exists(self.srcfilepath + expFile[:-4] + '_denoised_LfpForms' + '.npy'):
            lfpChId_last = np.load(self.srcfilepath + expFile[:-4] + '_denoised_LfpChIDs' + '.npy')
            lfpTimes_last = np.load(self.srcfilepath + expFile[:-4] + '_denoised_LfpTimes' + '.npy')
            LfpForms_last = np.load(self.srcfilepath + expFile[:-4] + '_denoised_LfpForms' + '.npy')
        else:
            ##################First highlight the duplicate cahnnels,empty cluster and Wrong name of cluster
            self.duplicate_channels(ChsGroups = ChsGroups)
            ##################Second lfpChId and lfpTimes denoise
            lfpChId_denoising, lfpTimes_denoising,tempWaveForms = self.denoising(lfpChId_raw=lfpChId, lfpTimes_raw=lfpTimes,
                                                                LfpForms=LfpForms, ChsGroups=ChsGroups,
                                                                stepVolt=stepVolt, signalInversion=signalInversion,threshold_low = threshold_low,threshold_high = threshold_high,samplingRate=samplingRate)
            #threshold: for the waveforms in the range(-threshold,+threshold) detected as noise or unactive channels
            ##################Third:Mean LFP rate threshold
            lfpChId_filter, lfpTimes_filter = self.MFR_denoise(lfpChId_raw=lfpChId_denoising, lfpTimes_raw=lfpTimes_denoising,recordingLength = recordingLength)
            #######################get the filter waveforms
            ID_filter = [i for i in range(len(lfpChId_denoising)) if lfpChId_denoising[i] in lfpChId_filter]
            LfpForms_filter = np.asarray([tempWaveForms[i] for i in ID_filter])
            ######################################hist filter
            lfpChId_last, lfpTimes_last ,LfpForms_last = self.histogram_filter(recordingLength=recordingLength, LFPTime=lfpTimes_filter,lfpChId=lfpChId_filter,LfpForms=LfpForms_filter)
            LfpForms_last = np.asarray(LfpForms_last).reshape(-1)
            ##############save the new results in .npy file
            np.save(self.srcfilepath + expFile[:-4] + '_denoised_LfpChIDs', lfpChId_last)
            np.save(self.srcfilepath + expFile[:-4] + '_denoised_LfpTimes', lfpTimes_last)
            np.save(self.srcfilepath + expFile[:-4] + '_denoised_LfpForms', LfpForms_last)
        print('Denoise is Done')
        return np.asarray(lfpChId_last), np.asarray(lfpTimes_last),np.asarray(LfpForms_last)